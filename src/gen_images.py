import time
import http.client
import os, re, json, base64, argparse
from pathlib import Path
from typing import Optional, Tuple, List, Dict, Any
from urllib.parse import urlparse


# ----------------------------
# Markdown split
# ----------------------------
def split_markdown(md: str) -> Tuple[str, List[Tuple[str, str]]]:
    parts = re.split(r'(?m)^#\s+', md)
    if len(parts) <= 1:
        raise ValueError("未找到一级标题（'# '）。请把4个章节写成一级标题。")

    intro = parts[0].strip()
    chapters = []
    for p in parts[1:]:
        lines = p.splitlines()
        title = (lines[0] if lines else "").strip()
        body = "\n".join(lines[1:]).strip() if len(lines) > 1 else ""
        chapters.append((title, body))
    return intro, chapters


# ----------------------------
# Run dir + tutorial md discovery
# ----------------------------
def resolve_run_dir(cli_run_dir: Optional[str]) -> Path:
    run_dir = cli_run_dir
    if not run_dir:
        run_dir = input(r"请输入run目录（例如 .\data\runs\20260109_153231 ）: ").strip().strip('"')
    p = Path(run_dir).expanduser().resolve()
    if not p.exists() or not p.is_dir():
        raise SystemExit(f"run目录不存在或不是目录：{p}")
    return p

def pick_tutorial_md(run_dir: Path) -> Path:
    mds = sorted(run_dir.glob("Final_Tutorial_*.md"))
    if not mds:
        raise SystemExit(f"在目录中未找到 Final_Tutorial_*.md：{run_dir}")

    if len(mds) == 1:
        return mds[0]

    print("找到多个教程文件，请选择：")
    for i, f in enumerate(mds, 1):
        print(f"{i}. {f.name}")
    idx = input("输入序号（回车默认1）: ").strip()
    if not idx:
        return mds[0]
    return mds[int(idx) - 1]


# ----------------------------
# Helpers
# ----------------------------
def get_message_text(resp: dict) -> str:
    choices = resp.get("choices") or []
    if not choices:
        return ""
    msg = (choices[0] or {}).get("message") or {}
    content = msg.get("content")

    if isinstance(content, str):
        return content.strip()

    # 兼容 content parts 数组
    if isinstance(content, list):
        texts = []
        for part in content:
            if isinstance(part, dict) and part.get("type") == "text":
                t = part.get("text")
                if isinstance(t, str):
                    texts.append(t)
        return "\n".join(texts).strip()

    return ""

def safe_name(s: str, max_len: int = 40) -> str:
    s = re.sub(r"[^0-9A-Za-z\u4e00-\u9fa5]+", "_", s).strip("_")
    return s[:max_len] if len(s) > max_len else s

def strip_code_fence(text: str) -> str:
    if not text:
        return text
    t = text.strip()
    t = re.sub(r"^\s*```[a-zA-Z0-9_-]*\s*", "", t)
    t = re.sub(r"\s*```\s*$", "", t)
    return t.strip()

def try_parse_json_object(raw: str) -> Optional[Dict[str, Any]]:
    if not raw:
        return None
    cleaned = strip_code_fence(raw)

    try:
        obj = json.loads(cleaned)
        return obj if isinstance(obj, dict) else None
    except:
        pass

    m = re.search(r"\{.*\}", cleaned, re.S)
    if not m:
        return None
    try:
        obj = json.loads(m.group(0))
        return obj if isinstance(obj, dict) else None
    except:
        return None


# ----------------------------
# VectorEngine chat(completions) - IMPORTANT: uses maxTokens (camelCase)
# ----------------------------
def ve_chat_completions(base_url, api_key, model, messages, max_tokens, temperature,
                        timeout=120, retries=3) -> Dict[str, Any]:
    u = urlparse(base_url)
    host = u.netloc
    base_path = u.path.rstrip("/")
    path = f"{base_path}/chat/completions"

    payload = {
        "model": model,
        "maxTokens": max_tokens,
        "temperature": temperature,
        "messages": messages,
    }

    last_exc = None
    for attempt in range(1, retries + 1):
        conn = None
        try:
            conn = http.client.HTTPSConnection(host, timeout=timeout)
            conn.request(
                "POST",
                path,
                body=json.dumps(payload),
                headers={
                    "Authorization": f"Bearer {api_key}",
                    "Content-Type": "application/json",
                    "Connection": "close",
                    "Accept-Encoding": "identity",
                },
            )
            resp = conn.getresponse()
            try:
                body = resp.read()
            except http.client.IncompleteRead as e:
                body = e.partial  # 关键：拿到已读到的部分

            data = body.decode("utf-8", errors="replace")

            if resp.status >= 400:
                raise RuntimeError(f"HTTP {resp.status} {resp.reason}\n{data}")

            return json.loads(data)

        except (http.client.IncompleteRead, ValueError, TimeoutError, OSError, ConnectionError) as e:
            last_exc = e
            if attempt == retries:
                raise
            time.sleep(1.2 * attempt)
        finally:
            if conn:
                try:
                    conn.close()
                except:
                    pass

    raise last_exc


def make_text_msg(role: str, text: str) -> Dict[str, Any]:
    # 按示例的“content数组 + type=text”更稳 [file:48]
    return {
        "role": role,
        "content": [{"type": "text", "text": text}]
    }


# ----------------------------
# Prompt extraction
# ----------------------------
def sanitize_md(snippet: str) -> str:
    # 去掉 fenced code block，避免模型进入“续写代码/教程”模式
    snippet = re.sub(r"```.*?```", "[代码块省略]", snippet, flags=re.S)
    # 弱化列表符号
    snippet = re.sub(r"(?m)^\s*[*\-+]\s+", "", snippet)
    snippet = snippet.replace("`", "")
    return snippet

def make_prompt(
    base_url: str,
    api_key: str,
    extract_model: str,
    title: str,
    text: str,
    kind: str,
    debug_dir: Path,
    idx: int
) -> Dict[str, str]:
    snippet = (text or "").strip()[:1200].replace("```", "`").replace("\r\n", "\n")
    snippet = sanitize_md(snippet)
    
    sys = (
    "你是科学教程的可视化设计师。\n"
    "你必须只输出一个JSON对象，不要代码块、不加解释。\n"
    "JSON结构固定为："
    "{\"title\": string, \"caption\": string, \"prompt\": string, \"negative_prompt\": string}\n"
    "风格硬性要求：白底、扁平矢量风、结构图/流程图/示意图；清晰线条与少量配色区分；不要3D、不要等距、不要拟物渲染。\n"
    "文字硬性要求：图中标签以简体中文为主（zh-cn）；非必要不要用英文；不要长段文字。\n"
    )
    user = f"""请为下面教程内容设计 1 张“串联知识点的结构图/流程图”，用于插入教程帮助读者理解。

    输出要求：
    - 只输出JSON对象：{{"title": "...", "caption":"...", "prompt":"...", "negative_prompt":"..."}}
    - prompt 要描述“图里有哪些模块/箭头/分组/输入输出”，而不是描述“海报/封面/氛围”
    - 图中标签以简体中文为主（zh-cn）；非必要不要出现英文
    - negative_prompt 请写清楚要避免：3D、照片、人物、复杂背景、大段文字

    【输入内容开始】
    标题：{title}
    内容摘录：
    {snippet}
    【输入内容结束】

    现在请只输出一个 JSON 对象（不要 Markdown、不要解释、不要复述输入内容）。
        """

    # 允许重试2次
    for attempt in [1, 2]:
        resp = ve_chat_completions(
            base_url=base_url,
            api_key=api_key,
            model=extract_model,
            messages=[make_text_msg("system", sys), make_text_msg("user", user)],
            max_tokens=1200,
            temperature=0.0,
        )
        raw = get_message_text(resp)

        raw = raw.strip()

        (debug_dir / f"debug_extract_raw_{idx:02d}_attempt{attempt}.txt").write_text(raw, encoding="utf-8")

        obj = try_parse_json_object(raw)
        if obj and isinstance(obj.get("prompt"), str):
            prompt = " ".join(obj["prompt"].splitlines()).strip()
            return {
                    "title": str(obj.get("title") or title),
                    "caption": str(obj.get("caption") or ""),
                    "prompt": prompt,
                    "negative_prompt": str(obj.get("negative_prompt") or "")
        }

        sys = (
            "严格只输出JSON对象，不要代码块。"
            "prompt 必须是单行，<=200字。"
            "不允许输出除JSON以外任何字符。"
        )

    raise ValueError(
        f"提取模型未能返回可解析JSON。请查看：{debug_dir}\\debug_extract_raw_{idx:02d}_attempt*.txt"
    )


# ----------------------------
# Image base64 extraction
# ----------------------------
def extract_base64_images_from_text(text: str) -> List[str]:
    b64s = []
    if not text:
        return b64s

    for m in re.finditer(r"data:image\/png;base64,([A-Za-z0-9+/=]+)", text):
        b64s.append(m.group(1))
    for m in re.finditer(r"data:image\/jpeg;base64,([A-Za-z0-9+/=]+)", text):
        b64s.append(m.group(1))

    # 如果内容本身是JSON
    try:
        obj = json.loads(text)
        if isinstance(obj, dict):
            for k in ["image_base64", "b64_json", "base64", "image"]:
                v = obj.get(k)
                if isinstance(v, str) and len(v) > 200:
                    b64s.append(v)
            if isinstance(obj.get("images"), list):
                for it in obj["images"]:
                    if isinstance(it, dict):
                        v = it.get("b64_json") or it.get("image_base64") or it.get("base64")
                        if isinstance(v, str) and len(v) > 200:
                            b64s.append(v)
    except:
        pass

    return b64s

def save_first_image(b64_list: List[str], out_path: Path) -> bool:
    if not b64_list:
        return False
    out_path.write_bytes(base64.b64decode(b64_list[0]))
    return True



# ----------------------------
# Main
# ----------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run_dir", default=None, help="运行目录；不传则交互输入")
    ap.add_argument("--base_url", default=os.getenv("VECTORENGINE_BASE_URL", "https://api.vectorengine.ai/v1"))
    ap.add_argument("--api_key", default=os.getenv("VECTORENGINE_API_KEY"))
    ap.add_argument("--extract_model", default="gemini-3-flash-preview")
    ap.add_argument("--image_model", default="gemini-3-pro-image-preview")
    ap.add_argument("--max_chapters", type=int, default=4)
    args = ap.parse_args()

    if not args.api_key:
        raise SystemExit("请设置环境变量 VECTORENGINE_API_KEY，或用 --api_key 传入。")

    run_dir = resolve_run_dir(args.run_dir)
    tutorial_path = pick_tutorial_md(run_dir)

    outdir = run_dir / "images"
    outdir.mkdir(parents=True, exist_ok=True)

    debug_dir = run_dir / "debug"
    debug_dir.mkdir(parents=True, exist_ok=True)

    md = tutorial_path.read_text(encoding="utf-8")
    intro, chapters = split_markdown(md)
    chapters = chapters[: args.max_chapters]

    items = [("全文开头", intro, "cover")]
    for i, (t, b) in enumerate(chapters, 1):
        items.append((f"第{i}章：{t}", b, "chapter"))

    # 1) 生成 prompts
    prompts = []
    for idx, (title, text, kind) in enumerate(items, 1):
        pj = make_prompt(
            base_url=args.base_url,
            api_key=args.api_key,
            extract_model=args.extract_model,
            title=title,
            text=text,
            kind=kind,
            debug_dir=debug_dir,
            idx=idx
        )
        prompts.append({"title": pj["title"], "prompt": pj["prompt"]})

    # 保存 prompts（自动生成，不需要你手工准备）
    prompt_json_path = run_dir / "image_generation_prompt.json"
    prompt_json_path.write_text(json.dumps({
        "run_dir": str(run_dir),
        "tutorial_md": tutorial_path.name,
        "extract_model": args.extract_model,
        "image_model": args.image_model,
        "items": prompts,
    }, ensure_ascii=False, indent=2), encoding="utf-8")

    # 2) 出图
    for idx, item in enumerate(prompts, 1):
        title = item["title"]
        prompt = (
            item["prompt"].rstrip()
            + "\n\n硬性要求：图中标签以简体中文为主（zh-cn）；白底扁平矢量结构图风格。"
        )

        neg = item.get("negative_prompt", "").strip()
        if neg:
            prompt += "\n\nAvoid: " + neg

        resp = ve_chat_completions(
            base_url=args.base_url,
            api_key=args.api_key,
            model=args.image_model,
            messages=[make_text_msg("system", "你是图片生成模型。请根据用户提示生成一张PNG图片。"),
                      make_text_msg("user", prompt)],
            max_tokens=4096,
            temperature=0.7,
        )

        raw = get_message_text(resp)

        b64s = extract_base64_images_from_text(raw)

        out_path = outdir / f"{idx:02d}_{safe_name(title)}.png"
        ok = save_first_image(b64s, out_path)

        if ok:
            print(f"[OK] {title} -> {out_path}")
        else:
            (outdir / f"{idx:02d}_raw.txt").write_text(str(raw), encoding="utf-8")
            print(f"[WARN] {title} 未解析到图片base64，已输出 raw.txt：{outdir}")

    print(f"\n完成。prompts已保存：{prompt_json_path}")
    print(f"debug目录：{debug_dir}")
    print(f"图片输出目录：{outdir}")


if __name__ == "__main__":
    main()
