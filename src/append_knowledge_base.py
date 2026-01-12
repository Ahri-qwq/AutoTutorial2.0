import os
import glob
import uuid
import time
import shutil
from typing import List

try:
    from pypdf import PdfReader
except ImportError:
    PdfReader = None


import chromadb
from chromadb import Documents, EmbeddingFunction, Embeddings

import dashscope
from dashscope import TextEmbedding
import yaml

try:
    from docx import Document
except ImportError:
    Document = None

# ================= 路径配置（与 build_knowledge_base.py 保持一致） =================

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(CURRENT_DIR)

CONFIG_PATH = os.path.join(ROOT_DIR, "config.yaml")
SOURCE_DIR = os.path.join(ROOT_DIR, "data", "knowledge_add")   # ★ 只读增量目录
DB_PATH = os.path.join(ROOT_DIR, "data", "chroma_db")
COLLECTION_NAME = "abacus_knowledge"

# ================= 读取 API Key（完全照抄 build_knowledge_base.py） =================

api_key = ""
try:
    if os.path.exists(CONFIG_PATH):
        with open(CONFIG_PATH, "r", encoding="utf-8") as f:
            config = yaml.safe_load(f)
        if config.get("llm"):
            api_key = config["llm"].get("api_key", "")
        else:
            api_key = config.get("api_key", "")
except Exception:
    pass

if not api_key:
    api_key = os.getenv("DASHSCOPE_API_KEY", "")

if api_key:
    dashscope.api_key = api_key
else:
    print("❌ 警告: 未找到 API Key，后续调用将失败。")

# ================= 核心组件（与 build_knowledge_base.py 保持一致） =================

class QwenEmbeddingFunction(EmbeddingFunction):
    """
    自定义的 ChromaDB 嵌入函数，使用阿里云 Qwen Embedding API
    """

    def __call__(self, input: Documents) -> Embeddings:
        # 阿里云限制 batch_size <= 10，这里保持和全量脚本一致
        batch_size = 5
        all_embeddings = []
        for i in range(0, len(input), batch_size):
            batch = input[i : i + batch_size]
            try:
                resp = TextEmbedding.call(
                    model=TextEmbedding.Models.text_embedding_v3,
                    input=batch
                )
                if resp.status_code == 200:
                    embeddings = [item["embedding"] for item in resp.output["embeddings"]]
                    all_embeddings.extend(embeddings)
                else:
                    print(f"⚠️ API Error: {resp.message}")
                    all_embeddings.extend([[0.0] * 1024 for _ in range(len(batch))])
            except Exception as e:
                print(f"⚠️ Network Exception: {e}")
                all_embeddings.extend([[0.0] * 1024 for _ in range(len(batch))])
            # 稍微限流
            time.sleep(0.2)
        return all_embeddings

# ================= 工具函数（与 build_knowledge_base.py 保持一致） =================

def read_file(filepath: str) -> str:
    ext = os.path.splitext(filepath)[1].lower()
    content = ""
    try:
        if ext in [".md", ".txt"]:
            with open(filepath, "r", encoding="utf-8") as f:
                content = f.read()
        elif ext == ".docx":
            if "Document" in globals() and Document is not None:
                doc = Document(filepath)
                content = "\n".join([p.text for p in doc.paragraphs])
        elif ext == ".pdf":
            if PdfReader is None:
                print(f"⚠️ 跳过 PDF（未安装 pypdf）: {filepath}")
                return ""
            reader = PdfReader(filepath)
            pages = []
            for page in reader.pages:
                t = page.extract_text() or ""
                pages.append(t)
            content = "\n".join(pages)
    except Exception as e:
        print(f"⚠️ 无法读取 {filepath}: {e}")
    return content


def split_text(text: str, chunk_size=800, overlap=100) -> List[str]:
    if not text:
        return []
    chunks: List[str] = []
    start = 0
    text_len = len(text)
    if text_len <= chunk_size:
        return [text]
    while start < text_len:
        end = min(start + chunk_size, text_len)
        chunk = text[start:end]
        chunks.append(chunk)
        if end == text_len:
            break
        start += (chunk_size - overlap)
    return chunks

# ================= 增量追加逻辑 =================

def append_db():
    print("==========================================")
    print(" 🧩 AutoTutorial 知识库增量追加 (Aliyun版) ")
    print("==========================================")
    print(f"💾 数据库路径: {DB_PATH}")
    print(f"📂 增量目录: {SOURCE_DIR}")

    client = chromadb.PersistentClient(path=DB_PATH)

    collection = client.get_or_create_collection(
        name=COLLECTION_NAME,
        embedding_function=QwenEmbeddingFunction()
    )

    if not os.path.exists(SOURCE_DIR):
        print(f"❌ 目录不存在: {SOURCE_DIR}")
        return

    files = (
        glob.glob(os.path.join(SOURCE_DIR, "**/*.md"), recursive=True)
        + glob.glob(os.path.join(SOURCE_DIR, "**/*.docx"), recursive=True)
        + glob.glob(os.path.join(SOURCE_DIR, "**/*.txt"), recursive=True)
        + glob.glob(os.path.join(SOURCE_DIR, "**/*.pdf"), recursive=True)
    )

    print(f"📦 找到 {len(files)} 个增量文件，开始处理...")
    if not files:
        return

    total_chunks = 0
    # 记录至少有一段成功写入的文件路径
    ingested_files = []

    for idx, filepath in enumerate(files):
        fname = os.path.basename(filepath)
        print(f"[{idx+1}/{len(files)}] 处理: {fname}...", end="", flush=True)

        content = read_file(filepath)
        if not content.strip():
            print(" [空]")
            continue

        chunks = split_text(content)
        if not chunks:
            print(" [无内容]")
            continue

        ids = [f"{fname}_{i}_{str(uuid.uuid4())[:8]}" for i in range(len(chunks))]
        metadatas = [
            {"source": fname, "chunk_index": i, "ingest_type": "append"}
            for i in range(len(chunks))
        ]

        try:
            collection.add(documents=chunks, metadatas=metadatas, ids=ids)
            print(f" ✅ {len(chunks)} 片段")
            total_chunks += len(chunks)
            ingested_files.append(filepath)  # 只要这一文件有一次成功写入，就标记归档
        except Exception as e:
            print(f" ❌ 写入失败: {e}")

    print(f"\n🎉 增量追加完成！本次新增片段: {total_chunks}")
    print(f"💾 数据库路径: {DB_PATH}")

    # ========== 自动归档到 knowledge_source ==========
    target_root = os.path.join(ROOT_DIR, "data", "knowledge_source")
    os.makedirs(target_root, exist_ok=True)

    moved = 0
    for src_path in ingested_files:
        # 保留相对目录结构：knowledge_add/... -> knowledge_source/...
        rel_path = os.path.relpath(src_path, SOURCE_DIR)
        dst_path = os.path.join(target_root, rel_path)
        dst_dir = os.path.dirname(dst_path)
        os.makedirs(dst_dir, exist_ok=True)
        try:
            shutil.move(src_path, dst_path)
            moved += 1
        except Exception as e:
            print(f"⚠️ 归档失败 {src_path} -> {dst_path}: {e}")

    print(f"📂 已自动归档 {moved} 个文件到 knowledge_source。")


if __name__ == "__main__":
    append_db()
