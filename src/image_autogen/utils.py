import os
import json
import glob
import re
import shutil
from datetime import datetime
from pathlib import Path
from typing import List, Dict, Any, Tuple, Optional

def get_timestamp_str() -> str:
    """生成当前时间戳字符串，用于 run_id"""
    return datetime.now().strftime("%Y%m%d_%H%M%S")

def load_json(path: str) -> Any:
    """读取 JSON 文件"""
    with open(path, 'r', encoding='utf-8') as f:
        return json.load(f)

def save_json(data: Any, path: str, indent: int = 2):
    """保存 JSON 文件"""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(data, f, ensure_ascii=False, indent=indent)

def find_input_md(run_dir: str) -> str:
    """在 run_dir 中查找唯一的 Final_Tutorial_*.md"""
    pattern = os.path.join(run_dir, "Final_Tutorial_*.md")
    files = glob.glob(pattern)
    
    if not files:
        raise FileNotFoundError(f"No file matching 'Final_Tutorial_*.md' found in {run_dir}")
    if len(files) > 1:
        raise ValueError(f"Multiple 'Final_Tutorial_*.md' files found in {run_dir}. Please ensure only one exists.")
    
    return files[0]

def prepare_output_dirs(run_dir: str, overwrite: bool = False) -> Dict[str, str]:
    """创建输出目录并返回路径字典"""
    base_out = os.path.join(run_dir, "images_autogen")
    images_out = os.path.join(base_out, "images_output")
    figures_dir = os.path.join(base_out, "figures") # 保留 figures 目录以防万一，虽然新计划主要用 images_output

    if overwrite and os.path.exists(base_out):
        print(f"⚠️ Overwriting output directory: {base_out}")
        shutil.rmtree(base_out)

    os.makedirs(base_out, exist_ok=True)
    os.makedirs(images_out, exist_ok=True)
    # os.makedirs(figures_dir, exist_ok=True) # 暂时不创建这个，除非有需求
    
    return {
        "base": base_out,
        "images": images_out,
        "chapter_plan": os.path.join(base_out, "chapter_plan.json"),
        "figurespec": os.path.join(base_out, "figurespec.json"),
        "figurespec_resolved": os.path.join(base_out, "figurespec.resolved.json"),
        "audit_report_md": os.path.join(base_out, "audit_report.md"),
        "audit_report_json": os.path.join(base_out, "audit_report.json"),
        "rework_queue": os.path.join(base_out, "rework_queue.json")
    }

class Chapter:
    def __init__(self, chapter_id: str, title: str, content: str, headings: List[str]):
        self.chapter_id = chapter_id
        self.title = title
        self.content = content
        self.headings = headings

    def to_dict(self):
        return {
            "chapter_id": self.chapter_id,
            "title": self.title,
            # content 太长通常不放在摘要里，但在 spec 生成时需要
            "headings_count": len(self.headings)
        }

def parse_chapters(md_text: str) -> List[Chapter]:
    """
    解析 Markdown 章节。
    策略：
    1. 查找所有以 # 开头的行（跳过代码块内部的行）。
    2. 判断主要层级（H1 还是 H2）。如果存在 H1，则按 H1 分章；否则按 H2 分章。
    """
    lines = md_text.split('\n')
    
    # 扫描所有 heading
    headings = []
    in_code_block = False
    for i, line in enumerate(lines):
        stripped = line.strip()
        
        # 处理代码块边界
        if stripped.startswith('```'):
            in_code_block = not in_code_block
            continue
            
        # 处于代码块中，跳过
        if in_code_block:
            continue

        # 识别标题行 (如 # 标题 或 ## 标题)
        match = re.match(r'^(#+)(.*)$', stripped)
        if match:
            hashes = match.group(1)
            title = match.group(2).strip()
            level = len(hashes)
            headings.append({'level': level, 'title': title, 'line_index': i, 'raw': stripped})
    
    if not headings:
        return []

    # 决定分章层级
    min_level = min(h['level'] for h in headings)
    chapter_level = min_level # H1 优先，如果没有 H1 则是 H2
    
    chapters = []
    current_chapter_start = 0
    current_chapter_title = "Intro"
    current_chapter_headings = []
    
    # 找到第一个章节标题的索引
    first_chapter_idx = -1
    for i, h in enumerate(headings):
        if h['level'] == chapter_level:
            first_chapter_idx = i
            break
            
    if first_chapter_idx == -1:
        # 理论上不会发生，除非全是 H2 但 min_level 计算错
        return [Chapter("ch00", "Full Document", md_text, [h['raw'] for h in headings])]

    # 处理第一个章节前的内容（如果有）
    if headings[first_chapter_idx]['line_index'] > 0:
        pre_content = '\n'.join(lines[:headings[first_chapter_idx]['line_index']])
        # 只有当前面有实质内容时才算作一章，或者直接忽略前言？
        # 为了简单，我们暂时忽略没有标题的前言，或者把它归为 Intro
        pass

    chapter_count = 0
    
    for i, h in enumerate(headings):
        if h['level'] == chapter_level:
            # 这是一个新章的开始
            
            # 如果之前有正在累积的章节，先保存
            if i > 0 and chapter_count > 0:
                # 获取上一章的结束行（当前章标题行的上一行）
                end_line = h['line_index']
                start_line = headings[headings.index(last_chapter_heading)]['line_index']
                content = '\n'.join(lines[start_line:end_line])
                
                # 收集上一章内部的子标题
                sub_headings = []
                for sub_h in headings:
                    if sub_h['line_index'] >= start_line and sub_h['line_index'] < end_line:
                        sub_headings.append(sub_h['raw'])

                chap_id = f"ch{chapter_count:02d}"
                chapters.append(Chapter(chap_id, last_chapter_heading['title'], content, sub_headings))
            
            last_chapter_heading = h
            chapter_count += 1
            
    # 保存最后一章
    if chapter_count > 0:
        start_line = last_chapter_heading['line_index']
        content = '\n'.join(lines[start_line:])
        sub_headings = [sh['raw'] for sh in headings if sh['line_index'] >= start_line]
        chap_id = f"ch{chapter_count:02d}"
        chapters.append(Chapter(chap_id, last_chapter_heading['title'], content, sub_headings))
        
    return chapters
        
    return chapters

def clean_json_string(text: str) -> str:
    """清理 LLM 返回的 JSON 字符串（去除 markdown code block）"""
    text = text.strip()
    # 移除 ```json ... ``` 或 ``` ... ```
    if text.startswith("```"):
        # 找到第一个换行符
        first_newline = text.find("\n")
        if first_newline != -1:
            # 去掉第一行 (```json)
            text = text[first_newline+1:]
            # 去掉结尾的 ```
            if text.endswith("```"):
                text = text[:-3]
    return text.strip()
