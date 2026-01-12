import os
import glob
import uuid
import time
from typing import List

import chromadb
from chromadb import Documents, EmbeddingFunction, Embeddings
import dashscope
from dashscope import TextEmbedding
import yaml

try:
    from docx import Document
except ImportError:
    pass

# ================= 路径配置 =================
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(CURRENT_DIR)
CONFIG_PATH = os.path.join(ROOT_DIR, "config.yaml")
SOURCE_DIR = os.path.join(ROOT_DIR, "data", "knowledge_source")
DB_PATH = os.path.join(ROOT_DIR, "data", "chroma_db")
COLLECTION_NAME = "abacus_knowledge"

# ================= 读取 API Key =================
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

# ================= 核心组件 =================

class QwenEmbeddingFunction(EmbeddingFunction):
    """
    自定义的 ChromaDB 嵌入函数，使用阿里云 Qwen Embedding API
    """
    def __call__(self, input: Documents) -> Embeddings:
        # 【关键修改】阿里云限制 batch_size <= 10，我们设为 5 以防万一
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
                    embeddings = [item['embedding'] for item in resp.output['embeddings']]
                    all_embeddings.extend(embeddings)
                else:
                    print(f"⚠️ API Error: {resp.message}")
                    all_embeddings.extend([[0.0]*1024 for _ in range(len(batch))])
            except Exception as e:
                print(f"⚠️ Network Exception: {e}")
                all_embeddings.extend([[0.0]*1024 for _ in range(len(batch))])
            
            # 稍微限流
            time.sleep(0.2)
            
        return all_embeddings

def read_file(filepath: str) -> str:
    ext = os.path.splitext(filepath)[1].lower()
    content = ""
    try:
        if ext in [".md", ".txt"]:
            with open(filepath, "r", encoding="utf-8") as f:
                content = f.read()
        elif ext == ".docx":
            if 'Document' in globals():
                doc = Document(filepath)
                content = "\n".join([p.text for p in doc.paragraphs])
    except Exception as e:
        print(f"⚠️ 无法读取 {filepath}: {e}")
    return content

def split_text(text: str, chunk_size=800, overlap=100) -> List[str]:
    if not text: return []
    chunks = []
    start = 0
    text_len = len(text)
    if text_len <= chunk_size: return [text]
    
    while start < text_len:
        end = min(start + chunk_size, text_len)
        chunk = text[start:end]
        chunks.append(chunk)
        if end == text_len: break
        start += (chunk_size - overlap)
    return chunks

def build_db():
    print(f"==========================================")
    print(f" 🚀 AutoTutorial 知识库构建 (Aliyun版) ")
    print(f"==========================================")
    print(f"💾 数据库路径: {DB_PATH}")
    
    # 1. 初始化
    client = chromadb.PersistentClient(path=DB_PATH)
    
    # 2. 清理旧集合
    try:
        client.delete_collection(COLLECTION_NAME)
        print(f"🧹 已清理旧集合")
    except Exception:
        pass

    # 3. 创建新集合
    print("🔌 连接阿里云 Embedding API...")
    collection = client.create_collection(
        name=COLLECTION_NAME,
        embedding_function=QwenEmbeddingFunction()
    )

    # 4. 扫描文件
    if not os.path.exists(SOURCE_DIR):
        print(f"❌ 目录不存在: {SOURCE_DIR}")
        return

    files = glob.glob(os.path.join(SOURCE_DIR, "**/*.md"), recursive=True) + \
            glob.glob(os.path.join(SOURCE_DIR, "**/*.docx"), recursive=True) + \
            glob.glob(os.path.join(SOURCE_DIR, "**/*.txt"), recursive=True)
    
    print(f"📦 找到 {len(files)} 个文件，开始处理...")
    total_chunks = 0
    
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
        metadatas = [{"source": fname, "chunk_index": i} for i in range(len(chunks))]
        
        # 写入
        try:
            collection.add(documents=chunks, metadatas=metadatas, ids=ids)
            print(f" ✅ {len(chunks)} 片段")
            total_chunks += len(chunks)
        except Exception as e:
            print(f" ❌ 写入失败: {e}")

    print(f"\n🎉 构建完成！总片段: {total_chunks}")
    print(f"💾 数据库已保存至: {DB_PATH}")

if __name__ == "__main__":
    build_db()
