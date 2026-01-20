import os
import chromadb
import dashscope
from dashscope import TextEmbedding
import yaml
import time

# 复用配置逻辑
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(CURRENT_DIR)
DB_PATH = os.path.join(ROOT_DIR, "data", "chroma_db")
CONFIG_PATH = os.path.join(ROOT_DIR, "config.yaml")

# 初始化 API Key (Qwen 的)
try:
    if os.path.exists(CONFIG_PATH):
        with open(CONFIG_PATH, "r", encoding="utf-8") as f:
            config = yaml.safe_load(f)
            # 注意：Retriever 这里只用阿里云的 Key 做 Embedding
            api_key = config.get("llm", {}).get("api_key") or config.get("api_key")
            if api_key:
                dashscope.api_key = api_key
except:
    pass
if not dashscope.api_key:
    dashscope.api_key = os.getenv("DASHSCOPE_API_KEY", "")

class LocalRetriever:
    def __init__(self):
        # 初始化 ChromaDB 客户端
        self.client = chromadb.PersistentClient(path=DB_PATH)
        
        # === 关键修改 ===
        # 不再传递 embedding_function 参数，避免与 knowledge_manager 创建时的默认配置冲突
        # 我们将在 search 时手动计算向量
        try:
            self.collection = self.client.get_collection(name="abacus_knowledge")
        except Exception:
            # 如果不存在，尝试创建（虽然后果可能还是不匹配，但至少防崩）
            self.collection = self.client.get_or_create_collection(name="abacus_knowledge")

    def _get_embedding(self, text: str):
        """
        手动调用阿里云 API 计算查询词的向量
        """
        try:
            resp = TextEmbedding.call(
                model=TextEmbedding.Models.text_embedding_v3,
                input=[text]
            )
            if resp.status_code == 200:
                return resp.output["embeddings"][0]["embedding"]
            else:
                print(f"[Retriever] Embedding Error: {resp.message}")
                return None
        except Exception as e:
            print(f"[Retriever] Network Error: {e}")
            return None

    def search(self, query: str, top_k: int = 5) -> str:
        """
        输入问题 -> 手动转向量 -> query_embeddings -> 返回文本
        """
        print(f"   [Retrieving] 搜索: {query} ...")
        
        # 1. 手动计算 Query 向量
        query_vec = self._get_embedding(query)
        if not query_vec:
            return ""

        # 2. 使用向量进行查询 (query_embeddings)
        results = self.collection.query(
            query_embeddings=[query_vec],  # <--- 这里改用向量查
            n_results=top_k
        )
        
        # 3. 拼接结果
        context_parts = []
        if results and results['documents']:
            docs = results['documents'][0]
            metas = results['metadatas'][0]
            for i, doc in enumerate(docs):
                source = metas[i].get('source', 'Unknown')
                context_parts.append(f"[资料{i+1} | Source: {source}]\n{doc}")
        
        return "\n\n".join(context_parts)

if __name__ == "__main__":
    # 测试代码
    r = LocalRetriever()
    print(r.search("DFT计算参数设置"))
