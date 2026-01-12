import os
import chromadb
from chromadb import Documents, EmbeddingFunction, Embeddings
import dashscope
from dashscope import TextEmbedding
import yaml
import time

# 复用配置逻辑
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(CURRENT_DIR)
DB_PATH = os.path.join(ROOT_DIR, "data", "chroma_db")
CONFIG_PATH = os.path.join(ROOT_DIR, "config.yaml")

# 初始化 API Key
try:
    if os.path.exists(CONFIG_PATH):
        with open(CONFIG_PATH, "r", encoding="utf-8") as f:
            config = yaml.safe_load(f)
            api_key = config.get("llm", {}).get("api_key") or config.get("api_key")
            if api_key:
                dashscope.api_key = api_key
except:
    pass
if not dashscope.api_key:
    dashscope.api_key = os.getenv("DASHSCOPE_API_KEY", "")

class QwenEmbeddingFunction(EmbeddingFunction):
    def __call__(self, input: Documents) -> Embeddings:
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
                    all_embeddings.extend([item['embedding'] for item in resp.output['embeddings']])
                else:
                    all_embeddings.extend([[0.0]*1024 for _ in batch])
            except:
                all_embeddings.extend([[0.0]*1024 for _ in batch])
        return all_embeddings

class LocalRetriever:
    def __init__(self):
        self.client = chromadb.PersistentClient(path=DB_PATH)
        self.collection = self.client.get_or_create_collection(
            name="abacus_knowledge",
            embedding_function=QwenEmbeddingFunction()
        )

    def search(self, query: str, top_k: int = 5) -> str:
        """
        输入问题，返回拼接好的参考资料文本
        """
        print(f"   [Retrieving] 搜索: {query} ...")
        results = self.collection.query(
            query_texts=[query],
            n_results=top_k
        )
        
        # 拼接结果
        context_parts = []
        if results['documents']:
            docs = results['documents'][0]
            metas = results['metadatas'][0]
            for doc, meta in zip(docs, metas):
                source = meta.get('source', 'unknown')
                context_parts.append(f"--- 来源: {source} ---\n{doc}")
        
        return "\n\n".join(context_parts)
