import os
import yaml
import time
from openai import OpenAI

class LLMClient:
    def __init__(self, config_path=None):
        # 1. 路径定位
        if config_path is None:
            base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            config_path = os.path.join(base_dir, "config.yaml")
        
        # 2. 设置默认值 (关键！先定义默认值，防止 AttributeError)
        self.api_key = ""
        self.base_url = "https://api.openai.com/v1" # 默认官方地址
        self.model_name = "gemini-1.5-pro"

        # 3. 读取配置文件
        if os.path.exists(config_path):
            with open(config_path, 'r', encoding='utf-8') as f:
                config = yaml.safe_load(f)
                llm_conf = config.get('llm', {})
                
                # 优先级逻辑：
                # Key: 优先 google_api_key > api_key
                self.api_key = llm_conf.get('google_api_key') or llm_conf.get('api_key') or config.get('api_key')
                
                # Base URL: 优先 google_base_url > base_url
                # 配置文件里写的是 google_base_url
                self.base_url = llm_conf.get('google_base_url') or llm_conf.get('base_url', self.base_url)
                
                # Model
                self.model_name = llm_conf.get('model', self.model_name)

        # 4. 环境变量兜底
        if not self.api_key:
            self.api_key = os.getenv("GOOGLE_API_KEY", "")

        if not self.api_key:
            print("⚠️ [LLMClient] Warning: No API Key found.")

        # 5. 打印调试信息 (方便你确认连对了地方)
        print(f"[LLM Init] Connecting to: {self.base_url}")
        print(f"[LLM Init] Model: {self.model_name}")

        # 6. 初始化 OpenAI 客户端
        # 注意：中转商的 base_url 通常需要以 /v1 结尾
        if not self.base_url.endswith("/v1"):
            self.base_url = self.base_url.rstrip("/") + "/v1"
            
        self.client = OpenAI(
            api_key=self.api_key, 
            base_url=self.base_url
        )

    def chat(self, prompt, history=None, max_retries=3, model_id=None, temperature=0.7):
        target_model = model_id if model_id else self.model_name
        
        # 构造消息
        messages = []
        if history:
            messages.extend(history)
        messages.append({"role": "user", "content": prompt})

        for attempt in range(max_retries):
            try:
                # 打印当前尝试的模型，方便调试
                if attempt == 0:
                    print(f"   [LLM] Calling {target_model}...")

                response = self.client.chat.completions.create(
                    model=target_model,
                    messages=messages,
                    temperature=temperature
                )

                content = response.choices[0].message.content
                if content:
                    return content
                else:
                    print(f"⚠️ [LLM] Empty response.")
                    time.sleep(1)

            except Exception as e:
                print(f"⚠️ [LLM] Exception (Attempt {attempt+1}/{max_retries}): {e}")
                time.sleep(2)
        
        print("❌ [LLM] Failed after retries.")
        return ""
