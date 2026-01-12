import os
import yaml
import time
from http import HTTPStatus
import dashscope
from dashscope import Generation

class LLMClient:
    def __init__(self, config_path=None):
        # 1. 加载配置
        if config_path is None:
            # 自动定位到根目录的 config.yaml
            base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            config_path = os.path.join(base_dir, "config.yaml")
        
        self.api_key = ""
        self.model_name = "qwen-max" # 默认模型

        if os.path.exists(config_path):
            with open(config_path, 'r', encoding='utf-8') as f:
                config = yaml.safe_load(f)
                # 兼容两种 config 格式
                llm_conf = config.get('llm', {})
                self.api_key = llm_conf.get('api_key') or config.get('api_key')
                # 如果 config 里写了 model 就用 config 的，否则默认 qwen-max
                self.model_name = llm_conf.get('model', 'qwen-max')

        # 环境变量兜底
        if not self.api_key:
            self.api_key = os.getenv("DASHSCOPE_API_KEY", "")

        if not self.api_key:
            print("⚠️ [LLMClient] Warning: No API Key found.")
        
        dashscope.api_key = self.api_key

    def chat(self, prompt, history=None, max_retries=3):
        """
        纯粹的对话接口，不包含任何检索逻辑
        """
        messages = [{'role': 'system', 'content': 'You are a helpful assistant.'}]
        if history:
            messages.extend(history)
        
        messages.append({'role': 'user', 'content': prompt})

        for attempt in range(max_retries):
            try:
                response = Generation.call(
                    model=self.model_name,
                    messages=messages,
                    result_format='message',  # 关键：设置为 message 格式
                    temperature=0.7
                )

                if response.status_code == HTTPStatus.OK:
                    # 提取内容
                    content = response.output.choices[0].message.content
                    return content
                else:
                    print(f"⚠️ [LLM] Request Failed: {response.code} - {response.message}")
                    time.sleep(1)
            
            except Exception as e:
                print(f"⚠️ [LLM] Exception: {e}")
                time.sleep(2)
        
        print("❌ [LLM] Failed after retries.")
        return ""
