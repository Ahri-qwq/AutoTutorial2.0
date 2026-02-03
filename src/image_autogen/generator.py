import os
import json
import dataclasses
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Any
from pathlib import Path

from src.llm_client import LLMClient
from src.image_autogen.utils import clean_json_string

@dataclass
class FigureSpec:
    figure_id: str
    title: str = ""
    caption: str = ""
    insert_after_heading: str = ""
    diagram_type: str = ""
    source_quotes: List[str] = field(default_factory=list)
    entities_whitelist: List[str] = field(default_factory=list)
    labels_required: List[str] = field(default_factory=list)
    labels_forbidden: List[str] = field(default_factory=list)
    assertions: List[str] = field(default_factory=list)
    prompt: str = ""
    negative_prompt: Optional[str] = ""

    @classmethod
    def from_dict(cls, data: Dict[str, Any]):
        # 安全地从字典创建对象，忽略多余字段，处理缺失字段
        # 这种方式会使用 dataclass 的默认值（如果定义了）
        fields = {f.name: f for f in dataclasses.fields(cls)}
        init_data = {}
        for name, f in fields.items():
            if name in data:
                init_data[name] = data[name]
            # 如果 data 中没有，dataclass 会在初始化时使用默认值
        return cls(**init_data)

@dataclass
class ChapterSpec:
    chapter_id: str
    chapter_title: str
    figures: List[FigureSpec]

class SpecGenerator:
    def __init__(self, llm_client: LLMClient):
        self.llm = llm_client
        self.template_path = os.path.join(
            os.path.dirname(__file__), 
            "prompts", 
            "step3_figurespec_zh.txt"
        )
        with open(self.template_path, 'r', encoding='utf-8') as f:
            self.prompt_template = f.read()

    def generate_chapter_spec(self, 
                              chapter_id: str, 
                              chapter_title: str, 
                              chapter_content: str, 
                              headings: List[str], 
                              num_figures: int,
                              model_name: Optional[str] = None) -> ChapterSpec:
        
        if num_figures <= 0:
            return ChapterSpec(chapter_id, chapter_title, [])

        # 构造 Prompt
        prompt = self.prompt_template.replace("CHAPTER_TITLE", chapter_title)
        # 截断过长的 chapter_content 以防超出 token 限制 (简单策略: 取前 20000 字符，通常够了)
        prompt = prompt.replace("CHAPTER_MD", chapter_content[:30000]) 
        prompt = prompt.replace("HEADINGS", json.dumps(headings, ensure_ascii=False))
        prompt = prompt.replace("NUM_FIGURES", str(num_figures))

        # 调用 LLM
        print(f"   [SpecGen] Generating spec for {chapter_id} ({num_figures} figures)...")
        response_text = self.llm.chat(prompt, model_id=model_name)
        
        # 解析 JSON
        try:
            cleaned_json = clean_json_string(response_text)
            data = json.loads(cleaned_json)
            
            figures_data = data.get("figures", [])
            figures = []
            for fig_data in figures_data:
                # 修正 figure_id 为全局唯一 ID (chXX-figYY)
                # 假设 LLM 生成的是 fig01, fig02...
                local_id = fig_data.get("figure_id", "fig00")
                # 简单的归一化处理
                if local_id.startswith("fig"):
                     # 尝试提取数字
                    try:
                        idx = int(local_id.replace("fig", ""))
                        fig_id = f"{chapter_id}-fig{idx:02d}"
                    except:
                        fig_id = f"{chapter_id}-{local_id}"
                else:
                    fig_id = f"{chapter_id}-{local_id}"
                
                fig_data['figure_id'] = fig_id
                figures.append(FigureSpec.from_dict(fig_data))
                
            return ChapterSpec(chapter_id, chapter_title, figures)
            
        except json.JSONDecodeError as e:
            print(f"❌ [SpecGen] JSON Parse Error: {e}")
            print(f"Raw response: {response_text[:500]}...")
            raise e
        except Exception as e:
            print(f"❌ [SpecGen] Error processing response: {e}")
            raise e

class ImageGenerator:
    def __init__(self, llm_client: LLMClient):
        self.llm = llm_client
        self.style_template_path = os.path.join(
            os.path.dirname(__file__), 
            "prompts", 
            "step4_image_style_zh.txt"
        )
        with open(self.style_template_path, 'r', encoding='utf-8') as f:
            self.style_prefix = f.read()

    def _generate_via_chat(self, spec: FigureSpec, prompt: str, output_path: str, model_name: str) -> Dict[str, Any]:
        """
        Special handler for models that generate images via Chat Completions (e.g. Gemini).
        Expects the model to return markdown image syntax or base64 data in the content.
        """
        import re
        import base64
        import time

        # Construct messages
        messages = [
            {"role": "system", "content": "你是图片生成模型。请根据用户提示生成一张PNG图片。"},
            {"role": "user", "content": prompt}
        ]

        # Use client.chat.completions.create directly
        # Note: Some proxies use 'maxTokens', others 'max_tokens'. 
        # OpenAI python client uses 'max_tokens'.
        try:
            print(f"   [ImageGen] Calling Chat API for {model_name}...")
            response = self.llm.client.chat.completions.create(
                model=model_name,
                messages=messages,
                max_tokens=4096, 
                temperature=0.7
            )
            
            content = response.choices[0].message.content
            
            # Extract base64
            # Look for: ![image](data:image/jpeg;base64,...) or just the base64 string
            b64_data = ""
            
            # Regex for markdown image with data URI
            match = re.search(r"data:image\/[a-zA-Z]+;base64,([A-Za-z0-9+/=]+)", content)
            if match:
                b64_data = match.group(1)
            else:
                # Try to find any long base64-like string
                # This is risky but might work if the model outputs raw base64
                # Filter for strings longer than 1000 chars that look like base64
                potential_b64s = re.findall(r"([A-Za-z0-9+/=]{1000,})", content)
                if potential_b64s:
                    b64_data = potential_b64s[0]

            if not b64_data:
                raise ValueError("No base64 image data found in response content.")

            # Decode and save
            image_bytes = base64.b64decode(b64_data)
            with open(output_path, 'wb') as f:
                f.write(image_bytes)

            return {
                "figure_id": spec.figure_id,
                "status": "success",
                "output_file": output_path,
                "final_prompt": prompt,
                "model": model_name,
                "revision": response.created # timestamp
            }

        except Exception as e:
            print(f"❌ [ImageGen] Chat-based generation failed for {spec.figure_id}: {e}")
            if 'content' in locals() and content:
                 print(f"    Response content start: {content[:200]}...")
            return {
                "figure_id": spec.figure_id,
                "status": "failed",
                "error": str(e)
            }

    def generate_image(self, 
                       spec: FigureSpec, 
                       output_path: str, 
                       model_name: str = "dall-e-3") -> Dict[str, Any]:
        
        # 组装 Prompt
        full_prompt = f"{self.style_prefix}\n\nSpecific Request:\n{spec.prompt}"
        
        if spec.negative_prompt:
             full_prompt += f"\n\nAvoid: {spec.negative_prompt}"
             
        # 标签要求通过 spec.prompt 内部控制，不再强制追加所有标签
        # DALL-E 3 中文生成能力较弱，强制追加中文标签会导致图片出现乱码
        # if spec.labels_required:
        #    labels_str = ", ".join(spec.labels_required)
        #    full_prompt += f"\n\nRequired Labels (Simplified Chinese): {labels_str}"

        print(f"   [ImageGen] Generating {spec.figure_id} using {model_name}...")
        
        try:
            # 兼容 gemini-3-pro-image-preview, nano banana2 等使用 Chat 接口生成图像的模型
            image_via_chat_models = ["gemini", "banana", "nano"]
            if any(m in model_name.lower() for m in image_via_chat_models):
                return self._generate_via_chat(spec, full_prompt, output_path, model_name)

            # 使用 LLMClient 的内部 client (OpenAI 兼容实例)
            response = self.llm.client.images.generate(
                model=model_name,
                prompt=full_prompt,
                n=1,
                size="1024x1024", 
                response_format="b64_json" 
            )
            
            # 保存图片
            import base64
            image_data = base64.b64decode(response.data[0].b64_json)
            
            with open(output_path, 'wb') as f:
                f.write(image_data)
                
            return {
                "figure_id": spec.figure_id,
                "status": "success",
                "output_file": output_path,
                "final_prompt": full_prompt,
                "model": model_name,
                "revision": response.created
            }
            
        except Exception as e:
            print(f"❌ [ImageGen] Generation failed for {spec.figure_id}: {e}")
            return {
                "figure_id": spec.figure_id,
                "status": "failed",
                "error": str(e)
            }
