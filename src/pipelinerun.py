import os
import re
import time
from src.llm_client import LLMClient
from src.retriever import LocalRetriever  # 新增引用

class AutoTutorialPipeline:
    def __init__(self, root_dir):
        self.root_dir = root_dir
        self.prompts_dir = os.path.join(root_dir, "prompts")
        self.llm = LLMClient()
        self.retriever = LocalRetriever()
        
        # === 生成基于时间戳的运行目录 ===
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        self.run_dir = os.path.join(root_dir, "data", "runs", f"{timestamp}")
        # 在该运行目录下创建 processed 子目录（用于存放过程文件）
        self.processed_dir = os.path.join(self.run_dir, "processed")
        os.makedirs(self.processed_dir, exist_ok=True)
        print(f"[Init] 本次运行的工作目录: {self.run_dir}")
        print(f"[Init] 中间文件将保存在: {self.processed_dir}")
        self.current_topic = ""
        self.current_special_instructions = ""

    def _load_prompt(self, filename):
        path = os.path.join(self.prompts_dir, filename)
        if os.path.exists(path):
            with open(path, 'r', encoding='utf-8') as f:
                return f.read()
        return ""

    def _fill_prompt(self, template, **kwargs):
        """智能填充 Prompt，如果模板里没有占位符，自动追加在末尾"""
        content = template
        context = kwargs.pop("CONTEXT", "")
        
        # 1. 替换标准占位符
        for k, v in kwargs.items():
            content = content.replace(f"{{{{{k}}}}}", str(v))
        
        # 2. 处理上下文 (CONTEXT)
        if "{{CONTEXT}}" in content:
            content = content.replace("{{CONTEXT}}", context)
        elif context:
            # 如果模板里没写 {{CONTEXT}}，我们手动加在最后
            content += f"\n\n# 参考资料库 (Knowledge Base)\n请基于以下检索到的资料进行回答，如果资料不足，请结合你的专业知识补充：\n{context}"
            
        return content

    def run_step1(self, topic):
        """Step 1: 知识预研"""
        self.current_topic = topic
        print(f"\n=== Step 1: 知识预研 (Topic: {topic}) ===")
        
        # 1. 检索
        retrieved_context = self.retriever.search(f"关于 {topic} 的物理原理、参数设置和计算流程", top_k=8)
        
        # 2. 加载 Prompt
        prompt_tmpl = self._load_prompt("step1_enrich.txt")
        if not prompt_tmpl: return

        # 3. 组装
        final_prompt = self._fill_prompt(prompt_tmpl, TOPIC=topic, CONTEXT=retrieved_context)

        # 4. 调用 LLM
        print(f"[Agent] 正在分析...")
        response = self.llm.chat(final_prompt)
        if not response: return

        # 5. 保存
        output_path = os.path.join(self.processed_dir, "step1_enrichment.md")
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(response)
        print(f"[Success] 知识预研完成: {output_path}")

    def run_step2(self):
        """Step 2: 大纲规划"""
        print("\n=== Step 2: 大纲规划 ===")
        
        step1_path = os.path.join(self.processed_dir, "step1_enrichment.md")
        if not os.path.exists(step1_path):
            print("[Error] 请先运行 Step 1")
            return
            
        with open(step1_path, 'r', encoding='utf-8') as f:
            background_info = f.read()

        prompt_tmpl = self._load_prompt("step2_outline.txt")
        if not prompt_tmpl: return

        # Step 2 不需要额外检索，因为它依赖 Step 1 的总结
        final_prompt = self._fill_prompt(
            prompt_tmpl, 
            BACKGROUND_INFO=background_info[:6000],
            TOPIC=self.current_topic
        )

        print("[Agent] 正在设计大纲...")
        response = self.llm.chat(final_prompt)
        if not response: return

        output_path = os.path.join(self.processed_dir, "step2_outline.md")
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(response)
        print(f"[Success] 大纲已生成: {output_path}")

    def run_step3(self):
        """Step 3: 分章撰写"""
        print("\n=== Step 3: 技术章节撰写 ===")
        outline_path = os.path.join(self.processed_dir, "step2_outline.md")
        if not os.path.exists(outline_path):
            print("[Error] 请先运行 Step 2")
            return

        with open(outline_path, 'r', encoding='utf-8') as f:
            outline = f.read()

        # === 新增：提取特殊指令 ===
        step1_path = os.path.join(self.processed_dir, "step1_enrichment.md")
        special_instructions = "无特殊指令。"
        if os.path.exists(step1_path):
            with open(step1_path, 'r', encoding='utf-8') as f:
                content = f.read()
                # 简单正则提取 "## 4. 教程编写特殊指令" 后的内容
                match = re.search(r"## 4\. 教程编写特殊指令(.*?)(##|\Z)", content, re.DOTALL)
                if match:
                    special_instructions = match.group(1).strip()
        # 将特殊指令存为类属性，供 _draft_chapter 使用
        self.current_special_instructions = special_instructions

        # 简单的切分逻辑：按HTML注释或二级标题切分
        parts = outline.split("<!-- CHAPTER_START -->")
        if len(parts) < 2:
            # 兼容没有注释的情况
            parts = outline.split("## ")
        
        # 跳过前言
        for idx, chunk in enumerate(parts):
            if idx == 0: continue # 跳过前言/Meta
            
            # 清理内容
            chunk = chunk.replace("CHAPTER_START -->", "").strip()
            if not chunk: continue
            if "##" not in chunk and not chunk.startswith("第"): 
                chunk = "## " + chunk # 补回标题标记

            # 提取标题
            lines = chunk.split('\n')
            title = lines[0].replace("#", "").strip()
            
            self._draft_chapter(title, chunk)

    def _draft_chapter(self, title, outline_context):
        print(f"\n[Drafting] 正在撰写: {title} ...")
        
        # 1. 针对本章标题进行精准检索
        # 比如标题是"能带计算"，我们会搜出具体的 KPT, INPUT 设置
        #retrieved_context = self.retriever.search(f"ABACUS教程: {title} 的具体参数、输入文件示例和注意事项", top_k=6)
        
        # （1）. 宽泛检索：搜原理
        theory_context = self.retriever.search(f"ABACUS原理 {title}", top_k=3)
        
        # （2）. 精准检索：搜参数 (增加 keyword 权重)
        # 强制带上 "INPUT" 或 "参数" 关键词
        param_context = self.retriever.search(f"ABACUS INPUT参数设置 {title} 示例", top_k=5)
        
        # 合并上下文
        full_context = f"【理论背景】\n{theory_context}\n\n【参数与实操】\n{param_context}"
        
        # 2. 加载 Prompt
        prompt_tmpl = self._load_prompt("step3_drafting.txt")
        if not prompt_tmpl: return

        # 3. 填充
        final_prompt = self._fill_prompt(
            prompt_tmpl, 
            CHAPTER_TITLE=title, 
            CHAPTER_DESCRIPTION=outline_context,
            #CONTEXT=retrieved_context
            CONTEXT=full_context,
            SPECIAL_INSTRUCTIONS=getattr(self, 'current_special_instructions', '')
        )

        # 4. 生成
        content = self.llm.chat(final_prompt)
        if not content: return
        # === [修复] 清理 LLM 可能重复生成的标题 ===
        content = re.sub(r'^#\s+.*?\n+', '', content.strip()).strip()
        # === 文件名增加数字前缀 ===
        # 尝试从标题中提取章节号 (例如 "第一章" -> "01", "3.1" -> "03_01")
        # 简单的映射逻辑：第X章 -> 0X
        prefix = "99" # 默认排序靠后
        chapter_map = {"一": "01", "二": "02", "三": "03", "四": "04", "五": "05", 
                       "六": "06", "七": "07", "八": "08", "九": "09", "十": "10"}
        
        match = re.search(r"第([一二三四五六七八九十]+)章", title)
        if match:
            cn_num = match.group(1)
            prefix = chapter_map.get(cn_num, "99")


        # 5. 保存
        safe_title = re.sub(r'[\\/*?:"<>|]', "", title).replace(" ", "_")
        # 新的文件名格式: 01_chapter_Title.md
        fname = f"{prefix}_chapter_{safe_title}.md"
        
        with open(os.path.join(self.processed_dir, fname), 'w', encoding='utf-8') as f:
            f.write(f"# {title}\n\n{content}")
        print(f" -> 已保存: {fname}")
        
        time.sleep(2)

    def run_step4(self):
        """Step 4: 智能组装全文 (增强版：带RAG的前言与附录)"""
        print("\n=== Step 4: 全文组装 (RAG Enhanced) ===")
        
        # 1. 读取章节
        chapter_files = sorted([f for f in os.listdir(self.processed_dir) if f.endswith(".md") and ("chapter_" in f or f.startswith("0"))])
        if not chapter_files:
            print(f"[Error] 在 {self.processed_dir} 没有找到章节文件。")
            return

        full_content = ""
        summaries = []
        for fname in chapter_files:
            path = os.path.join(self.processed_dir, fname)
            with open(path, 'r', encoding='utf-8') as f:
                content = f.read()
                full_content += content + "\n\n"
                summaries.append(f"[{fname}]: {content[:300]}...") # 增加摘要长度

        # 2. 宏观检索 (Macro Retrieval)
        # 我们去知识库里搜“学习路线”、“相关知识”、“前置基础”
        print("[Retrieving] 正在搜索宏观背景知识 (学习路线/相关概念)...")
        macro_query = f"{self.current_topic} 的前置知识、进阶学习路线、相关联的物理方法、推荐阅读文献"
        # 检索 top_k 设大一点，获取更多宽泛信息
        context = self.retriever.search(macro_query, top_k=6)

        # 3. 调用 Agent
        prompt_tmpl = self._load_prompt("step4_assembly.txt")
        if prompt_tmpl:
            print("[Agent] 正在撰写深度前言与附录...")
            
            # 这里我们利用 _fill_prompt 的自动追加功能，把 context 塞进去
            final_prompt = self._fill_prompt(
                prompt_tmpl,
                CHAPTER_SUMMARIES="\n".join(summaries),
                CONTEXT=context  # 注入宏观知识！
            )
            
            response = self.llm.chat(final_prompt)
            
            # JSON 解析逻辑 (保持不变)
            import json
            book_title = f"ABACUS 实战教程：{self.current_topic}"
            preface = ""
            appendix = ""
            
            try:
                clean_json = response.replace("```json", "").replace("```", "").strip()
                # 尝试找到第一个 { 和最后一个 }
                start = clean_json.find("{")
                end = clean_json.rfind("}") + 1
                if start != -1 and end != -1:
                    json_str = clean_json[start:end]
                    data = json.loads(json_str)
                    book_title = data.get("book_title", book_title)
                    preface = data.get("preface_markdown", "")
                    appendix = data.get("appendix_markdown", "")
                else:
                    print("[Warn] JSON 格式错误，尝试直接提取。")
            except Exception as e:
                print(f"[Warn] 解析失败 ({e})，将保留原始生成内容作为参考。")

        # 4. 最终拼接
        final_doc = f"# {book_title}\n\n"
        if preface: final_doc += f"{preface}\n\n---\n\n"
        final_doc += full_content
        if appendix: final_doc += f"\n---\n\n{appendix}"

        # === 核心修改：保存到本次运行的独立目录 ===
        # 文件名带上主题，更加清晰
        safe_topic = re.sub(r'[\\/*?:"<>|]', "", self.current_topic).replace(" ", "_")[:30]
        output_filename = f"Final_Tutorial_{safe_topic}.md"
        # 直接保存在 run_dir 下 (与 processed 同级)
        output_path = os.path.join(self.run_dir, output_filename)
        
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(final_doc)

        print(f"[Success] 全文生成完毕！")
        print(f"📂 本次运行归档位置: {self.run_dir}")
