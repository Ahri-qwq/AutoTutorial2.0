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
        #下面这行用于调试，修改最后一个参数为文件名去main注释掉已完成步骤即可
        #self.run_dir = os.path.join(root_dir, "data", "runs", "20260109_150115")
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
        """Step 2: 制定大纲 (增加流程检索)"""
        print("\n=== Step 2: 制定大纲 ===")
        
        # 1. 读取 Step 1 的 Metadata (基础)
        step1_path = os.path.join(self.processed_dir, "step1_enrichment.md")
        if not os.path.exists(step1_path):
            print("[Error] 请先运行 Step 1.")
            return

        with open(step1_path, 'r', encoding='utf-8') as f:
            step1_content = f.read()

        # 2. [新增] 检索计算流程 (Workflow Retrieval)
        # 目的：让大纲规划者知道"先做A，后做B"，避免逻辑断层
        workflow_query = f"{self.current_topic} 计算流程 步骤 脚本 工具 workflow"
        print(f"[Retrieving] 搜索大纲所需的流程信息: {workflow_query}")
        workflow_context = self.retriever.search(workflow_query, top_k=3)

        # 3. 混合 Context (Metadata + Workflow)
        full_context = f"{step1_content}\n\n【补充流程信息 (Workflow Context)】\n{workflow_context}"
        
        # 4. 生成大纲
        prompt_tmpl = self._load_prompt("step2_outline.txt")
        final_prompt = self._fill_prompt(prompt_tmpl, TOPIC=self.current_topic, BACKGROUND_INFO=full_context)
        
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
        """撰写单个章节 (增强检索版)"""
        print(f"\n[Drafting] 正在撰写: {title} ...")
        
        # 1. [优化] 构建更丰富的检索词
        # Base Query: 宏观原理
        base_query = f"ABACUS {title} 原理 参数"
        
        # Detail Query: 提取标题中的核心实体进行针对性检索
        # 例如: "Section 3.2: 应力-应变曲线拟合" -> 提取 "应力-应变曲线拟合"
        # 简单的清理逻辑：去掉 'Section', 'Chapter', 数字, 冒号
        clean_title = re.sub(r'(Chapter|Section|第\w+章|[\d\.:])', ' ', title).strip()
        detail_query = f"ABACUS {clean_title} INPUT参数设置 示例代码"
        
        print(f"[Retrieving] Queries: \n  1. {base_query}\n  2. {detail_query}")
        
        context_1 = self.retriever.search(base_query, top_k=3)
        context_2 = self.retriever.search(detail_query, top_k=3)
        
        # 合并检索结果，去重交给 LLM 处理
        full_context = f"{context_1}\n\n{context_2}"
        
        # 2. 填充 Prompt (注入 Special Instructions)
        prompt_tmpl = self._load_prompt("step3_drafting.txt")
        
        # 从类属性获取 Step 1 留下的锦囊妙计 (如果 run_step3 里没读取，这里用默认值)
        special_instr = getattr(self, 'current_special_instructions', '无特殊指令。')
        
        final_prompt = self._fill_prompt(
            prompt_tmpl, 
            CHAPTER_TITLE=title, 
            CHAPTER_DESCRIPTION=outline_context,
            CONTEXT=full_context,
            SPECIAL_INSTRUCTIONS=special_instr
        )
        
        # 3. 生成正文
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


        # 4. 保存
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

        full_content = []
        chapter_summaries = []
        
        for f in chapter_files:
            path = os.path.join(self.processed_dir, f)
            with open(path, 'r', encoding='utf-8') as file:
                content = file.read()
                full_content.append(content)
                # 提取摘要（简单取前300字）
                summary = content[:300].replace("\n", " ") + "..."
                chapter_summaries.append(f"- {f}: {summary}")

        # 2. RAG 检索：为了写好"前言"和"附录"，我们需要宏观知识
        print("[Retrieving] 正在搜索宏观背景知识 (学习路线/相关概念)...")
        rag_query = f"{self.current_topic} 的前置知识、进阶学习路线、相关联的物理方法、推荐阅读文献"
        background_knowledge = self.retriever.search(rag_query, top_k=6)

        # 3. 组装 Prompt
        prompt_tmpl = self._load_prompt("step4_assembly.txt")
        final_prompt = self._fill_prompt(
            prompt_tmpl,
            CHAPTER_SUMMARIES="\n".join(chapter_summaries),
            CONTEXT=background_knowledge 
        )
        print("[Agent] 正在撰写深度前言与附录...")
        # 强制要求 JSON 模式（部分模型支持，不支持则忽略）
        result = self.llm.chat(final_prompt)

        # 4. 解析 JSON (增强健壮性)
        try:
            # === 修复开始：更强的 JSON 清洗逻辑 ===
            import json
            import re
            
            # 1. 移除可能的 Markdown 代码块标记 ```json ... ```
            clean_result = re.sub(r'```json\s*', '', result)
            clean_result = re.sub(r'```\s*', '', clean_result)
            
            # 2. 移除可能导致解析失败的控制字符（关键步骤！）
            # 允许换行符 \n，但移除其他 ASCII 控制字符 (0-31, 127)
            # 注意：JSON 标准要求字符串内的换行必须是 \\n，而不是实际的换行符
            # 这里我们尝试简单粗暴地清理首尾空白
            clean_result = clean_result.strip()

            # 3. 尝试解析
            try:
                data = json.loads(clean_result)
            except json.JSONDecodeError:
                # 如果失败，尝试修复常见的 "非法换行" 问题
                # 很多 LLM 喜欢在 JSON value 里直接换行，这是非法的
                # 我们尝试把未转义的换行符替换为 \\n
                # 这是一个简化的修复，不能覆盖所有情况，但很有用
                clean_result = clean_result.replace('\n', '\\n')
                data = json.loads(clean_result)

            book_title = data.get("book_title", f"{self.current_topic} 实战教程")
            preface = data.get("preface_markdown", "")
            appendix = data.get("appendix_markdown", "")
            # === 修复结束 ===

        except Exception as e:
            print(f"[Warn] 解析失败 ({str(e)})，将使用默认格式。")
            book_title = f"{self.current_topic} 实战教程"
            preface = "## 前言\n（自动生成的前言解析失败，请手动补充）"
            appendix = "## 附录\n（自动生成的附录解析失败，请手动补充）"
            # 如果解析完全失败，保留原始输出来 debug
            with open(os.path.join(self.processed_dir, "step4_raw_output.txt"), "w", encoding="utf-8") as f:
                f.write(result)

        # 5. 最终拼接
        final_markdown = f"# {book_title}\n\n"
        final_markdown += f"{preface}\n\n"
        final_markdown += "---\n\n"
        final_markdown += "\n\n".join(full_content)
        final_markdown += "\n\n---\n\n"
        final_markdown += f"{appendix}\n"
        # 文件名带上主题，更加清晰
        safe_topic = re.sub(r'[\\/*?:"<>|]', "", self.current_topic).replace(" ", "_")[:30]
        output_path = os.path.join(self.run_dir, f"Final_Tutorial_{safe_topic}.md")
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(final_markdown)
            
        print(f"[Success] 全文生成完毕！已保存至: {output_path}")
