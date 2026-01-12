# AutoTutorial 2.0 — ABACUS 教程自动生成器

AutoTutorial 2.0 是一个“向量数据库（ChromaDB）+ RAG 检索 + LLM 写作”的自动化流水线，用于根据给定主题生成结构严谨的 ABACUS 实战教程。

## 功能概览

- Step 1：知识预研（从本地知识库检索，产出结构化元数据 `step1_enrichment.md`）。
- Step 2：定制大纲（读取 Step 1，同时补充检索“流程/工作流”信息，生成 `step2_outline.md`）。
- Step 3：分章撰写（针对每章进行更精准的检索，生成 `01_chapter_*.md`、`02_chapter_*.md`…）。
- Step 4：全文组装（基于各章节摘要生成前言与附录，并拼装最终教程）。

## 📁 项目结构
```
AutoTutorial/
├── main.py # 项目入口文件
├── config.yaml.example # 配置文件模板
├── requirements.txt # Python 依赖包列表
├── README.md # 项目说明文档
│
├── data/ # 数据目录
│ ├── chroma_db/ # 向量化数据库
│ ├── knowledge_source/ # 数据库源/数据库当前内容归档
│ ├── knowledge_add/ # 数据库待入库源
│ └── runs/ # 运行生成文件
│
├── prompts/ # LLM 提示词模板
│ ├── step1_knowledge_graph.txt
│ ├── step2_outline.txt
│ ├── step3_drafting.txt
│ └── step4_assembly.txt
│
└── src/ # 源代码模块
  ├── build_knowledge_base.py
  ├── append_knowledge_base.py
  ├── llm_client.py
  ├── pipeline.py
  └── retriever.py.
```

## 环境准备

### 1) Python 与依赖
建议使用 Python 3.9+ 并创建虚拟环境。

安装依赖（示例）：
```bash
pip install chromadb dashscope pyyaml python-docx
```

其中：
- `chromadb`：向量库持久化与检索。
- `dashscope`：调用 Qwen Embedding API 生成向量（用于入库与检索）。
- `python-docx`：可选，用于读取 `.docx` 知识源文件。

### 2) 配置 API Key
本项目的 Embedding 与（可能的）LLM 调用依赖 API Key。

支持两种方式之一：
- 在 `config.yaml` 中配置（推荐）。
- 或设置环境变量 `DASHSCOPE_API_KEY`。

`config.yaml` 示例（字段名以你本地实现为准）：
```yaml
llm:
  api_key: "YOUR_DASHSCOPE_API_KEY"
```

## 构建知识库（ChromaDB）

1. 将资料放入 `data/knowledge_source/`（支持 `md/txt/docx`）。
2. 执行构建脚本：

```bash
python build_knowledge_base.py
```

脚本会：
- 扫描 `data/knowledge_source/` 下所有文件。
- 按固定 chunk 规则切分文本并写入 ChromaDB。
- 将数据库保存到 `data/chroma_db/`。

> 注意：Embedding API 存在批量大小限制，脚本中已做了 batch 分片与简单限流处理。

## 生成教程（主流程）

运行：
```bash
python main.py
```

根据提示输入教程主题，例如：
- `基于ABACUS的弹性常数计算方法与实践`
- `基于 ATST-Tools 的过渡态搜索与验证（NEB/AutoNEB）`

`main.py` 会按顺序执行 Step 1–4，所有中间文件会保存在本次运行的 `data/runs/<timestamp>/processed/` 下（最终文件位置取决于 `pipeline.py` 的输出路径设置）。

## Prompt 说明（四阶段）

- `prompts/step1_enrich.txt`：让模型输出结构化元数据（概念、关键参数、接口、特殊写作指令、坑点）。
- `prompts/step2_outline.txt`：基于 Step 1 的素材生成“量体裁衣”的大纲（带 HTML 注释标记用于切分）。
- `prompts/step3_drafting.txt`：基于章节标题 + 检索到的上下文，输出章节正文（含示例与解释）。
- `prompts/step4_assembly.txt`：基于章节摘要生成书名、前言、附录（通常以 JSON 输出，便于程序解析）。