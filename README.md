# AutoTutorial 2.0

**AutoTutorial 2.0** 是一个专为 **ABACUS** 第一性原理计算软件设计的自动化教程生成系统。它结合了 **RAG (检索增强生成)** 技术与 **LLM (大语言模型)**，通过“知识预研 -> 大纲规划 -> 章节撰写 -> 全文组装”的四步流水线，自动产出高质量、包含实战参数配置的物理教程。

> **2026年1月 重要更新**：
> 本月版本实现了生成引擎的重大升级，生成部分从 Qwen 切换至 **Gemini (Pro/Flash)** 系列，并引入了“案例驱动生成”模式与分步温度控制策略，大幅提升了教程的逻辑严谨性与实战参考价值。

***

## 🌟 核心特性 (Key Features)

*   **🧠 混合双引擎架构 (Hybrid Engine)**
    *   **生成侧 (Right Brain)**：接入 **Google Gemini-3-Pro/Flash** (通过 OpenAI 兼容接口)，利用其强大的长文本逻辑与推理能力撰写教程。
    *   **检索侧 (Left Brain)**：保留 **Qwen Embedding (DashScope)**，利用其优秀的中文语义向量能力维护本地知识库。
*   **📂 案例驱动生成 (Case-Driven)**
    *   支持用户上传现成的 `.md` 或 `.docx` 案例文件。系统会自动提取案例中的核心参数、输入文件结构，并将其注入到大纲设计与章节撰写中，生成“高度复刻实战”的定制化教程。
*   **🌡️ 动态温度策略 (Step-wise Temperature)**
    *   针对不同生成阶段应用差异化温度，平衡严谨性与创造性：
        *   **Step 1 (0.4)**: 精准提取物理事实。
        *   **Step 2 (0.6)**: 创造性设计非模板化大纲。
        *   **Step 3 (0.4)**: 严格遵循参数设定，杜绝幻觉。
        *   **Step 4 (0.7)**: 润色文字，提升阅读体验。

***

## 🛠️ 快速开始 (Quick Start)

### 1. 环境配置
确保已安装 Python 3.8+，并安装依赖：
```bash
pip install openai dashscope chromadb pyyaml python-docx
```

### 2. 配置文件 (config.yaml)
在项目根目录创建或修改 `config.yaml`，配置混合双引擎：

```yaml
llm:
  # --- 检索侧 (Qwen Embedding) ---
  api_key: "sk-xxxxxxxxxxxx"  # 阿里云 DashScope Key
  
  # --- 生成侧 (Gemini via OpenAI Protocol) ---
  google_api_key: "sk-xxxxxxxxxxxx" # 中转站 Key 或 Google Key
  google_base_url: "https://api.vectorengine.ai/v1" # 中转站 API 地址
  model: "gemini-3-pro-preview-thinking" # 默认主力模型
  
temperatures:
  step1: 0.4
  step2: 0.6
  step3: 0.4
  step4: 0.7
```

***

## 📚 功能指南 (User Guide)

### 1. 知识库管理 (Knowledge Base)
在使用生成功能前，需确保本地 ChromaDB 向量库中有足够的 ABACUS 资料。

*   **资料准备**：将 PDF/Markdown/Txt 资料放入 `data/knowledge_source` 目录。
*   **执行入库**：
    ```bash
    python src/knowledge_manager.py
    ```
    *脚本会自动清洗数据、切片、调用 Qwen Embedding 计算向量并存入本地数据库。*

### 2. 标准教程生成
直接运行主程序，输入主题即可开始全自动生成。

```bash
python main.py
```
*   **交互提示**：`请输入教程主题: ` (例如：ABACUS 弹性常数计算)
*   **产出**：最终 Markdown 文件将保存在 `data/runs/<timestamp>/` 目录下。

### 3. 案例驱动生成 (Advanced) ✨
如果你手头有一个现成的案例（例如“硅的能带计算脚本”），希望生成的教程完全围绕这个案例展开：

1.  **准备案例**：将 `.md` 或 `.docx` 文件放入 `data/input/` 目录，用命令行启动项目。
2.  **运行生成**：
  例如主题为：ABACUS 弹性常数计算，案例文件为.\data\input\使用abacustest计算晶体的弹性性质.docx，在终端运行：
    ```bash
    python main.py "ABACUS 弹性常数计算" --input .\data\input\使用abacustest计算晶体的弹性性质.docx
    ```
3.  **效果**：系统会将该案例内容注入到 **Step 2 (大纲)** 和 **Step 3 (撰写)** 的上下文中，强制模型参考该案例的参数设置和文件结构进行讲解。

***

## ⚙️ 架构流水线 (Pipeline Workflow)

AutoTutorial 2.0 沿用经典的四步流水线，但在底层模型上进行了重构：

| 阶段 | 任务 | 模型配置 | 温度 | 作用 |
| :--- | :--- | :--- | :--- | :--- |
| **Step 1** | **知识预研** | Gemini-Flash | 0.4 | 从向量库检索 Top-8 资料，清洗并总结物理原理。 |
| **Step 2** | **大纲规划** | Gemini-Thinking | 0.6 | **(注入案例)** 结合预研知识与用户案例，设计逻辑连贯的实战大纲。 |
| **Step 3** | **分章撰写** | Gemini-Thinking | 0.4 | **(注入案例)** 对照案例参数，撰写具体的输入文件配置与结果分析。 |
| **Step 4** | **全文组装** | Gemini-Flash | 0.7 | 编写前言/附录，统一文风，生成最终 Markdown。 |

***

## 📊 成本估计

本月测试的推荐模型组合如下（价格基于 2026/01 中转站费率）：

*   **主力模型 (Step 1/2/3)**: `gemini-3-pro-preview-thinking`
    *   费率: Input ￥0.6/M, Output ￥3.6/M
*   **速推模型 (Step 4)**: `gemini-3-flash-preview`
    *   费率: Input ￥0.15/M, Output ￥0.9/M

**单次运行预估成本**: `~0.2 元` (包含约 40k Input tokens, 15k Output tokens)

***

## 📅 月度更新日志 (2026.01)

### 🚀 Major Changes
1.  **生成引擎迁移**: 完成从 Qwen 到 **Gemini** 的迁移，解决了长文本逻辑连贯性问题。
2.  **案例投喂功能**: 新增 `data/input` 接口，支持“带料加工”，解决了通用教程与特定算例不匹配的痛点。

### ⚡ Optimizations
1.  **温度分层**: 实施 Step-wise 温度控制 (0.4/0.6/0.4/0.7)，显著降低了技术细节的幻觉率。
2.  **代码重构**: 依据 Master Plan 重构了 `pipeline.py`，提升了代码的可维护性与扩展性。
3.  **文件名清洗**: 修复了 Windows 下文件名含换行符导致保存失败的 Bug。