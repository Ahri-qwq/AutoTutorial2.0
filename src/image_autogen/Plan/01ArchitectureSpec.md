## 1. 代码位置（硬约束）
所有新建代码必须放在：`.\src\image_autogen\`
严禁修改 `src/image_autogen/` 之外的任何代码。

## 2. 外部依赖
- `src.llm_client.LLMClient`: 用于所有 API 调用（Text & Image）。
  - `LLMClient` 已封装 `config.yaml` 读取逻辑，自动处理 `google_api_key` 和 `google_base_url`。

## 3. 模块划分（集成简化版）
为了减少文件数量并提高内聚性，建议采用以下结构：

- `main.py`
  - **CLI 入口**：解析参数 `--run_dir`, `--text_model`, `--image_model` 等。
  - **流程编排**：Step0 -> Step4 顺序执行。
  - **交互逻辑**：处理用户对每章图片数量的确认。

- `utils.py`
  - **IO 操作**：`find_input_md`, `prepare_output_dirs`。
  - **MD 解析**：`parse_chapters(md_text)`, `find_heading_positions`。
  - **通用工具**：JSON 读写、时间戳生成。

- `generator.py`
  - **数据模型**：`FigureSpec`, `ChapterSpec` (使用 dataclass 或 pydantic)。
  - **LLM 交互**：
    - `generate_spec(chapter_info, model_name)`: 调用 `LLMClient.chat()` 并解析 JSON。
    - `generate_image(prompt, model_name, output_path)`: 复用 `LLMClient.client.images.generate()`。
  - **校验逻辑**：检查 spec 字段完整性、heading 唯一性。

- `audit.py`
  - **报告生成**：生成 `audit_report.md`, `rework_queue.json`。
  - **记录追踪**：写入 `figurespec.resolved.json`。

## 4. 关键工程策略
- **API 复用**：
  - 在 `generator.py` 中：`from src.llm_client import LLMClient`。
  - 实例化：`llm = LLMClient()` (自动读取 config)。
  - 图像生成：利用 `llm.client` (OpenAI instance) 的 `images.generate` 接口。
- **输出规范**：
  - 图片路径：`images_autogen/images_output/{chapter_id}-{fig_index}.png`。
  - 不再执行 MD 插入操作。
- **JSON 解析增强**：
  - 文本模型可能返回 Markdown code block，需要在 `generator.py` 中增加 strip code block 的逻辑 (````json ... ````)。
