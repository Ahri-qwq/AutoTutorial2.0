## 1. 项目目标
对 `.\data\runs\<timestamp>\Final_Tutorial_*.md` 自动生成学术教程插图，并插回 Markdown，产出“带图版教程”。插图风格必须统一为**矢量扁平风**（白底、线条清晰、少量配色、少文字标签、标签为简体中文）。

## 2. 输入与输出约束
### 2.1 输入
- 输入目录：`.\data\runs\<timestamp>\`
- 输入文件：目录内**唯一**匹配 `Final_Tutorial_*.md` 的文件（前缀固定为 `Final_Tutorial_`）。

### 2.2 输出
在同一 run 目录下创建输出目录：`.\data\runs\<timestamp>\images_autogen\`，包含：
- `images_output/`：PNG 图片输出目录。
  - 命名格式：`{chapter_id}-{figure_id}.png`（例如 `第1章-第1图.png` 或 `ch01-fig01.png`）。
- `chapter_plan.json`：交互得到的“每章生成图片数量”计划。
- `figurespec.json`：Step3 输出的 FigureSpec（严格 JSON）。
- `figurespec.resolved.json`：Step4 落盘的“最终 prompt + 文件名 + revision + 模型名”等可追溯信息。
- `audit_report.json`、`audit_report.md`：人工抽查清单。
- `rework_queue.json`：返工队列（定位失败/校验失败/高风险项）。

> **注意**：本工具只负责生成图片，不负责插回 Markdown。图片生成后由人工审查并决定是否采纳。

## 3. 配置与模型调用
### 3.1 根目录 config.yaml & 代码复用
- 直接复用 `src.llm_client.LLMClient` 进行 API 调用。
- `config.yaml` 中相关配置：
  ```yaml
  llm:
    google_api_key: "..." 
    google_base_url: "https://api.vectorengine.ai" 
    model: "gemini-3-pro-preview-thinking"
  ```
- `LLMClient` 会自动处理 `google_api_key` 和 `google_base_url` 的读取。

### 3.2 调用方式
- **文本模型**：使用 `src.llm_client.LLMClient.chat()`。
- **图像模型**：复用 `LLMClient` 实例中的 `client` (OpenAI instance) 调用 `images.generate`。
- 脚本需支持 CLI 覆盖模型：
  - `--text_model`（默认取 config.yaml 的 model）
  - `--image_model`（默认写死一个值，但必须可覆盖）

## 4. 运行方式与交互
提供 CLI：
- `python -m src.image_autogen.main --run_dir .\data\runs\<timestamp>`

运行流程：
1. 解析输入 md 得到“章节列表”。  
2. 逐章询问“本章生成几张图？”默认值 2，允许 0/1/2/3…  
3. 确认后执行 Step3–Step4。

可选支持：
- `--non_interactive`：不交互，按 `--default_per_chapter` 统一生成。
- `--only_chapter ch02` 或 `--only_figure ch02-fig01`：用于返工重跑（作为 v2 目标，可先预留接口）。

## 5. 工作流（Step0–Step4）
### Step0：发现输入与建立输出目录
- 校验 run_dir 存在。
- glob 查找 `Final_Tutorial_*.md`，必须唯一，否则报错退出。
- 创建 `images_autogen/` 与 `images_autogen/images_output/`。
- 生成 `run_id`（时间戳字符串），写入各 JSON 产物中。

### Step1：解析 Markdown 章节结构
- 解析 headings（至少支持 `#`/`##`/`###`）。
- 章划分规则：
  - 优先 H1 作为章；
  - 若无 H1，则用 H2 作为章。
- 每章得到：`chapter_id`、`chapter_title`、`chapter_md`、`headings`（章内可用 heading 列表）。

### Step2：交互确定每章出图数量
- 对每章询问 `num_figures`（默认 2）。
- 将结果写入 `chapter_plan.json`。

### Step3：生成 FigureSpec（文本模型 → 严格 JSON）
- 对每章调用文本模型生成 `num_figures` 张图的 FigureSpec。
- 必须输出严格 JSON，顶层 `{ "figures": [ ... ] }`，禁止 markdown code block。 
- 每张图必须包含：插入定位字段（仅供参考）、审计字段（source_quotes、assertions 等），以便人工抽查与返工。

### Step4：门禁校验 → 图像生成 → 落盘
- 门禁校验（失败写入 `rework_queue.json`，跳过该图生成）：
  - JSON 解析成功、字段齐全。
  - `labels_required` 非空、`source_quotes` 非空（避免无依据图）。
- 图像生成：
  - 读取 Step4 风格前缀模板（你已放在 `.\src\image_autogen\prompts\`）。
  - `final_prompt = STYLE_PREFIX + "\n\n" + figure.prompt`
  - 若 negative_prompt 非空：追加 `"\n\nAvoid: " + negative_prompt`（避免出现空 `Avoid:`）。
  - 调用图像模型生成 PNG，保存到 `images_autogen/images_output/`，命名为 `{chapter_id}-{fig_index}.png`。
- 将最终 prompt、模型名、文件名、revision、时间戳写入 `figurespec.resolved.json`。


## 6. 验收标准（Acceptance Criteria）
- AC1：能正确定位并读取唯一输入 `Final_Tutorial_*.md`。  
- AC2：能解析章节并完成交互，生成 `chapter_plan.json`。  
- AC3：能生成 `figurespec.json` 且为严格 JSON 可解析。
- AC4：能生成 PNG 图片并保存至 `images_output`，命名符合规范。  
- AC5：能输出 `audit_report.*` 与 `rework_queue.json`，支持人工抽查与返工闭环。
