
## 1. ChapterPlan（images_autogen/chapter_plan.json）
```json
{
  "run_id": "20260127_1424xx",
  "input_md": "Final_Tutorial_xxx.md",
  "chapters": [
    {
      "chapter_id": "ch01",
      "chapter_title": "……",
      "num_figures": 2
    }
  ]
}
```

## 2. FigureSpec（images_autogen/figurespec.json）
说明：Step3 生成的结果需要被脚本消费，因此必须是严格 JSON 且结构稳定。

```json
{
  "run_id": "20260127_1424xx",
  "style": {
    "theme": "vector_flat",
    "label_lang": "zh-CN"
  },
  "chapters": [
    {
      "chapter_id": "ch01",
      "chapter_title": "……",
      "figures": [
        {
          "figure_id": "ch01-fig01",
          "title": "短图名",
          "caption": "1-3句图注，不引入原文没有的新结论。",
          "insert_after_heading": "## 原文中的某个小节标题",
          "diagram_type": "流程图",
          "source_quotes": [
            "原文摘录1",
            "原文摘录2"
          ],
          "entities_whitelist": [
            "术语A",
            "术语B"
          ],
          "labels_required": [
            "标签1",
            "标签2"
          ],
          "labels_forbidden": [
            "照片风格",
            "人物脸",
            "3D质感"
          ],
          "assertions": [
            "断言1：箭头表示……",
            "断言2：坐标轴含义是……"
          ],
          "prompt": "给图像模型的可变内容：画什么、有哪些元素与布局（不要长段文字）。",
          "negative_prompt": "可选：避免出现……"
        }
      ]
    }
  ]
}
```

## 3. ResolvedSpec（images_autogen/figurespec.resolved.json）
用于追溯与复现，建议结构：
```json
{
  "run_id": "20260127_1424xx",
  "generated_at": "2026-01-27T14:30:00",
  "items": [
    {
      "figure_id": "ch01-fig01",
      "revision": 1,
      "image_model": "gemini-3-pro-image-preview",
      "final_prompt": "STYLE_PREFIX + prompt + Avoid...",
      "output_file": "images_output/ch01-fig01.png"
    }
  ]
}
```

## 4. ReworkQueue（images_autogen/rework_queue.json）
```json
{
  "run_id": "20260127_1424xx",
  "items": [
    {
      "figure_id": "ch01-fig01",
      "chapter_id": "ch01",
      "reason_code": "HEADING_NOT_UNIQUE",
      "message": "insert_after_heading 在原文中命中多次，无法确定插入点。",
      "suggested_fix": "在 Step3 重新选择更唯一的 heading，或后续版本增加 heading_path/section_id。",
      "raw_context_ref": "images_autogen/figurespec.json"
    }
  ]
}
```
