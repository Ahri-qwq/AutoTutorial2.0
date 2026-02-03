import os
import json
from datetime import datetime
from typing import List, Dict, Any
from src.image_autogen.utils import save_json

class AuditManager:
    def __init__(self, output_paths: Dict[str, str], run_id: str):
        self.output_paths = output_paths
        self.run_id = run_id
        
        # 内存中维护的数据
        self.resolved_specs = []
        self.rework_items = []
        self.audit_summary = {
            "run_id": run_id,
            "generated_at": datetime.now().isoformat(),
            "total_figures_planned": 0,
            "success_count": 0,
            "failed_count": 0,
            "chapters": []
        }

    def record_plan(self, chapter_plan: Dict[str, Any]):
        """记录初始计划"""
        save_json(chapter_plan, self.output_paths['chapter_plan'])
        # 计算总计划数
        total = sum(ch.get('num_figures', 0) for ch in chapter_plan.get('chapters', []))
        self.audit_summary['total_figures_planned'] = total

    def record_full_spec(self, full_spec_data: Dict[str, Any]):
        """记录 Step 3 生成的完整 Spec"""
        save_json(full_spec_data, self.output_paths['figurespec'])

    def add_success(self, figure_id: str, chapter_id: str, result_data: Dict[str, Any]):
        """记录成功的生成"""
        self.resolved_specs.append({
            "figure_id": figure_id,
            "chapter_id": chapter_id,
            "timestamp": datetime.now().isoformat(),
            **result_data
        })
        self.audit_summary['success_count'] += 1

    def add_rework(self, 
                   figure_id: str, 
                   chapter_id: str, 
                   reason_code: str, 
                   message: str, 
                   suggested_fix: str,
                   raw_context: str = "figurespec.json"):
        """记录需要返工的项目"""
        self.rework_items.append({
            "figure_id": figure_id,
            "chapter_id": chapter_id,
            "reason_code": reason_code,
            "message": message,
            "suggested_fix": suggested_fix,
            "raw_context_ref": raw_context
        })
        self.audit_summary['failed_count'] += 1

    def finalize(self):
        """保存所有报告"""
        # 1. 保存 Resolved Spec
        resolved_data = {
            "run_id": self.run_id,
            "generated_at": datetime.now().isoformat(),
            "items": self.resolved_specs
        }
        save_json(resolved_data, self.output_paths['figurespec_resolved'])

        # 2. 保存 Rework Queue
        rework_data = {
            "run_id": self.run_id,
            "generated_at": datetime.now().isoformat(),
            "items": self.rework_items
        }
        save_json(rework_data, self.output_paths['rework_queue'])

        # 3. 保存 Audit Report (JSON)
        save_json(self.audit_summary, self.output_paths['audit_report_json'])

        # 4. 保存 Audit Report (Markdown)
        self._generate_md_report()

    def _generate_md_report(self):
        """生成人类可读的 MD 报告"""
        lines = [
            f"# Image Autogen Audit Report",
            f"- **Run ID**: {self.run_id}",
            f"- **Time**: {datetime.now().isoformat()}",
            f"- **Planned Figures**: {self.audit_summary['total_figures_planned']}",
            f"- **Success**: {self.audit_summary['success_count']}",
            f"- **Failed/Rework**: {self.audit_summary['failed_count']}",
            "",
            "## 1. Rework Queue (Failures)",
        ]

        if not self.rework_items:
            lines.append("No failures detected. All planned figures generated successfully.")
        else:
            for item in self.rework_items:
                lines.append(f"### {item['figure_id']} ({item['chapter_id']})")
                lines.append(f"- **Reason**: `{item['reason_code']}`")
                lines.append(f"- **Message**: {item['message']}")
                lines.append(f"- **Fix**: {item['suggested_fix']}")
                lines.append("")

        lines.append("")
        lines.append("## 2. Generated Figures")
        
        # 简单的表格展示
        lines.append("| Figure ID | Model | Output File |")
        lines.append("| --- | --- | --- |")
        for item in self.resolved_specs:
            # 简化路径显示
            rel_path = os.path.basename(item.get('output_file', ''))
            lines.append(f"| {item['figure_id']} | {item.get('model', 'N/A')} | {rel_path} |")

        with open(self.output_paths['audit_report_md'], 'w', encoding='utf-8') as f:
            f.write('\n'.join(lines))
