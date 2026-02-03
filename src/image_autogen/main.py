import os
import argparse
import sys
import dataclasses
from pathlib import Path

from src.image_autogen.utils import (
    get_timestamp_str, 
    find_input_md, 
    prepare_output_dirs, 
    parse_chapters,
    Chapter
)
from src.image_autogen.generator import SpecGenerator, ImageGenerator
from src.image_autogen.audit import AuditManager
from src.llm_client import LLMClient

def interactive_plan(chapters: list[Chapter], default_num: int, non_interactive: bool) -> dict:
    """交互式确定每章图片数量"""
    print("\n=== Step 2: Plan Confirmation ===")
    planned_chapters = []
    
    for ch in chapters:
        if non_interactive:
            num = default_num
            print(f"Chapter: {ch.chapter_id} - {ch.title} -> {num} figures (Auto)")
        else:
            while True:
                user_input = input(f"Chapter: {ch.chapter_id} - {ch.title} [Default: {default_num}] > ")
                if not user_input.strip():
                    num = default_num
                    break
                try:
                    num = int(user_input)
                    if num < 0:
                        raise ValueError
                    break
                except ValueError:
                    print("Please enter a valid non-negative integer.")
        
        planned_chapters.append({
            "chapter_id": ch.chapter_id,
            "chapter_title": ch.title,
            "num_figures": num,
            # 传递引用，后续生成 spec 用
            "_chapter_obj": ch 
        })
    
    return {"chapters": planned_chapters}

def main():
    parser = argparse.ArgumentParser(description="Auto-generate figures for academic tutorials.")
    parser.add_argument("--run_dir", required=True, help="Directory containing the Final_Tutorial_*.md")
    parser.add_argument("--text_model", help="Override text model name (default from config)")
    parser.add_argument("--image_model", default="gemini-3-pro-image-preview", help="Override image model name (default: gemini-3-pro-image-preview)")
    parser.add_argument("--non_interactive", action="store_true", help="Skip interactive planning")
    parser.add_argument("--default_per_chapter", type=int, default=2, help="Default figures per chapter")
    parser.add_argument("--force", action="store_true", help="Overwrite existing output directory")
    
    args = parser.parse_args()
    
    print("\n🚀 Starting Image Autogen Workflow...")
    
    # === Step 0: Setup ===
    try:
        input_md_path = find_input_md(args.run_dir)
        print(f"✅ Found input: {os.path.basename(input_md_path)}")
    except Exception as e:
        print(f"❌ Error: {e}")
        sys.exit(1)
        
    out_paths = prepare_output_dirs(args.run_dir, overwrite=args.force)
    run_id = get_timestamp_str()
    print(f"✅ Output directories prepared. Run ID: {run_id}")
    
    audit_mgr = AuditManager(out_paths, run_id)
    
    # Init LLM Client
    try:
        llm_client = LLMClient()
        print("✅ LLM Client initialized.")
    except Exception as e:
        print(f"❌ LLM Client init failed: {e}")
        sys.exit(1)
        
    spec_gen = SpecGenerator(llm_client)
    img_gen = ImageGenerator(llm_client)
    
    # === Step 1: Parse MD ===
    print("\n=== Step 1: Parsing Markdown ===")
    with open(input_md_path, 'r', encoding='utf-8') as f:
        md_text = f.read()
    
    chapters = parse_chapters(md_text)
    print(f"✅ Parsed {len(chapters)} chapters.")
    
    # === Step 2: Interactive Plan ===
    plan_data = interactive_plan(chapters, args.default_per_chapter, args.non_interactive)
    plan_data['run_id'] = run_id
    plan_data['input_md'] = os.path.basename(input_md_path)
    
    # 移除内部使用的 _chapter_obj 对象以便序列化保存
    plan_for_save = {
        **plan_data, 
        "chapters": [{k:v for k,v in ch.items() if not k.startswith('_')} for ch in plan_data['chapters']]
    }
    audit_mgr.record_plan(plan_for_save)
    print("✅ Plan saved.")
    
    # === Step 3: Generate Specs ===
    print("\n=== Step 3: Generating Figure Specs ===")
    full_specs = {
        "run_id": run_id,
        "style": {"theme": "vector_flat", "label_lang": "zh-CN"},
        "chapters": []
    }
    
    for ch_plan in plan_data['chapters']:
        num = ch_plan['num_figures']
        if num == 0:
            continue
            
        ch_obj: Chapter = ch_plan['_chapter_obj']
        
        try:
            ch_spec = spec_gen.generate_chapter_spec(
                ch_obj.chapter_id,
                ch_obj.title,
                ch_obj.content,
                ch_obj.headings,
                num,
                model_name=args.text_model
            )
            
            # Convert to dict for storage
            # dataclass to dict
            figures_dict = [dataclasses.asdict(fig) for fig in ch_spec.figures]
            
            ch_spec_dict = {
                "chapter_id": ch_spec.chapter_id,
                "chapter_title": ch_spec.chapter_title,
                "figures": figures_dict
            }
            full_specs['chapters'].append(ch_spec_dict)
            
        except Exception as e:
            print(f"⚠️ Failed to generate spec for {ch_obj.chapter_id}: {e}")
            audit_mgr.add_rework(
                f"{ch_obj.chapter_id}-ALL", 
                ch_obj.chapter_id, 
                "SPEC_GEN_FAILED", 
                str(e), 
                "Manual check required"
            )
            
    audit_mgr.record_full_spec(full_specs)
    print("✅ Figure Specs generated and saved.")
    
    # === Step 4: Generate Images ===
    print("\n=== Step 4: Generating Images ===")
    
    for ch in full_specs['chapters']:
        chapter_id = ch['chapter_id']
        for fig in ch['figures']:
            fig_id = fig['figure_id']
            
            # Simple Gatekeeping
            if not fig.get('labels_required'):
                 print(f"⚠️ Skipping {fig_id}: Missing labels_required")
                 audit_mgr.add_rework(fig_id, chapter_id, "MISSING_LABELS", "labels_required is empty", "Add labels in spec")
                 continue
                 
            # 检查 heading 唯一性 (可选，仅作警告)
            # 这里如果不重新解析全文很难检查全局唯一性，generator 里是针对章内上下文。
            # 暂且跳过严格唯一性检查，只生成图片。
            
            # 生成图片
            output_filename = f"{fig_id}.png"
            output_path = os.path.join(out_paths['images'], output_filename)
            
            # Convert dict back to FigureSpec object if needed, or just pass dict fields
            # ImageGenerator expects FigureSpec object
            from src.image_autogen.generator import FigureSpec
            fig_obj = FigureSpec.from_dict(fig)
            
            result = img_gen.generate_image(
                fig_obj,
                output_path,
                model_name=args.image_model
            )
            
            if result['status'] == 'success':
                print(f"✅ Generated: {output_filename}")
                audit_mgr.add_success(fig_id, chapter_id, result)
            else:
                print(f"❌ Failed: {output_filename}")
                audit_mgr.add_rework(fig_id, chapter_id, "IMAGE_GEN_FAILED", result.get('error', 'Unknown'), "Check prompt or API quota")

    # === Finalize ===
    audit_mgr.finalize()
    print(f"\n🎉 Workflow Completed! Reports saved to {out_paths['base']}")

if __name__ == "__main__":
    main()
