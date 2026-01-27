import os
import sys
import argparse  # 新增
from src.pipeline import AutoTutorialPipeline

def main():
    # 获取项目根目录
    root_dir = os.path.dirname(os.path.abspath(__file__))

    # === 修改 1: 使用 argparse 处理命令行参数 ===
    parser = argparse.ArgumentParser(description="AutoTutorial 2.0 - ABACUS 教程生成器")
    parser.add_argument("topic", nargs="?", help="教程主题") # 使其可选，支持交互式输入
    parser.add_argument("-i", "--input", help="重要案例输入文件路径 (docx/md)", default=None)
    args = parser.parse_args()

    # 初始化流水线
    pipeline = AutoTutorialPipeline(root_dir)

    print("==========================================")
    print("  AutoTutorial 2.0 - ABACUS 教程生成器   ")
    print("==========================================")

    # 获取主题 (优先使用命令行参数，否则交互输入)
    topic = args.topic
    if not topic:
        topic = input("请输入教程主题: ").strip()
    
    if not topic:
        print("主题不能为空，程序退出。")
        return

    # === 修改 2: 处理输入案例 ===
    if args.input:
        input_path = os.path.abspath(args.input)
        if os.path.exists(input_path):
            print(f"[Info] 检测到输入案例: {input_path}")
            pipeline.set_case_study(input_path)
        else:
            print(f"[Warn] 输入文件不存在: {input_path}，将忽略该参数。")

    # 按顺序执行步骤
    try:
        # Step 1: 知识预研
        pipeline.run_step1(topic)
        # Step 2: 大纲规划
        pipeline.run_step2()
        # Step 3: 分章撰写
        pipeline.run_step3()
        # Step 4: 组装
        pipeline.run_step4()

    except KeyboardInterrupt:
        print("\n[Warn] 用户中断了程序。")
    except Exception as e:
        print(f"\n[Error] 程序运行出错: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
