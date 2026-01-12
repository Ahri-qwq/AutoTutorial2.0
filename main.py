import os
import sys
from src.pipeline import AutoTutorialPipeline

def main():
    # 获取项目根目录
    root_dir = os.path.dirname(os.path.abspath(__file__))
    
    # 初始化流水线
    pipeline = AutoTutorialPipeline(root_dir)
    
    print("==========================================")
    print("   AutoTutorial 2.0 - ABACUS 教程生成器   ")
    print("==========================================")
    
    # 获取用户输入的主题
    topic = input("请输入教程主题: ").strip()
    if not topic:
        print("主题不能为空，程序退出。")
        return

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
