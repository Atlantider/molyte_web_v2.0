#!/usr/bin/env python3
"""
使用 polling worker 的原有函数来修复任务
直接调用 _handle_job_completion() 来处理结果
"""

import os
import sys
from pathlib import Path

# 添加deployment目录到路径，以便导入 polling worker
current_dir = Path(__file__).parent
deployment_dir = current_dir / "deployment"
if str(deployment_dir) not in sys.path:
    sys.path.insert(0, str(deployment_dir))

# 需要修复的任务列表
TASKS_TO_FIX = [
    {
        'job_id': 179,
        'work_dir': '/public/home/xiaoji/molyte_web/data/md_work/EL-20251205-0002-Li-PF6-EC-DEC-PC-MD1-298K',
        'job_name': 'EL-20251205-0002-Li-PF6-EC-DEC-PC-MD1-298K',
        'type': 'md',
        'slurm_job_id': '22092'
    },
    {
        'job_id': 178,
        'work_dir': '/public/home/xiaoji/molyte_web/data/md_work/EL-20251205-0001-Li-PF6-EC-DEC-PC-MD1-298K',
        'job_name': 'EL-20251205-0001-Li-PF6-EC-DEC-PC-MD1-298K',
        'type': 'md',
        'slurm_job_id': '22091'
    }
]

def create_polling_worker():
    """创建polling worker实例"""
    try:
        # 导入polling worker
        from polling_worker import PollingWorker

        # 创建worker实例，指定正确的配置文件
        config_path = "deployment/polling_worker_config_tencent.yaml"
        worker = PollingWorker(config_path=config_path)
        print(f"✅ Polling worker 实例创建成功")
        print(f"   配置文件: {config_path}")
        return worker

    except Exception as e:
        print(f"❌ 创建polling worker失败: {e}")
        import traceback
        traceback.print_exc()
        return None

def check_work_directory(work_dir):
    """检查工作目录"""
    try:
        work_path = Path(work_dir)
        
        if not work_path.exists():
            print(f"❌ 工作目录不存在: {work_dir}")
            return False
        
        print(f"📁 检查工作目录: {work_dir}")
        
        # 统计总文件数
        all_files = list(work_path.iterdir())
        print(f"   📋 工作目录共有 {len(all_files)} 个文件")
        
        return True
        
    except Exception as e:
        print(f"❌ 检查工作目录失败: {e}")
        return False

def simulate_job_completion(worker, task_info):
    """模拟任务完成，调用polling worker的处理方法"""
    try:
        job_id = task_info['job_id']
        work_dir = Path(task_info['work_dir'])
        
        print(f"🔄 模拟任务 {job_id} 完成，调用 _handle_job_completion...")
        
        # 构造job_info字典，模拟从API获取的任务信息
        api_job_info = {
            'id': job_id,
            'type': task_info['type'],
            'work_dir': str(work_dir),
            'slurm_job_id': task_info.get('slurm_job_id'),
            'status': 'RUNNING'  # 模拟当前状态
        }
        
        # 调用polling worker的任务完成处理方法
        # 这会执行：
        # 1. _upload_results_to_oss() - 上传文件到COS
        # 2. _parse_md_results() - 解析MD结果
        # 3. _upload_md_results() - 上传解析结果到API
        # 4. _update_job_status() - 更新任务状态
        worker._handle_job_completion(job_id, api_job_info)
        
        print(f"✅ 任务 {job_id} 完成处理成功")
        return True
        
    except Exception as e:
        print(f"❌ 处理任务完成失败: {e}")
        import traceback
        traceback.print_exc()
        return False

def fix_single_task(worker, task_info):
    """修复单个任务"""
    print(f"\n{'='*60}")
    print(f"🔧 开始修复任务 {task_info['job_id']}: {task_info['job_name']}")
    print(f"   类型: {task_info['type'].upper()}")
    print(f"   工作目录: {task_info['work_dir']}")
    
    # 1. 检查工作目录
    if not check_work_directory(task_info['work_dir']):
        print(f"❌ 任务 {task_info['job_id']} 修复失败: 工作目录不存在")
        return False
    
    # 2. 调用polling worker的任务完成处理方法
    print(f"\n🔄 调用 polling worker 的完成处理流程...")
    success = simulate_job_completion(worker, task_info)
    
    if success:
        print(f"✅ 任务 {task_info['job_id']} 修复完成!")
        print(f"\n🎉 Polling worker 已完成以下操作:")
        print(f"   📤 上传结果文件到COS")
        print(f"   📊 解析MD结果 (RDF、MSD等)")
        print(f"   📡 上传解析结果到腾讯云API")
        print(f"   🔄 更新任务状态为COMPLETED")
        print(f"\n💾 解析后的结果数据已存储到数据库!")
    else:
        print(f"❌ 任务 {task_info['job_id']} 修复失败!")
    
    return success

def main():
    """主函数 - 批量修复任务"""
    print("🔧 开始使用 polling worker 修复任务...")
    print(f"📋 共有 {len(TASKS_TO_FIX)} 个任务需要修复")
    
    # 1. 创建polling worker实例
    worker = create_polling_worker()
    if not worker:
        return False
    
    # 2. 批量修复任务
    success_count = 0
    failed_count = 0
    
    for task_info in TASKS_TO_FIX:
        try:
            if fix_single_task(worker, task_info):
                success_count += 1
            else:
                failed_count += 1
        except Exception as e:
            print(f"❌ 任务 {task_info['job_id']} 修复时发生异常: {e}")
            failed_count += 1
    
    # 3. 总结
    print(f"\n{'='*60}")
    print(f"🎉 批量修复完成!")
    print(f"   ✅ 成功修复: {success_count} 个任务")
    print(f"   ❌ 修复失败: {failed_count} 个任务")
    print(f"   📊 总计: {len(TASKS_TO_FIX)} 个任务")
    
    if success_count > 0:
        print(f"\n🎉 成功修复的任务:")
        print(f"   📤 文件已上传到COS")
        print(f"   📊 MD结果已解析 (RDF、MSD等)")
        print(f"   💾 解析后的数据已存储到数据库")
        print(f"   🔄 任务状态已更新为COMPLETED")
        print(f"\n🌐 前端应该可以看到这些任务的完整结果数据了!")
    
    return success_count > 0

if __name__ == "__main__":
    try:
        success = main()
        if success:
            print(f"\n🎉 脚本执行成功!")
            sys.exit(0)
        else:
            print(f"\n💥 脚本执行失败!")
            sys.exit(1)
    except KeyboardInterrupt:
        print(f"\n⚠️ 用户中断执行")
        sys.exit(1)
    except Exception as e:
        print(f"\n💥 脚本执行时发生异常: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
