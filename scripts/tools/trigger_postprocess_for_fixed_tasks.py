#!/usr/bin/env python3
"""
为已修复的任务手动触发后处理任务
确保RDF、MSD等数据被正确解析并存储到数据库
"""

import sys
import os
from pathlib import Path

# 添加后端模块路径
backend_path = Path(__file__).parent / "backend"
if str(backend_path) not in sys.path:
    sys.path.insert(0, str(backend_path))

# 已修复的任务列表
FIXED_TASKS = [
    {
        'id': 179,
        'type': 'md',
        'name': 'EL-20251205-0002-Li-PF6-EC-DEC-PC-MD1-298K'
    },
    {
        'id': 178,
        'type': 'md', 
        'name': 'EL-20251205-0001-Li-PF6-EC-DEC-PC-MD1-298K'
    }
]

def load_config():
    """加载polling worker配置"""
    try:
        config_path = Path('/public/home/xiaoji/molyte_web/deployment/polling_worker_config_tencent.yaml')
        if not config_path.exists():
            print(f"❌ 配置文件不存在: {config_path}")
            return None
        
        with open(config_path, 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)
        
        return config
        
    except Exception as e:
        print(f"❌ 加载配置文件失败: {e}")
        return None

def trigger_postprocess(config, job_id, job_name):
    """触发单个任务的后处理"""
    try:
        api_base_url = config['api']['base_url']
        api_headers = {
            'Content-Type': 'application/json',
            'Authorization': f"Bearer {config['api']['worker_token']}"
        }
        
        # 使用后端API的触发后处理端点
        endpoint = f"{api_base_url}/jobs/{job_id}/trigger_postprocess"
        
        print(f"🔄 触发任务 {job_id} ({job_name}) 的后处理...")
        print(f"   端点: {endpoint}")
        
        response = requests.post(
            endpoint,
            headers=api_headers,
            timeout=config['api']['timeout']
        )
        
        if response.status_code == 200:
            result = response.json()
            print(f"✅ 后处理任务已触发")
            print(f"   任务ID: {result.get('task_id', 'N/A')}")
            print(f"   消息: {result.get('message', 'N/A')}")
            return True
        else:
            print(f"❌ 触发后处理失败: HTTP {response.status_code}")
            print(f"   响应: {response.text}")
            return False
            
    except Exception as e:
        print(f"❌ 触发后处理失败: {e}")
        return False

def check_task_status(config, job_id):
    """检查任务状态"""
    try:
        api_base_url = config['api']['base_url']
        api_headers = {
            'Content-Type': 'application/json',
            'Authorization': f"Bearer {config['api']['worker_token']}"
        }
        
        endpoint = f"{api_base_url}/jobs/{job_id}"
        
        response = requests.get(
            endpoint,
            headers=api_headers,
            timeout=config['api']['timeout']
        )
        
        if response.status_code == 200:
            job_data = response.json()
            status = job_data.get('status', 'UNKNOWN')
            print(f"📊 任务 {job_id} 当前状态: {status}")
            return status
        else:
            print(f"⚠️ 无法获取任务 {job_id} 状态: HTTP {response.status_code}")
            return None
            
    except Exception as e:
        print(f"❌ 检查任务状态失败: {e}")
        return None

def main():
    """主函数"""
    print("🔧 开始为已修复的MD任务触发后处理...")
    print(f"📋 共有 {len(FIXED_TASKS)} 个MD任务需要触发后处理")
    
    # 1. 加载配置
    config = load_config()
    if not config:
        return False
    
    # 2. 为每个MD任务触发后处理
    success_count = 0
    failed_count = 0
    
    for task in FIXED_TASKS:
        if task['type'] != 'md':
            print(f"⏭️ 跳过非MD任务: {task['id']} ({task['type']})")
            continue
            
        print(f"\n{'='*60}")
        print(f"🔄 处理任务 {task['id']}: {task['name']}")
        
        # 检查当前状态
        current_status = check_task_status(config, task['id'])
        
        if current_status == 'COMPLETED':
            # 触发后处理
            if trigger_postprocess(config, task['id'], task['name']):
                success_count += 1
                print(f"✅ 任务 {task['id']} 后处理触发成功")
            else:
                failed_count += 1
                print(f"❌ 任务 {task['id']} 后处理触发失败")
        else:
            print(f"⚠️ 任务 {task['id']} 状态不是COMPLETED ({current_status})，跳过")
            failed_count += 1
    
    # 3. 总结
    print(f"\n{'='*60}")
    print(f"🎉 后处理触发完成!")
    print(f"   ✅ 成功触发: {success_count} 个任务")
    print(f"   ❌ 触发失败: {failed_count} 个任务")
    print(f"   📊 总计: {len([t for t in FIXED_TASKS if t['type'] == 'md'])} 个MD任务")
    
    if success_count > 0:
        print(f"\n📊 后处理任务已在后台执行，将自动解析RDF、MSD等数据并存储到数据库")
        print(f"💡 你可以在前端查看任务详情页面来确认后处理是否完成")
    
    return success_count > 0

if __name__ == "__main__":
    try:
        success = main()
        if success:
            print(f"\n🎉 脚本执行成功!")
        else:
            print(f"\n💥 脚本执行失败!")
    except KeyboardInterrupt:
        print(f"\n⚠️ 用户中断执行")
    except Exception as e:
        print(f"\n💥 脚本执行时发生异常: {e}")
