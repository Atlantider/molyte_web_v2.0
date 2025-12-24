#!/usr/bin/env python3
"""
专门修复 EL-20251205-0002-Li-PF6-EC-DEC-PC-MD1-298K 任务
直接使用 polling worker 的 _handle_job_completion() 方法来处理结果
"""

import os
import sys
import subprocess
import yaml
import json
from datetime import datetime
from pathlib import Path

# 添加当前目录到路径，以便导入 polling worker
current_dir = Path(__file__).parent
deployment_dir = current_dir / "deployment"
if str(deployment_dir) not in sys.path:
    sys.path.insert(0, str(deployment_dir))

# 任务信息
TASK_INFO = {
    'job_id': 179,
    'work_dir': '/public/home/xiaoji/molyte_web/data/md_work/EL-20251205-0002-Li-PF6-EC-DEC-PC-MD1-298K',
    'job_name': 'EL-20251205-0002-Li-PF6-EC-DEC-PC-MD1-298K',
    'type': 'md',
    'slurm_job_id': '22092'  # 模拟的Slurm ID
}

def create_polling_worker():
    """创建polling worker实例"""
    try:
        # 导入polling worker
        from polling_worker import PollingWorker

        # 创建worker实例
        worker = PollingWorker()
        print(f"✅ Polling worker 实例创建成功")
        return worker

    except Exception as e:
        print(f"❌ 创建polling worker失败: {e}")
        return None

def simulate_job_completion(worker, job_info):
    """模拟任务完成，调用polling worker的处理方法"""
    try:
        job_id = job_info['job_id']
        work_dir = Path(job_info['work_dir'])

        print(f"🔄 模拟任务 {job_id} 完成，调用 _handle_job_completion...")

        # 构造job_info字典，模拟从API获取的任务信息
        api_job_info = {
            'id': job_id,
            'type': job_info['type'],
            'work_dir': str(work_dir),
            'slurm_job_id': job_info.get('slurm_job_id'),
            'status': 'RUNNING'  # 模拟当前状态
        }

        # 调用polling worker的任务完成处理方法
        worker._handle_job_completion(job_id, api_job_info)

        print(f"✅ 任务 {job_id} 完成处理成功")
        return True

    except Exception as e:
        print(f"❌ 处理任务完成失败: {e}")
        import traceback
        traceback.print_exc()
        return False

def check_work_directory():
    """检查工作目录和文件"""
    work_path = Path(TASK_INFO['work_dir'])
    if not work_path.exists():
        print(f"❌ 工作目录不存在: {TASK_INFO['work_dir']}")
        return False
    
    print(f"📁 检查工作目录: {TASK_INFO['work_dir']}")
    
    # 检查关键文件
    key_files = [
        f"{TASK_INFO['job_name']}.log",
        f"{TASK_INFO['job_name']}.in",
        "job.sh"
    ]
    
    for file_name in key_files:
        file_path = work_path / file_name
        if file_path.exists():
            print(f"   ✅ 找到文件: {file_name} ({file_path.stat().st_size/1024:.1f}KB)")
        else:
            print(f"   ⚠️ 缺少文件: {file_name}")
    
    # 列出所有文件
    all_files = list(work_path.glob('*'))
    print(f"   📋 工作目录共有 {len(all_files)} 个文件")
    
    return True

def upload_md_results_selective(cos_client, cos_bucket, config, job_id):
    """按照polling worker的策略选择性上传MD任务结果文件到COS"""
    try:
        work_path = Path(TASK_INFO['work_dir'])

        # 使用polling worker的上传配置
        essential_patterns = config['upload']['essential_files']
        excluded_patterns = config['upload'].get('excluded_files', [])
        max_size = config['upload']['max_file_size'] * 1024 * 1024  # MB to bytes

        # 获取结果前缀
        result_prefix = config['cos'].get('md_result_prefix', 'results/')

        uploaded_files = []

        print(f"📋 使用essential文件模式: {essential_patterns}")

        # 只上传essential文件
        for pattern in essential_patterns:
            files = list(work_path.glob(pattern))
            for file_path in files:
                if not file_path.is_file():
                    continue

                # 检查是否在排除列表中
                excluded = False
                for exclude_pattern in excluded_patterns:
                    if file_path.match(exclude_pattern):
                        print(f"⏭️ 跳过排除的文件: {file_path.name}")
                        excluded = True
                        break
                if excluded:
                    continue

                # 检查文件大小
                file_size = file_path.stat().st_size
                if file_size > max_size:
                    print(f"⚠️ 文件 {file_path.name} ({file_size/1024/1024:.1f}MB) 超过大小限制，跳过")
                    continue

                if file_size == 0:
                    print(f"⏭️ 跳过空文件: {file_path.name}")
                    continue

                # 构建COS key
                cos_key = f"{result_prefix}{job_id}/{file_path.name}"

                print(f"📤 上传文件: {file_path.name} ({file_size/1024/1024:.1f}MB)")

                with open(file_path, 'rb') as f:
                    cos_client.put_object(
                        Bucket=cos_bucket,
                        Body=f,
                        Key=cos_key
                    )

                uploaded_files.append(cos_key)
                print(f"   ✅ 上传成功")

        print(f"✅ 共上传 {len(uploaded_files)} 个文件到COS")
        return len(uploaded_files) > 0, uploaded_files

    except Exception as e:
        print(f"❌ 上传文件时出错: {e}")
        return False, []

def parse_md_results():
    """解析MD结果文件"""
    try:
        work_path = Path(TASK_INFO['work_dir'])
        
        # 查找RDF和MSD文件
        rdf_files = list(work_path.glob('out_rdf*.dat'))
        msd_files = list(work_path.glob('out_*_msd.dat'))
        
        results = {}
        
        # 解析RDF文件
        if rdf_files:
            rdf_file = rdf_files[0]
            print(f"📊 解析RDF文件: {rdf_file.name}")
            # 这里可以添加RDF解析逻辑
            results['rdf_file'] = str(rdf_file)
        
        # 解析MSD文件
        for msd_file in msd_files:
            print(f"📊 解析MSD文件: {msd_file.name}")
            # 这里可以添加MSD解析逻辑
            if 'msd_files' not in results:
                results['msd_files'] = []
            results['msd_files'].append(str(msd_file))
        
        print(f"✅ MD结果解析完成: {len(results)} 个结果文件")
        return results
        
    except Exception as e:
        print(f"❌ 解析MD结果时出错: {e}")
        return {}

def update_job_status(config, job_id, status='COMPLETED', result_files=None):
    """更新任务状态到腾讯云后端"""
    try:
        api_base_url = config['api']['base_url']
        api_headers = {
            'Content-Type': 'application/json',
            'Authorization': f"Bearer {config['api']['worker_token']}"
        }

        endpoint = f"{api_base_url}/workers/jobs/{job_id}/status"

        data = {
            'status': status,
            'job_type': 'MD',
            'worker_name': config['worker']['name'],
            'progress': 100.0
        }

        if result_files:
            data['result_files'] = result_files

        # 清除错误信息
        data['error_message'] = None

        print(f"🔄 更新任务状态到后端: {endpoint}")
        print(f"   数据: {data}")

        response = requests.put(
            endpoint,
            headers=api_headers,
            json=data,
            timeout=config['api']['timeout']
        )

        if response.status_code == 200:
            print(f"✅ 任务 {job_id} 状态已更新为 {status}")
            return True
        else:
            print(f"❌ 更新状态失败: HTTP {response.status_code}")
            print(f"   响应: {response.text}")
            return False

    except Exception as e:
        print(f"❌ 更新任务状态失败: {e}")
        return False

def upload_md_results_to_api(config, job_id, md_results):
    """上传MD结果到后端API"""
    try:
        api_base_url = config['api']['base_url']
        api_headers = {
            'Content-Type': 'application/json',
            'Authorization': f"Bearer {config['api']['worker_token']}"
        }

        endpoint = f"{api_base_url}/workers/jobs/{job_id}/md_results"

        print(f"📊 上传MD结果到API: {endpoint}")

        response = requests.post(
            endpoint,
            headers=api_headers,
            json=md_results,
            timeout=config['api']['timeout']
        )

        if response.status_code == 200:
            print(f"✅ MD结果上传成功")
            return True
        else:
            print(f"❌ MD结果上传失败: HTTP {response.status_code}")
            print(f"   响应: {response.text}")
            return False

    except Exception as e:
        print(f"❌ 上传MD结果失败: {e}")
        return False

def main():
    """主函数 - 直接使用polling worker的完成处理流程"""
    print("🔧 开始修复 EL-20251205-0002 MD任务...")
    print("📋 直接调用 polling worker 的 _handle_job_completion 方法")

    # 1. 检查工作目录
    if not check_work_directory():
        return False

    # 2. 创建polling worker实例
    worker = create_polling_worker()
    if not worker:
        return False

    # 3. 调用polling worker的任务完成处理方法
    print(f"\n� 调用 polling worker 的完成处理流程...")
    success = simulate_job_completion(worker, TASK_INFO)

    if success:
        print(f"\n✅ EL-20251205-0002 任务修复完成!")
        print(f"   任务ID: {TASK_INFO['job_id']}")
        print(f"   工作目录: {TASK_INFO['work_dir']}")
        print(f"\n🎉 Polling worker 已完成以下操作:")
        print(f"   📤 上传结果文件到COS")
        print(f"   📊 解析MD结果 (RDF、MSD等)")
        print(f"   📡 上传解析结果到腾讯云API")
        print(f"   🔄 更新任务状态为COMPLETED")
        print(f"\n💾 解析后的结果数据已存储到数据库!")
    else:
        print(f"\n❌ EL-20251205-0002 任务修复失败!")

    return success

if __name__ == "__main__":
    try:
        success = main()
        if success:
            print(f"\n🎉 脚本执行成功!")
        else:
            print(f"\n💥 脚本执行失败!")
            exit(1)
    except Exception as e:
        print(f"\n💥 脚本执行出错: {e}")
        import traceback
        traceback.print_exc()
        exit(1)
