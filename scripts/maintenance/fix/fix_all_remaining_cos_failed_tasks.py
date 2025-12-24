#!/usr/bin/env python3
"""
批量修复所有因COS问题而失败的任务
模拟 polling worker 的完成处理流程
"""

import os
import sys
import subprocess
import yaml
import json
import requests
from datetime import datetime
from pathlib import Path

# 导入腾讯云COS SDK
try:
    from qcloud_cos import CosConfig, CosS3Client
except ImportError:
    print("请安装腾讯云 COS SDK: pip install cos-python-sdk-v5")
    exit(1)

# 需要修复的任务列表 (只包含工作目录存在的任务)
FAILED_TASKS = [
    {
        'id': 178,
        'type': 'md',
        'work_dir': '/public/home/xiaoji/molyte_web/data/md_work/EL-20251205-0001-Li-PF6-EC-DEC-PC-MD1-298K',
        'name': 'EL-20251205-0001-Li-PF6-EC-DEC-PC-MD1-298K'
    },
    {
        'id': 1836,
        'type': 'qc',
        'work_dir': '/public/home/xiaoji/molyte_web/data/qc_work/QC-1836-MD176_dimer_PF6_1507',
        'name': 'QC-1836-MD176_dimer_PF6_1507'
    },
    {
        'id': 1838,
        'type': 'qc',
        'work_dir': '/public/home/xiaoji/molyte_web/data/qc_work/QC-1838-MD176_intermediate_DEC_1_PF6_1_1507',
        'name': 'QC-1838-MD176_intermediate_DEC_1_PF6_1_1507'
    },
    {
        'id': 1839,
        'type': 'qc',
        'work_dir': '/public/home/xiaoji/molyte_web/data/qc_work/QC-1839-MD176_intermediate_DEC_2_1507',
        'name': 'QC-1839-MD176_intermediate_DEC_2_1507'
    },
    {
        'id': 1840,
        'type': 'qc',
        'work_dir': '/public/home/xiaoji/molyte_web/data/qc_work/QC-1840-MD176_intermediate_DEC_2_PF6_1_1507',
        'name': 'QC-1840-MD176_intermediate_DEC_2_PF6_1_1507'
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

def init_cos_client(config):
    """初始化COS客户端"""
    try:
        if 'cos' not in config:
            print("❌ 配置文件中未找到COS配置")
            return None, None
        
        cos_config = config['cos']
        cos_cfg = CosConfig(
            Region=cos_config['region'],
            SecretId=cos_config['secret_id'],
            SecretKey=cos_config['secret_key'],
            Scheme='https'
        )
        
        cos_client = CosS3Client(cos_cfg)
        cos_bucket = cos_config['bucket']
        
        print(f"✅ COS客户端初始化成功")
        return cos_client, cos_bucket
        
    except Exception as e:
        print(f"❌ 初始化COS客户端失败: {e}")
        return None, None

def check_work_directory(work_dir):
    """检查工作目录"""
    try:
        work_path = Path(work_dir)
        if not work_path.exists():
            print(f"❌ 工作目录不存在: {work_dir}")
            return False
        
        files = list(work_path.iterdir())
        print(f"📁 工作目录: {work_dir}")
        print(f"   📋 共有 {len(files)} 个文件")
        return True
        
    except Exception as e:
        print(f"❌ 检查工作目录失败: {e}")
        return False

def upload_results_selective(cos_client, cos_bucket, config, job_id, work_dir, job_type):
    """按照polling worker的策略选择性上传结果文件到COS"""
    try:
        work_path = Path(work_dir)
        
        # 使用polling worker的上传配置
        essential_patterns = config['upload']['essential_files']
        excluded_patterns = config['upload'].get('excluded_files', [])
        max_size = config['upload']['max_file_size'] * 1024 * 1024  # MB to bytes
        
        # 获取结果前缀
        if job_type.lower() == 'md':
            result_prefix = config['cos'].get('md_result_prefix', 'MD_results/')
        else:
            result_prefix = config['cos'].get('qc_result_prefix', 'QC_results/')
        
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

def parse_md_results(work_dir):
    """解析MD结果"""
    try:
        work_path = Path(work_dir)
        results = {}
        
        # 查找RDF文件
        rdf_files = list(work_path.glob("out_rdf.dat"))
        if rdf_files:
            print(f"📊 解析RDF文件: {rdf_files[0].name}")
            results['rdf'] = str(rdf_files[0])
        
        # 查找MSD文件
        msd_files = list(work_path.glob("out_*_msd.dat"))
        for msd_file in msd_files:
            print(f"📊 解析MSD文件: {msd_file.name}")
            results[f'msd_{msd_file.stem}'] = str(msd_file)
        
        print(f"✅ MD结果解析完成: {len(results)} 个结果文件")
        return results
        
    except Exception as e:
        print(f"❌ 解析MD结果失败: {e}")
        return {}

def parse_qc_results(work_dir):
    """解析QC结果"""
    try:
        work_path = Path(work_dir)
        results = {}
        
        # 查找输出文件
        out_files = list(work_path.glob("*.out"))
        log_files = list(work_path.glob("*.log"))
        fchk_files = list(work_path.glob("*.fchk"))
        
        if out_files:
            results['output_file'] = str(out_files[0])
        if log_files:
            results['log_file'] = str(log_files[0])
        if fchk_files:
            results['fchk_file'] = str(fchk_files[0])
        
        print(f"✅ QC结果解析完成: {len(results)} 个结果文件")
        return results
        
    except Exception as e:
        print(f"❌ 解析QC结果失败: {e}")
        return {}

def upload_results_to_api(config, job_id, job_type, results):
    """上传结果到后端API"""
    try:
        api_base_url = config['api']['base_url']
        api_headers = {
            'Content-Type': 'application/json',
            'Authorization': f"Bearer {config['api']['worker_token']}"
        }
        
        if job_type.lower() == 'md':
            endpoint = f"{api_base_url}/workers/jobs/{job_id}/md_results"
        else:
            endpoint = f"{api_base_url}/workers/jobs/{job_id}/qc_results"
        
        print(f"📊 上传{job_type.upper()}结果到API: {endpoint}")
        
        response = requests.post(
            endpoint,
            headers=api_headers,
            json=results,
            timeout=config['api']['timeout']
        )
        
        if response.status_code == 200:
            print(f"✅ {job_type.upper()}结果上传成功")
            return True
        else:
            print(f"❌ {job_type.upper()}结果上传失败: HTTP {response.status_code}")
            print(f"   响应: {response.text}")
            return False
            
    except Exception as e:
        print(f"❌ 上传{job_type.upper()}结果失败: {e}")
        return False

def update_job_status(config, job_id, job_type, status='COMPLETED', result_files=None):
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
            'job_type': job_type.upper(),
            'worker_name': config['worker']['name'],
            'progress': 100.0,
            'error_message': None
        }
        
        if result_files:
            data['result_files'] = result_files
        
        print(f"🔄 更新任务状态到后端: {endpoint}")
        
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

def fix_single_task(config, cos_client, cos_bucket, task):
    """修复单个任务"""
    print(f"\n{'='*60}")
    print(f"🔧 开始修复任务 {task['id']}: {task['name']}")
    print(f"   类型: {task['type'].upper()}")
    print(f"   工作目录: {task['work_dir']}")
    
    # 1. 检查工作目录
    if not check_work_directory(task['work_dir']):
        print(f"❌ 任务 {task['id']} 修复失败: 工作目录不存在")
        return False
    
    # 2. 上传结果文件到COS
    print(f"\n📤 步骤1: 上传结果文件到COS...")
    upload_success, uploaded_files = upload_results_selective(
        cos_client, cos_bucket, config, task['id'], task['work_dir'], task['type']
    )
    
    if not upload_success:
        print(f"❌ 任务 {task['id']} 修复失败: 文件上传失败")
        return False
    
    # 3. 解析结果
    print(f"\n📊 步骤2: 解析{task['type'].upper()}结果...")
    if task['type'].lower() == 'md':
        results = parse_md_results(task['work_dir'])
    else:
        results = parse_qc_results(task['work_dir'])
    
    # 4. 上传结果到API
    if results:
        print(f"\n📡 步骤3: 上传{task['type'].upper()}结果到API...")
        api_success = upload_results_to_api(config, task['id'], task['type'], results)
        if not api_success:
            print(f"⚠️ {task['type'].upper()}结果上传到API失败，但继续执行...")
    else:
        print(f"⚠️ 未能解析{task['type'].upper()}结果数据")
    
    # 5. 更新任务状态为COMPLETED
    print(f"\n🔄 步骤4: 更新任务状态...")
    status_success = update_job_status(config, task['id'], task['type'], 'COMPLETED', uploaded_files)
    
    if not status_success:
        print(f"❌ 任务 {task['id']} 修复失败: 状态更新失败")
        return False
    
    print(f"✅ 任务 {task['id']} 修复完成!")
    print(f"   上传文件数: {len(uploaded_files)}")
    print(f"   解析结果数: {len(results)}")
    print(f"   任务状态: COMPLETED")
    
    return True

def main():
    """主函数 - 批量修复所有失败的任务"""
    print("🔧 开始批量修复因COS问题而失败的任务...")
    print(f"📋 共有 {len(FAILED_TASKS)} 个任务需要修复")
    
    # 1. 加载配置
    config = load_config()
    if not config:
        return False
    
    # 2. 初始化COS客户端
    cos_client, cos_bucket = init_cos_client(config)
    if not cos_client:
        return False
    
    # 3. 批量修复任务
    success_count = 0
    failed_count = 0
    
    for task in FAILED_TASKS:
        try:
            if fix_single_task(config, cos_client, cos_bucket, task):
                success_count += 1
            else:
                failed_count += 1
        except Exception as e:
            print(f"❌ 任务 {task['id']} 修复时发生异常: {e}")
            failed_count += 1
    
    # 4. 总结
    print(f"\n{'='*60}")
    print(f"🎉 批量修复完成!")
    print(f"   ✅ 成功修复: {success_count} 个任务")
    print(f"   ❌ 修复失败: {failed_count} 个任务")
    print(f"   📊 总计: {len(FAILED_TASKS)} 个任务")
    
    if success_count > 0:
        print(f"\n🎉 成功修复的任务前端状态应该已更新!")
    
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
        sys.exit(1)
