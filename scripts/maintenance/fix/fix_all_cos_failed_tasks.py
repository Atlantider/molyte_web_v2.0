#!/usr/bin/env python3
"""
修复所有因COS账户欠费而失败的任务
包括MD任务和QC任务
"""

import os
import subprocess
import yaml
from datetime import datetime
from pathlib import Path

# 导入腾讯云COS SDK
try:
    from qcloud_cos import CosConfig, CosS3Client
except ImportError:
    print("请安装腾讯云 COS SDK: pip install cos-python-sdk-v5")
    exit(1)

# 失败的任务信息
FAILED_TASKS = [
    {
        'id': 178,
        'type': 'md',
        'slurm_id': 22091,
        'work_dir': '/public/home/xiaoji/molyte_web/data/md_work/EL-20251205-0001-Li-PF6-EC-DEC-PC-MD1-298K',
        'failed_time': '2025-12-15 01:06:48'
    },
    {
        'id': 179,  # EL-20251205-0002 任务 (推测的ID)
        'type': 'md',
        'slurm_id': 22092,  # 推测的Slurm ID
        'work_dir': '/public/home/xiaoji/molyte_web/data/md_work/EL-20251205-0002-Li-PF6-EC-DEC-PC-MD1-298K',
        'failed_time': '2025-12-15 01:06:48'
    },
    {
        'id': 1574,
        'type': 'qc',
        'slurm_id': 7153,  # 从日志中找到的Slurm ID
        'work_dir': '/public/home/xiaoji/molyte_web/data/qc_work/QC-1574-MD133_cluster_minus_PF6_2_1466',
        'failed_time': '2025-12-08 17:49:09'
    },
    {
        'id': 1584,
        'type': 'qc',
        'slurm_id': 7154,  # 推测的Slurm ID
        'work_dir': '/public/home/xiaoji/molyte_web/data/qc_work/QC-1584-MD133_cluster_minus_EC_1_DEC_1_1466',
        'failed_time': '2025-12-08 17:56:25'
    },
    {
        'id': 1836,
        'type': 'qc',
        'slurm_id': 22085,  # 最后一次成功的Slurm ID
        'work_dir': '/public/home/xiaoji/molyte_web/data/qc_work/QC-1836-MD176_dimer_PF6_1507',
        'failed_time': '2025-12-14 08:29:38'
    },
    {
        'id': 1838,
        'type': 'qc',
        'slurm_id': 21837,
        'work_dir': '/public/home/xiaoji/molyte_web/data/qc_work/QC-1838-MD176_intermediate_DEC_1_PF6_1_1507',
        'failed_time': '2025-12-14 08:29:38'
    },
    {
        'id': 1839,
        'type': 'qc',
        'slurm_id': 21838,
        'work_dir': '/public/home/xiaoji/molyte_web/data/qc_work/QC-1839-MD176_intermediate_DEC_2_1507',
        'failed_time': '2025-12-15 00:08:23'
    },
    {
        'id': 1840,
        'type': 'qc',
        'slurm_id': 21839,
        'work_dir': '/public/home/xiaoji/molyte_web/data/qc_work/QC-1840-MD176_intermediate_DEC_2_PF6_1_1507',
        'failed_time': '2025-12-14 17:01:01'
    }
]

def init_cos_client():
    """初始化COS客户端"""
    try:
        config_path = Path('/public/home/xiaoji/molyte_web/deployment/polling_worker_config_tencent.yaml')
        if not config_path.exists():
            print(f"❌ 配置文件不存在: {config_path}")
            return None, None
        
        with open(config_path, 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)
        
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
        print(f"   Bucket: {cos_bucket}")
        print(f"   Region: {cos_config['region']}")
        return cos_client, cos_bucket
        
    except Exception as e:
        print(f"❌ 初始化COS客户端失败: {e}")
        return None, None

def get_slurm_job_info(slurm_job_id):
    """获取Slurm作业的详细信息"""
    try:
        cmd = f"sacct -j {slurm_job_id} --format=JobID,State,CPUTime,Elapsed,Start,End --parsable2 --noheader"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if result.returncode == 0:
            lines = result.stdout.strip().split('\n')
            for line in lines:
                if line and not line.endswith('.batch') and not line.endswith('.extern'):
                    parts = line.split('|')
                    if len(parts) >= 6:
                        job_id, state, cpu_time, elapsed_time, start_time, end_time = parts[:6]
                        return {
                            'job_id': job_id,
                            'state': state,
                            'cpu_time': cpu_time,
                            'elapsed_time': elapsed_time,
                            'start_time': start_time,
                            'end_time': end_time
                        }
        return None
    except Exception as e:
        print(f"❌ 获取Slurm信息失败: {e}")
        return None

def calculate_cpu_hours(cpu_time_str):
    """计算CPU小时数"""
    try:
        if not cpu_time_str or cpu_time_str == 'Unknown':
            return 0.0
        
        if '-' in cpu_time_str:
            days, time_part = cpu_time_str.split('-')
            days = int(days)
        else:
            days = 0
            time_part = cpu_time_str
        
        time_parts = time_part.split(':')
        if len(time_parts) == 3:
            hours = int(time_parts[0])
            minutes = int(time_parts[1])
            seconds = int(time_parts[2])
            
            total_hours = days * 24 + hours + minutes / 60.0 + seconds / 3600.0
            return round(total_hours, 2)
        
        return 0.0
    except Exception as e:
        print(f"❌ 计算CPU小时数失败: {e}")
        return 0.0

def upload_task_results(task_info, cos_client, cos_bucket):
    """上传任务结果文件到COS"""
    try:
        job_id = task_info['id']
        job_type = task_info['type']
        work_dir = task_info['work_dir']
        
        work_path = Path(work_dir)
        if not work_path.exists():
            print(f"❌ 工作目录不存在: {work_dir}")
            return False
        
        print(f"📁 检查工作目录: {work_dir}")
        
        # 根据任务类型确定上传前缀和文件模式
        if job_type == 'md':
            result_prefix = 'results'  # MD任务使用 results 前缀
            file_patterns = ['*.data', '*.log', '*.out', '*.err', '*.xyz', '*.pdb', '*.lammpstrj', 'out_*.dat']
        else:  # qc
            result_prefix = 'results'  # QC任务也使用 results 前缀
            file_patterns = ['*.log', '*.out', '*.err', '*.xyz', '*.pdb', '*.gjf', '*.chk', '*.fchk', '*.png']
        
        uploaded_files = []
        
        # 上传文件
        for pattern in file_patterns:
            files = list(work_path.glob(pattern))
            for file_path in files:
                if file_path.stat().st_size > 0:  # 只上传非空文件
                    cos_key = f"{result_prefix}/{job_id}/{file_path.name}"
                    file_size = file_path.stat().st_size
                    
                    print(f"📤 上传文件: {file_path.name} ({file_size/1024:.1f}KB)")
                    
                    with open(file_path, 'rb') as f:
                        cos_client.put_object(
                            Bucket=cos_bucket,
                            Body=f,
                            Key=cos_key
                        )
                    
                    uploaded_files.append(cos_key)
                    print(f"   ✅ 上传成功")
        
        print(f"✅ 任务{job_id}共上传 {len(uploaded_files)} 个文件到COS")
        return len(uploaded_files) > 0
            
    except Exception as e:
        print(f"❌ 上传任务{task_info['id']}文件时出错: {e}")
        return False

def fix_task(task_info):
    """修复单个任务"""
    print(f"\n🔧 开始修复任务{task_info['id']} ({task_info['type'].upper()})...")
    
    # 获取Slurm信息
    print(f"🔍 获取Slurm作业信息...")
    slurm_info = get_slurm_job_info(task_info['slurm_id'])
    
    if slurm_info:
        print(f"   ✅ Slurm状态: {slurm_info['state']}")
        print(f"   CPU时间: {slurm_info['cpu_time']}")
        print(f"   运行时间: {slurm_info['elapsed_time']}")
        print(f"   开始时间: {slurm_info['start_time']}")
        print(f"   结束时间: {slurm_info['end_time']}")
        
        cpu_hours = calculate_cpu_hours(slurm_info['cpu_time'])
        print(f"   计算的CPU小时: {cpu_hours}")
    else:
        print("   ⚠️ 无法获取Slurm信息")
        cpu_hours = 0.0
        slurm_info = {}
    
    # 初始化COS客户端
    cos_client, cos_bucket = init_cos_client()
    if not cos_client:
        return False
    
    # 重新上传结果文件
    print(f"📤 重新上传结果文件...")
    upload_success = upload_task_results(task_info, cos_client, cos_bucket)
    
    if upload_success:
        print(f"✅ 任务{task_info['id']}修复成功!")
        
        # 生成数据库更新SQL
        table_name = 'md_jobs' if task_info['type'] == 'md' else 'qc_jobs'
        
        print(f"\n💡 数据库更新SQL:")
        print(f"UPDATE {table_name} SET")
        print(f"  status = 'COMPLETED',")
        print(f"  error_message = NULL,")
        print(f"  cpu_hours = {cpu_hours},")
        if slurm_info.get('end_time') and slurm_info['end_time'] != 'Unknown':
            print(f"  finished_at = '{slurm_info['end_time']}'")
        else:
            print(f"  finished_at = NOW()")
        print(f"WHERE id = {task_info['id']};")
        
        return True
    else:
        print(f"❌ 任务{task_info['id']}修复失败!")
        return False

def main():
    """主函数"""
    print("🔧 开始修复所有因COS账户欠费而失败的任务...")
    print(f"📋 共发现 {len(FAILED_TASKS)} 个失败任务")
    
    success_count = 0
    
    for i, task_info in enumerate(FAILED_TASKS, 1):
        print(f"\n{'='*60}")
        print(f"修复进度: {i}/{len(FAILED_TASKS)}")
        
        if fix_task(task_info):
            success_count += 1
    
    print(f"\n{'='*60}")
    print(f"🎉 修复完成!")
    print(f"   成功修复: {success_count}/{len(FAILED_TASKS)} 个任务")
    
    if success_count < len(FAILED_TASKS):
        print(f"   失败任务: {len(FAILED_TASKS) - success_count} 个")
    
    print(f"\n📝 接下来需要在腾讯云服务器上执行上述SQL命令来更新数据库")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"\n💥 脚本执行出错: {e}")
        import traceback
        traceback.print_exc()
        exit(1)
