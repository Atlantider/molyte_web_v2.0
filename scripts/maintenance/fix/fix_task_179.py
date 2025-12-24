#!/usr/bin/env python3
"""
修复任务179的脚本
- 清除COS错误信息
- 重新获取Slurm核时信息
- 重新上传结果到COS
- 更新任务状态为COMPLETED
"""

import sys
import os
import subprocess
import yaml
from datetime import datetime
from pathlib import Path

# 添加项目路径
sys.path.insert(0, '/public/home/xiaoji/molyte_web/backend')

from app.database import SessionLocal
from app.models.job import MDJob
from app.core.config import settings

# 导入腾讯云COS SDK
try:
    from qcloud_cos import CosConfig, CosS3Client
except ImportError:
    print("请安装腾讯云 COS SDK: pip install cos-python-sdk-v5")
    sys.exit(1)

def get_slurm_job_info(slurm_job_id):
    """获取Slurm作业的详细信息"""
    try:
        # 使用sacct命令获取作业信息
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
        print(f"获取Slurm信息失败: {e}")
        return None

def calculate_cpu_hours(cpu_time_str):
    """计算CPU小时数"""
    try:
        if not cpu_time_str or cpu_time_str == 'Unknown':
            return 0.0
        
        # 解析格式如 "01:23:45" 或 "1-02:34:56"
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
        print(f"计算CPU小时数失败: {e}")
        return 0.0

def init_cos_client():
    """初始化COS客户端"""
    try:
        # 加载配置文件
        config_path = Path('/public/home/xiaoji/molyte_web/deployment/polling_worker_config_tencent.yaml')
        if not config_path.exists():
            print(f"配置文件不存在: {config_path}")
            return None, None

        with open(config_path, 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)

        if 'cos' not in config:
            print("配置文件中未找到COS配置")
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

        print(f"✅ COS客户端初始化成功 (Bucket: {cos_bucket})")
        return cos_client, cos_bucket

    except Exception as e:
        print(f"初始化COS客户端失败: {e}")
        return None, None

def upload_result_file(job_id, work_dir):
    """重新上传结果文件到COS"""
    try:
        # 初始化COS客户端
        cos_client, cos_bucket = init_cos_client()
        if not cos_client:
            return False

        # 查找结果文件
        data_file = None
        work_path = Path(work_dir)
        for file_path in work_path.iterdir():
            if file_path.suffix == '.data':
                data_file = file_path
                break

        if not data_file or not data_file.exists():
            print(f"未找到结果文件在目录: {work_dir}")
            return False

        # 上传到COS
        cos_key = f"MD_results/{job_id}/{data_file.name}"
        print(f"上传文件: {data_file} -> {cos_key}")

        file_size = data_file.stat().st_size
        print(f"文件大小: {file_size/1024/1024:.1f}MB")

        # 执行上传
        with open(data_file, 'rb') as f:
            cos_client.put_object(
                Bucket=cos_bucket,
                Body=f,
                Key=cos_key
            )

        print(f"✅ 文件上传成功: {cos_key}")
        return True

    except Exception as e:
        print(f"上传文件时出错: {e}")
        return False

def fix_task_179():
    """修复任务179"""
    db = SessionLocal()
    
    try:
        # 1. 查询任务179
        job = db.query(MDJob).filter(MDJob.id == 179).first()
        
        if not job:
            print("❌ 未找到任务179")
            return False
        
        print(f"📋 当前任务状态:")
        print(f"   ID: {job.id}")
        print(f"   状态: {job.status}")
        print(f"   Slurm ID: {job.slurm_job_id}")
        print(f"   工作目录: {job.work_dir}")
        print(f"   错误信息: {job.error_message}")
        print(f"   CPU小时: {job.cpu_hours}")
        
        # 2. 获取Slurm作业信息
        print(f"\n🔍 获取Slurm作业信息...")
        slurm_info = get_slurm_job_info(job.slurm_job_id)
        
        if slurm_info:
            print(f"   Slurm状态: {slurm_info['state']}")
            print(f"   CPU时间: {slurm_info['cpu_time']}")
            print(f"   运行时间: {slurm_info['elapsed_time']}")
            print(f"   开始时间: {slurm_info['start_time']}")
            print(f"   结束时间: {slurm_info['end_time']}")
            
            # 计算CPU小时数
            cpu_hours = calculate_cpu_hours(slurm_info['cpu_time'])
            print(f"   计算的CPU小时: {cpu_hours}")
        else:
            print("   ⚠️ 无法获取Slurm信息")
            cpu_hours = 0.0
        
        # 3. 重新上传结果文件
        print(f"\n📤 重新上传结果文件...")
        upload_success = upload_result_file(job.id, job.work_dir)
        
        if not upload_success:
            print("❌ 文件上传失败，请检查COS账户状态")
            return False
        
        # 4. 更新数据库
        print(f"\n💾 更新数据库...")
        job.status = 'COMPLETED'
        job.error_message = None
        job.cpu_hours = cpu_hours
        job.finished_at = datetime.now()
        
        if slurm_info and slurm_info['end_time'] != 'Unknown':
            try:
                # 解析Slurm的结束时间
                job.finished_at = datetime.strptime(slurm_info['end_time'], '%Y-%m-%dT%H:%M:%S')
            except:
                pass
        
        db.commit()
        
        print(f"✅ 任务179修复完成!")
        print(f"   新状态: {job.status}")
        print(f"   CPU小时: {job.cpu_hours}")
        print(f"   完成时间: {job.finished_at}")
        
        return True
        
    except Exception as e:
        print(f"❌ 修复过程中出错: {e}")
        db.rollback()
        return False
    finally:
        db.close()

if __name__ == "__main__":
    print("🔧 开始修复任务179...")
    try:
        success = fix_task_179()

        if success:
            print("\n🎉 任务179修复成功!")
        else:
            print("\n💥 任务179修复失败!")
            sys.exit(1)
    except Exception as e:
        print(f"\n💥 脚本执行出错: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
