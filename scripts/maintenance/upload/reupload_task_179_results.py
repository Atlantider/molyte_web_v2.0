#!/usr/bin/env python3
"""
重新上传任务179的结果文件到COS
由于数据库在腾讯云上，校园网无法直接访问，所以这个脚本只负责重新上传文件
数据库的更新需要通过API或者在腾讯云服务器上执行
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

def init_cos_client():
    """初始化COS客户端"""
    try:
        # 加载配置文件
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
        print(f"❌ 获取Slurm信息失败: {e}")
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
        print(f"❌ 计算CPU小时数失败: {e}")
        return 0.0

def upload_result_files(job_id, work_dir):
    """重新上传结果文件到COS"""
    try:
        # 初始化COS客户端
        cos_client, cos_bucket = init_cos_client()
        if not cos_client:
            return False
        
        work_path = Path(work_dir)
        if not work_path.exists():
            print(f"❌ 工作目录不存在: {work_dir}")
            return False
        
        print(f"\n📁 检查工作目录: {work_dir}")
        
        # 查找结果文件
        uploaded_files = []
        
        # 查找.data文件
        data_files = list(work_path.glob('*.data'))
        if data_files:
            for data_file in data_files:
                cos_key = f"MD_results/{job_id}/{data_file.name}"
                file_size = data_file.stat().st_size
                
                print(f"📤 上传文件: {data_file.name} ({file_size/1024/1024:.1f}MB)")
                print(f"   目标: {cos_key}")
                
                # 执行上传
                with open(data_file, 'rb') as f:
                    cos_client.put_object(
                        Bucket=cos_bucket,
                        Body=f,
                        Key=cos_key
                    )
                
                uploaded_files.append(cos_key)
                print(f"   ✅ 上传成功")
        else:
            print("⚠️  未找到.data文件")
        
        # 查找其他可能的结果文件
        other_patterns = ['*.log', '*.out', '*.err', '*.xyz', '*.pdb']
        for pattern in other_patterns:
            files = list(work_path.glob(pattern))
            for file_path in files:
                if file_path.stat().st_size > 0:  # 只上传非空文件
                    cos_key = f"MD_results/{job_id}/{file_path.name}"
                    file_size = file_path.stat().st_size
                    
                    print(f"📤 上传额外文件: {file_path.name} ({file_size/1024:.1f}KB)")
                    
                    with open(file_path, 'rb') as f:
                        cos_client.put_object(
                            Bucket=cos_bucket,
                            Body=f,
                            Key=cos_key
                        )
                    
                    uploaded_files.append(cos_key)
                    print(f"   ✅ 上传成功")
        
        print(f"\n✅ 共上传 {len(uploaded_files)} 个文件到COS")
        for file_key in uploaded_files:
            print(f"   - {file_key}")
        
        return len(uploaded_files) > 0
            
    except Exception as e:
        print(f"❌ 上传文件时出错: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """主函数"""
    print("🔧 开始重新上传任务179的结果文件...")
    
    # 任务信息
    job_id = 179
    slurm_job_id = 22092
    work_dir = "/public/home/xiaoji/molyte_web/data/md_work/EL-20251205-0002-Li-PF6-EC-DEC-PC-MD1-298K"
    
    print(f"\n📋 任务信息:")
    print(f"   任务ID: {job_id}")
    print(f"   Slurm ID: {slurm_job_id}")
    print(f"   工作目录: {work_dir}")
    
    # 1. 获取Slurm作业信息
    print(f"\n🔍 获取Slurm作业信息...")
    slurm_info = get_slurm_job_info(slurm_job_id)
    
    if slurm_info:
        print(f"   ✅ Slurm状态: {slurm_info['state']}")
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
    
    # 2. 重新上传结果文件
    print(f"\n📤 重新上传结果文件...")
    upload_success = upload_result_files(job_id, work_dir)
    
    if upload_success:
        print(f"\n🎉 任务179结果文件重新上传成功!")
        print(f"\n📝 接下来需要在腾讯云服务器上更新数据库:")
        print(f"   1. 清除 error_message 字段")
        print(f"   2. 设置 cpu_hours = {cpu_hours}")
        print(f"   3. 设置 status = 'COMPLETED'")
        if slurm_info and slurm_info['end_time'] != 'Unknown':
            print(f"   4. 设置 finished_at = '{slurm_info['end_time']}'")
        else:
            print(f"   4. 设置 finished_at = 当前时间")
        
        print(f"\n💡 建议的SQL命令:")
        print(f"UPDATE md_jobs SET")
        print(f"  status = 'COMPLETED',")
        print(f"  error_message = NULL,")
        print(f"  cpu_hours = {cpu_hours},")
        if slurm_info and slurm_info['end_time'] != 'Unknown':
            print(f"  finished_at = '{slurm_info['end_time']}'")
        else:
            print(f"  finished_at = NOW()")
        print(f"WHERE id = {job_id};")
        
        return True
    else:
        print(f"\n💥 任务179结果文件上传失败!")
        return False

if __name__ == "__main__":
    try:
        success = main()
        exit(0 if success else 1)
    except Exception as e:
        print(f"\n💥 脚本执行出错: {e}")
        import traceback
        traceback.print_exc()
        exit(1)
