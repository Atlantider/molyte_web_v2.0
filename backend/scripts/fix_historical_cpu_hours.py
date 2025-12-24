#!/usr/bin/env python3
"""
修复历史任务的CPU核时数据

这个脚本用于修复那些在系统实现CPU核时追踪之前完成的任务。
这些任务的actual_cpu_hours为0，但根据运行时间和配置可以估算出实际的核时消耗。
"""

import sys
import os
sys.path.append('/opt/molyte_web_v1.0/backend')
os.environ['DATABASE_URL'] = 'postgresql://molyte:molyte2025@127.0.0.1:5432/molyte_db'

from app.database import SessionLocal
from app.models.job import MDJob, JobStatus
from app.models.billing import QuotaTransaction
from app.models.user import User

def fix_historical_cpu_hours():
    """修复历史任务的CPU核时数据"""
    db = SessionLocal()
    
    try:
        # 查询所有核时为0的已完成任务
        zero_cpu_jobs = db.query(MDJob).filter(
            MDJob.status == JobStatus.COMPLETED,
            MDJob.actual_cpu_hours == 0.0,
            MDJob.started_at.isnot(None),
            MDJob.finished_at.isnot(None)
        ).all()
        
        print(f"Found {len(zero_cpu_jobs)} completed jobs with 0 CPU hours\n")
        
        if not zero_cpu_jobs:
            print("No jobs to fix")
            return
        
        total_fixed = 0
        total_cpu_hours = 0.0
        
        for job in zero_cpu_jobs:
            # 计算运行时间
            duration_seconds = (job.finished_at - job.started_at).total_seconds()
            duration_hours = duration_seconds / 3600.0
            
            # 从配置中获取CPU核心数
            config = job.config or {}
            slurm_ntasks = config.get('slurm_ntasks', 1)
            slurm_cpus_per_task = config.get('slurm_cpus_per_task', 1)
            total_cpus = slurm_ntasks * slurm_cpus_per_task
            
            # 估算核时 = 运行时间 * CPU核心数
            estimated_cpu_hours = duration_hours * total_cpus
            
            # 更新任务的actual_cpu_hours
            job.actual_cpu_hours = estimated_cpu_hours
            
            print(f"Job {job.id}:")
            print(f"  Duration: {duration_seconds:.0f}s ({duration_hours:.2f}h)")
            print(f"  CPUs: {slurm_ntasks} tasks × {slurm_cpus_per_task} cpus/task = {total_cpus} total")
            print(f"  Estimated CPU hours: {estimated_cpu_hours:.2f}h")
            
            total_fixed += 1
            total_cpu_hours += estimated_cpu_hours
            print()
        
        # 提交更改
        db.commit()
        
        print(f"✅ Fixed {total_fixed} jobs")
        print(f"✅ Total CPU hours added: {total_cpu_hours:.2f}h")
        
        # 显示受影响的用户
        print("\nAffected users:")
        affected_users = db.query(User).filter(
            User.id.in_([job.user_id for job in zero_cpu_jobs])
        ).all()
        
        for user in affected_users:
            user_jobs = [j for j in zero_cpu_jobs if j.user_id == user.id]
            user_cpu_hours = sum(j.actual_cpu_hours for j in user_jobs)
            print(f"  {user.username}: {len(user_jobs)} jobs, {user_cpu_hours:.2f}h")
        
    except Exception as e:
        print(f"❌ Error: {e}")
        db.rollback()
    finally:
        db.close()

if __name__ == '__main__':
    fix_historical_cpu_hours()

