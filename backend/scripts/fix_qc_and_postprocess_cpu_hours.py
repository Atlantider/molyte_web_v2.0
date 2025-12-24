#!/usr/bin/env python3
"""
修复QC任务和后处理任务的CPU核时数据

这个脚本用于修复那些在系统实现CPU核时追踪之前完成的QC和后处理任务。
这些任务的actual_cpu_hours为0，但根据运行时间和配置可以估算出实际的核时消耗。
"""

import sys
import os
sys.path.append('/opt/molyte_web_v1.0/backend')
os.environ['DATABASE_URL'] = 'postgresql://molyte:molyte2025@127.0.0.1:5432/molyte_db'

from app.database import SessionLocal
from app.models.qc import QCJob, QCJobStatus
from app.models.job import PostprocessJob, JobStatus
from sqlalchemy import func

def fix_qc_cpu_hours():
    """修复QC任务的CPU核时数据"""
    db = SessionLocal()
    
    try:
        # 查询所有核时为0的已完成QC任务
        zero_cpu_jobs = db.query(QCJob).filter(
            QCJob.status == QCJobStatus.COMPLETED,
            QCJob.actual_cpu_hours == 0.0,
            QCJob.started_at.isnot(None),
            QCJob.finished_at.isnot(None)
        ).all()
        
        print(f"Found {len(zero_cpu_jobs)} completed QC jobs with 0 CPU hours\n")
        
        total_fixed = 0
        total_cpu_hours = 0.0
        
        for job in zero_cpu_jobs:
            # 计算运行时间
            duration_seconds = (job.finished_at - job.started_at).total_seconds()
            duration_hours = duration_seconds / 3600.0
            
            # 从配置中获取CPU核心数
            slurm_cpus = job.slurm_cpus or 16  # 默认16核
            
            # 估算核时 = 运行时间 * CPU核心数
            estimated_cpu_hours = duration_hours * slurm_cpus
            
            # 更新任务的actual_cpu_hours
            job.actual_cpu_hours = estimated_cpu_hours
            
            print(f"QC Job {job.id}:")
            print(f"  Duration: {duration_seconds:.0f}s ({duration_hours:.2f}h)")
            print(f"  CPUs: {slurm_cpus}")
            print(f"  Estimated CPU hours: {estimated_cpu_hours:.2f}h")
            
            total_fixed += 1
            total_cpu_hours += estimated_cpu_hours
            print()
        
        # 提交更改
        db.commit()
        
        print(f"✅ Fixed {total_fixed} QC jobs")
        print(f"✅ Total CPU hours added: {total_cpu_hours:.2f}h\n")
        
        return total_cpu_hours
        
    except Exception as e:
        print(f"❌ Error fixing QC jobs: {e}")
        db.rollback()
        return 0.0
    finally:
        db.close()


def fix_postprocess_cpu_hours():
    """修复后处理任务的CPU核时数据"""
    db = SessionLocal()
    
    try:
        # 查询所有核时为0的已完成后处理任务
        zero_cpu_jobs = db.query(PostprocessJob).filter(
            PostprocessJob.status == JobStatus.COMPLETED,
            PostprocessJob.actual_cpu_hours == 0.0,
            PostprocessJob.started_at.isnot(None),
            PostprocessJob.finished_at.isnot(None)
        ).all()
        
        print(f"Found {len(zero_cpu_jobs)} completed postprocess jobs with 0 CPU hours\n")
        
        total_fixed = 0
        total_cpu_hours = 0.0
        
        for job in zero_cpu_jobs:
            # 计算运行时间
            duration_seconds = (job.finished_at - job.started_at).total_seconds()
            duration_hours = duration_seconds / 3600.0
            
            # 后处理任务通常使用较少的CPU核心（假设4核）
            cpu_cores = 4
            
            # 估算核时 = 运行时间 * CPU核心数
            estimated_cpu_hours = duration_hours * cpu_cores
            
            # 更新任务的actual_cpu_hours
            job.actual_cpu_hours = estimated_cpu_hours
            
            print(f"Postprocess Job {job.id}:")
            print(f"  Duration: {duration_seconds:.0f}s ({duration_hours:.2f}h)")
            print(f"  CPUs: {cpu_cores}")
            print(f"  Estimated CPU hours: {estimated_cpu_hours:.2f}h")
            
            total_fixed += 1
            total_cpu_hours += estimated_cpu_hours
            print()
        
        # 提交更改
        db.commit()
        
        print(f"✅ Fixed {total_fixed} postprocess jobs")
        print(f"✅ Total CPU hours added: {total_cpu_hours:.2f}h\n")
        
        return total_cpu_hours
        
    except Exception as e:
        print(f"❌ Error fixing postprocess jobs: {e}")
        db.rollback()
        return 0.0
    finally:
        db.close()


if __name__ == '__main__':
    print("=" * 60)
    print("Fixing QC and Postprocess Job CPU Hours")
    print("=" * 60)
    print()
    
    qc_hours = fix_qc_cpu_hours()
    postprocess_hours = fix_postprocess_cpu_hours()
    
    print("=" * 60)
    print("Summary")
    print("=" * 60)
    print(f"QC jobs CPU hours: {qc_hours:.2f}h")
    print(f"Postprocess jobs CPU hours: {postprocess_hours:.2f}h")
    print(f"Total CPU hours: {qc_hours + postprocess_hours:.2f}h")

