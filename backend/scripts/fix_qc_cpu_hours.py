#!/usr/bin/env python3
"""
修复 QC 任务的核时记录

对于已完成但 actual_cpu_hours = 0 的 QC 任务，
尝试从 Slurm 获取实际的 CPU 核时数据
"""

import sys
import os
import subprocess
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from app.database import SessionLocal
from app.models.qc import QCJob, QCJobStatus

def get_slurm_cpu_hours(slurm_job_id: str) -> float:
    """从 Slurm 获取任务的 CPU 核时"""
    try:
        # sacct -j JOB_ID -o CPUTimeRAW,NCPUS -n -X
        result = subprocess.run(
            ['sacct', '-j', slurm_job_id, '-o', 'CPUTimeRAW,NCPUS,Elapsed', '-n', '-X'],
            capture_output=True,
            text=True,
            timeout=10
        )

        if result.returncode == 0 and result.stdout.strip():
            parts = result.stdout.strip().split()
            if len(parts) >= 2:
                cpu_time_seconds = int(parts[0])
                cpu_hours = cpu_time_seconds / 3600.0
                return cpu_hours

        return 0.0

    except Exception as e:
        print(f"    ❌ 获取核时失败: {e}")
        return 0.0

def main():
    db = SessionLocal()
    
    try:
        print("=== 修复 QC 任务核时记录 ===\n")
        
        # 查找所有已完成但核时为0的QC任务
        qc_jobs = db.query(QCJob).filter(
            QCJob.status == QCJobStatus.COMPLETED,
            QCJob.actual_cpu_hours == 0.0,
            QCJob.slurm_job_id.isnot(None)
        ).all()
        
        print(f"找到 {len(qc_jobs)} 个需要修复的 QC 任务")
        
        if not qc_jobs:
            print("没有需要修复的任务")
            return
        
        fixed_count = 0
        total_hours = 0.0
        
        for qc in qc_jobs:
            print(f"\n处理 QC 任务 {qc.id} (Slurm: {qc.slurm_job_id}):")
            
            # 从 Slurm 获取核时
            cpu_hours = get_slurm_cpu_hours(qc.slurm_job_id)
            
            if cpu_hours > 0:
                qc.actual_cpu_hours = cpu_hours
                db.commit()
                print(f"    ✅ 更新核时: {cpu_hours:.2f}h")
                fixed_count += 1
                total_hours += cpu_hours
            else:
                print(f"    ⚠️  无法获取核时数据")
        
        print(f"\n=== 修复完成 ===")
        print(f"成功修复任务数: {fixed_count}")
        print(f"总核时: {total_hours:.2f}h")
        
        # 验证结果
        print(f"\n=== 验证结果 ===")
        remaining_jobs = db.query(QCJob).filter(
            QCJob.status == QCJobStatus.COMPLETED,
            QCJob.actual_cpu_hours == 0.0,
            QCJob.slurm_job_id.isnot(None)
        ).count()
        print(f"剩余未修复任务: {remaining_jobs}")
        
        # 重新检查 cluster analysis 任务
        print(f"\n=== 重新检查 Cluster Analysis 任务 ===")
        from app.models.job import AdvancedClusterJob
        
        for job_id in [29, 27]:
            job = db.query(AdvancedClusterJob).filter(AdvancedClusterJob.id == job_id).first()
            if job:
                qc_jobs = db.query(QCJob).filter(QCJob.cluster_analysis_job_id == job_id).all()
                total_qc_hours = sum(qc.actual_cpu_hours for qc in qc_jobs)
                print(f"任务 {job_id}: {len(qc_jobs)} 个QC任务, 总核时 {total_qc_hours:.2f}h")
        
    except Exception as e:
        print(f"❌ 错误: {e}")
        import traceback
        traceback.print_exc()
    finally:
        db.close()

if __name__ == "__main__":
    main()
