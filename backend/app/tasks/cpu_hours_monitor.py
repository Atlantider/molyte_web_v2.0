"""
CPU核时监控任务

定期更新运行中任务的实时核时数据，使用户能看到任务的实时消耗
"""

import logging
import subprocess
from typing import Dict, Any
from datetime import datetime

from sqlalchemy.orm import Session
from app.celery_app import celery_app
from app.database import SessionLocal
from app.models.job import MDJob, JobStatus
from app.models.qc import QCJob, QCJobStatus

logger = logging.getLogger(__name__)


def get_running_job_cpu_hours(slurm_job_id: str) -> float:
    """
    获取运行中任务的实时CPU核时数
    
    使用sstat获取运行中任务的实时核时数据
    
    Args:
        slurm_job_id: Slurm任务ID
        
    Returns:
        CPU核时数（小时）
    """
    try:
        # 使用sstat获取运行中任务的实时核时
        result = subprocess.run(
            ['sstat', '-j', slurm_job_id, '-o', 'CPUTimeRAW', '-n', '-X'],
            capture_output=True,
            text=True,
            timeout=10
        )
        
        if result.returncode == 0 and result.stdout.strip():
            try:
                cpu_time_seconds = int(result.stdout.strip().split()[0])
                if cpu_time_seconds > 0:
                    cpu_hours = cpu_time_seconds / 3600.0
                    logger.debug(f"Job {slurm_job_id}: CPUTimeRAW={cpu_time_seconds}s, CPU hours={cpu_hours:.2f}h (from sstat)")
                    return cpu_hours
            except (ValueError, IndexError) as e:
                logger.warning(f"Failed to parse CPUTimeRAW for job {slurm_job_id}: {e}")
        
        return 0.0
        
    except subprocess.TimeoutExpired:
        logger.warning(f"sstat timeout for job {slurm_job_id}")
        return 0.0
    except Exception as e:
        logger.error(f"Failed to get CPU hours for job {slurm_job_id}: {e}")
        return 0.0


@celery_app.task(
    name="app.tasks.cpu_hours_monitor.update_running_jobs_cpu_hours",
    bind=True,
)
def update_running_jobs_cpu_hours(self) -> Dict[str, Any]:
    """
    定期更新所有运行中任务的实时CPU核时数据
    
    这个任务应该被定期调度（如每分钟一次）
    """
    db = SessionLocal()
    
    try:
        logger.info("[CPU Hours Monitor] Starting to update running jobs CPU hours...")
        
        # 查询所有运行中的MD任务
        running_md_jobs = db.query(MDJob).filter(
            MDJob.status == JobStatus.RUNNING
        ).all()
        
        # 查询所有运行中的QC任务
        running_qc_jobs = db.query(QCJob).filter(
            QCJob.status == QCJobStatus.RUNNING
        ).all()
        
        md_updated = 0
        qc_updated = 0
        
        # 更新MD任务的核时
        for job in running_md_jobs:
            if not job.slurm_job_id:
                continue
            
            try:
                cpu_hours = get_running_job_cpu_hours(job.slurm_job_id)
                if cpu_hours > 0 and cpu_hours != job.actual_cpu_hours:
                    job.actual_cpu_hours = cpu_hours
                    md_updated += 1
                    logger.debug(f"Updated MD Job {job.id} CPU hours to {cpu_hours:.2f}h")
            except Exception as e:
                logger.error(f"Failed to update CPU hours for MD Job {job.id}: {e}")
        
        # 更新QC任务的核时
        for job in running_qc_jobs:
            if not job.slurm_job_id:
                continue
            
            try:
                cpu_hours = get_running_job_cpu_hours(job.slurm_job_id)
                if cpu_hours > 0 and cpu_hours != job.actual_cpu_hours:
                    job.actual_cpu_hours = cpu_hours
                    qc_updated += 1
                    logger.debug(f"Updated QC Job {job.id} CPU hours to {cpu_hours:.2f}h")
            except Exception as e:
                logger.error(f"Failed to update CPU hours for QC Job {job.id}: {e}")
        
        # 提交更改
        if md_updated > 0 or qc_updated > 0:
            db.commit()
            logger.info(f"[CPU Hours Monitor] Updated {md_updated} MD jobs and {qc_updated} QC jobs")
        
        return {
            "status": "ok",
            "md_jobs_updated": md_updated,
            "qc_jobs_updated": qc_updated,
            "total_updated": md_updated + qc_updated
        }
        
    except Exception as e:
        logger.error(f"[CPU Hours Monitor] Error updating CPU hours: {e}", exc_info=True)
        return {
            "status": "error",
            "error": str(e)
        }
    finally:
        db.close()

