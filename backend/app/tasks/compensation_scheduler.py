"""
Scheduled tasks for compensation and CPU hours expiration
"""
import logging
from datetime import datetime, timedelta
from sqlalchemy.orm import Session

from app.celery_app import celery_app
from app.core.database import SessionLocal
from app.models import User, CompensationRecord, CompensationStatus, MDJob, JobStatus
from app.services.compensation import CompensationService
from app.services.quota_service import QuotaService

logger = logging.getLogger(__name__)


class DatabaseTask:
    """Base task with database session"""
    def __init__(self):
        self.db = None

    def __call__(self, *args, **kwargs):
        self.db = SessionLocal()
        try:
            return self.run(*args, **kwargs)
        finally:
            if self.db:
                self.db.close()

    def run(self, *args, **kwargs):
        raise NotImplementedError


@celery_app.task(
    bind=True,
    name="app.tasks.compensation_scheduler.check_cpu_hours_expiration",
    default_retry_delay=3600,  # 1 hour
    max_retries=3
)
def check_cpu_hours_expiration(self):
    """
    定期检查核时过期
    应该每天执行一次
    """
    db = SessionLocal()
    try:
        logger.info("Starting CPU hours expiration check")
        
        # 获取所有活跃用户
        users = db.query(User).filter(User.is_active == True).all()
        
        expired_count = 0
        for user in users:
            success, message = CompensationService.check_cpu_hours_expiration(db, user.id)
            if success and "已过期" in message:
                expired_count += 1
                logger.info(f"User {user.id}: {message}")
        
        logger.info(f"CPU hours expiration check completed. {expired_count} users had expired hours")
        return {
            "success": True,
            "checked_users": len(users),
            "expired_users": expired_count
        }
    except Exception as e:
        logger.error(f"Error checking CPU hours expiration: {e}", exc_info=True)
        raise self.retry(exc=e)
    finally:
        db.close()


@celery_app.task(
    bind=True,
    name="app.tasks.compensation_scheduler.process_pending_compensations",
    default_retry_delay=1800,  # 30 minutes
    max_retries=3
)
def process_pending_compensations(self):
    """
    定期处理待处理的补偿记录
    应该每30分钟执行一次
    """
    db = SessionLocal()
    try:
        logger.info("Starting pending compensations processing")
        
        # 获取所有待处理的补偿记录
        pending_records = db.query(CompensationRecord).filter(
            CompensationRecord.status == CompensationStatus.PENDING
        ).all()
        
        logger.info(f"Found {len(pending_records)} pending compensation records")
        
        # 这里可以添加自动批准逻辑
        # 例如：自动批准小于某个金额的补偿
        auto_approved = 0
        for record in pending_records:
            # 自动批准小于等于 10 核时的补偿
            if record.amount <= 10:
                success, message = CompensationService.approve_compensation(
                    db=db,
                    record_id=record.id,
                    approved_by=None,  # 系统自动批准
                    approval_reason="自动批准（金额小于10核时）"
                )
                if success:
                    auto_approved += 1
                    logger.info(f"Auto-approved compensation record {record.id}: {message}")
        
        logger.info(f"Pending compensations processing completed. {auto_approved} records auto-approved")
        return {
            "success": True,
            "pending_records": len(pending_records),
            "auto_approved": auto_approved
        }
    except Exception as e:
        logger.error(f"Error processing pending compensations: {e}", exc_info=True)
        raise self.retry(exc=e)
    finally:
        db.close()


@celery_app.task(
    bind=True,
    name="app.tasks.compensation_scheduler.cleanup_old_records",
    default_retry_delay=86400,  # 24 hours
    max_retries=3
)
def cleanup_old_records(self):
    """
    定期清理旧的补偿记录
    应该每天执行一次
    """
    db = SessionLocal()
    try:
        logger.info("Starting old records cleanup")
        
        # 清理超过 1 年的已完成补偿记录
        cutoff_date = datetime.now() - timedelta(days=365)
        
        old_records = db.query(CompensationRecord).filter(
            CompensationRecord.status == CompensationStatus.COMPLETED,
            CompensationRecord.completed_at < cutoff_date
        ).all()
        
        logger.info(f"Found {len(old_records)} old compensation records to clean up")
        
        # 这里可以选择删除或归档这些记录
        # 目前只是记录日志
        
        return {
            "success": True,
            "old_records_found": len(old_records)
        }
    except Exception as e:
        logger.error(f"Error cleaning up old records: {e}", exc_info=True)
        raise self.retry(exc=e)
    finally:
        db.close()


@celery_app.task(
    bind=True,
    name="app.tasks.compensation_scheduler.process_failed_job_refunds",
    default_retry_delay=1800,  # 30 minutes
    max_retries=3
)
def process_failed_job_refunds(self):
    """
    定期处理失败任务的退款
    应该每30分钟执行一次

    当任务失败时，自动退款用户消费的配额
    """
    db = SessionLocal()
    try:
        logger.info("Starting failed job refunds processing")

        # 获取所有失败但还未处理退款的任务
        # 检查 config 中是否有 quota_consumed 标记
        failed_jobs = db.query(MDJob).filter(
            MDJob.status == JobStatus.FAILED
        ).all()

        refunded_count = 0
        total_refunded = 0.0

        for job in failed_jobs:
            try:
                # 检查是否已经处理过退款
                if job.config and job.config.get("quota_refunded"):
                    continue

                # 检查是否消费过配额
                if not (job.config and job.config.get("quota_consumed")):
                    continue

                # 获取消费的配额
                quota_consumed = job.config.get("quota_consumed", 0.0)
                if quota_consumed <= 0:
                    continue

                # 执行退款
                user = db.query(User).filter(User.id == job.user_id).first()
                if user:
                    success, message = QuotaService.refund_quota(user, quota_consumed, db)
                else:
                    success = False
                if success:
                    # 标记为已退款
                    if job.config is None:
                        job.config = {}
                    job.config["quota_refunded"] = True
                    job.config["refunded_at"] = datetime.now().isoformat()
                    job.config["refund_amount"] = quota_consumed
                    db.commit()

                    refunded_count += 1
                    total_refunded += quota_consumed
                    logger.info(f"Refunded {quota_consumed} hours for failed job {job.id}: {message}")
                else:
                    logger.warning(f"Failed to refund quota for job {job.id}: {message}")

            except Exception as e:
                logger.error(f"Error processing refund for job {job.id}: {e}", exc_info=True)
                continue

        logger.info(f"Failed job refunds processing completed. {refunded_count} jobs refunded, total {total_refunded:.2f} hours")
        return {
            "success": True,
            "refunded_jobs": refunded_count,
            "total_refunded": total_refunded
        }

    except Exception as e:
        logger.error(f"Error processing failed job refunds: {e}", exc_info=True)
        raise self.retry(exc=e)
    finally:
        db.close()

