"""
Compensation service for automatic CPU hours refund and expiration
"""
import logging
import json
from typing import Tuple, List, Optional
from datetime import datetime, timedelta
from sqlalchemy.orm import Session

from app.models import (
    CompensationRule, CompensationRecord, CPUHoursExpiration,
    CompensationRuleType, CompensationStatus, User, MDJob, QuotaTransaction
)
from app.services.billing import BillingService

logger = logging.getLogger(__name__)


class CompensationService:
    """Service for managing compensation rules and records"""

    @staticmethod
    def create_compensation_rule(
        db: Session,
        name: str,
        rule_type: CompensationRuleType,
        config: dict,
        description: Optional[str] = None
    ) -> Tuple[bool, str, Optional[CompensationRule]]:
        """创建补偿规则"""
        try:
            rule = CompensationRule(
                name=name,
                description=description,
                rule_type=rule_type,
                config=config,
                is_active=True
            )
            db.add(rule)
            db.commit()
            logger.info(f"Created compensation rule: {name}")
            return True, f"补偿规则 '{name}' 创建成功", rule
        except Exception as e:
            db.rollback()
            logger.error(f"Failed to create compensation rule: {e}")
            return False, f"创建补偿规则失败: {str(e)}", None

    @staticmethod
    def process_job_failure_compensation(
        db: Session,
        job: MDJob
    ) -> Tuple[bool, str]:
        """处理任务失败补偿"""
        # 获取任务失败补偿规则
        rule = db.query(CompensationRule).filter(
            CompensationRule.rule_type == CompensationRuleType.JOB_FAILURE,
            CompensationRule.is_active == True
        ).first()

        if not rule:
            return True, "没有激活的任务失败补偿规则"

        config = rule.config
        failure_types = config.get("failure_types", ["FAILED", "TIMEOUT"])
        refund_percentage = config.get("refund_percentage", 100)
        min_cpu_hours = config.get("min_cpu_hours", 0.1)

        # 检查是否应该补偿
        if job.status.value not in failure_types:
            return True, f"任务状态 {job.status.value} 不在补偿范围内"

        # 计算应该退款的核时
        total_cpu_hours = BillingService.calculate_job_cpu_hours(job)
        refund_amount = total_cpu_hours * (refund_percentage / 100)

        if refund_amount < min_cpu_hours:
            return True, f"退款核时 {refund_amount:.2f} 小于最小值 {min_cpu_hours}"

        try:
            user = db.query(User).filter(User.id == job.user_id).first()
            if not user:
                return False, "用户不存在"

            # 创建补偿记录
            record = CompensationRecord(
                rule_id=rule.id,
                user_id=job.user_id,
                status=CompensationStatus.APPROVED,
                reason=f"任务 #{job.id} 失败自动退款",
                amount=refund_amount,
                reference_id=job.id,
                reference_type="job",
                approved_by=None,  # 自动批准
                approved_at=datetime.now(),
                completed_at=datetime.now()
            )
            db.add(record)

            # 退款给用户
            balance_before = user.balance_cpu_hours
            user.balance_cpu_hours += refund_amount

            # 记录交易
            trans = QuotaTransaction(
                user_id=user.id,
                type="refund",
                amount=refund_amount,
                balance_before=balance_before,
                balance_after=user.balance_cpu_hours,
                reference_id=job.id,
                reference_type="job",
                description=f"任务 #{job.id} 失败自动退款 {refund_amount:.2f} 核时"
            )
            db.add(trans)
            db.commit()

            logger.info(f"Job {job.id} failure compensation: refunded {refund_amount:.2f} hours to user {job.user_id}")
            return True, f"已自动退款 {refund_amount:.2f} 核时"
        except Exception as e:
            db.rollback()
            logger.error(f"Failed to process job failure compensation: {e}")
            return False, f"处理补偿失败: {str(e)}"

    @staticmethod
    def check_cpu_hours_expiration(
        db: Session,
        user_id: int
    ) -> Tuple[bool, str]:
        """检查并处理核时过期"""
        rule = db.query(CompensationRule).filter(
            CompensationRule.rule_type == CompensationRuleType.EXPIRATION,
            CompensationRule.is_active == True
        ).first()

        if not rule:
            return True, "没有激活的核时过期规则"

        config = rule.config
        expiration_days = config.get("expiration_days", 365)

        try:
            user = db.query(User).filter(User.id == user_id).first()
            if not user:
                return False, "用户不存在"

            # 计算过期时间
            expiration_date = datetime.now() - timedelta(days=expiration_days)

            # 查询用户的最后一次消费时间
            last_transaction = db.query(QuotaTransaction).filter(
                QuotaTransaction.user_id == user_id,
                QuotaTransaction.type == "consume"
            ).order_by(QuotaTransaction.created_at.desc()).first()

            if not last_transaction or last_transaction.created_at < expiration_date:
                # 核时已过期，记录过期
                if user.balance_cpu_hours > 0:
                    expired_amount = user.balance_cpu_hours
                    
                    expiration_record = CPUHoursExpiration(
                        user_id=user_id,
                        amount=expired_amount,
                        reason=f"核时未使用超过 {expiration_days} 天自动过期",
                        expired_at=datetime.now()
                    )
                    db.add(expiration_record)

                    # 清空用户余额
                    user.balance_cpu_hours = 0
                    db.commit()

                    logger.info(f"User {user_id} CPU hours expired: {expired_amount:.2f} hours")
                    return True, f"已过期 {expired_amount:.2f} 核时"
                else:
                    return True, "用户没有可过期的核时"
            else:
                return True, "用户核时未过期"
        except Exception as e:
            db.rollback()
            logger.error(f"Failed to check CPU hours expiration: {e}")
            return False, f"检查过期失败: {str(e)}"

    @staticmethod
    def approve_compensation(
        db: Session,
        record_id: int,
        approved_by: int,
        approval_reason: Optional[str] = None
    ) -> Tuple[bool, str]:
        """批准补偿记录"""
        try:
            record = db.query(CompensationRecord).filter(CompensationRecord.id == record_id).first()
            if not record:
                return False, "补偿记录不存在"

            if record.status != CompensationStatus.PENDING:
                return False, f"补偿记录状态为 {record.status.value}，无法批准"

            user = db.query(User).filter(User.id == record.user_id).first()
            if not user:
                return False, "用户不存在"

            # 批准补偿
            record.status = CompensationStatus.APPROVED
            record.approved_by = approved_by
            record.approval_reason = approval_reason
            record.approved_at = datetime.now()

            # 退款给用户
            balance_before = user.balance_cpu_hours
            user.balance_cpu_hours += record.amount

            # 记录交易
            trans = QuotaTransaction(
                user_id=user.id,
                type="refund",
                amount=record.amount,
                balance_before=balance_before,
                balance_after=user.balance_cpu_hours,
                reference_id=record.id,
                reference_type="compensation",
                description=f"补偿记录 #{record.id} 批准，退款 {record.amount:.2f} 核时"
            )
            db.add(trans)

            # 标记为已完成
            record.status = CompensationStatus.COMPLETED
            record.completed_at = datetime.now()

            db.commit()
            logger.info(f"Compensation record {record_id} approved and completed")
            return True, f"补偿已批准并完成，退款 {record.amount:.2f} 核时"
        except Exception as e:
            db.rollback()
            logger.error(f"Failed to approve compensation: {e}")
            return False, f"批准补偿失败: {str(e)}"

    @staticmethod
    def reject_compensation(
        db: Session,
        record_id: int,
        rejection_reason: str
    ) -> Tuple[bool, str]:
        """拒绝补偿记录"""
        try:
            record = db.query(CompensationRecord).filter(CompensationRecord.id == record_id).first()
            if not record:
                return False, "补偿记录不存在"

            if record.status != CompensationStatus.PENDING:
                return False, f"补偿记录状态为 {record.status.value}，无法拒绝"

            record.status = CompensationStatus.REJECTED
            record.approval_reason = rejection_reason
            db.commit()

            logger.info(f"Compensation record {record_id} rejected")
            return True, "补偿已拒绝"
        except Exception as e:
            db.rollback()
            logger.error(f"Failed to reject compensation: {e}")
            return False, f"拒绝补偿失败: {str(e)}"

