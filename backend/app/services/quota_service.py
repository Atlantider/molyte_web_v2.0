"""
Unified quota management service
Scheme B refactoring - Phase 2

Handles quota checking and consumption for all account types
"""

from typing import Tuple, Optional, Dict, Any
from sqlalchemy.orm import Session
from app.models.user import User, AccountType
from app.models.organization_v2 import MasterAccount, SubAccount
from app.models.billing import QuotaTransaction, TransactionType
from datetime import datetime
import logging

logger = logging.getLogger(__name__)


class QuotaService:
    """Service for unified quota management"""
    
    @staticmethod
    def get_available_quota(user: User, db: Session) -> float:
        """
        Get total available quota for a user

        Calculation:
        - Personal: user.balance_cpu_hours
        - Master account: user.balance_cpu_hours
        - Sub account: balance_cpu_hours + allocated_quota（两个池的总和）

        Args:
            user: User object
            db: Database session

        Returns:
            Total available CPU hours
        """
        account_type = user.account_type

        if account_type == AccountType.PERSONAL.value:
            return user.balance_cpu_hours

        elif account_type == AccountType.MASTER_ACCOUNT.value:
            # 主账号的配额来自 User 表
            return user.balance_cpu_hours

        elif account_type == AccountType.SUB_ACCOUNT.value:
            sub_account = db.query(SubAccount).filter(
                SubAccount.user_id == user.id
            ).first()
            if not sub_account:
                return 0.0

            # 子账号的总可用 = 个人充值余额 + 主账号分配配额（两个池的总和）
            personal_balance = user.balance_cpu_hours
            allocated_quota = sub_account.allocated_quota or 0.0
            return personal_balance + allocated_quota

        return 0.0
    
    @staticmethod
    def check_quota(user: User, required_hours: float, db: Session) -> Tuple[bool, str]:
        """
        Check if user has sufficient quota
        
        Args:
            user: User object
            required_hours: Required CPU hours
            db: Database session
            
        Returns:
            Tuple of (has_quota: bool, message: str)
        """
        available = QuotaService.get_available_quota(user, db)
        
        if available >= required_hours:
            return True, f"Sufficient quota available: {available} hours"
        else:
            return False, f"Insufficient quota. Required: {required_hours}, Available: {available}"
    
    @staticmethod
    def consume_quota(
        user: User,
        hours: float,
        db: Session,
        reason: str = "Job execution",
        job_id: Optional[int] = None
    ) -> Tuple[bool, str]:
        """
        Consume quota from user's account

        Consumption logic:
        - Personal users: deduct from user.balance_cpu_hours
        - Master account: deduct from user.balance_cpu_hours
        - Sub account: deduct from user.balance_cpu_hours (also check allocated_quota limit)

        Args:
            user: User object
            hours: Hours to consume
            db: Database session
            reason: Reason for consumption
            job_id: Related job ID

        Returns:
            Tuple of (success: bool, message: str)
        """
        account_type = user.account_type

        try:
            if account_type == AccountType.PERSONAL.value:
                if user.balance_cpu_hours < hours:
                    return False, "Insufficient personal quota"
                user.balance_cpu_hours -= hours
                user.frozen_cpu_hours += hours

            elif account_type == AccountType.MASTER_ACCOUNT.value:
                # 主账号的配额来自 User 表
                if user.balance_cpu_hours < hours:
                    return False, "Insufficient master account quota"
                user.balance_cpu_hours -= hours
                user.frozen_cpu_hours += hours

            elif account_type == AccountType.SUB_ACCOUNT.value:
                sub_account = db.query(SubAccount).filter(
                    SubAccount.user_id == user.id
                ).first()
                if not sub_account:
                    return False, "Sub account record not found"

                # 子账号的总可用 = 个人充值余额 + 主账号分配配额（两个池的总和）
                personal_balance = user.balance_cpu_hours
                allocated_quota = sub_account.allocated_quota or 0.0
                total_available = personal_balance + allocated_quota

                if total_available < hours:
                    return False, f"Insufficient quota. Personal: {personal_balance:.2f}, Allocated: {allocated_quota:.2f}, Total: {total_available:.2f}"

                # 优先从个人充值池消费，然后从主账号分配池消费
                personal_consumed = 0.0
                allocated_consumed = 0.0
                source = "personal"

                if personal_balance >= hours:
                    # 个人充值池足够
                    personal_consumed = hours
                    user.balance_cpu_hours -= hours
                else:
                    # 个人充值池不足，需要使用主账号分配池
                    personal_consumed = personal_balance
                    allocated_consumed = hours - personal_balance
                    user.balance_cpu_hours = 0.0
                    sub_account.allocated_quota -= allocated_consumed
                    source = "mixed" if personal_consumed > 0 else "allocated"

                user.frozen_cpu_hours += hours

                # 记录消费来源信息供后续释放使用
                # 这些信息会在创建 QuotaTransaction 时使用
                user._quota_source = source
                user._personal_consumed = personal_consumed
                user._allocated_consumed = allocated_consumed

            # Record transaction with source information
            source = getattr(user, '_quota_source', 'personal')
            personal_consumed = getattr(user, '_personal_consumed', hours if account_type != AccountType.SUB_ACCOUNT.value else 0.0)
            allocated_consumed = getattr(user, '_allocated_consumed', 0.0)

            transaction = QuotaTransaction(
                user_id=user.id,
                type=TransactionType.CONSUME.value,
                amount=hours,
                balance_before=0.0,  # 这里应该记录消费前的余额，但需要在消费前获取
                balance_after=user.balance_cpu_hours,
                reference_id=job_id,
                reference_type="job",
                source=source,
                personal_consumed=personal_consumed,
                allocated_consumed=allocated_consumed,
                description=reason
            )
            db.add(transaction)
            db.commit()

            logger.info(f"Consumed {hours} hours from user {user.id} ({user.username}), source: {source}")
            return True, f"Successfully consumed {hours} hours"

        except Exception as e:
            db.rollback()
            logger.error(f"Error consuming quota: {e}")
            return False, f"Error consuming quota: {str(e)}"
    
    @staticmethod
    def release_quota(
        user: User,
        hours: float,
        db: Session,
        reason: str = "Job cancellation",
        job_id: Optional[int] = None
    ) -> Tuple[bool, str]:
        """
        Release frozen quota back to available quota

        For sub-accounts, restore to the same source where it was consumed from.

        Args:
            user: User object
            hours: Hours to release
            db: Database session
            reason: Reason for release
            job_id: Related job ID to find the original transaction

        Returns:
            Tuple of (success: bool, message: str)
        """
        account_type = user.account_type

        try:
            if account_type == AccountType.PERSONAL.value:
                user.frozen_cpu_hours -= hours
                user.balance_cpu_hours += hours

            elif account_type == AccountType.MASTER_ACCOUNT.value:
                # 主账号的配额来自 User 表
                user.frozen_cpu_hours -= hours
                user.balance_cpu_hours += hours

            elif account_type == AccountType.SUB_ACCOUNT.value:
                # 子账号：根据消费来源恢复配额
                sub_account = db.query(SubAccount).filter(
                    SubAccount.user_id == user.id
                ).first()

                if not sub_account:
                    return False, "Sub account record not found"

                # 查找原始消费记录以获取消费来源
                source = "personal"
                personal_consumed = 0.0
                allocated_consumed = 0.0

                if job_id:
                    transaction = db.query(QuotaTransaction).filter(
                        QuotaTransaction.user_id == user.id,
                        QuotaTransaction.reference_id == job_id,
                        QuotaTransaction.type == TransactionType.CONSUMPTION.value
                    ).first()

                    if transaction:
                        source = transaction.source
                        personal_consumed = transaction.personal_consumed
                        allocated_consumed = transaction.allocated_consumed

                # 根据消费来源恢复配额
                if source == "personal":
                    # 全部从个人池消费，恢复到个人池
                    user.balance_cpu_hours += hours
                elif source == "allocated":
                    # 全部从分配池消费，恢复到分配池
                    sub_account.allocated_quota += hours
                elif source == "mixed":
                    # 从两个池消费，分别恢复
                    user.balance_cpu_hours += personal_consumed
                    sub_account.allocated_quota += allocated_consumed

                user.frozen_cpu_hours -= hours

            db.commit()
            logger.info(f"Released {hours} hours for user {user.id}")
            return True, f"Successfully released {hours} hours"

        except Exception as e:
            db.rollback()
            logger.error(f"Error releasing quota: {e}")
            return False, f"Error releasing quota: {str(e)}"

