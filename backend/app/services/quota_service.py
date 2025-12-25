"""
Unified quota management service
Scheme B refactoring - Phase 2

Handles quota checking and consumption for all account types
"""

from typing import Tuple, Optional, Dict, Any
from sqlalchemy.orm import Session
from sqlalchemy import func
from app.models.user import User, AccountType, UserRole
from app.models.organization_v2 import MasterAccount, SubAccount  
from app.models.billing import QuotaTransaction, TransactionType
from app.models.job import MDJob, JobStatus
from datetime import datetime, date
import logging

logger = logging.getLogger(__name__)


class QuotaService:
    """Service for unified quota management"""
    
    @staticmethod
    def check_can_submit(
        user: User,
        estimated_hours: float,
        db: Session,
        job_type: str = "md"
    ) -> Tuple[bool, str, Dict[str, Any]]:
        """
        检查用户是否可以提交任务（唯一权威入口）
        
        检查顺序（任何一项失败立即返回）：
        1. 账户状态 (is_active, 子账号状态)
        2. 余额充足 (balance_cpu_hours >= estimated_hours) ← 核心
        3. 资源保护 (并发限制, 每日限制)
        
        Args:
            user: 用户对象
            estimated_hours: 预估核时
            db: 数据库会话
            job_type: 任务类型（用于统计）
            
        Returns:
            (can_submit: bool, message: str, info: dict)
            
            info = {
                "balance": float,          # 当前余额
                "available": float,        # 可用余额（余额-冻结）
                "running_jobs": int,       # 运行中任务数
                "today_jobs": int,        # 今日任务数
                "warnings": {              # 软限制警告
                    "low_balance": bool,
                    "storage_high": bool
                }
            }
        """
        # ============ 1. 账户状态检查 ============
        
        # 账户激活状态
        if not user.is_active:
            return False, "账户已被禁用，请联系管理员", {}
        
        # 子账号特殊检查
        if user.account_type == AccountType.SUB_ACCOUNT.value:
            sub_account = db.query(SubAccount).filter(
                SubAccount.user_id == user.id
            ).first()
            if sub_account and not sub_account.is_active:
                return False, "子账号已被主账号禁用", {}
        
        # ============ 2. 余额检查 (核心) ============
        
        # 管理员豁免余额检查
        if user.role == UserRole.ADMIN:
            logger.info(f"Admin user {user.id} bypassed balance check")
            available_balance = float('inf')  # 无限余额
        else:
            # 获取可用配额（考虑子账号分配池）
            available_balance = QuotaService.get_available_quota(user, db)
            
            # 余额不足 - 立即拒绝
            if available_balance < estimated_hours:
                return False, (
                    f"余额不足。需要 {estimated_hours:.2f} 核时，"
                    f"当前可用 {available_balance:.2f} 核时。请充值。"
                ), {
                    "balance": user.balance_cpu_hours,
                    "available": available_balance,
                    "required": estimated_hours
                }
        
        # ============ 3. 资源保护检查 ============
        
        warnings = {}
        
        # 检查是否启用资源限制
        if user.enable_resource_limits:
            
            # 3.1 并发任务限制 (硬)
            if user.concurrent_job_limit > 0:
                running_count = QuotaService._get_running_job_count(user.id, db)
                
                if running_count >= user.concurrent_job_limit:
                    return False, (
                        f"并发任务数已达上限。"
                        f"当前运行 {running_count} 个任务，"
                        f"限制 {user.concurrent_job_limit} 个。"
                        f"请等待任务完成后再提交。"
                    ), {
                        "running_jobs": running_count,
                        "limit": user.concurrent_job_limit
                    }
            
            # 3.2 每日任务限制 (硬)
            if user.daily_job_limit > 0:
                today_count = QuotaService._get_today_job_count(user.id, db)
                
                if today_count >= user.daily_job_limit:
                    return False, (
                        f"今日提交任务数已达上限。"
                        f"今日已提交 {today_count} 个任务，"
                        f"限制 {user.daily_job_limit} 个。"
                        f"明日 00:00 后可继续提交。"
                    ), {
                        "today_jobs": today_count,
                        "limit": user.daily_job_limit
                    }
            
            # 3.3 存储配额检查 (软 - 仅警告)
            if user.storage_quota_gb > 0:
                used_storage = QuotaService._get_user_storage_gb(user.id, db)
                if used_storage > user.storage_quota_gb * 0.9:  # 90%警告
                    warnings["storage_high"] = True
        
        # ============ 4. 余额预警 (软) ============
        
        # 余额低于100核时预警
        if available_balance < 100 and user.role != UserRole.ADMIN:
            warnings["low_balance"] = True
        
        # ============ 成功 ============
        
        running_count = QuotaService._get_running_job_count(user.id, db)
        today_count = QuotaService._get_today_job_count(user.id, db)
        
        return True, "检查通过，可以提交任务", {
            "balance": user.balance_cpu_hours,
            "available": available_balance if available_balance != float('inf') else 999999,
            "running_jobs": running_count,
            "today_jobs": today_count,
            "warnings": warnings
        }
    
    @staticmethod
    def _get_running_job_count(user_id: int, db: Session) -> int:
        """获取运行中任务数"""
        # 导入QCJob如果存在
        try:
            from app.models.qc_job import QCJob, QCJobStatus
            qc_count = db.query(QCJob).filter(
                QCJob.user_id == user_id,
                QCJob.status.in_(['queued', 'running'])
            ).count()
        except ImportError:
            qc_count = 0
        
        md_count = db.query(MDJob).filter(
            MDJob.user_id == user_id,
            MDJob.status.in_([JobStatus.QUEUED, JobStatus.RUNNING])
        ).count()
        
        return md_count + qc_count
    
    @staticmethod
    def _get_today_job_count(user_id: int, db: Session) -> int:
        """获取今日提交任务数"""
        try:
            from app.models.qc_job import QCJob
            qc_exists = True
        except ImportError:
            qc_exists = False
        
        today = date.today()
        
        md_count = db.query(MDJob).filter(
            MDJob.user_id == user_id,
            func.date(MDJob.created_at) == today
        ).count()
        
        if qc_exists:
            from app.models.qc_job import QCJob
            qc_count = db.query(QCJob).filter(
                QCJob.user_id == user_id,
                func.date(QCJob.created_at) == today
            ).count()
        else:
            qc_count = 0
        
        return md_count + qc_count
    
    @staticmethod
    def _get_user_storage_gb(user_id: int, db: Session) -> float:
        """获取用户已用存储（GB）"""
        # TODO: 实现存储统计逻辑
        # 可以从job work_dir统计文件大小
        return 0.0
    
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

