"""
计费服务 - 配额管理、充值、扣费等核心逻辑
"""
from datetime import datetime, timedelta
from typing import Optional, Dict, Any, Tuple
from sqlalchemy.orm import Session
from sqlalchemy import func
import uuid

from app.models.user import User
from app.models.job import MDJob, JobStatus, AdvancedClusterJob
from app.models.billing import (
    SystemConfig, RechargeOrder, QuotaTransaction,
    PaymentMethod, PaymentStatus, TransactionType
)
from app.models.user_stats import UserUsageStats
from app.core.logger import logger
from datetime import date


class BillingService:
    """计费服务类"""
    
    # 默认配置
    DEFAULT_CPU_HOUR_PRICE = 0.1  # 元/核时
    DEFAULT_MIN_RECHARGE = 10.0   # 最低充值金额
    DEFAULT_MAX_DEBT = 100.0      # 最大欠费机时
    
    @staticmethod
    def get_config(db: Session, key: str, default: str = None) -> Optional[str]:
        """获取系统配置"""
        config = db.query(SystemConfig).filter(SystemConfig.key == key).first()
        return config.value if config else default
    
    @staticmethod
    def set_config(db: Session, key: str, value: str, description: str = None, updated_by: int = None):
        """设置系统配置"""
        config = db.query(SystemConfig).filter(SystemConfig.key == key).first()
        if config:
            config.value = value
            config.updated_at = datetime.now()
            config.updated_by = updated_by
            if description:
                config.description = description
        else:
            config = SystemConfig(
                key=key,
                value=value,
                description=description,
                updated_by=updated_by
            )
            db.add(config)
        db.commit()
        return config
    
    @staticmethod
    def get_cpu_hour_price(db: Session, user_id: int = None) -> float:
        """
        获取机时单价

        优先级：
        1. 用户自定义价格（如果设置了）
        2. 全局系统配置价格
        3. 默认价格

        Args:
            db: 数据库会话
            user_id: 用户 ID（可选）

        Returns:
            float: 机时单价（元/核时）
        """
        # 如果指定了用户 ID，先检查用户自定义价格
        if user_id:
            user = db.query(User).filter(User.id == user_id).first()
            if user and user.custom_cpu_hour_price is not None:
                return float(user.custom_cpu_hour_price)

        # 使用全局系统配置
        price = BillingService.get_config(db, 'cpu_hour_price')
        return float(price) if price else BillingService.DEFAULT_CPU_HOUR_PRICE
    
    @staticmethod
    def set_user_cpu_hour_price(db: Session, user_id: int, price: float, admin_id: int) -> Tuple[bool, str]:
        """
        设置用户的自定义核时单价

        Args:
            db: 数据库会话
            user_id: 用户 ID
            price: 新的核时单价（元/核时），为 None 时删除自定义价格
            admin_id: 管理员 ID

        Returns:
            Tuple[bool, str]: (成功标志, 消息)
        """
        user = db.query(User).filter(User.id == user_id).first()
        if not user:
            return False, f"用户 {user_id} 不存在"

        old_price = user.custom_cpu_hour_price
        user.custom_cpu_hour_price = price
        user.price_updated_at = datetime.now()
        user.price_updated_by = admin_id

        db.commit()

        if price is None:
            message = f"已删除用户 {user.username} 的自定义价格，将使用全局定价"
        else:
            message = f"已设置用户 {user.username} 的核时单价为 ¥{price:.4f}/核时（原价：¥{old_price or '全局定价'}/核时）"

        logger.info(f"Admin {admin_id} set user {user_id} price to {price}: {message}")
        return True, message

    @staticmethod
    def get_user_balance(db: Session, user_id: int) -> Dict[str, float]:
        """获取用户余额信息"""
        user = db.query(User).filter(User.id == user_id).first()
        if not user:
            return None

        # 计算欠费：balance_cpu_hours 为负数时表示欠费
        debt = max(0, -user.balance_cpu_hours)
        balance = max(0, user.balance_cpu_hours)

        return {
            "balance": balance,
            "frozen": user.frozen_cpu_hours,
            "debt": debt,
            "available": max(0, user.balance_cpu_hours - user.frozen_cpu_hours),
        }
    
    @staticmethod
    def can_submit_job(db: Session, user: User) -> Tuple[bool, str]:
        """
        检查用户是否可以提交任务
        
        DEPRECATED: 此方法已弃用，请使用 QuotaService.check_can_submit()
        保留此方法仅为向后兼容，将在未来版本中移除
        
        新的统一配额检查应使用：
        from app.services.quota_service import QuotaService
        can_submit, message, info = QuotaService.check_can_submit(
            user=user, 
            estimated_hours=hours,
            db=db
        )
        """
        import warnings
        warnings.warn(
            "BillingService.can_submit_job() is deprecated. "
            "Use QuotaService.check_can_submit() instead.",
            DeprecationWarning,
            stacklevel=2
        )
        
        # 重定向到新的统一接口
        from app.services.quota_service import QuotaService
        can_submit, message, info = QuotaService.check_can_submit(
            user=user,
            estimated_hours=1.0,  # 默认估算
            db=db
        )
        
        return can_submit, message
    
    @staticmethod
    def generate_order_no() -> str:
        """生成订单号"""
        now = datetime.now()
        return f"RC{now.strftime('%Y%m%d%H%M%S')}{uuid.uuid4().hex[:8].upper()}"
    
    @staticmethod
    def create_recharge_order(
        db: Session,
        user_id: int,
        amount: float,
        payment_method: PaymentMethod = PaymentMethod.SIMULATED,
        remark: str = None
    ) -> RechargeOrder:
        """创建充值订单"""
        # 使用用户级别的定价（如果有的话）
        price_per_hour = BillingService.get_cpu_hour_price(db, user_id)
        cpu_hours = amount / price_per_hour

        order = RechargeOrder(
            order_no=BillingService.generate_order_no(),
            user_id=user_id,
            amount=amount,
            cpu_hours=cpu_hours,
            price_per_hour=price_per_hour,
            payment_method=payment_method.value,  # 使用枚举值
            payment_status=PaymentStatus.PENDING.value,  # 使用枚举值
            expired_at=datetime.now() + timedelta(hours=2),
            remark=remark
        )
        db.add(order)
        db.commit()
        db.refresh(order)

        logger.info(f"Created recharge order: {order.order_no} for user {user_id}, amount={amount}")
        return order

    @staticmethod
    def complete_payment(db: Session, order_id: int, transaction_id: str = None) -> Tuple[bool, str]:
        """
        完成支付（模拟支付或回调确认）
        1. 增加用户余额
        2. 如有欠费，自动偿还
        3. 解锁被锁定的任务结果
        """
        order = db.query(RechargeOrder).filter(RechargeOrder.id == order_id).first()
        if not order:
            return False, "订单不存在"

        if order.payment_status != PaymentStatus.PENDING:
            return False, f"订单状态异常: {order.payment_status.value}"

        user = db.query(User).filter(User.id == order.user_id).first()
        if not user:
            return False, "用户不存在"

        # 更新订单状态
        order.payment_status = PaymentStatus.PAID.value  # 使用枚举值
        order.paid_at = datetime.now()
        order.transaction_id = transaction_id or f"SIM_{uuid.uuid4().hex[:12].upper()}"

        balance_before = user.balance_cpu_hours

        # 1. 统一的核时系统：直接增加 balance_cpu_hours
        # 如果 balance < 0（欠费），充值会先还债，然后增加可用余额
        user.balance_cpu_hours += order.cpu_hours

        # 更新充值核时统计（用于核时来源追踪）
        user.recharge_cpu_hours += order.cpu_hours

        # 记录充值流水
        recharge_trans = QuotaTransaction(
            user_id=user.id,
            type=TransactionType.RECHARGE.value,
            amount=order.cpu_hours,
            balance_before=balance_before,
            balance_after=user.balance_cpu_hours,
            reference_id=order.id,
            reference_type="order",
            description=f"充值 ¥{order.amount:.2f}，获得 {order.cpu_hours:.2f} 核时"
        )
        db.add(recharge_trans)

        # 2. 解锁该用户所有被锁定的任务结果（只要有余额就解锁）
        if user.balance_cpu_hours >= 0:
            locked_jobs = db.query(MDJob).filter(
                MDJob.user_id == user.id,
                MDJob.result_locked == True
            ).all()
            for job in locked_jobs:
                job.result_locked = False
                job.locked_reason = None
            logger.info(f"Unlocked {len(locked_jobs)} jobs for user {user.id}")

        db.commit()
        logger.info(f"Payment completed for order {order.order_no}, user {user.id} balance: {user.balance_cpu_hours}")

        return True, f"支付成功，获得 {order.cpu_hours:.2f} 核时"

    @staticmethod
    def calculate_job_cpu_hours(job: MDJob) -> float:
        """
        计算任务实际消耗的机时

        包括 MD 计算核时和 RESP 电荷计算核时的总和。
        使用从 Slurm 获取的实际核时（CPUTimeRAW），而不是时间差。
        这确保只计算真正在集群上运行的时间，不包括排队时间。

        核时计算公式：
        总核时 = MD 核时 (actual_cpu_hours) + RESP 核时 (resp_cpu_hours)
        """
        # 计算 MD 核时
        md_hours = job.actual_cpu_hours if job.actual_cpu_hours and job.actual_cpu_hours > 0 else 0.0

        # 计算 RESP 核时
        resp_hours = job.resp_cpu_hours if job.resp_cpu_hours and job.resp_cpu_hours > 0 else 0.0

        # 返回总核时
        total_hours = md_hours + resp_hours
        return total_hours

    @staticmethod
    def settle_job(db: Session, job: MDJob) -> Tuple[bool, str]:
        """
        任务完成后结算
        1. 计算实际消耗机时
        2. 从用户余额扣除
        3. 余额不足则记录欠费并锁定结果
        4. 更新用户使用统计
        """
        if job.billed:
            return True, "任务已结算"

        user = db.query(User).filter(User.id == job.user_id).first()
        if not user:
            return False, "用户不存在"

        # 管理员不扣费
        if user.role.value == "ADMIN":
            job.billed = True
            db.commit()
            return True, "管理员任务免费"

        # 计算实际机时（包括 MD 核时和 RESP 核时）
        total_cpu_hours = BillingService.calculate_job_cpu_hours(job)

        # 注意：不要覆盖 job.actual_cpu_hours 和 job.resp_cpu_hours
        # 这两个字段由 Worker 从 Slurm 获取并上报，应该保持原值
        # 这里只是用于计费的临时变量

        balance_before = user.balance_cpu_hours

        # 统一的核时系统：直接从 balance_cpu_hours 扣费
        # 如果 balance < 0，表示欠费
        user.balance_cpu_hours -= total_cpu_hours

        if user.balance_cpu_hours >= 0:
            # 余额足够
            job.result_locked = False
        else:
            # 余额不足，产生欠费（balance < 0）
            debt = abs(user.balance_cpu_hours)
            job.result_locked = True
            job.locked_reason = f"余额不足，欠费 {debt:.2f} 核时"
            logger.warning(f"Job {job.id} completed with debt: {debt:.2f} hours")

        # 记录消费流水
        trans = QuotaTransaction(
            user_id=user.id,
            type=TransactionType.CONSUME.value,  # 使用枚举值
            amount=-total_cpu_hours,
            balance_before=balance_before,
            balance_after=user.balance_cpu_hours,
            reference_id=job.id,
            reference_type="job",
            description=f"任务 #{job.id} 消耗 {total_cpu_hours:.2f} 核时 (MD: {job.actual_cpu_hours:.2f}h + RESP: {job.resp_cpu_hours:.2f}h)"
        )
        db.add(trans)

        # 更新用户使用统计
        from app.models.user_stats import UserUsageStats
        from datetime import date

        today = date.today()
        stats = db.query(UserUsageStats).filter(
            UserUsageStats.user_id == user.id,
            UserUsageStats.date == today
        ).first()

        if not stats:
            stats = UserUsageStats(
                user_id=user.id,
                date=today,
                jobs_submitted=0,
                jobs_completed=0,
                jobs_failed=0,
                jobs_cancelled=0,
                cpu_hours_used=0.0,
                cluster_analysis_cpu_hours=0.0,
                cluster_analysis_task_count=0,
                storage_used_gb=0.0,
                max_concurrent_jobs=0
            )
            db.add(stats)

        # 更新统计数据
        stats.jobs_completed += 1
        stats.cpu_hours_used += total_cpu_hours

        job.billed = True
        db.commit()

        return True, f"结算完成，消耗 {total_cpu_hours:.2f} 核时 (MD: {job.actual_cpu_hours:.2f}h + RESP: {job.resp_cpu_hours:.2f}h)"

    @staticmethod
    def get_transactions(db: Session, user_id: int, skip: int = 0, limit: int = 20):
        """获取用户配额变更记录"""
        return db.query(QuotaTransaction).filter(
            QuotaTransaction.user_id == user_id
        ).order_by(QuotaTransaction.created_at.desc()).offset(skip).limit(limit).all()

    @staticmethod
    def get_orders(db: Session, user_id: int, skip: int = 0, limit: int = 20):
        """获取用户充值订单"""
        return db.query(RechargeOrder).filter(
            RechargeOrder.user_id == user_id
        ).order_by(RechargeOrder.created_at.desc()).offset(skip).limit(limit).all()

    @staticmethod
    def admin_adjust_balance(
        db: Session,
        admin_id: int,
        user_id: int,
        amount: float,
        reason: str
    ) -> Tuple[bool, str]:
        """管理员手动调整用户余额"""
        user = db.query(User).filter(User.id == user_id).first()
        if not user:
            return False, "用户不存在"

        balance_before = user.balance_cpu_hours
        user.balance_cpu_hours += amount

        # 防止负数
        if user.balance_cpu_hours < 0:
            user.balance_cpu_hours = 0

        trans = QuotaTransaction(
            user_id=user.id,
            type=TransactionType.ADMIN_ADJUST.value,  # 使用枚举值
            amount=amount,
            balance_before=balance_before,
            balance_after=user.balance_cpu_hours,
            reference_id=admin_id,
            reference_type="admin",
            description=f"管理员调整: {reason}"
        )
        db.add(trans)

        # 如果调整后用户有余额，自动解锁所有被锁定的任务
        if user.balance_cpu_hours >= 0:
            locked_jobs = db.query(MDJob).filter(
                MDJob.user_id == user.id,
                MDJob.result_locked == True
            ).all()
            for job in locked_jobs:
                job.result_locked = False
                job.locked_reason = None
            logger.info(f"Unlocked {len(locked_jobs)} jobs for user {user_id} after admin adjust")

        db.commit()

        logger.info(f"Admin {admin_id} adjusted user {user_id} balance by {amount}: {reason}")
        return True, f"调整成功，当前余额 {user.balance_cpu_hours:.2f} 核时"

    @staticmethod
    def admin_grant_cpu_hours(
        db: Session,
        admin_id: int,
        user_id: int,
        amount: float,
        reason: str
    ) -> Tuple[bool, str]:
        """管理员赠送核时给用户（区别于调整余额）"""
        user = db.query(User).filter(User.id == user_id).first()
        if not user:
            return False, "用户不存在"

        if amount <= 0:
            return False, "赠送核时必须大于 0"

        balance_before = user.balance_cpu_hours

        # 同时更新余额和管理员赠送核时统计
        user.balance_cpu_hours += amount
        user.admin_granted_cpu_hours += amount

        trans = QuotaTransaction(
            user_id=user.id,
            type="admin_grant",  # 新的交易类型
            amount=amount,
            balance_before=balance_before,
            balance_after=user.balance_cpu_hours,
            reference_id=admin_id,
            reference_type="admin",
            description=f"管理员赠送: {reason}"
        )
        db.add(trans)

        # 如果赠送后用户有余额，自动解锁所有被锁定的任务
        if user.balance_cpu_hours >= 0:
            locked_jobs = db.query(MDJob).filter(
                MDJob.user_id == user.id,
                MDJob.result_locked == True
            ).all()
            for job in locked_jobs:
                job.result_locked = False
                job.locked_reason = None
            logger.info(f"Unlocked {len(locked_jobs)} jobs for user {user_id} after admin grant")

        db.commit()

        logger.info(f"Admin {admin_id} granted {amount} cpu hours to user {user_id}: {reason}")
        return True, f"赠送成功，用户 {user.username} 获得 {amount:.2f} 核时，当前余额 {user.balance_cpu_hours:.2f} 核时"

    @staticmethod
    def settle_cluster_analysis_job(db: Session, job: AdvancedClusterJob) -> Tuple[bool, str]:
        """
        后处理分析任务完成后结算
        1. 计算实际消耗核时
        2. 计算任务计数
        3. 从用户余额扣除
        4. 更新用户统计
        """
        if job.status.value != 'COMPLETED':
            return False, "任务未完成"

        user = db.query(User).filter(User.id == job.user_id).first()
        if not user:
            return False, "用户不存在"

        # 计算实际核时：所有关联 QC 任务的真实 Slurm CPU 核时总和
        from app.models.qc import QCJob
        qc_jobs = db.query(QCJob).filter(
            QCJob.cluster_analysis_job_id == job.id
        ).all()

        # 核时定义：所有 QC 任务的实际 Slurm CPUTimeRAW 总和（单位：小时）
        cpu_hours = sum(qc.actual_cpu_hours for qc in qc_jobs)

        job.cpu_hours_used = cpu_hours  # 记录实际使用的核时（用于显示）
        job.task_count = job.calculate_task_count()

        # 管理员不扣费，但要记录实际使用量
        if user.role.value == "ADMIN":
            db.commit()
            return True, f"管理员任务免费，实际使用 {cpu_hours:.2f} 核时"

        balance_before = user.balance_cpu_hours

        # 统一的核时系统：直接从 balance_cpu_hours 扣费
        user.balance_cpu_hours -= cpu_hours

        if user.balance_cpu_hours < 0:
            # 余额不足，产生欠费（balance < 0）
            debt = abs(user.balance_cpu_hours)
            logger.warning(f"Cluster analysis job {job.id} completed with debt: {debt:.2f} hours")

        # 记录消费流水
        trans = QuotaTransaction(
            user_id=user.id,
            type=TransactionType.CONSUME.value,
            amount=-cpu_hours,
            balance_before=balance_before,
            balance_after=user.balance_cpu_hours,
            reference_id=job.id,
            reference_type="cluster_analysis_job",
            description=f"后处理分析任务 #{job.id} 消耗 {cpu_hours:.2f} 核时"
        )
        db.add(trans)

        # 更新用户统计
        today = date.today()
        stats = db.query(UserUsageStats).filter(
            UserUsageStats.user_id == user.id,
            UserUsageStats.date == today
        ).first()

        if not stats:
            stats = UserUsageStats(
                user_id=user.id,
                date=today,
                jobs_submitted=0,
                jobs_completed=0,
                jobs_failed=0,
                jobs_cancelled=0,
                cpu_hours_used=0.0,
                cluster_analysis_cpu_hours=0.0,
                cluster_analysis_task_count=0,
                storage_used_gb=0.0,
                max_concurrent_jobs=0
            )
            db.add(stats)

        stats.jobs_completed += 1
        stats.cpu_hours_used += cpu_hours
        stats.cluster_analysis_cpu_hours += cpu_hours
        stats.cluster_analysis_task_count += job.task_count

        db.commit()

        return True, f"结算完成，消耗 {cpu_hours:.2f} 核时，任务计数 {job.task_count}"

    @staticmethod
    def settle_qc_job(db: Session, qc_job) -> Tuple[bool, str]:
        """
        QC任务完成后结算
        
        处理独立提交的QC计算任务的计费
        
        Args:
            db: 数据库会话
            qc_job: QC任务对象
            
        Returns:
            Tuple[bool, str]: (成功标志, 消息)
        """
        # 导入放在函数内避免循环依赖
        from app.models.qc import QCJob
        
        if qc_job.billed:
            return True, "任务已结算"
        
        user = db.query(User).filter(User.id == qc_job.user_id).first()
        if not user:
            return False, "用户不存在"
        
        # 管理员不扣费
        if user.role.value == "ADMIN":
            qc_job.billed = True
            db.commit()
            return True, "管理员任务免费"
        
        # 计算实际核时
        cpu_hours = qc_job.actual_cpu_hours if qc_job.actual_cpu_hours and qc_job.actual_cpu_hours > 0 else 0.0
        
        if cpu_hours == 0:
            # 没有核时记录，可能是复用任务或失败任务
            qc_job.billed = True
            db.commit()
            return True, "无需扣费（核时为0）"
        
        balance_before = user.balance_cpu_hours
        
        # 统一的核时系统：直接从 balance_cpu_hours 扣费
        user.balance_cpu_hours -= cpu_hours
        
        # 记录消费流水
        trans = QuotaTransaction(
            user_id=user.id,
            type=TransactionType.CONSUME.value,
            amount=-cpu_hours,
            balance_before=balance_before,
            balance_after=user.balance_cpu_hours,
            reference_id=qc_job.id,
            reference_type="qc_job",
            description=f"QC任务 #{qc_job.id} ({qc_job.molecule_name}) 消耗 {cpu_hours:.2f} 核时"
        )
        db.add(trans)
        
        # 更新用户统计
        today = date.today()
        stats = db.query(UserUsageStats).filter(
            UserUsageStats.user_id == user.id,
            UserUsageStats.date == today
        ).first()
        
        if not stats:
            stats = UserUsageStats(
                user_id=user.id,
                date=today,
                jobs_submitted=0,
                jobs_completed=0,
                jobs_failed=0,
                jobs_cancelled=0,
                cpu_hours_used=0.0,
                cluster_analysis_cpu_hours=0.0,
                cluster_analysis_task_count=0,
                storage_used_gb=0.0,
                max_concurrent_jobs=0
            )
            db.add(stats)
        
        stats.jobs_completed += 1
        stats.cpu_hours_used += cpu_hours
        
        qc_job.billed = True
        db.commit()
        
        logger.info(f"QC job {qc_job.id} settled: {cpu_hours:.2f} hours, balance: {user.balance_cpu_hours:.2f}")
        return True, f"结算完成，消耗 {cpu_hours:.2f} 核时"

