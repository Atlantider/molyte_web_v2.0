"""
User model
"""
from sqlalchemy import Column, Integer, String, DateTime, Enum, Boolean, Float, JSON, Text
from sqlalchemy.orm import relationship
from sqlalchemy.sql import func
from app.database import Base
import enum


class UserRole(str, enum.Enum):
    """User role enumeration"""
    ADMIN = "ADMIN"
    PREMIUM = "PREMIUM"
    USER = "USER"
    GUEST = "GUEST"


class UserType(str, enum.Enum):
    """User type enumeration - 用户类型"""
    STUDENT = "STUDENT"        # 学生
    RESEARCHER = "RESEARCHER"  # 研究者
    COMPANY = "COMPANY"        # 企业用户


class AccountType(str, enum.Enum):
    """Account type enumeration - 账号类型（互斥）"""
    PERSONAL = "personal"              # 个人用户
    MASTER_ACCOUNT = "master_account"  # 主账号
    SUB_ACCOUNT = "sub_account"        # 子账号


class BillingMode(str, enum.Enum):
    """Billing mode enumeration - 计费模式"""
    CORE_HOUR = "CORE_HOUR"    # 按核时计费（默认）
    TASK_TYPE = "TASK_TYPE"    # 按任务类型计费


class User(Base):
    """User model"""
    __tablename__ = "users"

    id = Column(Integer, primary_key=True, index=True)
    email = Column(String, unique=True, index=True, nullable=False)
    username = Column(String, unique=True, index=True, nullable=False)
    password_hash = Column(String, nullable=False)
    role = Column(Enum(UserRole), default=UserRole.USER, nullable=False)

    # 手机号（可选，用于验证）
    phone = Column(String(20), unique=True, index=True, nullable=True)
    phone_verified = Column(Boolean, default=False, nullable=False)

    # 用户类型和单位信息（用户注册时填写的个人信息）
    user_type = Column(Enum(UserType), default=UserType.STUDENT, nullable=False)  # 用户类型
    organization = Column(String(200), nullable=True)  # 单位名称
    department = Column(String(100), nullable=True)    # 部门（可选）
    real_name = Column(String(50), nullable=True)  # 真实姓名

    # 账号类型（互斥）- 用于区分个人用户、主账号、子账号
    # 临时使用 String 类型避免 enum 缓存问题
    account_type = Column(String(50), default=AccountType.PERSONAL.value, nullable=False, index=True)

    # 邮箱验证
    email_verified = Column(Boolean, default=False, nullable=False)  # 邮箱是否验证
    verification_token = Column(String(100), nullable=True)          # 验证令牌
    verification_expires = Column(DateTime(timezone=True), nullable=True)  # 令牌过期时间

    # Status
    is_active = Column(Boolean, default=True, nullable=False)
    last_login_at = Column(DateTime(timezone=True), nullable=True)

    # ============================================================================
    # 核心配额系统 (经济控制 - Single Source of Truth)
    # ============================================================================
    balance_cpu_hours = Column(Float, default=100.0, nullable=False,
        comment="可用核时余额：正数=可用，负数=欠费，0=无余额（唯一真实来源）")
    frozen_cpu_hours = Column(Float, default=0.0, nullable=False,
        comment="冻结核时：运行中任务占用的核时")

    # 核时来源追踪（用于统计和分析，不影响实际可用核时）
    free_cpu_hours_granted = Column(Float, default=100.0, nullable=False,
        comment="统计：初始赠送的免费核时")
    recharge_cpu_hours = Column(Float, default=0.0, nullable=False,
        comment="统计：充值获得的核时")
    admin_granted_cpu_hours = Column(Float, default=0.0, nullable=False,
        comment="统计：管理员赠送的核时")

    # ============================================================================
    # 资源保护系统 (集群保护 - Resource Protection)
    # ============================================================================
    # 防止资源滥用和集群过载，与余额配额互补
    concurrent_job_limit = Column(Integer, default=5, nullable=False,
        comment="并发任务硬限制：同时运行的任务数上限，0=不限制（仅ADMIN）")
    daily_job_limit = Column(Integer, default=20, nullable=False,
        comment="每日任务硬限制：每日提交任务数上限，0=不限制（仅ADMIN）")
    storage_quota_gb = Column(Float, default=50.0, nullable=False,
        comment="存储配额软限制：超出90%时警告但不阻止，0=不限制")
    
    # 资源保护开关
    enable_resource_limits = Column(Boolean, default=True, nullable=False,
        comment="是否启用资源保护限制（管理员可关闭，普通用户启用）")

    # Queue/Partition permissions (JSON array of allowed partition names)
    allowed_partitions = Column(JSON, nullable=True)
    allowed_modules = Column(JSON, nullable=True)  # 允许访问的模块列表
    
    # QC Engine permissions
    can_use_gaussian = Column(Boolean, default=False, nullable=False)  # 是否允许使用Gaussian(需license)
    
    custom_cpu_hour_price = Column(Float, nullable=True)  # 自定义核时单价
    billing_mode = Column(String(20), default=BillingMode.CORE_HOUR.value, nullable=False,
        comment="计费模式: CORE_HOUR=按核时, TASK_TYPE=按任务类型")
    custom_task_prices = Column(JSON, nullable=True,
        comment="自定义任务价格: {MD: 10.0, QC: 50.0, ...}")
    price_updated_at = Column(DateTime(timezone=True), nullable=True)  # 定价最后更新时间
    price_updated_by = Column(Integer, nullable=True)  # 更新定价的管理员 ID

    # Timestamps
    created_at = Column(DateTime(timezone=True), server_default=func.now(), nullable=False)
    updated_at = Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now(), nullable=False)

    # Relationships
    projects = relationship("Project", back_populates="owner", cascade="all, delete-orphan")
    md_jobs = relationship("MDJob", back_populates="user", cascade="all, delete-orphan",
                          foreign_keys="MDJob.user_id")
    qc_jobs = relationship("QCJob", back_populates="user", cascade="all, delete-orphan",
                          foreign_keys="QCJob.user_id")
    resp_jobs = relationship("RESPJob", back_populates="user", cascade="all, delete-orphan",
                            foreign_keys="RESPJob.user_id")
    usage_stats = relationship("UserUsageStats", back_populates="user", cascade="all, delete-orphan")
    recharge_orders = relationship("RechargeOrder", back_populates="user", cascade="all, delete-orphan")
    quota_transactions = relationship("QuotaTransaction", back_populates="user", cascade="all, delete-orphan")

    # User preferences
    solvent_combinations = relationship("UserSolventCombination", back_populates="user", cascade="all, delete-orphan")
    ion_combinations = relationship("UserIonCombination", back_populates="user", cascade="all, delete-orphan")

    # Anion generation jobs
    anion_generation_jobs = relationship("AnionGenerationJob", back_populates="user", cascade="all, delete-orphan",
                                        foreign_keys="AnionGenerationJob.user_id")
    
    # Reaction network jobs
    reaction_network_jobs = relationship("ReactionNetworkJob", back_populates="user", cascade="all, delete-orphan",
                                        foreign_keys="ReactionNetworkJob.user_id")


    # 贡献统计和配额使用
    public_data_count = Column(Integer, default=0, nullable=False)  # 公开数据数量
    contribution_points = Column(Float, default=0.0, nullable=False)  # 贡献积分
    private_quota_used = Column(Integer, default=0, nullable=False)  # 私有配额已使用
    private_quota_limit = Column(Integer, default=0, nullable=False)  # 私有配额限制

    def __repr__(self):
        return f"<User(id={self.id}, username={self.username}, email={self.email}, role={self.role}, account_type={self.account_type}, user_type={self.user_type})>"


# 用户类型配额定义
USER_TYPE_QUOTAS = {
    UserType.STUDENT: {
        "max_delay_years": 1,
        "default_cpu_hours": 100.0,
        "daily_job_limit": 10,
        "concurrent_job_limit": 3,
        "storage_quota_gb": 10.0,
    },
    UserType.RESEARCHER: {
        "max_delay_years": 3,
        "default_cpu_hours": 500.0,
        "daily_job_limit": 50,
        "concurrent_job_limit": 10,
        "storage_quota_gb": 100.0,
    },
    UserType.COMPANY: {
        "max_delay_years": 5,
        "default_cpu_hours": 5000.0,
        "daily_job_limit": 500,
        "concurrent_job_limit": 100,
        "storage_quota_gb": 1000.0,
    },
}
