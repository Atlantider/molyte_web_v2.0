"""
定价配置模型
Pricing Configuration Models

支持灵活的计费系统,包括:
- 全局计费模式配置
- 任务类型核时单价
- 任务包定价
- 超限计费策略
- 用户折扣配置
"""
from sqlalchemy import Column, Integer, String, Float, DateTime, Boolean, ForeignKey, Text, Enum as SQLEnum
from sqlalchemy.dialects.postgresql import JSONB
from sqlalchemy.orm import relationship
from sqlalchemy.sql import func
from app.database import Base
import enum


class BillingMode(str, enum.Enum):
    """计费模式"""
    CORE_HOUR = "CORE_HOUR"  # 按核时计费
    PACKAGE = "PACKAGE"      # 按任务包计费


class TaskType(str, enum.Enum):
    """任务类型"""
    MD = "MD"
    QC = "QC"
    POSTPROCESS = "POSTPROCESS"
    REACTION_NETWORK = "REACTION_NETWORK"
    FORCEFIELD = "FORCEFIELD"


class PackageType(str, enum.Enum):
    """任务包类型"""
    # MD包
    FAST = "FAST"
    STANDARD = "STANDARD"
    PRECISE = "PRECISE"
    CUSTOM = "CUSTOM"
    
    # QC包
    QC_FAST = "QC_FAST"
    QC_PRECISE = "QC_PRECISE"
    
    # Postprocess包
    LIGAND = "LIGAND"
    DIMER = "DIMER"
    CLUSTER = "CLUSTER"
    
    # 反应网络包
    SMALL = "SMALL"
    MEDIUM = "MEDIUM"
    LARGE = "LARGE"


class OverflowStrategy(str, enum.Enum):
    """超限处理策略"""
    PROPORTIONAL = "PROPORTIONAL"              # 按比例计费
    SWITCH_TO_CORE_HOUR = "SWITCH_TO_CORE_HOUR"  # 切换到按核时计费


class GlobalBillingConfig(Base):
    """全局计费模式配置"""
    __tablename__ = "global_billing_configs"
    
    id = Column(Integer, primary_key=True, index=True)
    task_type = Column(String(50), unique=True, nullable=False, index=True, comment="任务类型")
    billing_mode = Column(String(20), nullable=False, comment="计费模式: CORE_HOUR/PACKAGE")
    is_active = Column(Boolean, default=True, nullable=False, comment="是否启用")
    
    # 审计字段
    updated_by = Column(Integer, ForeignKey("users.id", ondelete="SET NULL"), comment="更新者ID")
    updated_at = Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now(), comment="更新时间")
    created_at = Column(DateTime(timezone=True), server_default=func.now(), comment="创建时间")
    
    def __repr__(self):
        return f"<GlobalBillingConfig(task_type={self.task_type}, mode={self.billing_mode})>"


class TaskTypePricing(Base):
    """任务类型核时单价表"""
    __tablename__ = "task_type_pricing"
    
    id = Column(Integer, primary_key=True, index=True)
    task_type = Column(String(50), unique=True, nullable=False, index=True, comment="任务类型")
    price_per_hour = Column(Float, nullable=False, comment="核时单价(元/核时)")
    is_active = Column(Boolean, default=True, nullable=False, comment="是否启用")
    
    # 审计字段
    updated_by = Column(Integer, ForeignKey("users.id", ondelete="SET NULL"), comment="更新者ID")
    updated_at = Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now(), comment="更新时间")
    created_at = Column(DateTime(timezone=True), server_default=func.now(), comment="创建时间")
    
    def __repr__(self):
        return f"<TaskTypePricing(task_type={self.task_type}, price={self.price_per_hour})>"


class TaskPackage(Base):
    """任务包定价表"""
    __tablename__ = "task_packages"
    
    id = Column(Integer, primary_key=True, index=True)
    task_type = Column(String(50), nullable=False, index=True, comment="任务类型")
    package_type = Column(String(50), nullable=False, index=True, comment="包类型")
    package_name = Column(String(100), nullable=False, comment="包名称")
    base_price = Column(Float, nullable=False, comment="基础价格(元)")
    
    # 配置参数 (JSONB)
    # MD: {npt_time_ns: 2.0, nvt_time_ns: 5.0}
    # 反应网络: {max_generations: 3, max_species: 50}
    # QC: {method: "B3LYP", basis_set: "6-31G*"}
    config_params = Column(JSONB, default={}, comment="配置参数")
    
    is_active = Column(Boolean, default=True, nullable=False, comment="是否启用")
    description = Column(Text, comment="描述")
    
    # 审计字段
    updated_by = Column(Integer, ForeignKey("users.id", ondelete="SET NULL"), comment="更新者ID")
    updated_at = Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now(), comment="更新时间")
    created_at = Column(DateTime(timezone=True), server_default=func.now(), comment="创建时间")
    
    def __repr__(self):
        return f"<TaskPackage(task_type={self.task_type}, package={self.package_type}, price={self.base_price})>"


class OverflowPricingConfig(Base):
    """超限计费配置"""
    __tablename__ = "overflow_pricing_configs"
    
    id = Column(Integer, primary_key=True, index=True)
    task_type = Column(String(50), unique=True, nullable=False, index=True, comment="任务类型")
    strategy = Column(String(50), nullable=False, comment="超限策略: PROPORTIONAL/SWITCH_TO_CORE_HOUR")
    base_package_type = Column(String(50), nullable=False, comment="用于比例计算的基准包")
    
    # 审计字段
    updated_by = Column(Integer, ForeignKey("users.id", ondelete="SET NULL"), comment="更新者ID")
    updated_at = Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now(), comment="更新时间")
    created_at = Column(DateTime(timezone=True), server_default=func.now(), comment="创建时间")
    
    def __repr__(self):
        return f"<OverflowPricingConfig(task_type={self.task_type}, strategy={self.strategy})>"


class MDCustomPackageConfig(Base):
    """MD自定义包配置"""
    __tablename__ = "md_custom_package_configs"
    
    id = Column(Integer, primary_key=True, index=True)
    base_package_type = Column(String(50), nullable=False, comment="基准包类型: FAST/STANDARD/PRECISE")
    calculation_formula = Column(Text, comment="计算公式说明")
    
    # 审计字段
    updated_by = Column(Integer, ForeignKey("users.id", ondelete="SET NULL"), comment="更新者ID")
    updated_at = Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now(), comment="更新时间")
    created_at = Column(DateTime(timezone=True), server_default=func.now(), comment="创建时间")
    
    def __repr__(self):
        return f"<MDCustomPackageConfig(base_package={self.base_package_type})>"


class UserDiscountConfig(Base):
    """用户折扣配置"""
    __tablename__ = "user_discount_configs"
    
    id = Column(Integer, primary_key=True, index=True)
    
    # 用户类型折扣 (user_type不为空时)
    user_type = Column(String(50), index=True, comment="用户类型: STUDENT/RESEARCHER/COMPANY/PREMIUM")
    
    # 特定用户折扣 (custom_user_id不为空时)
    custom_user_id = Column(Integer, ForeignKey("users.id", ondelete="CASCADE"), index=True, comment="特定用户ID")
    
    discount_rate = Column(Float, nullable=False, comment="折扣率 (0.7表示7折)")
    is_active = Column(Boolean, default=True, nullable=False, comment="是否启用")
    
    # 审计字段
    updated_by = Column(Integer, ForeignKey("users.id", ondelete="SET NULL"), comment="更新者ID")
    updated_at = Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now(), comment="更新时间")
    created_at = Column(DateTime(timezone=True), server_default=func.now(), comment="创建时间")
    
    # 关系
    custom_user = relationship("User", foreign_keys=[custom_user_id])
    
    def __repr__(self):
        if self.user_type:
            return f"<UserDiscountConfig(user_type={self.user_type}, rate={self.discount_rate})>"
        else:
            return f"<UserDiscountConfig(user_id={self.custom_user_id}, rate={self.discount_rate})>"
