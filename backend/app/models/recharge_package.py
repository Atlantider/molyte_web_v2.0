"""
充值套餐模型
"""
import enum
from datetime import datetime
from sqlalchemy import Column, Integer, String, Float, Boolean, DateTime, Text, Enum, func
from app.database import Base


class UserTypeEnum(str, enum.Enum):
    """用户类型"""
    STUDENT = "STUDENT"           # 学生
    RESEARCHER = "RESEARCHER"     # 研究人员
    COMPANY = "COMPANY"           # 企业
    TEAM = "TEAM"                 # 合作团队


class RechargePackage(Base):
    """充值套餐表"""
    __tablename__ = "recharge_packages"

    id = Column(Integer, primary_key=True, index=True)
    
    # 基本信息
    name = Column(String(100), nullable=False)                    # 套餐名称，如"学生基础套餐"
    description = Column(Text, nullable=True)                     # 套餐描述
    user_type = Column(String(50), nullable=False, index=True)    # 目标用户类型 (STUDENT, RESEARCHER, COMPANY, TEAM)
    
    # 价格和配额
    price = Column(Float, nullable=False)                         # 价格（元）
    cpu_hours = Column(Float, nullable=False)                     # 对应的核时数
    
    # 显示信息
    display_order = Column(Integer, default=0)                    # 显示顺序
    badge = Column(String(50), nullable=True)                     # 徽章，如"推荐"、"热销"
    color = Column(String(20), default="#1890ff")                 # 卡片颜色
    icon = Column(String(50), nullable=True)                      # 图标
    
    # 状态
    is_active = Column(Boolean, default=True, index=True)         # 是否启用
    
    # 时间戳
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    updated_at = Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now())
    
    def __repr__(self):
        return f"<RechargePackage(id={self.id}, name={self.name}, user_type={self.user_type}, price={self.price})>"


class PackageDiscount(Base):
    """套餐折扣表（用于批量购买或特殊活动）"""
    __tablename__ = "package_discounts"

    id = Column(Integer, primary_key=True, index=True)
    
    # 关联信息
    package_id = Column(Integer, nullable=False, index=True)      # 套餐ID
    
    # 折扣信息
    name = Column(String(100), nullable=False)                    # 折扣名称，如"买3送1"
    discount_type = Column(String(20), nullable=False)            # 折扣类型: percentage(百分比), fixed(固定金额), bonus(赠送)
    discount_value = Column(Float, nullable=False)                # 折扣值
    
    # 条件
    min_quantity = Column(Integer, default=1)                     # 最少购买数量
    max_quantity = Column(Integer, nullable=True)                 # 最多购买数量（NULL表示无限制）
    
    # 状态
    is_active = Column(Boolean, default=True, index=True)
    
    # 时间戳
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    updated_at = Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now())
    
    def __repr__(self):
        return f"<PackageDiscount(id={self.id}, package_id={self.package_id}, name={self.name})>"

