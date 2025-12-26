"""
Improved organization models with hierarchy, master/sub accounts, and approval system
"""
from sqlalchemy import Column, Integer, String, Float, Boolean, DateTime, ForeignKey, Text, JSON, func, Enum
from sqlalchemy.orm import relationship
from datetime import datetime
import enum

from app.database import Base


class MasterAccount(Base):
    """主账号模型 - 用于管理主账号下的子账号"""
    __tablename__ = "master_accounts"

    id = Column(Integer, primary_key=True, index=True)

    # 主账号信息
    user_id = Column(Integer, ForeignKey("users.id", ondelete="CASCADE"), nullable=False, unique=True, index=True)
    # 注意：移除了 organization_id 外键，因为组织表不存在，主账号直接关联用户

    # 注意：主账号的配额来自 User 表（balance_cpu_hours, frozen_cpu_hours 等）
    # 不在此表中重复存储，避免数据不一致

    # 子账号限制
    max_sub_accounts = Column(Integer, default=10, nullable=False)
    current_sub_accounts = Column(Integer, default=0, nullable=False)

    # 状态
    is_active = Column(Boolean, default=True, nullable=False, index=True)

    # 时间戳
    created_at = Column(DateTime(timezone=True), server_default=func.now(), nullable=False)
    updated_at = Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now(), nullable=False)

    # 关系
    sub_accounts = relationship("SubAccount", back_populates="master_account", cascade="all, delete-orphan")

    def __repr__(self):
        return f"<MasterAccount(user_id={self.user_id})>"


class SubAccount(Base):
    """子账号模型 - 用于管理主账号下的子账号"""
    __tablename__ = "sub_accounts"

    id = Column(Integer, primary_key=True, index=True)

    # 账号关系
    master_account_id = Column(Integer, ForeignKey("master_accounts.id", ondelete="CASCADE"), nullable=False, index=True)
    user_id = Column(Integer, ForeignKey("users.id", ondelete="CASCADE"), nullable=False, unique=True, index=True)

    # 配额信息
    # allocated_quota: 主账号分配给子账号的配额（子账号实际可用 = min(User.balance_cpu_hours, allocated_quota)）
    allocated_quota = Column(Float, default=0.0, nullable=False)  # 主账号分配的配额

    # 状态
    is_active = Column(Boolean, default=True, nullable=False, index=True)

    # 时间戳
    created_at = Column(DateTime(timezone=True), server_default=func.now(), nullable=False)
    updated_at = Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now(), nullable=False)

    # 关系
    master_account = relationship("MasterAccount", back_populates="sub_accounts")

    def __repr__(self):
        return f"<SubAccount(user_id={self.user_id}, allocated_quota={self.allocated_quota})>"


