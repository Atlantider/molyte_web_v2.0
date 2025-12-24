"""
Compensation rules models for automatic CPU hours refund and expiration
"""
from sqlalchemy import Column, Integer, String, Float, Boolean, DateTime, ForeignKey, Text, JSON, func, Enum
from sqlalchemy.orm import relationship
from datetime import datetime
import enum

from app.database import Base


class CompensationRuleType(str, enum.Enum):
    """补偿规则类型"""
    JOB_FAILURE = "job_failure"  # 任务失败时退款
    EXPIRATION = "expiration"  # 核时过期
    MANUAL = "manual"  # 手动补偿


class CompensationStatus(str, enum.Enum):
    """补偿状态"""
    PENDING = "pending"  # 待处理
    APPROVED = "approved"  # 已批准
    REJECTED = "rejected"  # 已拒绝
    COMPLETED = "completed"  # 已完成


class CompensationRule(Base):
    """补偿规则配置"""
    __tablename__ = "compensation_rules"

    id = Column(Integer, primary_key=True, index=True)
    name = Column(String(200), nullable=False)  # 规则名称
    description = Column(Text, nullable=True)  # 规则描述
    
    # 规则类型和配置
    rule_type = Column(Enum(CompensationRuleType), nullable=False)  # 规则类型
    is_active = Column(Boolean, default=True, nullable=False)  # 是否激活
    
    # 规则配置（JSON格式）
    config = Column(JSON, default={}, nullable=False)
    # 对于 JOB_FAILURE: {
    #   "failure_types": ["FAILED", "TIMEOUT"],  # 哪些失败类型触发补偿
    #   "refund_percentage": 100,  # 退款百分比（0-100）
    #   "min_cpu_hours": 0.1  # 最小退款核时
    # }
    # 对于 EXPIRATION: {
    #   "expiration_days": 365,  # 多少天未使用后过期
    #   "warning_days": 30  # 提前多少天发送警告
    # }
    
    # 时间戳
    created_at = Column(DateTime(timezone=True), server_default=func.now(), nullable=False)
    updated_at = Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now(), nullable=False)
    
    # 关系
    compensations = relationship("CompensationRecord", back_populates="rule", cascade="all, delete-orphan")

    def __repr__(self):
        return f"<CompensationRule(name={self.name}, type={self.rule_type})>"


class CompensationRecord(Base):
    """补偿记录"""
    __tablename__ = "compensation_records"

    id = Column(Integer, primary_key=True, index=True)
    rule_id = Column(Integer, ForeignKey("compensation_rules.id", ondelete="SET NULL"), nullable=True)
    user_id = Column(Integer, ForeignKey("users.id", ondelete="CASCADE"), nullable=False, index=True)
    
    # 补偿信息
    status = Column(Enum(CompensationStatus), default=CompensationStatus.PENDING, nullable=False)
    reason = Column(String(500), nullable=False)  # 补偿原因
    amount = Column(Float, nullable=False)  # 补偿核时数
    
    # 关联信息
    reference_id = Column(Integer, nullable=True)  # 关联的任务ID或其他ID
    reference_type = Column(String(50), nullable=True)  # job / order / manual
    
    # 审批信息
    approved_by = Column(Integer, ForeignKey("users.id", ondelete="SET NULL"), nullable=True)  # 审批人
    approval_reason = Column(String(500), nullable=True)  # 审批意见
    
    # 时间戳
    created_at = Column(DateTime(timezone=True), server_default=func.now(), nullable=False)
    approved_at = Column(DateTime(timezone=True), nullable=True)
    completed_at = Column(DateTime(timezone=True), nullable=True)
    
    # 关系
    rule = relationship("CompensationRule", back_populates="compensations")
    user = relationship("User", foreign_keys=[user_id])
    approver = relationship("User", foreign_keys=[approved_by])

    def __repr__(self):
        return f"<CompensationRecord(user_id={self.user_id}, amount={self.amount}, status={self.status})>"


class CPUHoursExpiration(Base):
    """核时过期记录"""
    __tablename__ = "cpu_hours_expirations"

    id = Column(Integer, primary_key=True, index=True)
    user_id = Column(Integer, ForeignKey("users.id", ondelete="CASCADE"), nullable=False, index=True)
    
    # 过期信息
    amount = Column(Float, nullable=False)  # 过期的核时数
    reason = Column(String(500), nullable=False)  # 过期原因
    
    # 时间戳
    created_at = Column(DateTime(timezone=True), server_default=func.now(), nullable=False)
    expired_at = Column(DateTime(timezone=True), nullable=False)  # 过期时间
    
    # 关系
    user = relationship("User", foreign_keys=[user_id])

    def __repr__(self):
        return f"<CPUHoursExpiration(user_id={self.user_id}, amount={self.amount})>"

