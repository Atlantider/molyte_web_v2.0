"""
消息/通知模型
"""
import enum
from datetime import datetime
from sqlalchemy import Column, Integer, String, Text, Boolean, DateTime, ForeignKey, Index
from sqlalchemy.orm import relationship
from app.database import Base


class NotificationTypeEnum(str, enum.Enum):
    """消息类型"""
    SYSTEM = "SYSTEM"              # 系统通知
    DEBT_WARNING = "DEBT_WARNING"  # 欠费警告
    DEBT_CLEARED = "DEBT_CLEARED"  # 欠费已清
    JOB_COMPLETED = "JOB_COMPLETED"  # 任务完成
    JOB_FAILED = "JOB_FAILED"      # 任务失败
    QUOTA_LOW = "QUOTA_LOW"        # 配额不足
    QUOTA_RECHARGED = "QUOTA_RECHARGED"  # 配额已充值
    ACCOUNT_CREATED = "ACCOUNT_CREATED"  # 账户创建
    SUB_ACCOUNT_CREATED = "SUB_ACCOUNT_CREATED"  # 子账户创建
    ADMIN_ALERT = "ADMIN_ALERT"    # 管理员告警


class NotificationPriorityEnum(str, enum.Enum):
    """消息优先级"""
    LOW = "LOW"
    NORMAL = "NORMAL"
    HIGH = "HIGH"
    CRITICAL = "CRITICAL"


class Notification(Base):
    """消息/通知表"""
    __tablename__ = "notifications"

    id = Column(Integer, primary_key=True, index=True)
    user_id = Column(Integer, ForeignKey("users.id"), nullable=False, index=True)
    
    # 消息内容
    type = Column(String(50), nullable=False, index=True)  # NotificationTypeEnum
    title = Column(String(200), nullable=False)
    message = Column(Text, nullable=False)
    priority = Column(String(20), default="NORMAL")  # NotificationPriorityEnum
    
    # 消息状态
    is_read = Column(Boolean, default=False, index=True)
    read_at = Column(DateTime, nullable=True)
    
    # 关联数据
    related_id = Column(String(100), nullable=True)  # 关联的资源ID（如job_id、order_id等）
    related_type = Column(String(50), nullable=True)  # 关联的资源类型（如job、order等）
    
    # 时间戳
    created_at = Column(DateTime, default=datetime.utcnow, index=True)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)
    
    # 索引
    __table_args__ = (
        Index('idx_user_is_read', 'user_id', 'is_read'),
        Index('idx_user_created', 'user_id', 'created_at'),
    )

    def mark_as_read(self):
        """标记为已读"""
        self.is_read = True
        self.read_at = datetime.utcnow()
        self.updated_at = datetime.utcnow()

    def mark_as_unread(self):
        """标记为未读"""
        self.is_read = False
        self.read_at = None
        self.updated_at = datetime.utcnow()

