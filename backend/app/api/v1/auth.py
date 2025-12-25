"""
Authentication API endpoints
"""
from typing import Dict, Any,Optional
from fastapi import APIRouter, Depends, status
from fastapi.security import OAuth2PasswordRequestForm
from sqlalchemy.orm import Session
from sqlalchemy import func
from datetime import timedelta, datetime, date
from app.database import get_db
from app.dependencies import get_current_user
from app.models.user import User, UserRole, UserType
from app.models.job import MDJob, JobStatus
from app.schemas.user import UserCreate, User as UserSchema, ChangePassword
from app.schemas.token import Token
from app.core.security import verify_password, get_password_hash, create_access_token
from app.config import settings
from app.core.logger import logger
# Import error handling
from app.core.errors import (
    ErrorCode,
    invalid_input_error,
    unauthorized_error,
    not_found_error,
    already_exists_error
)


router = APIRouter()


# Helper functions with type hints
def get_running_job_count(user_id: int, db: Session) -> int:
    """
    Get count of running jobs for a user
    
    Args:
        user_id: User ID
        db: Database session
        
    Returns:
        Number of running jobs
    """
    return db.query(func.count(MDJob.id)).filter(
        MDJob.user_id == user_id,
        MDJob.status == JobStatus.RUNNING
    ).scalar() or 0


def calculate_user_cpu_hours(user_id: int, db: Session) -> float:
    """
    Calculate total CPU hours used by a user
    
    Args:
        user_id: User ID
        db: Database session
        
    Returns:
        Total CPU hours consumed
    """
    from app.models.billing import QuotaTransaction
    total = db.query(func.sum(QuotaTransaction.amount)).filter(
        QuotaTransaction.user_id == user_id,
        QuotaTransaction.type == 'consume'
    ).scalar() or 0.0
    return abs(total)


def get_today_job_count(user_id: int, db: Session) -> int:
    """
    Get count of jobs created today by a user
    
    Args:
        user_id: User ID
        db: Database session
        
    Returns:
        Number of jobs created today
    """
    today = date.today()
    return db.query(func.count(MDJob.id)).filter(
        MDJob.user_id == user_id,
        func.date(MDJob.created_at) == today
    ).scalar() or 0



@router.post("/register", response_model=UserSchema, status_code=status.HTTP_201_CREATED)
def register(
    user_data: UserCreate,
    db: Session = Depends(get_db)
) -> User:
    """
    Register a new user with initial 100 free CPU hours
    
    Args:
        user_data: User registration data including email, username, password
        db: Database session
        
    Returns:
        Created user object with default quotas
        
    Raises:
        APIError 400: If email, username, or phone already exists
        APIError 400: If phone verification code is invalid
    """
    from app.services.sms import sms_service

    # Check if email exists
    if db.query(User).filter(User.email == user_data.email).first():
        raise already_exists_error(
            resource="邮箱",
            details={"email": user_data.email}
        )

    # Check if username exists
    if db.query(User).filter(User.username == user_data.username).first():
        raise already_exists_error(
            resource="用户名",
            details={"username": user_data.username}
        )

    # 检查手机号（如果提供）
    phone_verified: bool = False
    if user_data.phone:
        # 检查手机号是否已被注册
        if db.query(User).filter(User.phone == user_data.phone).first():
            raise already_exists_error(
                resource="手机号",
                details={"phone": user_data.phone}
            )

        # 验证手机验证码
        if not user_data.phone_code:
            raise invalid_input_error(
                message="请输入手机验证码",
                field="phone_code"
            )

        success, message = sms_service.verify_code(
            phone=user_data.phone,
            code=user_data.phone_code,
            purpose="register"
        )
        if not success:
            raise invalid_input_error(
                message=message,
                field="phone_code"
            )
        phone_verified = True

    # Create new user with default quotas
    hashed_password: str = get_password_hash(user_data.password)
    db_user = User(
        email=user_data.email,
        username=user_data.username,
        password_hash=hashed_password,
        role=UserRole.USER,  # 强制设置为普通用户，不接受前端传入的role
        phone=user_data.phone,
        phone_verified=phone_verified,
        user_type=user_data.user_type or UserType.STUDENT,
        organization=user_data.organization,
        department=user_data.department,
        real_name=user_data.real_name,
        # 默认赠送100核时
        balance_cpu_hours=100.0,
    )

    db.add(db_user)
    db.commit()
    db.refresh(db_user)

    logger.info(f"New user registered: {db_user.username}")
    return db_user




@router.post("/login", response_model=Token)
def login(
    form_data: OAuth2PasswordRequestForm = Depends(),
    db: Session = Depends(get_db)
) -> Dict[str, str]:
    """
    Login and get access token
    
    Args:
        form_data: Login form data (username/email and password)
        db: Database session
        
    Returns:
        Access token and token type
        
    Raises:
        APIError 401: If credentials are invalid
    """
    # Try to get user by username first, then by email
    user: Optional[User] = db.query(User).filter(User.username == form_data.username).first()

    # If not found by username, try email
    if not user:
        user = db.query(User).filter(User.email == form_data.username).first()

    # Verify user and password
    if not user or not verify_password(form_data.password, user.password_hash):
        raise unauthorized_error(
            message="用户名/邮箱或密码错误",
            details={"username": form_data.username}
        )

    # Create access token
    access_token_expires = timedelta(minutes=settings.ACCESS_TOKEN_EXPIRE_MINUTES)
    access_token: str = create_access_token(
        data={"sub": user.username},
        expires_delta=access_token_expires
    )

    logger.info(f"User logged in: {user.username}")
    return {"access_token": access_token, "token_type": "bearer"}


@router.get("/me")
def get_current_user_info(
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
) -> Dict[str, Any]:
    """
    Get current user information with computed fields
    
    Args:
        current_user: Current authenticated user
        db: Database session
        
    Returns:
        User data dictionary with quota statistics and job counts
    """
    # 计算任务统计
    total_jobs: int = db.query(func.count(MDJob.id)).filter(MDJob.user_id == current_user.id).scalar() or 0
    completed_jobs: int = db.query(func.count(MDJob.id)).filter(
        MDJob.user_id == current_user.id,
        MDJob.status == JobStatus.COMPLETED
    ).scalar() or 0
    running_jobs: int = get_running_job_count(current_user.id, db)
    failed_jobs: int = db.query(func.count(MDJob.id)).filter(
        MDJob.user_id == current_user.id,
        MDJob.status == JobStatus.FAILED
    ).scalar() or 0

    # 计算配额使用情况
    used_cpu_hours: float = calculate_user_cpu_hours(current_user.id, db)
    today_jobs: int = get_today_job_count(current_user.id, db)

    # 获取余额系统数据（统一的核时系统）
    balance: float = current_user.balance_cpu_hours or 0.0
    frozen: float = current_user.frozen_cpu_hours or 0.0
    free_granted: float = current_user.free_cpu_hours_granted or 100.0
    recharge_hours: float = current_user.recharge_cpu_hours or 0.0
    admin_granted: float = current_user.admin_granted_cpu_hours or 0.0

    # 计算总消耗核时（从交易记录统计）
    from app.models.billing import QuotaTransaction
    total_consumed: float = db.query(func.sum(QuotaTransaction.amount)).filter(
        QuotaTransaction.user_id == current_user.id,
        QuotaTransaction.type == 'consume'
    ).scalar() or 0.0
    total_consumed = abs(total_consumed)  # consume是负数

    # 对于管理员，如果没有消费记录，计算实际使用的核时（用于显示）
    if current_user.role.value == "ADMIN" and total_consumed == 0:
        # 计算管理员所有已完成任务的实际核时
        admin_cpu_hours: float = db.query(
            func.sum(MDJob.actual_cpu_hours + func.coalesce(MDJob.resp_cpu_hours, 0))
        ).filter(
            MDJob.user_id == current_user.id,
            MDJob.status == JobStatus.COMPLETED
        ).scalar() or 0.0
        total_consumed = admin_cpu_hours

    # 可用余额 = 账户余额 - 冻结
    available_balance: float = max(0, balance - frozen)

    # 总核时 = 初始赠送 + 充值 + 管理员赠送
    total_cpu_hours: float = free_granted + recharge_hours + admin_granted

    # 欠费 = 当 balance < 0 时的绝对值
    debt: float = abs(balance) if balance < 0 else 0.0

    return {
        # 基本信息
        "id": current_user.id,
        "username": current_user.username,
        "email": current_user.email,
        "role": current_user.role.value if current_user.role else "USER",
        "user_type": current_user.user_type.value if current_user.user_type else "STUDENT",
        "organization": current_user.organization,
        "department": current_user.department,
        "email_verified": getattr(current_user, 'email_verified', False),
        "is_active": current_user.is_active,
        "created_at": current_user.created_at.isoformat() if current_user.created_at else None,
        "updated_at": current_user.updated_at.isoformat() if current_user.updated_at else None,

        # 配额配置
        "daily_job_limit": current_user.daily_job_limit,
        "concurrent_job_limit": current_user.concurrent_job_limit,
        "storage_quota_gb": current_user.storage_quota_gb,
        "allowed_partitions": current_user.allowed_partitions,
        "allowed_modules": current_user.allowed_modules,  # 模块权限

        # 核时余额系统
        "balance_cpu_hours": round(balance, 2),           # 账户余额
        "frozen_cpu_hours": round(frozen, 2),             # 冻结核时（运行中的任务）
        "available_cpu_hours": round(available_balance, 2), # 可用余额
        "debt_cpu_hours": round(debt, 2),                 # 欠费核时

        # 核时统计
        "free_cpu_hours_granted": round(free_granted, 2),   # 初始赠送
        "total_recharged": round(recharge_hours, 2),       # 总充值
        "total_consumed": round(total_consumed, 2),         # 总消耗
        "total_cpu_hours": round(total_cpu_hours, 2),       # 总核时额度

        # 兼容旧字段（给前端使用）
        "used_cpu_hours": round(total_consumed, 2),
        "remaining_cpu_hours": round(available_balance, 2),

        # 贡献统计
        "public_data_count": current_user.public_data_count or 0,
        "contribution_points": round(current_user.contribution_points or 0.0, 1),
        "private_quota_used": current_user.private_quota_used or 0,
        "private_quota_limit": getattr(current_user, 'private_quota_limit', 0),

        # 使用情况
        "today_jobs": today_jobs,
        "running_jobs": running_jobs,

        # 任务统计
        "total_jobs": total_jobs,
        "completed_jobs": completed_jobs,
        "failed_jobs": failed_jobs,

        # 计算使用率（基于余额系统）
        "cpu_hours_usage_percent": round(min(100, total_consumed / total_cpu_hours * 100), 1) if total_cpu_hours > 0 else 0,
        "daily_jobs_usage_percent": round(min(100, today_jobs / current_user.daily_job_limit * 100), 1) if current_user.daily_job_limit > 0 else 0,
    }


@router.post("/change-password", status_code=status.HTTP_200_OK)
def change_password(
    password_data: ChangePassword,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
) -> Dict[str, str]:
    """
    Change user password
    
    Args:
        password_data: Old and new password
        current_user: Current authenticated user
        db: Database session
        
    Returns:
        Success message
        
    Raises:
        APIError 400: If old password is incorrect
    """
    # Verify old password
    if not verify_password(password_data.old_password, current_user.password_hash):
        raise invalid_input_error(
            message="原密码不正确",
            field="old_password"
        )

    # Update password
    current_user.password_hash = get_password_hash(password_data.new_password)
    db.commit()

    logger.info(f"Password changed for user: {current_user.username}")
    return {"message": "密码修改成功"}


@router.get("/me/profile")
def get_user_profile(
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    获取当前用户完整资料（包含配额使用情况和贡献统计）

    Returns:
        dict: 用户完整资料
    """
    # 计算任务统计
    total_jobs = db.query(func.count(MDJob.id)).filter(MDJob.user_id == current_user.id).scalar() or 0
    completed_jobs = db.query(func.count(MDJob.id)).filter(
        MDJob.user_id == current_user.id,
        MDJob.status == JobStatus.COMPLETED
    ).scalar() or 0
    running_jobs = get_running_job_count(current_user.id, db)
    failed_jobs = db.query(func.count(MDJob.id)).filter(
        MDJob.user_id == current_user.id,
        MDJob.status == JobStatus.FAILED
    ).scalar() or 0

    # 计算配额使用情况
    used_cpu_hours = calculate_user_cpu_hours(current_user.id, db)
    today_jobs = get_today_job_count(current_user.id, db)

    # 获取余额系统数据
    balance = current_user.balance_cpu_hours or 0.0
    frozen = current_user.frozen_cpu_hours or 0.0
    free_granted = current_user.free_cpu_hours_granted or 100.0
    recharge_hours = current_user.recharge_cpu_hours or 0.0
    admin_granted = current_user.admin_granted_cpu_hours or 0.0

    # 计算总消耗核时（从交易记录统计）
    from app.models.billing import QuotaTransaction
    total_consumed = db.query(func.sum(QuotaTransaction.amount)).filter(
        QuotaTransaction.user_id == current_user.id,
        QuotaTransaction.type == 'consume'
    ).scalar() or 0.0
    total_consumed = abs(total_consumed)  # consume是负数

    # 对于管理员，如果没有消费记录，计算实际使用的核时（用于显示）
    if current_user.role.value == "ADMIN" and total_consumed == 0:
        # 计算管理员所有已完成任务的实际核时
        admin_cpu_hours = db.query(
            func.sum(MDJob.actual_cpu_hours + func.coalesce(MDJob.resp_cpu_hours, 0))
        ).filter(
            MDJob.user_id == current_user.id,
            MDJob.status == JobStatus.COMPLETED
        ).scalar() or 0.0
        total_consumed = admin_cpu_hours

    # 可用余额 = 账户余额 - 冻结
    available_balance = max(0, balance - frozen)

    # 总核时 = 初始赠送 + 充值 + 管理员赠送
    total_cpu_hours = free_granted + recharge_hours + admin_granted

    # 欠费 = 当 balance < 0 时的绝对值
    debt = abs(balance) if balance < 0 else 0.0

    return {
        # 基本信息
        "id": current_user.id,
        "username": current_user.username,
        "email": current_user.email,
        "role": current_user.role.value if current_user.role else "USER",
        "user_type": current_user.user_type.value if current_user.user_type else "STUDENT",
        "organization": current_user.organization,
        "department": current_user.department,
        "email_verified": getattr(current_user, 'email_verified', False),
        "is_active": current_user.is_active,
        "created_at": current_user.created_at.isoformat() if current_user.created_at else None,
        "last_login_at": current_user.last_login_at.isoformat() if current_user.last_login_at else None,

        # 配额配置
        "daily_job_limit": current_user.daily_job_limit,
        "concurrent_job_limit": current_user.concurrent_job_limit,
        "storage_quota_gb": current_user.storage_quota_gb,
        "allowed_partitions": current_user.allowed_partitions,
        "allowed_modules": current_user.allowed_modules,  # 模块权限

        # 核时余额系统
        "balance_cpu_hours": round(balance, 2),           # 账户余额
        "frozen_cpu_hours": round(frozen, 2),             # 冻结核时（运行中的任务）
        "available_cpu_hours": round(available_balance, 2), # 可用余额
        "debt_cpu_hours": round(debt, 2),                 # 欠费核时

        # 核时统计
        "free_cpu_hours_granted": round(free_granted, 2),   # 初始赠送
        "total_recharged": round(recharge_hours, 2),       # 总充值
        "total_consumed": round(total_consumed, 2),         # 总消耗
        "total_cpu_hours": round(total_cpu_hours, 2),       # 总核时额度

        # 兼容旧字段（给前端使用）
        "used_cpu_hours": round(total_consumed, 2),
        "remaining_cpu_hours": round(available_balance, 2),

        # 贡献统计
        "public_data_count": current_user.public_data_count or 0,
        "contribution_points": round(current_user.contribution_points or 0.0, 1),
        "private_quota_used": current_user.private_quota_used or 0,
        "private_quota_limit": getattr(current_user, 'private_quota_limit', 0),

        # 使用情况
        "today_jobs": today_jobs,
        "running_jobs": running_jobs,

        # 任务统计
        "total_jobs": total_jobs,
        "completed_jobs": completed_jobs,
        "failed_jobs": failed_jobs,

        # 计算使用率（基于余额系统）
        "cpu_hours_usage_percent": round(min(100, total_consumed / total_cpu_hours * 100), 1) if total_cpu_hours > 0 else 0,
        "daily_jobs_usage_percent": round(min(100, today_jobs / current_user.daily_job_limit * 100), 1) if current_user.daily_job_limit > 0 else 0,
    }

