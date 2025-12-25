"""
User management API endpoints
"""
from fastapi import APIRouter, Depends, status
from sqlalchemy.orm import Session
from typing import List, Optional, Dict, Any
from datetime import datetime, date, timedelta
from app.database import get_db
from app.models.user import User, UserRole
from app.models.user_stats import UserUsageStats
from app.schemas.user import User as UserSchema, UserUpdate
from app.dependencies import get_current_active_user, get_current_admin_user
from app.core.security import get_password_hash
from app.core.logger import logger
from app.services.quota_service import QuotaService
# Import error handling
from app.core.errors import (
    not_found_error,
    forbidden_error,
    invalid_input_error
)

router = APIRouter()


@router.get("/", response_model=List[UserSchema])
def list_users(
    skip: int = 0,
    limit: int = 100,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_admin_user)
):
    """
    List all users (admin only)
    
    Args:
        skip: Number of records to skip
        limit: Maximum number of records to return
        db: Database session
        current_user: Current admin user
        
    Returns:
        List[User]: List of users
    """
    users = db.query(User).offset(skip).limit(limit).all()
    return users



@router.get("/{user_id}", response_model=UserSchema)
def get_user(
    user_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
) -> User:
    """
    Get user by ID
    
    Args:
        user_id: User ID
        db: Database session
        current_user: Current authenticated user
        
    Returns:
        User: User data
        
    Raises:
        APIError 403: If no permission to view user
        APIError 404: If user not found
    """
    # Users can only view their own profile unless they are admin
    if current_user.id != user_id and current_user.role != UserRole.ADMIN:
        raise forbidden_error("您无权查看其他用户的信息")
    
    user: Optional[User] = db.query(User).filter(User.id == user_id).first()
    if not user:
        raise not_found_error("用户", user_id)
    
    return user


@router.put("/{user_id}", response_model=UserSchema)
def update_user(
    user_id: int,
    user_update: UserUpdate,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
) -> User:
    """
    Update user
    
    Args:
        user_id: User ID
        user_update: User update data
        db: Database session
        current_user: Current authenticated user
        
    Returns:
        User: Updated user
        
    Raises:
        APIError 403: If no permission to update user
        APIError 404: If user not found
    """
    # Users can only update their own profile unless they are admin
    if current_user.id != user_id and current_user.role != UserRole.ADMIN:
        raise forbidden_error("您无权修改其他用户的信息")
    
    user: Optional[User] = db.query(User).filter(User.id == user_id).first()
    if not user:
        raise not_found_error("用户", user_id)
    
    # Update fields
    update_data = user_update.model_dump(exclude_unset=True)
    
    # Hash password if provided
    if "password" in update_data:
        update_data["password_hash"] = get_password_hash(update_data.pop("password"))
    
    # Only admin can change role
    if "role" in update_data and current_user.role != UserRole.ADMIN:
        del update_data["role"]
    
    for field, value in update_data.items():
        setattr(user, field, value)
    
    db.commit()
    db.refresh(user)
    
    logger.info(f"User updated: {user.username}")
    return user


@router.delete("/{user_id}", status_code=status.HTTP_204_NO_CONTENT)
def delete_user(
    user_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_admin_user)
):
    """
    Delete user (admin only)
    
    Args:
        user_id: User ID
        db: Database session
        current_user: Current admin user
        
    Raises:
        HTTPException: If user not found
    """
    user = db.query(User).filter(User.id == user_id).first()
    if not user:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="User not found"
        )
    
    db.delete(user)
    db.commit()

    logger.info(f"User deleted: {user.username}")
    return None


@router.get("/me/daily-stats")
def get_my_daily_stats(
    days: int = 7,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    Get current user's daily usage statistics

    Args:
        days: Number of days to retrieve (default: 7)
        db: Database session
        current_user: Current user

    Returns:
        List of daily statistics
    """
    end_date = date.today()
    start_date = end_date - timedelta(days=days-1)

    stats = db.query(UserUsageStats).filter(
        UserUsageStats.user_id == current_user.id,
        UserUsageStats.date >= start_date,
        UserUsageStats.date <= end_date
    ).order_by(UserUsageStats.date.desc()).all()

    return [
        {
            "date": s.date,
            "jobs_submitted": s.jobs_submitted,
            "jobs_completed": s.jobs_completed,
            "jobs_failed": s.jobs_failed,
            "jobs_cancelled": s.jobs_cancelled,
            "cpu_hours_used": s.cpu_hours_used,
            "cluster_analysis_cpu_hours": s.cluster_analysis_cpu_hours,
            "cluster_analysis_task_count": s.cluster_analysis_task_count,
            "max_concurrent_jobs": s.max_concurrent_jobs,
            "storage_used_gb": s.storage_used_gb,
        }
        for s in stats
    ]


@router.get("/{user_id}/daily-stats")
def get_user_daily_stats(
    user_id: int,
    days: int = 7,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    Get user's daily usage statistics (admin can view any user, regular users can only view themselves)

    Args:
        user_id: User ID
        days: Number of days to retrieve (default: 7)
        db: Database session
        current_user: Current user

    Returns:
        List of daily statistics
    """
    # Check permission
    if user_id != current_user.id and current_user.role.value != 'ADMIN':
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="You can only view your own statistics"
        )

    # Check if user exists
    user = db.query(User).filter(User.id == user_id).first()
    if not user:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="User not found"
        )

    end_date = date.today()
    start_date = end_date - timedelta(days=days-1)

    stats = db.query(UserUsageStats).filter(
        UserUsageStats.user_id == user_id,
        UserUsageStats.date >= start_date,
        UserUsageStats.date <= end_date
    ).order_by(UserUsageStats.date.desc()).all()

    return [
        {
            "date": s.date,
            "jobs_submitted": s.jobs_submitted,
            "jobs_completed": s.jobs_completed,
            "jobs_failed": s.jobs_failed,
            "jobs_cancelled": s.jobs_cancelled,
            "cpu_hours_used": s.cpu_hours_used,
            "cluster_analysis_cpu_hours": s.cluster_analysis_cpu_hours,
            "cluster_analysis_task_count": s.cluster_analysis_task_count,
            "max_concurrent_jobs": s.max_concurrent_jobs,
            "storage_used_gb": s.storage_used_gb,
        }
        for s in stats
    ]


# ============================================================================
# Phase 2: Scheme B 重构 - 新增端点
# ============================================================================

@router.get("/me/account-info")
def get_account_info(
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db)
):
    """
    获取当前用户的账号信息

    包括：
    - 账号类型（个人/主账号/子账号）
    - 账号详细信息
    - 用户权限

    Args:
        current_user: 当前认证用户
        db: 数据库会话

    Returns:
        dict: 账号信息和权限
    """
    try:
        # 基本账号信息
        # 计算总核时（核时来源之和）
        total_cpu_hours = (
            current_user.free_cpu_hours_granted +
            current_user.recharge_cpu_hours +
            current_user.admin_granted_cpu_hours
        )

        account_info = {
            "user_id": current_user.id,
            "username": current_user.username,
            "email": current_user.email,
            "organization": current_user.organization,
            "department": current_user.department,
            "user_type": current_user.user_type.value if current_user.user_type else None,
            "balance_cpu_hours": current_user.balance_cpu_hours,
            "frozen_cpu_hours": current_user.frozen_cpu_hours,
            "total_cpu_hours": total_cpu_hours,
        }

        # 权限信息
        permissions = {
            "can_create_jobs": current_user.is_active,
            "can_manage_master_account": current_user.account_type == "master_account",
            "can_manage_sub_accounts": current_user.account_type == "master_account",
            "can_access_admin": current_user.role.value == "ADMIN",
        }

        return {
            "account_type": current_user.account_type,
            "account_info": account_info,
            "permissions": permissions
        }
    except Exception as e:
        logger.error(f"Error getting account info for user {current_user.id}: {str(e)}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to get account info"
        )


@router.get("/me", response_model=schemas.UserProfile)
def get_current_user_profile(
    current_user: User = Depends(get_current_active_user)
):
    """Get current user's profile"""
    return current_user


@router.get("/me/quota-status")
def get_quota_status(
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db)
):
    """
    获取用户配额状态(用于前端显示和权限判断)
    
    Returns:
        配额状态信息,包括余额、冻结、可用、欠费等
    """
    from app.services.quota_service import QuotaService
    
    # 获取可用配额
    available_quota = QuotaService.get_available_quota(current_user, db)
    
    # 判断是否欠费
    has_debt = current_user.balance_cpu_hours < 0
    debt_amount = abs(current_user.balance_cpu_hours) if has_debt else 0
    
    # 计算实际余额(正数部分)
    balance = max(0, current_user.balance_cpu_hours)
    
    # 低余额阈值
    warning_threshold = 10.0
    is_low_balance = available_quota < warning_threshold
    
    # 是否可以提交任务和查看结果
    can_submit_jobs = available_quota >= 1.0 and not has_debt
    can_view_results = not has_debt
    
    return {
        'balance': balance,
        'frozen': current_user.frozen_cpu_hours,
        'available': available_quota,
        'has_debt': has_debt,
        'debt_amount': debt_amount,
        'can_submit_jobs': can_submit_jobs,
        'can_view_results': can_view_results,
        'warning_threshold': warning_threshold,
        'is_low_balance': is_low_balance,
        'account_type': current_user.account_type
    }


@router.get("/me/quota")
def get_user_quota(
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db)
):
    """
    获取当前用户的配额信息（统一的核时系统）

    核时系统说明：
    - balance_cpu_hours 是唯一的真实来源
    - 正数 = 可用核时
    - 负数 = 欠费核时
    - 0 = 无可用核时

    核时来源追踪（仅用于统计）：
    - free_cpu_hours_granted: 初始赠送
    - recharge_cpu_hours: 充值获得
    - admin_granted_cpu_hours: 管理员赠送
    - 总核时 = 上述三项之和 - 已消费

    Args:
        current_user: 当前认证用户
        db: 数据库会话

    Returns:
        dict: 配额信息
    """
    try:
        from app.models.billing import QuotaTransaction
        from sqlalchemy import func

        # 计算已消费的核时（从 QuotaTransaction 表统计）
        total_consumed = db.query(func.sum(QuotaTransaction.amount)).filter(
            QuotaTransaction.user_id == current_user.id,
            QuotaTransaction.type == 'consume'
        ).scalar() or 0.0
        used_cpu_hours = abs(total_consumed)  # consume 是负数

        # 计算总核时（核时来源之和）
        total_cpu_hours = (
            current_user.free_cpu_hours_granted +
            current_user.recharge_cpu_hours +
            current_user.admin_granted_cpu_hours
        )

        return {
            # 核时系统（唯一真实来源）
            "balance_cpu_hours": current_user.balance_cpu_hours,  # 可用余额（正数）或欠费（负数）
            "frozen_cpu_hours": current_user.frozen_cpu_hours,    # 冻结核时

            # 核时来源追踪（用于统计和分析）
            "quota_sources": {
                "free_granted": current_user.free_cpu_hours_granted,
                "recharge": current_user.recharge_cpu_hours,
                "admin_granted": current_user.admin_granted_cpu_hours,
            },

            # 使用统计
            "used_cpu_hours": used_cpu_hours,      # 已消费核时
            "total_cpu_hours": total_cpu_hours,    # 总核时（来源之和）

            # 账户信息
            "account_type": current_user.account_type,

            # 兼容字段（保留用于前端兼容）
            "available_quota": current_user.balance_cpu_hours,
            "account_details": {
                "balance_cpu_hours": current_user.balance_cpu_hours,
                "frozen_cpu_hours": current_user.frozen_cpu_hours,
                "used_cpu_hours": used_cpu_hours,
                "total_cpu_hours": total_cpu_hours,
                "free_granted": current_user.free_cpu_hours_granted,
                "recharge": current_user.recharge_cpu_hours,
                "admin_granted": current_user.admin_granted_cpu_hours,
            }
        }
    except Exception as e:
        logger.error(f"Error getting quota for user {current_user.id}: {str(e)}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to get quota info"
        )

