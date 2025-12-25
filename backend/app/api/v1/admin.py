"""
Admin API endpoints
"""
from typing import List, Optional
from datetime import datetime, date, timedelta
from fastapi import APIRouter, Depends, HTTPException, status, Request
from sqlalchemy.orm import Session
from sqlalchemy import func, desc
import logging

from app.database import get_db
from app.dependencies import get_current_admin_user
from app.models.user import User, UserRole
from app.models.job import MDJob, JobStatus
from app.models.user_stats import UserUsageStats, AuditLog
from app.schemas.admin import (
    UserListItem, UserDetail, UserUpdate, UserCreate,
    GlobalStats, UserUsageStatsItem, UserRanking,
    TrendDataPoint, StatisticsSummary, AuditLogItem,
    QuotaCheckResponse
)
from app.utils.audit import (
    log_user_update, log_quota_update, log_user_status_change,
    log_job_cancel, create_audit_log
)
from app.core.security import get_password_hash
from app.services.slurm import list_partitions
from app.api.v1.auth import calculate_user_cpu_hours, get_running_job_count, get_today_job_count

logger = logging.getLogger(__name__)
router = APIRouter(prefix="/admin", tags=["admin"])


# ============ User Management Endpoints ============

@router.get("/users", response_model=List[UserListItem])
async def get_all_users(
    skip: int = 0,
    limit: int = 100,
    role: Optional[UserRole] = None,
    user_type: Optional[str] = None,  # UserType枚举值
    is_active: Optional[bool] = None,
    search: Optional[str] = None,  # 搜索用户名、邮箱、单位
    sort_by: str = "created_at",  # created_at, username, email, last_login_at
    sort_order: str = "desc",  # asc, desc
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """
    Get all users with enhanced filtering (admin only)

    Args:
        skip: Number of records to skip
        limit: Maximum number of records to return
        role: Filter by role (ADMIN, PREMIUM, USER, GUEST)
        user_type: Filter by user type (STUDENT, RESEARCHER, COMPANY)
        is_active: Filter by active status
        search: Search in username, email, organization
        sort_by: Sort field (created_at, username, email, last_login_at)
        sort_order: Sort order (asc, desc)
        db: Database session
        admin: Current admin user

    Returns:
        List of users
    """
    # Admin permission is already checked by get_current_admin_user dependency

    query = db.query(User)

    # 角色筛选
    if role:
        query = query.filter(User.role == role)

    # 用户类型筛选
    if user_type:
        query = query.filter(User.user_type == user_type)

    # 激活状态筛选
    if is_active is not None:
        query = query.filter(User.is_active == is_active)

    # 搜索功能（用户名、邮箱、单位）
    if search:
        search_pattern = f"%{search}%"
        query = query.filter(
            (User.username.ilike(search_pattern)) |
            (User.email.ilike(search_pattern)) |
            (User.organization.ilike(search_pattern))
        )

    # 排序
    sort_column = {
        "created_at": User.created_at,
        "username": User.username,
        "email": User.email,
        "last_login_at": User.last_login_at,
    }.get(sort_by, User.created_at)

    if sort_order == "asc":
        query = query.order_by(sort_column.asc())
    else:
        query = query.order_by(sort_column.desc())

    users = query.offset(skip).limit(limit).all()
    return users


@router.get("/users/count/total")
async def get_users_count(
    role: Optional[UserRole] = None,
    user_type: Optional[str] = None,
    is_active: Optional[bool] = None,
    search: Optional[str] = None,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """
    Get total count of users with filters (admin only)

    Args:
        role: Filter by role
        user_type: Filter by user type
        is_active: Filter by active status
        search: Search in username, email, organization
        db: Database session
        admin: Current admin user

    Returns:
        Total count of users matching filters
    """
    from app.models.user import UserType

    query = db.query(User)

    if role:
        query = query.filter(User.role == role)

    if user_type:
        query = query.filter(User.user_type == user_type)

    if is_active is not None:
        query = query.filter(User.is_active == is_active)

    if search:
        search_pattern = f"%{search}%"
        query = query.filter(
            (User.username.ilike(search_pattern)) |
            (User.email.ilike(search_pattern)) |
            (User.organization.ilike(search_pattern))
        )

    total = query.count()
    return {"total": total}


@router.get("/users/{user_id}", response_model=UserDetail)
async def get_user_detail(
    user_id: int,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """
    Get detailed user information (admin only)
    
    Args:
        user_id: User ID
        db: Database session
        admin: Current admin user
        
    Returns:
        Detailed user information
    """
    user = db.query(User).filter(User.id == user_id).first()
    if not user:
        raise HTTPException(status_code=404, detail="User not found")
    
    # Calculate usage statistics
    used_cpu_hours = calculate_user_cpu_hours(user_id, db)
    today_jobs = get_today_job_count(user_id, db)
    running_jobs = get_running_job_count(user_id, db)
    
    total_jobs = db.query(MDJob).filter(MDJob.user_id == user_id).count()
    completed_jobs = db.query(MDJob).filter(
        MDJob.user_id == user_id,
        MDJob.status == JobStatus.COMPLETED
    ).count()
    failed_jobs = db.query(MDJob).filter(
        MDJob.user_id == user_id,
        MDJob.status == JobStatus.FAILED
    ).count()
    
    # Convert to dict and add statistics
    user_dict = {
        "id": user.id,
        "username": user.username,
        "email": user.email,
        "role": user.role,
        "user_type": user.user_type,
        "organization": user.organization,
        "department": user.department,
        "is_active": user.is_active,
        "balance_cpu_hours": user.balance_cpu_hours,
        "frozen_cpu_hours": user.frozen_cpu_hours,
        "free_cpu_hours_granted": user.free_cpu_hours_granted,
        "recharge_cpu_hours": user.recharge_cpu_hours,
        "admin_granted_cpu_hours": user.admin_granted_cpu_hours,
        "daily_job_limit": user.daily_job_limit,
        "concurrent_job_limit": user.concurrent_job_limit,
        "storage_quota_gb": user.storage_quota_gb,
        "allowed_partitions": user.allowed_partitions,
        "allowed_modules": user.allowed_modules,
        "custom_cpu_hour_price": user.custom_cpu_hour_price,
        "last_login_at": user.last_login_at,
        "created_at": user.created_at,
        "updated_at": user.updated_at,
        "used_cpu_hours": used_cpu_hours,
        "today_jobs": today_jobs,
        "running_jobs": running_jobs,
        "total_jobs": total_jobs,
        "completed_jobs": completed_jobs,
        "failed_jobs": failed_jobs,
    }

    return UserDetail(**user_dict)


@router.put("/users/{user_id}", response_model=UserDetail)
async def update_user(
    user_id: int,
    user_update: UserUpdate,
    request: Request,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """
    Update user information and quotas (admin only)
    
    Args:
        user_id: User ID
        user_update: User update data
        request: FastAPI request
        db: Database session
        admin: Current admin user
        
    Returns:
        Updated user information
    """
    user = db.query(User).filter(User.id == user_id).first()
    if not user:
        raise HTTPException(status_code=404, detail="User not found")
    
    # Track changes for audit log
    changes = {}
    old_quotas = {}
    new_quotas = {}
    
    # Update fields
    update_data = user_update.model_dump(exclude_unset=True)
    for field, value in update_data.items():
        if hasattr(user, field):
            old_value = getattr(user, field)
            if old_value != value:
                changes[field] = {"old": str(old_value), "new": str(value)}
                setattr(user, field, value)

                # Track quota changes
                if field in ["balance_cpu_hours", "free_cpu_hours_granted", "recharge_cpu_hours",
                            "admin_granted_cpu_hours", "daily_job_limit", "concurrent_job_limit", "storage_quota_gb"]:
                    old_quotas[field] = old_value
                    new_quotas[field] = value
    
    db.commit()
    db.refresh(user)
    
    # Log the update
    if changes:
        log_user_update(db, admin, user, changes, request)
        
        if old_quotas:
            log_quota_update(db, admin, user, old_quotas, new_quotas, request)
    
    # Return detailed user info
    return await get_user_detail(user_id, db, admin)


@router.post("/users", response_model=UserDetail, status_code=status.HTTP_201_CREATED)
async def create_user(
    user_create: UserCreate,
    request: Request,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """
    Create a new user (admin only)

    Args:
        user_create: User creation data
        request: FastAPI request
        db: Database session
        admin: Current admin user

    Returns:
        Created user information
    """
    # Check if username already exists
    existing_user = db.query(User).filter(User.username == user_create.username).first()
    if existing_user:
        raise HTTPException(status_code=400, detail="Username already exists")

    # Check if email already exists
    existing_email = db.query(User).filter(User.email == user_create.email).first()
    if existing_email:
        raise HTTPException(status_code=400, detail="Email already exists")

    # Create new user
    new_user = User(
        username=user_create.username,
        email=user_create.email,
        password_hash=get_password_hash(user_create.password),
        role=user_create.role,
        is_active=True,
        balance_cpu_hours=user_create.balance_cpu_hours,
        free_cpu_hours_granted=user_create.free_cpu_hours_granted,
        daily_job_limit=user_create.daily_job_limit,
        concurrent_job_limit=user_create.concurrent_job_limit,
        storage_quota_gb=user_create.storage_quota_gb,
        allowed_partitions=user_create.allowed_partitions
    )

    db.add(new_user)
    db.commit()
    db.refresh(new_user)

    # Log the creation
    create_audit_log(
        db=db,
        user=admin,
        action="create_user",
        resource_type="user",
        resource_id=new_user.id,
        details={
            "username": new_user.username,
            "email": new_user.email,
            "role": new_user.role.value
        },
        request=request
    )

    return await get_user_detail(new_user.id, db, admin)


@router.delete("/users/{user_id}", status_code=status.HTTP_204_NO_CONTENT)
async def delete_user(
    user_id: int,
    request: Request,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """
    Delete a user (admin only)

    Args:
        user_id: User ID
        request: FastAPI request
        db: Database session
        admin: Current admin user
    """
    user = db.query(User).filter(User.id == user_id).first()
    if not user:
        raise HTTPException(status_code=404, detail="User not found")

    # Prevent deleting yourself
    if user.id == admin.id:
        raise HTTPException(status_code=400, detail="Cannot delete yourself")

    # Log the deletion
    create_audit_log(
        db=db,
        user=admin,
        action="delete_user",
        resource_type="user",
        resource_id=user.id,
        details={
            "username": user.username,
            "email": user.email
        },
        request=request
    )

    db.delete(user)
    db.commit()


@router.put("/users/{user_id}/status")
async def update_user_status(
    user_id: int,
    is_active: bool,
    request: Request,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """
    Enable or disable a user (admin only)

    Args:
        user_id: User ID
        is_active: New active status
        request: FastAPI request
        db: Database session
        admin: Current admin user

    Returns:
        Success message
    """
    user = db.query(User).filter(User.id == user_id).first()
    if not user:
        raise HTTPException(status_code=404, detail="User not found")

    # Prevent disabling yourself
    if user.id == admin.id and not is_active:
        raise HTTPException(status_code=400, detail="Cannot disable yourself")

    old_status = user.is_active
    user.is_active = is_active
    db.commit()

    # Log the status change
    log_user_status_change(db, admin, user, old_status, is_active, request)

    return {
        "message": f"User {'enabled' if is_active else 'disabled'} successfully",
        "user_id": user_id,
        "is_active": is_active
    }


@router.get("/users/{user_id}/quota", response_model=QuotaCheckResponse)
async def check_user_quota_endpoint(
    user_id: int,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """
    Check user's quota status (admin only)

    Args:
        user_id: User ID
        db: Database session
        admin: Current admin user

    Returns:
        Quota check result
    """
    user = db.query(User).filter(User.id == user_id).first()
    if not user:
        raise HTTPException(status_code=404, detail="User not found")

    result = check_user_quota(user, db)
    return QuotaCheckResponse(**result)


# ============ Statistics and Monitoring Endpoints ============

@router.get("/stats/overview", response_model=GlobalStats)
async def get_global_stats(
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """
    Get global system statistics (admin only)

    Args:
        db: Database session
        admin: Current admin user

    Returns:
        Global statistics
    """
    # User statistics
    total_users = db.query(User).count()
    active_users = db.query(User).filter(User.is_active == True).count()

    # Job statistics
    total_jobs = db.query(MDJob).count()
    running_jobs = db.query(MDJob).filter(MDJob.status == JobStatus.RUNNING).count()
    queued_jobs = db.query(MDJob).filter(MDJob.status == JobStatus.QUEUED).count()
    completed_jobs = db.query(MDJob).filter(MDJob.status == JobStatus.COMPLETED).count()
    failed_jobs = db.query(MDJob).filter(MDJob.status == JobStatus.FAILED).count()

    # Resource statistics
    # 计算所有用户的总配额（来自所有来源）
    total_cpu_hours_allocated = db.query(
        func.sum(User.free_cpu_hours_granted + User.recharge_cpu_hours + User.admin_granted_cpu_hours)
    ).scalar() or 0.0
    total_storage_allocated_gb = db.query(func.sum(User.storage_quota_gb)).scalar() or 0.0

    # Calculate total CPU hours used - 优化版本：使用交易记录统计
    from app.models.billing import QuotaTransaction
    total_consumed_from_transactions = db.query(func.sum(QuotaTransaction.amount)).filter(
        QuotaTransaction.type == 'consume'
    ).scalar() or 0.0
    total_cpu_hours_used = abs(total_consumed_from_transactions)  # consume是负数

    # 如果没有交易记录（旧系统），回退到慢速计算（但只计算前100个用户）
    if total_cpu_hours_used == 0:
        logger.warning("No billing transactions found, using slow calculation method (limited to 100 users)")
        sample_users = db.query(User).limit(100).all()
        # 直接计算而不是调用外部函数
        for user in sample_users:
            user_consumed = db.query(func.sum(QuotaTransaction.amount)).filter(
                QuotaTransaction.user_id == user.id,
                QuotaTransaction.type == 'consume'
            ).scalar() or 0.0
            total_cpu_hours_used += abs(user_consumed)

    # Storage used (placeholder - would need actual implementation)
    total_storage_used_gb = 0.0

    return GlobalStats(
        total_users=total_users,
        active_users=active_users,
        total_jobs=total_jobs,
        running_jobs=running_jobs,
        queued_jobs=queued_jobs,
        completed_jobs=completed_jobs,
        failed_jobs=failed_jobs,
        total_cpu_hours_used=total_cpu_hours_used,
        total_cpu_hours_allocated=total_cpu_hours_allocated,
        total_storage_used_gb=total_storage_used_gb,
        total_storage_allocated_gb=total_storage_allocated_gb
    )


@router.get("/stats/users", response_model=List[UserUsageStatsItem])
async def get_user_usage_stats(
    skip: int = 0,
    limit: int = 100,
    sort_by: str = "used_cpu_hours",  # used_cpu_hours, total_jobs, usage_percentage
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """
    Get usage statistics for all users (admin only)

    Args:
        skip: Number of records to skip
        limit: Maximum number of records to return
        sort_by: Sort field
        db: Database session
        admin: Current admin user

    Returns:
        List of user usage statistics
    """
    try:
        users = db.query(User).offset(skip).limit(limit).all()

        stats_list = []
        for user in users:
            try:
                used_cpu_hours = calculate_user_cpu_hours(user.id, db)
            except Exception as e:
                logger.error(f"Failed to calculate CPU hours for user {user.id}: {e}")
                used_cpu_hours = 0.0

            try:
                total_jobs = db.query(MDJob).filter(MDJob.user_id == user.id).count()
                running_jobs = get_running_job_count(user.id, db)
                completed_jobs = db.query(MDJob).filter(
                    MDJob.user_id == user.id,
                    MDJob.status == JobStatus.COMPLETED
                ).count()
                failed_jobs = db.query(MDJob).filter(
                    MDJob.user_id == user.id,
                    MDJob.status == JobStatus.FAILED
                ).count()

                # Get last job time
                last_job = db.query(MDJob).filter(
                    MDJob.user_id == user.id
                ).order_by(desc(MDJob.created_at)).first()
                last_job_at = last_job.created_at if last_job else None
            except Exception as e:
                logger.error(f"Failed to get job stats for user {user.id}: {e}")
                total_jobs = running_jobs = completed_jobs = failed_jobs = 0
                last_job_at = None

            # Calculate total CPU hours from sources
            total_cpu_hours = (user.free_cpu_hours_granted or 0) + (user.recharge_cpu_hours or 0) + (user.admin_granted_cpu_hours or 0)
            usage_percentage = (used_cpu_hours / total_cpu_hours * 100) if total_cpu_hours > 0 else 0

            stats_list.append(UserUsageStatsItem(
                user_id=user.id,
                username=user.username,
                email=user.email,
                role=user.role,
                used_cpu_hours=used_cpu_hours,
                total_cpu_hours=total_cpu_hours,
                usage_percentage=usage_percentage,
                total_jobs=total_jobs,
                running_jobs=running_jobs,
                completed_jobs=completed_jobs,
                failed_jobs=failed_jobs,
                last_job_at=last_job_at
            ))

        # Sort the results
        if sort_by == "used_cpu_hours":
            stats_list.sort(key=lambda x: x.used_cpu_hours, reverse=True)
        elif sort_by == "total_jobs":
            stats_list.sort(key=lambda x: x.total_jobs, reverse=True)
        elif sort_by == "usage_percentage":
            stats_list.sort(key=lambda x: x.usage_percentage, reverse=True)

        return stats_list
    except Exception as e:
        logger.error(f"Failed to get user usage stats: {e}")
        raise HTTPException(status_code=500, detail=f"Failed to get user usage stats: {str(e)}")


@router.get("/stats/ranking/cpu", response_model=List[UserRanking])
async def get_cpu_usage_ranking(
    limit: int = 10,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """
    Get top users by CPU hours usage (admin only)

    Args:
        limit: Number of top users to return
        db: Database session
        admin: Current admin user

    Returns:
        List of top users by CPU usage
    """
    try:
        from app.models.billing import QuotaTransaction

        users = db.query(User).all()

        user_cpu_hours = []
        for user in users:
            try:
                # 优先从交易记录计算总消耗
                total_consumed = db.query(func.sum(QuotaTransaction.amount)).filter(
                    QuotaTransaction.user_id == user.id,
                    QuotaTransaction.type == 'consume'
                ).scalar() or 0.0
                total_consumed = abs(total_consumed)  # consume是负数

                # 如果没有交易记录，使用旧的计算方法（这种情况很少见）
                # 直接返回0，因为没有交易记录就说明没有消费
                if total_consumed == 0:
                    total_consumed = 0.0

                user_cpu_hours.append({
                    "user_id": user.id,
                    "username": user.username,
                    "email": user.email,
                    "cpu_hours": total_consumed
                })
            except Exception as e:
                logger.error(f"Failed to calculate CPU hours for user {user.id}: {e}")
                user_cpu_hours.append({
                    "user_id": user.id,
                    "username": user.username,
                    "email": user.email,
                    "cpu_hours": 0.0
                })

        # Sort by CPU hours
        user_cpu_hours.sort(key=lambda x: x["cpu_hours"], reverse=True)

        # Create ranking
        ranking = []
        for rank, item in enumerate(user_cpu_hours[:limit], start=1):
            ranking.append(UserRanking(
                user_id=item["user_id"],
                username=item["username"],
                email=item["email"],
                metric_value=item["cpu_hours"],
                rank=rank
            ))

        return ranking
    except Exception as e:
        logger.error(f"Failed to get CPU usage ranking: {e}")
        raise HTTPException(status_code=500, detail=f"Failed to get CPU usage ranking: {str(e)}")


@router.get("/stats/ranking/jobs", response_model=List[UserRanking])
async def get_job_count_ranking(
    limit: int = 10,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """
    Get top users by job count (admin only)

    Args:
        limit: Number of top users to return
        db: Database session
        admin: Current admin user

    Returns:
        List of top users by job count
    """
    # Query job counts grouped by user
    job_counts = db.query(
        MDJob.user_id,
        func.count(MDJob.id).label("job_count")
    ).group_by(MDJob.user_id).order_by(desc("job_count")).limit(limit).all()

    ranking = []
    for rank, (user_id, job_count) in enumerate(job_counts, start=1):
        user = db.query(User).filter(User.id == user_id).first()
        if user:
            ranking.append(UserRanking(
                user_id=user.id,
                username=user.username,
                email=user.email,
                metric_value=float(job_count),
                rank=rank
            ))

    return ranking


# ============ Job Management Endpoints ============

@router.get("/jobs/all")
async def get_all_jobs(
    skip: int = 0,
    limit: int = 100,
    status: Optional[JobStatus] = None,
    user_id: Optional[int] = None,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """
    Get all jobs from all users (admin only)

    Args:
        skip: Number of records to skip
        limit: Maximum number of records to return
        status: Filter by job status
        user_id: Filter by user ID
        db: Database session
        admin: Current admin user

    Returns:
        List of jobs
    """
    query = db.query(MDJob)

    if status:
        query = query.filter(MDJob.status == status)
    if user_id:
        query = query.filter(MDJob.user_id == user_id)

    jobs = query.order_by(desc(MDJob.created_at)).offset(skip).limit(limit).all()

    # Convert to dict and add user info
    result = []
    for job in jobs:
        job_dict = {
            "id": job.id,
            "user_id": job.user_id,
            "system_id": job.system_id,
            "status": job.status,
            "slurm_job_id": job.slurm_job_id,
            "created_at": job.created_at,
            "started_at": job.started_at,
            "finished_at": job.finished_at,
            "config": job.config
        }

        # Add user info
        user = db.query(User).filter(User.id == job.user_id).first()
        if user:
            job_dict["username"] = user.username
            job_dict["user_email"] = user.email

        result.append(job_dict)

    return result


@router.post("/jobs/{job_id}/cancel")
async def admin_cancel_job(
    job_id: int,
    reason: Optional[str] = None,
    request: Request = None,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """
    Cancel a job (admin only)

    Args:
        job_id: Job ID
        reason: Reason for cancellation
        request: FastAPI request
        db: Database session
        admin: Current admin user

    Returns:
        Success message
    """
    job = db.query(MDJob).filter(MDJob.id == job_id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    # Check if job can be cancelled
    if job.status not in [JobStatus.QUEUED, JobStatus.RUNNING, JobStatus.POSTPROCESSING]:
        raise HTTPException(
            status_code=400,
            detail=f"Cannot cancel job with status {job.status}"
        )

    # Update job status
    job.status = JobStatus.CANCELLED
    job.finished_at = datetime.now()
    db.commit()

    # Log the cancellation
    log_job_cancel(db, admin, job_id, reason, request)

    return {
        "message": "Job cancelled successfully",
        "job_id": job_id,
        "reason": reason
    }


# ============ Audit Log Endpoints ============

@router.get("/logs", response_model=List[AuditLogItem])
async def get_audit_logs(
    skip: int = 0,
    limit: int = 100,
    action: Optional[str] = None,
    user_id: Optional[int] = None,
    resource_type: Optional[str] = None,
    start_date: Optional[datetime] = None,
    end_date: Optional[datetime] = None,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """
    Get audit logs (admin only)

    Args:
        skip: Number of records to skip
        limit: Maximum number of records to return
        action: Filter by action
        user_id: Filter by user ID
        resource_type: Filter by resource type
        start_date: Filter by start date
        end_date: Filter by end date
        db: Database session
        admin: Current admin user

    Returns:
        List of audit logs
    """
    query = db.query(AuditLog)

    if action:
        query = query.filter(AuditLog.action == action)
    if user_id:
        query = query.filter(AuditLog.user_id == user_id)
    if resource_type:
        query = query.filter(AuditLog.resource_type == resource_type)
    if start_date:
        query = query.filter(AuditLog.created_at >= start_date)
    if end_date:
        query = query.filter(AuditLog.created_at <= end_date)

    logs = query.order_by(desc(AuditLog.created_at)).offset(skip).limit(limit).all()

    # Add username to each log
    result = []
    for log in logs:
        log_dict = {
            "id": log.id,
            "user_id": log.user_id,
            "action": log.action,
            "resource_type": log.resource_type,
            "resource_id": log.resource_id,
            "details": log.details,
            "ip_address": log.ip_address,
            "created_at": log.created_at,
            "username": None
        }

        if log.user_id:
            user = db.query(User).filter(User.id == log.user_id).first()
            if user:
                log_dict["username"] = user.username

        result.append(AuditLogItem(**log_dict))

    return result


# ============ Statistics Report Endpoints ============

@router.get("/statistics/summary", response_model=StatisticsSummary)
async def get_statistics_summary(
    period: str = "7days",  # today, 7days, 30days
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """
    Get statistics summary for a period (admin only)

    Args:
        period: Time period (today, 7days, 30days)
        db: Database session
        admin: Current admin user

    Returns:
        Statistics summary
    """
    # Calculate date range
    end_date = datetime.now()
    if period == "today":
        start_date = datetime.now().replace(hour=0, minute=0, second=0, microsecond=0)
    elif period == "7days":
        start_date = end_date - timedelta(days=7)
    elif period == "30days":
        start_date = end_date - timedelta(days=30)
    else:
        raise HTTPException(status_code=400, detail="Invalid period")

    # Query jobs in the period
    jobs = db.query(MDJob).filter(
        MDJob.created_at >= start_date,
        MDJob.created_at <= end_date
    ).all()

    jobs_submitted = len(jobs)
    jobs_completed = sum(1 for j in jobs if j.status == JobStatus.COMPLETED)
    jobs_failed = sum(1 for j in jobs if j.status == JobStatus.FAILED)

    # Calculate CPU hours and average duration
    total_cpu_hours = sum(calculate_user_cpu_hours(j.user_id, db) for j in jobs if j.status == JobStatus.COMPLETED)

    completed_durations = []
    for job in jobs:
        if job.status == JobStatus.COMPLETED and job.started_at and job.finished_at:
            duration = (job.finished_at - job.started_at).total_seconds() / 3600.0
            completed_durations.append(duration)

    avg_duration = sum(completed_durations) / len(completed_durations) if completed_durations else 0.0

    # Peak concurrent jobs (simplified - would need time-series data for accurate calculation)
    peak_concurrent = db.query(MDJob).filter(
        MDJob.status.in_([JobStatus.RUNNING, JobStatus.QUEUED])
    ).count()

    return StatisticsSummary(
        period=period,
        jobs_submitted=jobs_submitted,
        jobs_completed=jobs_completed,
        jobs_failed=jobs_failed,
        cpu_hours_used=total_cpu_hours,
        avg_job_duration_hours=avg_duration,
        peak_concurrent_jobs=peak_concurrent
    )

# ============ Master Account Management ============

@router.get("/master-accounts", response_model=List[dict])
async def get_all_master_accounts(
    skip: int = 0,
    limit: int = 100,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """
    Get all master accounts (admin only)

    Args:
        skip: Number of records to skip
        limit: Maximum number of records to return
        db: Database session
        admin: Current admin user

    Returns:
        List of master accounts
    """
    try:
        from app.models.organization_v2 import MasterAccount

        master_accounts = db.query(MasterAccount).offset(skip).limit(limit).all()

        result = []
        for ma in master_accounts:
            user = db.query(User).filter(User.id == ma.user_id).first()
            if user:
                result.append({
                    "id": ma.id,
                    "user_id": ma.user_id,
                    "username": user.username,
                    "email": user.email,
                    "organization": user.organization,
                    "custom_cpu_hour_price": user.custom_cpu_hour_price,
                    # 主账号的配额来自 User 表
                    "balance_cpu_hours": user.balance_cpu_hours,
                    "frozen_cpu_hours": user.frozen_cpu_hours,
                    "free_cpu_hours_granted": user.free_cpu_hours_granted,
                    "recharge_cpu_hours": user.recharge_cpu_hours,
                    "admin_granted_cpu_hours": user.admin_granted_cpu_hours,
                    "current_sub_accounts": ma.current_sub_accounts,
                    "max_sub_accounts": ma.max_sub_accounts,
                    "is_active": ma.is_active,
                    "created_at": ma.created_at,
                    "updated_at": ma.updated_at,
                })

        return result
    except Exception as e:
        logger.error(f"Failed to get master accounts: {e}")
        raise HTTPException(status_code=500, detail=f"Failed to get master accounts: {str(e)}")


@router.post("/master-accounts", response_model=dict, status_code=status.HTTP_201_CREATED)
async def create_master_account(
    request: Request,
    user_id: Optional[int] = None,
    username: Optional[str] = None,
    email: Optional[str] = None,
    password: Optional[str] = None,
    organization: Optional[str] = None,
    max_sub_accounts: int = 10,
    initial_cpu_hours: float = 0.0,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """
    Create a new master account (admin only)

    支持两种方式：
    1. 提供 user_id：将已有用户升级为主账号
    2. 提供 username, email, password：创建新用户并设为主账号

    Args:
        user_id: Existing user ID to upgrade (方式1)
        username: Username for new account (方式2)
        email: Email for new account (方式2)
        password: Password for new account (方式2)
        organization: Organization name (optional)
        max_sub_accounts: Maximum number of sub-accounts allowed
        initial_cpu_hours: Initial CPU hours to grant (仅方式2有效)

    Returns:
        Created master account info
    """
    try:
        from app.models.organization_v2 import MasterAccount
        from app.models.user import AccountType

        # 方式1: 将已有用户升级为主账号
        if user_id:
            existing_user = db.query(User).filter(User.id == user_id).first()
            if not existing_user:
                raise HTTPException(status_code=404, detail="用户不存在")

            # 检查是否已经是主账号
            existing_master = db.query(MasterAccount).filter(MasterAccount.user_id == user_id).first()
            if existing_master:
                raise HTTPException(status_code=400, detail="该用户已经是主账号")

            # 检查是否是子账号
            if existing_user.account_type == AccountType.SUB_ACCOUNT.value:
                raise HTTPException(status_code=400, detail="子账号不能升级为主账号，请先解除子账号关系")

            # 升级为主账号
            existing_user.account_type = AccountType.MASTER_ACCOUNT.value
            if organization:
                existing_user.organization = organization

            # 创建主账号记录
            master_account = MasterAccount(
                user_id=existing_user.id,
                max_sub_accounts=max_sub_accounts,
                current_sub_accounts=0,
                is_active=True,
            )
            db.add(master_account)
            db.commit()

            # Create audit log
            create_audit_log(
                db=db,
                user=admin,
                action="upgrade_to_master_account",
                resource_type="master_account",
                resource_id=master_account.id,
                details={
                    "user_id": existing_user.id,
                    "username": existing_user.username,
                    "max_sub_accounts": max_sub_accounts,
                    "admin": admin.username
                },
                request=request
            )

            logger.info(f"User {existing_user.username} upgraded to master account by admin")

            return {
                "id": master_account.id,
                "user_id": existing_user.id,
                "username": existing_user.username,
                "email": existing_user.email,
                "organization": existing_user.organization,
                "balance_cpu_hours": existing_user.balance_cpu_hours,
                "frozen_cpu_hours": existing_user.frozen_cpu_hours,
                "free_cpu_hours_granted": existing_user.free_cpu_hours_granted,
                "max_sub_accounts": master_account.max_sub_accounts,
                "current_sub_accounts": master_account.current_sub_accounts,
                "is_active": master_account.is_active,
                "created_at": master_account.created_at,
                "updated_at": master_account.updated_at,
            }

        # 方式2: 创建新用户并设为主账号
        if not username or not email or not password:
            raise HTTPException(
                status_code=400,
                detail="请提供 user_id（升级已有用户）或提供 username, email, password（创建新用户）"
            )

        # Check if username already exists
        existing_user = db.query(User).filter(User.username == username).first()
        if existing_user:
            raise HTTPException(status_code=400, detail="用户名已存在")

        # Check if email already exists
        existing_email = db.query(User).filter(User.email == email).first()
        if existing_email:
            raise HTTPException(status_code=400, detail="邮箱已存在")

        # Create new user
        new_user = User(
            username=username,
            email=email,
            password_hash=get_password_hash(password),
            account_type=AccountType.MASTER_ACCOUNT.value,
            organization=organization or "",
            is_active=True,
            email_verified=True,
            balance_cpu_hours=initial_cpu_hours,
            free_cpu_hours_granted=initial_cpu_hours,
        )
        db.add(new_user)
        db.flush()

        # Create master account record
        master_account = MasterAccount(
            user_id=new_user.id,
            max_sub_accounts=max_sub_accounts,
            current_sub_accounts=0,
            is_active=True,
        )
        db.add(master_account)
        db.commit()

        # Create audit log
        create_audit_log(
            db=db,
            user=admin,
            action="create_master_account",
            resource_type="master_account",
            resource_id=master_account.id,
            details={
                "username": username,
                "email": email,
                "organization": organization,
                "max_sub_accounts": max_sub_accounts,
                "initial_cpu_hours": initial_cpu_hours,
                "admin": admin.username
            },
            request=request
        )

        logger.info(f"Master account created by admin: {username}")

        return {
            "id": master_account.id,
            "user_id": new_user.id,
            "username": new_user.username,
            "email": new_user.email,
            "organization": new_user.organization,
            "balance_cpu_hours": new_user.balance_cpu_hours,
            "frozen_cpu_hours": new_user.frozen_cpu_hours,
            "free_cpu_hours_granted": new_user.free_cpu_hours_granted,
            "max_sub_accounts": master_account.max_sub_accounts,
            "current_sub_accounts": master_account.current_sub_accounts,
            "is_active": master_account.is_active,
            "created_at": master_account.created_at,
            "updated_at": master_account.updated_at,
        }
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Failed to create master account: {e}")
        raise HTTPException(status_code=500, detail=f"创建主账号失败: {str(e)}")


@router.put("/master-accounts/{master_id}", response_model=dict)
async def update_master_account(
    master_id: int,
    request: Request,
    max_sub_accounts: Optional[int] = None,
    is_active: Optional[bool] = None,
    balance_cpu_hours: Optional[float] = None,
    custom_cpu_hour_price: Optional[float] = None,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """
    Update master account (admin only)

    注意：修改主账号的 balance_cpu_hours 只会修改 User 表中的余额，
    不会影响该用户作为个人账号时的余额（因为是同一个用户）。

    Args:
        master_id: Master account ID
        max_sub_accounts: New max sub-accounts limit
        is_active: New active status
        balance_cpu_hours: New balance CPU hours
        custom_cpu_hour_price: Custom CPU hour price for this user

    Returns:
        Updated master account info
    """
    try:
        from app.models.organization_v2 import MasterAccount

        ma = db.query(MasterAccount).filter(MasterAccount.id == master_id).first()
        if not ma:
            raise HTTPException(status_code=404, detail="主账号不存在")

        user = db.query(User).filter(User.id == ma.user_id).first()
        if not user:
            raise HTTPException(status_code=404, detail="用户不存在")

        changes = {}

        if max_sub_accounts is not None:
            changes["max_sub_accounts"] = f"{ma.max_sub_accounts} -> {max_sub_accounts}"
            ma.max_sub_accounts = max_sub_accounts

        if is_active is not None:
            changes["is_active"] = f"{ma.is_active} -> {is_active}"
            ma.is_active = is_active

        if balance_cpu_hours is not None:
            changes["balance_cpu_hours"] = f"{user.balance_cpu_hours} -> {balance_cpu_hours}"
            # 计算差额，更新 admin_granted_cpu_hours
            diff = balance_cpu_hours - user.balance_cpu_hours
            user.balance_cpu_hours = balance_cpu_hours
            user.admin_granted_cpu_hours = (user.admin_granted_cpu_hours or 0) + diff

        if custom_cpu_hour_price is not None:
            changes["custom_cpu_hour_price"] = f"{user.custom_cpu_hour_price} -> {custom_cpu_hour_price}"
            user.custom_cpu_hour_price = custom_cpu_hour_price

        db.commit()

        # Create audit log
        create_audit_log(
            db=db,
            user=admin,
            action="update_master_account",
            resource_type="master_account",
            resource_id=master_id,
            details={
                "username": user.username,
                "changes": changes,
                "admin": admin.username
            },
            request=request
        )

        logger.info(f"Master account {master_id} updated: {changes}")

        return {
            "id": ma.id,
            "user_id": ma.user_id,
            "username": user.username,
            "email": user.email,
            "organization": user.organization,
            "custom_cpu_hour_price": user.custom_cpu_hour_price,
            "balance_cpu_hours": user.balance_cpu_hours,
            "frozen_cpu_hours": user.frozen_cpu_hours,
            "free_cpu_hours_granted": user.free_cpu_hours_granted,
            "admin_granted_cpu_hours": user.admin_granted_cpu_hours,
            "max_sub_accounts": ma.max_sub_accounts,
            "current_sub_accounts": ma.current_sub_accounts,
            "is_active": ma.is_active,
            "created_at": ma.created_at,
            "updated_at": ma.updated_at,
        }
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Failed to update master account: {e}")
        raise HTTPException(status_code=500, detail=f"更新主账号失败: {str(e)}")


@router.get("/master-accounts/{master_id}", response_model=dict)
async def get_master_account_detail(
    master_id: int,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """
    Get master account detail (admin only)

    Args:
        master_id: Master account ID
        db: Database session
        admin: Current admin user

    Returns:
        Master account detail
    """
    try:
        from app.models.organization_v2 import MasterAccount, SubAccount

        ma = db.query(MasterAccount).filter(MasterAccount.id == master_id).first()
        if not ma:
            raise HTTPException(status_code=404, detail="Master account not found")

        user = db.query(User).filter(User.id == ma.user_id).first()
        if not user:
            raise HTTPException(status_code=404, detail="User not found")

        # Get sub-accounts
        sub_accounts = db.query(SubAccount).filter(SubAccount.master_account_id == master_id).all()

        return {
            "id": ma.id,
            "user_id": ma.user_id,
            "username": user.username,
            "email": user.email,
            "organization": user.organization,
            # 主账号的配额来自 User 表
            "balance_cpu_hours": user.balance_cpu_hours,
            "frozen_cpu_hours": user.frozen_cpu_hours,
            "free_cpu_hours_granted": user.free_cpu_hours_granted,
            "recharge_cpu_hours": user.recharge_cpu_hours,
            "admin_granted_cpu_hours": user.admin_granted_cpu_hours,
            "current_sub_accounts": ma.current_sub_accounts,
            "max_sub_accounts": ma.max_sub_accounts,
            "is_active": ma.is_active,
            "created_at": ma.created_at,
            "updated_at": ma.updated_at,
            "sub_accounts": [
                {
                    "id": sa.id,
                    "user_id": sa.user_id,
                    "username": db.query(User).filter(User.id == sa.user_id).first().username if sa.user_id else None,
                    "balance_cpu_hours": db.query(User).filter(User.id == sa.user_id).first().balance_cpu_hours if sa.user_id else 0,
                    "frozen_cpu_hours": db.query(User).filter(User.id == sa.user_id).first().frozen_cpu_hours if sa.user_id else 0,
                    "allocated_quota": sa.allocated_quota,
                    "is_active": sa.is_active,
                    "created_at": sa.created_at,
                }
                for sa in sub_accounts
            ]
        }
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Failed to get master account detail: {e}")
        raise HTTPException(status_code=500, detail=f"Failed to get master account detail: {str(e)}")


@router.get("/master-accounts/{master_id}/sub-accounts", response_model=List[dict])
async def get_master_account_sub_accounts(
    master_id: int,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """
    Get sub-accounts of a master account (admin only)

    Args:
        master_id: Master account ID
        db: Database session
        admin: Current admin user

    Returns:
        List of sub-accounts
    """
    try:
        from app.models.organization_v2 import MasterAccount, SubAccount

        ma = db.query(MasterAccount).filter(MasterAccount.id == master_id).first()
        if not ma:
            raise HTTPException(status_code=404, detail="Master account not found")

        sub_accounts = db.query(SubAccount).filter(SubAccount.master_account_id == master_id).all()

        result = []
        for sa in sub_accounts:
            user = db.query(User).filter(User.id == sa.user_id).first()
            if user:
                result.append({
                    "id": sa.id,
                    "user_id": sa.user_id,
                    "username": user.username,
                    "email": user.email,
                    # 子账号的配额来自两个来源
                    "balance_cpu_hours": user.balance_cpu_hours,  # 子账号自己充值的余额
                    "frozen_cpu_hours": user.frozen_cpu_hours,    # 冻结核时
                    "allocated_quota": sa.allocated_quota,        # 主账号分配的配额
                    "is_active": sa.is_active,
                    "created_at": sa.created_at,
                    "updated_at": sa.updated_at,
                })

        return result
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Failed to get sub-accounts: {e}")
        raise HTTPException(status_code=500, detail=f"Failed to get sub-accounts: {str(e)}")


@router.post("/master-accounts/{master_id}/sub-accounts", response_model=dict, status_code=status.HTTP_201_CREATED)
async def create_sub_account_admin(
    master_id: int,
    request: Request,
    username: str,
    email: str,
    password: str,
    allocated_quota: Optional[float] = None,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """
    Create sub-account (admin only)

    Args:
        master_id: Master account ID
        request: FastAPI request
        username: Sub-account username
        email: Sub-account email
        password: Sub-account password
        allocated_quota: Quota allocated by master account (optional)
        db: Database session
        admin: Current admin user

    Returns:
        Created sub-account info
    """
    try:
        from app.models.organization_v2 import MasterAccount, SubAccount

        # Verify master account exists
        ma = db.query(MasterAccount).filter(MasterAccount.id == master_id).first()
        if not ma:
            raise HTTPException(status_code=404, detail="Master account not found")

        # Check if username already exists
        existing_user = db.query(User).filter(User.username == username).first()
        if existing_user:
            raise HTTPException(status_code=400, detail="Username already exists")

        # Check if email already exists
        existing_email = db.query(User).filter(User.email == email).first()
        if existing_email:
            raise HTTPException(status_code=400, detail="Email already exists")

        # Create new user
        new_user = User(
            username=username,
            email=email,
            password_hash=get_password_hash(password),
            account_type=AccountType.SUB_ACCOUNT.value,
            is_active=True,
            email_verified=True,
        )
        db.add(new_user)
        db.flush()

        # Create sub-account record
        sub_account = SubAccount(
            master_account_id=master_id,
            user_id=new_user.id,
            allocated_quota=allocated_quota or 0.0,  # 主账号分配的配额
            is_active=True,
        )
        db.add(sub_account)

        # Update master account sub-account count
        ma.current_sub_accounts += 1

        db.commit()

        # Create audit log
        create_audit_log(
            db=db,
            user=admin,
            action="create_sub_account",
            resource_type="sub_account",
            resource_id=sub_account.id,
            details={
                "master_account_id": master_id,
                "username": username,
                "email": email,
                "allocated_quota": allocated_quota or 0.0,
                "admin": admin.username
            },
            request=request
        )

        logger.info(f"Sub-account created by admin: {username} under master account {master_id}")

        return {
            "id": sub_account.id,
            "master_account_id": sub_account.master_account_id,
            "user_id": new_user.id,
            "username": new_user.username,
            "email": new_user.email,
            "balance_cpu_hours": new_user.balance_cpu_hours,  # 子账号自己充值的余额
            "frozen_cpu_hours": new_user.frozen_cpu_hours,    # 冻结核时
            "allocated_quota": sub_account.allocated_quota,   # 主账号分配的配额
            "is_active": sub_account.is_active,
            "created_at": sub_account.created_at,
            "updated_at": sub_account.updated_at,
        }
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Failed to create sub-account: {e}")
        raise HTTPException(status_code=500, detail=f"Failed to create sub-account: {str(e)}")


@router.delete("/master-accounts/{master_id}", status_code=status.HTTP_204_NO_CONTENT)
async def delete_master_account(
    master_id: int,
    request: Request,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """
    Delete master account (admin only)

    将主账号转换为个人账号，所有子账号也转换为个人账号

    Args:
        master_id: Master account ID
        request: FastAPI request
        db: Database session
        admin: Current admin user
    """
    try:
        from app.models.organization_v2 import MasterAccount, SubAccount
        from app.models.user import AccountType

        # Verify master account exists
        ma = db.query(MasterAccount).filter(MasterAccount.id == master_id).first()
        if not ma:
            raise HTTPException(status_code=404, detail="Master account not found")

        # Get master account user
        master_user = db.query(User).filter(User.id == ma.user_id).first()
        if not master_user:
            raise HTTPException(status_code=404, detail="Master account user not found")

        # Prevent deleting yourself
        if master_user.id == admin.id:
            raise HTTPException(status_code=400, detail="Cannot delete your own master account")

        # Get all sub-accounts
        sub_accounts = db.query(SubAccount).filter(SubAccount.master_account_id == master_id).all()

        # Convert all sub-accounts to personal accounts
        for sub_account in sub_accounts:
            sub_user = db.query(User).filter(User.id == sub_account.user_id).first()
            if sub_user:
                sub_user.account_type = AccountType.PERSONAL.value
                logger.info(f"Sub-account {sub_account.id} (user {sub_user.username}) converted to personal account")

        # Convert master account to personal account
        master_user.account_type = AccountType.PERSONAL.value

        # Delete master account record
        db.delete(ma)

        # Log the action
        create_audit_log(
            db=db,
            user=admin,
            action="delete_master_account",
            resource_type="master_account",
            resource_id=master_id,
            details={
                "username": master_user.username,
                "email": master_user.email,
                "sub_accounts_count": len(sub_accounts)
            },
            request=request
        )

        db.commit()
        logger.info(f"Master account {master_id} deleted and converted to personal account")

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Failed to delete master account: {e}")
        raise HTTPException(status_code=500, detail=f"Failed to delete master account: {str(e)}")


@router.put("/master-accounts/{master_id}/sub-accounts/{sub_account_id}")
async def update_sub_account(
    master_id: int,
    sub_account_id: int,
    request: Request,
    allocated_quota: Optional[float] = None,
    is_active: Optional[bool] = None,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """
    Update sub-account (admin only)

    Args:
        master_id: Master account ID
        sub_account_id: Sub-account ID
        request: FastAPI request
        allocated_quota: New allocated quota from master account (optional)
        is_active: New active status (optional)
        db: Database session
        admin: Current admin user

    Returns:
        Updated sub-account info
    """
    try:
        from app.models.organization_v2 import MasterAccount, SubAccount

        # Verify master account exists
        ma = db.query(MasterAccount).filter(MasterAccount.id == master_id).first()
        if not ma:
            raise HTTPException(status_code=404, detail="Master account not found")

        # Get sub-account
        sub_account = db.query(SubAccount).filter(
            SubAccount.id == sub_account_id,
            SubAccount.master_account_id == master_id
        ).first()

        if not sub_account:
            raise HTTPException(status_code=404, detail="Sub-account not found")

        # Get user info
        user = db.query(User).filter(User.id == sub_account.user_id).first()
        if not user:
            raise HTTPException(status_code=404, detail="User not found")

        # Track changes for audit log
        changes = {}
        old_allocated_quota = sub_account.allocated_quota
        old_active = sub_account.is_active

        # Update fields
        if allocated_quota is not None:
            sub_account.allocated_quota = allocated_quota
            changes["allocated_quota"] = {"old": old_allocated_quota, "new": allocated_quota}
        if is_active is not None:
            sub_account.is_active = is_active
            user.is_active = is_active
            changes["is_active"] = {"old": old_active, "new": is_active}

        db.commit()
        db.refresh(sub_account)

        # Create audit log
        if changes:
            create_audit_log(
                db=db,
                user=admin,
                action="update_sub_account",
                resource_type="sub_account",
                resource_id=sub_account_id,
                details={
                    "master_account_id": master_id,
                    "sub_account_user": user.username,
                    "changes": changes,
                    "admin": admin.username
                },
                request=request
            )

        return {
            "id": sub_account.id,
            "user_id": sub_account.user_id,
            "username": user.username,
            "email": user.email,
            # 子账号的配额来自两个来源
            "balance_cpu_hours": user.balance_cpu_hours,  # 子账号自己充值的余额
            "frozen_cpu_hours": user.frozen_cpu_hours,    # 冻结核时
            "allocated_quota": sub_account.allocated_quota,  # 主账号分配的配额
            "is_active": sub_account.is_active,
            "created_at": sub_account.created_at,
            "updated_at": sub_account.updated_at,
        }
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Failed to update sub-account: {e}")
        raise HTTPException(status_code=500, detail=f"Failed to update sub-account: {str(e)}")


# ============ Slurm Partition Management ============

@router.get("/partitions")
async def get_partitions(
    admin: User = Depends(get_current_admin_user)
):
    """
    Get available Slurm partitions (admin only)

    Args:
        admin: Current admin user

    Returns:
        List of available partitions
    """
    try:
        partitions = list_partitions()
        return {"partitions": partitions}
    except Exception as e:
        logger.error(f"Failed to list partitions: {e}")
        raise HTTPException(status_code=500, detail="Failed to list partitions")


# ============ Permission Management Endpoints ============

@router.put("/users/{user_id}/gaussian-permission")
async def set_gaussian_permission(
    user_id: int,
    can_use: bool,
    request: Request,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """
    设置用户Gaussian使用权限 (admin only)
    
    Gaussian引擎需要商业许可证，只有授权用户才能使用。
    
    Args:
        user_id: 用户ID
        can_use: 是否允许使用Gaussian
        request: FastAPI request
        db: 数据库会话
        admin: 当前管理员
        
    Returns:
        更新结果
    """
    user = db.query(User).filter(User.id == user_id).first()
    if not user:
        raise HTTPException(status_code=404, detail="User not found")
    
    old_permission = user.can_use_gaussian
    user.can_use_gaussian = can_use
    db.commit()
    
    # 记录审计日志
    create_audit_log(
        db=db,
        user=admin,
        action="update_gaussian_permission",
        resource_type="user",
        resource_id=user.id,
        details={
            "username": user.username,
            "old_permission": old_permission,
            "new_permission": can_use
        },
        request=request
    )
    
    logger.info(f"Admin {admin.id} set user {user_id} Gaussian permission: {can_use}")
    
    return {
        "message": f"Gaussian permission {'granted' if can_use else 'revoked'} successfully",
        "user_id": user_id,
        "username": user.username,
        "can_use_gaussian": can_use
    }


@router.get("/users/gaussian-users")
async def get_gaussian_users(
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """
    获取所有有Gaussian权限的用户列表 (admin only)
    
    Args:
        db: 数据库会话
        admin: 当前管理员
        
    Returns:
        有Gaussian权限的用户列表
    """
    users = db.query(User).filter(User.can_use_gaussian == True).all()
    
    return {
        "total": len(users),
        "users": [
            {
                "id": u.id,
                "username": u.username,
                "email": u.email,
                "role": u.role.value,
                "organization": u.organization
            }
                "available_nodes": 3,
                "total_cpus": 128,
                "available_cpus": 96,
                "max_time": "3-00:00:00",
            },
            {
                "name": "debug",
                "state": "up",
                "total_nodes": 2,
                "available_nodes": 2,
                "total_cpus": 64,
                "available_cpus": 64,
                "max_time": "01:00:00",
            },
        ]


