"""
账户管理 API - 主账号和子账号管理
"""
from fastapi import APIRouter, Depends, HTTPException, status
from pydantic import BaseModel
from sqlalchemy.orm import Session
from typing import List, Optional
from datetime import datetime

from app.database import get_db
from app.models.user import User, AccountType
from app.models.organization_v2 import MasterAccount, SubAccount
from app.dependencies import get_current_active_user
from app.core.security import get_password_hash
from app.core.logger import logger

router = APIRouter(prefix="/accounts", tags=["Accounts"])


# ============ 请求模型 ============

class CreateSubAccountRequest(BaseModel):
    username: str
    email: str
    password: str
    allocated_quota: Optional[float] = None  # 主账号分配的配额


class UpdateSubAccountRequest(BaseModel):
    allocated_quota: Optional[float] = None  # 主账号分配的配额
    is_active: Optional[bool] = None


class AddExistingUserAsSubAccountRequest(BaseModel):
    """添加现有用户为子账号的请求"""
    username_or_email: str  # 用户名或邮箱
    allocated_quota: Optional[float] = 0.0  # 分配的配额


# ============ 子账号管理 ============

@router.get("/my-sub-accounts", response_model=List[dict])
def get_my_sub_accounts(
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db)
):
    """
    获取当前用户的子账号列表（仅主账号可用）
    """
    if current_user.account_type != AccountType.MASTER_ACCOUNT.value:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Only master accounts can manage sub-accounts"
        )
    
    master_account = db.query(MasterAccount).filter(
        MasterAccount.user_id == current_user.id
    ).first()
    
    if not master_account:
        return []
    
    sub_accounts = db.query(SubAccount).filter(
        SubAccount.master_account_id == master_account.id
    ).all()
    
    result = []
    for sub in sub_accounts:
        user = db.query(User).filter(User.id == sub.user_id).first()
        result.append({
            "id": sub.id,
            "master_account_id": sub.master_account_id,
            "user_id": sub.user_id,
            "username": user.username if user else None,
            "email": user.email if user else None,
            # 子账号的配额来自两个来源
            "balance_cpu_hours": user.balance_cpu_hours if user else 0.0,  # 子账号自己充值的余额
            "frozen_cpu_hours": user.frozen_cpu_hours if user else 0.0,    # 冻结核时
            "allocated_quota": sub.allocated_quota,                         # 主账号分配的配额
            "is_active": sub.is_active,
            "created_at": sub.created_at,
            "updated_at": sub.updated_at,
        })

    return result


@router.post("/my-sub-accounts", response_model=dict, status_code=status.HTTP_201_CREATED)
def create_sub_account(
    request: CreateSubAccountRequest,
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db)
):
    """
    创建子账号（仅主账号可用）
    """
    if current_user.account_type != AccountType.MASTER_ACCOUNT.value:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Only master accounts can create sub-accounts"
        )

    # 检查用户名是否已存在
    existing_user = db.query(User).filter(User.username == request.username).first()
    if existing_user:
        raise HTTPException(status_code=400, detail="Username already exists")

    # 获取主账号
    master_account = db.query(MasterAccount).filter(
        MasterAccount.user_id == current_user.id
    ).first()

    if not master_account:
        raise HTTPException(status_code=400, detail="Master account not found")

    # 创建新用户
    new_user = User(
        username=request.username,
        email=request.email,
        password_hash=get_password_hash(request.password),
        account_type=AccountType.SUB_ACCOUNT.value,
        is_active=True,
        email_verified=True,
    )
    db.add(new_user)
    db.flush()

    # 创建子账号记录
    sub_account = SubAccount(
        master_account_id=master_account.id,
        user_id=new_user.id,
        allocated_quota=request.allocated_quota or 0.0,  # 主账号分配的配额
        is_active=True,
    )
    db.add(sub_account)

    # 更新主账号的子账号计数
    master_account.current_sub_accounts += 1

    db.commit()

    logger.info(f"Sub-account created: {request.username} by {current_user.username}")

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


@router.post("/my-sub-accounts/add-existing", response_model=dict, status_code=status.HTTP_201_CREATED)
def add_existing_user_as_sub_account(
    request: AddExistingUserAsSubAccountRequest,
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db)
):
    """
    将现有个人用户添加为子账号（仅主账号可用）

    只能添加 account_type 为 personal 的用户
    """
    if current_user.account_type != AccountType.MASTER_ACCOUNT.value:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Only master accounts can add sub-accounts"
        )

    # 获取主账号
    master_account = db.query(MasterAccount).filter(
        MasterAccount.user_id == current_user.id
    ).first()

    if not master_account:
        raise HTTPException(status_code=400, detail="Master account not found")

    # 检查子账号数量限制
    if master_account.current_sub_accounts >= master_account.max_sub_accounts:
        raise HTTPException(
            status_code=400,
            detail=f"已达到子账号数量上限 ({master_account.max_sub_accounts})"
        )

    # 查找目标用户（按用户名或邮箱）
    target_user = db.query(User).filter(
        (User.username == request.username_or_email) |
        (User.email == request.username_or_email)
    ).first()

    if not target_user:
        raise HTTPException(status_code=404, detail="用户不存在")

    # 检查目标用户是否是个人账户
    if target_user.account_type != AccountType.PERSONAL.value:
        if target_user.account_type == AccountType.MASTER_ACCOUNT.value:
            raise HTTPException(status_code=400, detail="不能将主账号添加为子账号")
        elif target_user.account_type == AccountType.SUB_ACCOUNT.value:
            raise HTTPException(status_code=400, detail="该用户已经是其他主账号的子账号")
        else:
            raise HTTPException(status_code=400, detail="只能将个人账户添加为子账号")

    # 检查是否已经是子账号
    existing_sub = db.query(SubAccount).filter(SubAccount.user_id == target_user.id).first()
    if existing_sub:
        raise HTTPException(status_code=400, detail="该用户已经是子账号")

    # 更新用户的 account_type
    target_user.account_type = AccountType.SUB_ACCOUNT.value

    # 创建子账号记录
    sub_account = SubAccount(
        master_account_id=master_account.id,
        user_id=target_user.id,
        allocated_quota=request.allocated_quota or 0.0,
        is_active=True,
    )
    db.add(sub_account)

    # 更新主账号的子账号计数
    master_account.current_sub_accounts += 1

    db.commit()
    db.refresh(sub_account)

    logger.info(f"Existing user {target_user.username} added as sub-account by {current_user.username}")

    return {
        "id": sub_account.id,
        "master_account_id": sub_account.master_account_id,
        "user_id": target_user.id,
        "username": target_user.username,
        "email": target_user.email,
        "balance_cpu_hours": target_user.balance_cpu_hours,
        "frozen_cpu_hours": target_user.frozen_cpu_hours,
        "allocated_quota": sub_account.allocated_quota,
        "is_active": sub_account.is_active,
        "created_at": sub_account.created_at,
        "updated_at": sub_account.updated_at,
    }


@router.put("/my-sub-accounts/{sub_account_id}", response_model=dict)
def update_sub_account(
    sub_account_id: int,
    request: UpdateSubAccountRequest,
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db)
):
    """
    更新子账号配额和状态
    """
    if current_user.account_type != AccountType.MASTER_ACCOUNT.value:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Only master accounts can manage sub-accounts"
        )
    
    master_account = db.query(MasterAccount).filter(
        MasterAccount.user_id == current_user.id
    ).first()
    
    sub_account = db.query(SubAccount).filter(
        SubAccount.id == sub_account_id,
        SubAccount.master_account_id == master_account.id
    ).first()
    
    if not sub_account:
        raise HTTPException(status_code=404, detail="Sub-account not found")

    if request.allocated_quota is not None:
        sub_account.allocated_quota = request.allocated_quota
    if request.is_active is not None:
        sub_account.is_active = request.is_active

    sub_account.updated_at = datetime.utcnow()
    db.commit()

    logger.info(f"Sub-account updated: {sub_account_id}")

    user = db.query(User).filter(User.id == sub_account.user_id).first()
    return {
        "id": sub_account.id,
        "master_account_id": sub_account.master_account_id,
        "user_id": sub_account.user_id,
        "username": user.username if user else None,
        "email": user.email if user else None,
        # 子账号的配额来自两个来源
        "balance_cpu_hours": user.balance_cpu_hours if user else 0.0,  # 子账号自己充值的余额
        "frozen_cpu_hours": user.frozen_cpu_hours if user else 0.0,    # 冻结核时
        "allocated_quota": sub_account.allocated_quota,                 # 主账号分配的配额
        "is_active": sub_account.is_active,
        "created_at": sub_account.created_at,
        "updated_at": sub_account.updated_at,
    }


@router.delete("/my-sub-accounts/{sub_account_id}", status_code=status.HTTP_204_NO_CONTENT)
def delete_sub_account(
    sub_account_id: int,
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db)
):
    """
    删除子账号

    删除子账号时，将该子账号用户的 account_type 改为 PERSONAL_ACCOUNT
    """
    if current_user.account_type != AccountType.MASTER_ACCOUNT.value:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Only master accounts can manage sub-accounts"
        )

    master_account = db.query(MasterAccount).filter(
        MasterAccount.user_id == current_user.id
    ).first()

    sub_account = db.query(SubAccount).filter(
        SubAccount.id == sub_account_id,
        SubAccount.master_account_id == master_account.id
    ).first()

    if not sub_account:
        raise HTTPException(status_code=404, detail="Sub-account not found")

    # 获取子账号关联的用户
    sub_account_user = db.query(User).filter(User.id == sub_account.user_id).first()

    # 删除子账号前，将用户的 account_type 改为 PERSONAL_ACCOUNT
    if sub_account_user:
        sub_account_user.account_type = AccountType.PERSONAL_ACCOUNT.value

    # 删除子账号记录
    db.delete(sub_account)

    # 更新主账号的子账号计数
    if master_account:
        master_account.current_sub_accounts = max(0, master_account.current_sub_accounts - 1)

    db.commit()

    logger.info(f"Sub-account deleted: {sub_account_id}, user {sub_account.user_id} converted to personal account")


# ============ 主账号信息 ============

@router.get("/my-master-account", response_model=dict)
def get_my_master_account(
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db)
):
    """
    获取当前用户的主账号信息
    """
    if current_user.account_type != AccountType.MASTER_ACCOUNT.value:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Only master accounts can access this endpoint"
        )
    
    master_account = db.query(MasterAccount).filter(
        MasterAccount.user_id == current_user.id
    ).first()
    
    if not master_account:
        raise HTTPException(status_code=404, detail="Master account not found")
    
    return {
        "id": master_account.id,
        "user_id": master_account.user_id,
        "username": current_user.username,
        "email": current_user.email,
        "organization": current_user.organization,
        # 主账号的配额来自 User 表
        "balance_cpu_hours": current_user.balance_cpu_hours,
        "frozen_cpu_hours": current_user.frozen_cpu_hours,
        "free_cpu_hours_granted": current_user.free_cpu_hours_granted,
        "recharge_cpu_hours": current_user.recharge_cpu_hours,
        "admin_granted_cpu_hours": current_user.admin_granted_cpu_hours,
        "current_sub_accounts": master_account.current_sub_accounts,
        "max_sub_accounts": master_account.max_sub_accounts,
        "is_active": master_account.is_active,
        "created_at": master_account.created_at,
        "updated_at": master_account.updated_at,
    }


@router.get("/my-sub-account", response_model=dict)
def get_my_sub_account(
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db)
):
    """
    获取当前用户的子账号信息（仅子账号可用）
    """
    if current_user.account_type != AccountType.SUB_ACCOUNT.value:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Only sub-accounts can access this endpoint"
        )
    
    sub_account = db.query(SubAccount).filter(
        SubAccount.user_id == current_user.id
    ).first()
    
    if not sub_account:
        raise HTTPException(status_code=404, detail="Sub-account not found")
    
    master_user = db.query(User).filter(
        User.id == sub_account.master_account.user_id
    ).first()

    return {
        "id": sub_account.id,
        "master_account_id": sub_account.master_account_id,
        "user_id": sub_account.user_id,
        "username": current_user.username,
        "email": current_user.email,
        "master_username": master_user.username if master_user else None,
        "master_email": master_user.email if master_user else None,
        # 子账号的配额来自两个来源
        "balance_cpu_hours": current_user.balance_cpu_hours,  # 子账号自己充值的余额
        "frozen_cpu_hours": current_user.frozen_cpu_hours,    # 冻结核时
        "allocated_quota": sub_account.allocated_quota,       # 主账号分配的配额
        "is_active": sub_account.is_active,
        "created_at": sub_account.created_at,
        "updated_at": sub_account.updated_at,
    }


# ============ 子账号任务查看 ============

@router.get("/my-sub-accounts/{sub_account_id}/jobs", response_model=dict)
def get_sub_account_jobs(
    sub_account_id: int,
    skip: int = 0,
    limit: int = 20,
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db)
):
    """
    获取子账号的任务列表（仅主账号可用）
    
    包括MD和QC任务
    """
    from app.models.job import MDJob, JobStatus
    from app.models.qc import QCJob, QCJobStatus

    if current_user.account_type != AccountType.MASTER_ACCOUNT.value:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Only master accounts can view sub-account jobs"
        )

    master_account = db.query(MasterAccount).filter(
        MasterAccount.user_id == current_user.id
    ).first()

    sub_account = db.query(SubAccount).filter(
        SubAccount.id == sub_account_id,
        SubAccount.master_account_id == master_account.id
    ).first()

    if not sub_account:
        raise HTTPException(status_code=404, detail="Sub-account not found")

    # 获取子账号用户信息
    sub_user = db.query(User).filter(User.id == sub_account.user_id).first()

    # 查询该子账号的MD和QC任务
    md_jobs = db.query(MDJob).filter(
        MDJob.user_id == sub_account.user_id,
        MDJob.is_deleted == False
    ).all()

    qc_jobs = db.query(QCJob).filter(
        QCJob.user_id == sub_account.user_id
    ).all()

    # 合并任务并按时间排序
    all_jobs = []
    
    for job in md_jobs:
        all_jobs.append({
            "id": job.id,
            "type": "md",
            "system_id": job.system_id,
            "status": job.status,
            "created_at": job.created_at,
            "started_at": job.started_at,
            "finished_at": job.finished_at,
            "actual_cpu_hours": job.actual_cpu_hours,
            "estimated_cpu_hours": job.estimated_cpu_hours,
        })
    
    for job in qc_jobs:
        all_jobs.append({
            "id": job.id,
            "type": "qc",
            "molecule_name": job.molecule_name,
            "status": job.status,
            "created_at": job.created_at,
            "started_at": job.started_at,
            "finished_at": job.finished_at,
            "actual_cpu_hours": job.actual_cpu_hours,
            "estimated_cpu_hours": getattr(job, 'estimated_cpu_hours', None),
        })
    
    # 按创建时间降序排序
    all_jobs.sort(key=lambda x: x["created_at"], reverse=True)
    
    total = len(all_jobs)
    paginated_jobs = all_jobs[skip:skip+limit]

    return {
        "sub_account_id": sub_account_id,
        "username": sub_user.username if sub_user else None,
        "total": total,
        "md_total": len(md_jobs),
        "qc_total": len(qc_jobs),
        "skip": skip,
        "limit": limit,
        "jobs": paginated_jobs
    }


@router.get("/my-sub-accounts/all-jobs", response_model=dict)
def get_all_sub_accounts_jobs(
    skip: int = 0,
    limit: int = 50,
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db)
):
    """
    获取所有子账号的任务汇总（仅主账号可用）
    
    包括MD和QC任务
    """
    from app.models.job import MDJob, JobStatus
    from app.models.qc import QCJob, QCJobStatus

    if current_user.account_type != AccountType.MASTER_ACCOUNT.value:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Only master accounts can view sub-account jobs"
        )

    master_account = db.query(MasterAccount).filter(
        MasterAccount.user_id == current_user.id
    ).first()

    if not master_account:
        raise HTTPException(status_code=400, detail="Master account not found")

    # 获取所有子账号的 user_id
    sub_accounts = db.query(SubAccount).filter(
        SubAccount.master_account_id == master_account.id
    ).all()

    sub_user_ids = [sa.user_id for sa in sub_accounts]

    if not sub_user_ids:
        return {
            "total": 0,
            "skip": skip,
            "limit": limit,
            "jobs": []
        }

    # 查询所有子账号的任务
    total = db.query(MDJob).filter(
        MDJob.user_id.in_(sub_user_ids),
        MDJob.is_deleted == False
    ).count()

    jobs = db.query(MDJob, User).join(User, MDJob.user_id == User.id).filter(
        MDJob.user_id.in_(sub_user_ids),
        MDJob.is_deleted == False
    ).order_by(MDJob.created_at.desc()).offset(skip).limit(limit).all()

    return {
        "total": total,
        "skip": skip,
        "limit": limit,
        "jobs": [
            {
                "id": job.id,
                "system_id": job.system_id,
                "status": job.status,
                "username": user.username,
                "user_id": user.id,
                "created_at": job.created_at,
                "started_at": job.started_at,
                "finished_at": job.finished_at,
                "actual_cpu_hours": job.actual_cpu_hours,
                "estimated_cpu_hours": job.estimated_cpu_hours,
            }
            for job, user in jobs
        ]
    }


# ============ 增强功能: 使用统计和概览 ============

@router.get("/my-sub-accounts/{sub_account_id}/usage-stats", response_model=dict)
def get_sub_account_usage_stats(
    sub_account_id: int,
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db)
):
    """
    获取子账号的使用统计（仅主账号可用）
    
    包括配额使用、任务统计、CPU hours分类等
    """
    from app.models.job import MDJob, JobStatus
    from app.models.qc import QCJob, QCJobStatus
    from app.models.billing import QuotaTransaction
    from sqlalchemy import func

    if current_user.account_type != AccountType.MASTER_ACCOUNT.value:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Only master accounts can view sub-account stats"
        )

    master_account = db.query(MasterAccount).filter(
        MasterAccount.user_id == current_user.id
    ).first()

    sub_account = db.query(SubAccount).filter(
        SubAccount.id == sub_account_id,
        SubAccount.master_account_id == master_account.id
    ).first()

    if not sub_account:
        raise HTTPException(status_code=404, detail="Sub-account not found")

    # 获取子账号用户信息
    sub_user = db.query(User).filter(User.id == sub_account.user_id).first()

    # 配额统计
    personal_quota = sub_user.balance_cpu_hours
    allocated_quota = sub_account.allocated_quota
    frozen_quota = sub_user.frozen_cpu_hours
    total_available = personal_quota + allocated_quota - frozen_quota

    # 计算已使用的核时（从交易记录）
    total_consumed = db.query(func.sum(QuotaTransaction.amount)).filter(
        QuotaTransaction.user_id == sub_user.id,
        QuotaTransaction.type == 'consume'
    ).scalar() or 0.0
    total_used = abs(total_consumed)

    # 总配额 = 个人 + 分配
    total_quota = personal_quota + allocated_quota
    usage_percentage = (total_used / total_quota * 100) if total_quota > 0 else 0

    # MD任务统计
    md_total = db.query(MDJob).filter(
        MDJob.user_id == sub_user.id,
        MDJob.is_deleted == False
    ).count()
    
    md_completed = db.query(MDJob).filter(
        MDJob.user_id == sub_user.id,
        MDJob.status == JobStatus.COMPLETED,
        MDJob.is_deleted == False
    ).count()
    
    md_running = db.query(MDJob).filter(
        MDJob.user_id == sub_user.id,
        MDJob.status.in_([JobStatus.RUNNING, JobStatus.QUEUED, JobStatus.SUBMITTED]),
        MDJob.is_deleted == False
    ).count()
    
    md_failed = db.query(MDJob).filter(
        MDJob.user_id == sub_user.id,
        MDJob.status == JobStatus.FAILED,
        MDJob.is_deleted == False
    ).count()

    # MD任务CPU hours
    md_cpu_hours = db.query(func.sum(MDJob.actual_cpu_hours)).filter(
        MDJob.user_id == sub_user.id,
        MDJob.actual_cpu_hours.isnot(None),
        MDJob.is_deleted == False
    ).scalar() or 0.0

    # QC任务统计
    qc_total = db.query(QCJob).filter(QCJob.user_id == sub_user.id).count()
    
    qc_completed = db.query(QCJob).filter(
        QCJob.user_id == sub_user.id,
        QCJob.status == QCJobStatus.COMPLETED
    ).count()
    
    qc_running = db.query(QCJob).filter(
        QCJob.user_id == sub_user.id,
        QCJob.status.in_([QCJobStatus.RUNNING, QCJobStatus.QUEUED, QCJobStatus.SUBMITTED, QCJobStatus.RETRYING])
    ).count()
    
    qc_failed = db.query(QCJob).filter(
        QCJob.user_id == sub_user.id,
        QCJob.status == QCJobStatus.FAILED
    ).count()

    # QC任务CPU hours
    qc_cpu_hours = db.query(func.sum(QCJob.actual_cpu_hours)).filter(
        QCJob.user_id == sub_user.id,
        QCJob.actual_cpu_hours.isnot(None)
    ).scalar() or 0.0

    # 总任务统计
    total_jobs = md_total + qc_total
    total_completed = md_completed + qc_completed
    total_running = md_running + qc_running
    total_failed = md_failed + qc_failed
    success_rate = (total_completed / total_jobs * 100) if total_jobs > 0 else 0

    return {
        "sub_account_id": sub_account_id,
        "username": sub_user.username,
        "email": sub_user.email,
        "is_active": sub_account.is_active,
        "quota": {
            "personal": round(personal_quota, 2),
            "allocated": round(allocated_quota, 2),
            "frozen": round(frozen_quota, 2),
            "total_available": round(total_available, 2),
            "total_used": round(total_used, 2),
            "total_quota": round(total_quota, 2),
            "usage_percentage": round(usage_percentage, 1)
        },
        "jobs": {
            "total": total_jobs,
            "completed": total_completed,
            "running": total_running,
            "failed": total_failed,
            "success_rate": round(success_rate, 1),
            "md_jobs": md_total,
            "qc_jobs": qc_total
        },
        "cpu_hours_breakdown": {
            "total": round(total_used, 2),
            "md_jobs": round(md_cpu_hours, 2),
            "qc_jobs": round(qc_cpu_hours, 2)
        }
    }


@router.get("/my-sub-accounts/overview", response_model=dict)
def get_sub_accounts_overview(
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db)
):
    """
    获取所有子账号的概览（仅主账号可用）
    
    提供快速监控仪表板，显示所有子账号的关键指标
    """
    from app.models.job import MDJob, JobStatus
    from app.models.qc import QCJob, QCJobStatus
    from app.models.billing import QuotaTransaction
    from sqlalchemy import func

    if current_user.account_type != AccountType.MASTER_ACCOUNT.value:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Only master accounts can view overview"
        )

    master_account = db.query(MasterAccount).filter(
        MasterAccount.user_id == current_user.id
    ).first()

    if not master_account:
        raise HTTPException(status_code=400, detail="Master account not found")

    # 获取所有子账号
    sub_accounts = db.query(SubAccount).filter(
        SubAccount.master_account_id == master_account.id
    ).all()

    # 主账号信息
    master_info = {
        "username": current_user.username,
        "balance": round(current_user.balance_cpu_hours, 2),
        "total_sub_accounts": len(sub_accounts),
        "active_sub_accounts": sum(1 for sa in sub_accounts if sa.is_active),
        "max_sub_accounts": master_account.max_sub_accounts
    }

    # 计算总分配配额
    total_allocated = sum(sa.allocated_quota for sa in sub_accounts)
    master_info["total_allocated"] = round(total_allocated, 2)

    # 子账号概览
    sub_accounts_summary = []
    for sub_account in sub_accounts:
        sub_user = db.query(User).filter(User.id == sub_account.user_id).first()
        if not sub_user:
            continue

        # 计算已使用核时
        total_consumed = db.query(func.sum(QuotaTransaction.amount)).filter(
            QuotaTransaction.user_id == sub_user.id,
            QuotaTransaction.type == 'consume'
        ).scalar() or 0.0
        used_cpu_hours = abs(total_consumed)

        # 总可用配额
        total_quota = sub_user.balance_cpu_hours + sub_account.allocated_quota
        usage_percentage = (used_cpu_hours / total_quota * 100) if total_quota > 0 else 0

        # 运行中的任务
        md_active = db.query(MDJob).filter(
            MDJob.user_id == sub_user.id,
            MDJob.status.in_([JobStatus.RUNNING, JobStatus.QUEUED, JobStatus.SUBMITTED]),
            MDJob.is_deleted == False
        ).count()

        qc_active = db.query(QCJob).filter(
            QCJob.user_id == sub_user.id,
            QCJob.status.in_([QCJobStatus.RUNNING, QCJobStatus.QUEUED, QCJobStatus.SUBMITTED, QCJobStatus.RETRYING])
        ).count()

        active_jobs = md_active + qc_active

        sub_accounts_summary.append({
            "sub_account_id": sub_account.id,
            "user_id": sub_user.id,
            "username": sub_user.username,
            "email": sub_user.email,
            "is_active": sub_account.is_active,
            "allocated_quota": round(sub_account.allocated_quota, 2),
            "personal_quota": round(sub_user.balance_cpu_hours, 2),
            "total_quota": round(total_quota, 2),
            "used_cpu_hours": round(used_cpu_hours, 2),
            "usage_percentage": round(usage_percentage, 1),
            "active_jobs": active_jobs,
            "frozen_cpu_hours": round(sub_user.frozen_cpu_hours, 2)
        })

    # 按使用率降序排序
    sub_accounts_summary.sort(key=lambda x: x["usage_percentage"], reverse=True)

    return {
        "master_account": master_info,
        "sub_accounts_summary": sub_accounts_summary
    }
