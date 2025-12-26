"""
Compensation API endpoints
"""
import logging
from typing import List, Optional
from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy.orm import Session
from pydantic import BaseModel, Field

from app.database import get_db
from app.dependencies import get_current_admin_user, get_current_active_user
from app.models import User, CompensationRecord, CompensationStatus
from app.services.compensation import CompensationService

logger = logging.getLogger(__name__)
router = APIRouter(prefix="/compensation", tags=["compensation"])


# ============ Schemas ============

class CompensationRecordResponse(BaseModel):
    """补偿记录响应"""
    id: int
    user_id: int
    status: str
    reason: str
    amount: float
    reference_type: Optional[str]
    created_at: str
    approved_at: Optional[str]
    completed_at: Optional[str]

    class Config:
        from_attributes = True


class ApproveCompensationRequest(BaseModel):
    """批准补偿请求"""
    approval_reason: Optional[str] = Field(None, max_length=500)


class RejectCompensationRequest(BaseModel):
    """拒绝补偿请求"""
    rejection_reason: str = Field(..., max_length=500)


# ============ Admin Endpoints ============

@router.get("/admin/pending", response_model=List[CompensationRecordResponse])
async def get_pending_compensations(
    skip: int = 0,
    limit: int = 50,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """获取待处理的补偿记录（管理员）"""
    records = db.query(CompensationRecord).filter(
        CompensationRecord.status == CompensationStatus.PENDING
    ).order_by(CompensationRecord.created_at.desc()).offset(skip).limit(limit).all()

    return [CompensationRecordResponse(
        id=r.id,
        user_id=r.user_id,
        status=r.status.value,
        reason=r.reason,
        amount=r.amount,
        reference_type=r.reference_type,
        created_at=r.created_at.isoformat(),
        approved_at=r.approved_at.isoformat() if r.approved_at else None,
        completed_at=r.completed_at.isoformat() if r.completed_at else None
    ) for r in records]


@router.get("/admin/user/{user_id}", response_model=List[CompensationRecordResponse])
async def get_user_compensations(
    user_id: int,
    skip: int = 0,
    limit: int = 50,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """获取用户的补偿记录（管理员）"""
    user = db.query(User).filter(User.id == user_id).first()
    if not user:
        raise HTTPException(status_code=404, detail="用户不存在")

    records = db.query(CompensationRecord).filter(
        CompensationRecord.user_id == user_id
    ).order_by(CompensationRecord.created_at.desc()).offset(skip).limit(limit).all()

    return [CompensationRecordResponse(
        id=r.id,
        user_id=r.user_id,
        status=r.status.value,
        reason=r.reason,
        amount=r.amount,
        reference_type=r.reference_type,
        created_at=r.created_at.isoformat(),
        approved_at=r.approved_at.isoformat() if r.approved_at else None,
        completed_at=r.completed_at.isoformat() if r.completed_at else None
    ) for r in records]


@router.post("/admin/{record_id}/approve", response_model=dict)
async def approve_compensation(
    record_id: int,
    request: ApproveCompensationRequest,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """批准补偿记录（管理员）"""
    success, message = CompensationService.approve_compensation(
        db=db,
        record_id=record_id,
        approved_by=admin.id,
        approval_reason=request.approval_reason
    )

    if not success:
        raise HTTPException(status_code=400, detail=message)

    return {"success": True, "message": message}


@router.post("/admin/{record_id}/reject", response_model=dict)
async def reject_compensation(
    record_id: int,
    request: RejectCompensationRequest,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """拒绝补偿记录（管理员）"""
    success, message = CompensationService.reject_compensation(
        db=db,
        record_id=record_id,
        rejection_reason=request.rejection_reason
    )

    if not success:
        raise HTTPException(status_code=400, detail=message)

    return {"success": True, "message": message}


@router.get("/admin/stats", response_model=dict)
async def get_compensation_stats(
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """获取补偿统计信息（管理员）"""
    pending = db.query(CompensationRecord).filter(
        CompensationRecord.status == CompensationStatus.PENDING
    ).count()
    
    approved = db.query(CompensationRecord).filter(
        CompensationRecord.status == CompensationStatus.APPROVED
    ).count()
    
    rejected = db.query(CompensationRecord).filter(
        CompensationRecord.status == CompensationStatus.REJECTED
    ).count()
    
    completed = db.query(CompensationRecord).filter(
        CompensationRecord.status == CompensationStatus.COMPLETED
    ).count()

    total_amount = db.query(CompensationRecord).filter(
        CompensationRecord.status == CompensationStatus.COMPLETED
    ).with_entities(
        db.func.sum(CompensationRecord.amount)
    ).scalar() or 0.0

    return {
        "pending": pending,
        "approved": approved,
        "rejected": rejected,
        "completed": completed,
        "total_compensated_amount": total_amount
    }

