"""
Permission dependencies for billing and quota control

提供权限控制的依赖函数
"""
from fastapi import HTTPException, Depends
from sqlalchemy.orm import Session
from app.models.user import User
from app.api.deps import get_current_user, get_db


async def require_no_debt(
    current_user: User = Depends(get_current_user)
) -> User:
    """
    要求用户无欠费
    
    用于需要查看结果的API端点
    
    Raises:
        HTTPException: 如果用户有欠费
    
    Returns:
        User对象
    """
    if current_user.balance_cpu_hours < 0:
        debt_amount = abs(current_user.balance_cpu_hours)
        raise HTTPException(
            status_code=403,
            detail={
                'error': 'account_in_debt',
                'message': f'账户欠费{debt_amount:.2f}核时，请充值后查看结果',
                'debt': debt_amount,
                'balance': current_user.balance_cpu_hours
            }
        )
    return current_user


async def check_quota_for_submission(
    estimated_hours: float,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
) -> User:
    """
    检查用户配额是否足够提交任务
    
    Args:
        estimated_hours: 预估核时消耗
    
    Raises:
        HTTPException: 如果配额不足
    
    Returns:
        User对象
    """
    from app.services.quota_check import check_submission_quota
    
    quota_check = check_submission_quota(
        user=current_user,
        db=db,
        estimated_hours=estimated_hours,
        min_balance_required=1.0
    )
    
    if not quota_check['can_submit']:
        raise HTTPException(
            status_code=403,
            detail={
                'error': 'insufficient_quota',
                'message': quota_check['reason'],
                'available_quota': quota_check['available_quota'],
                'required_quota': quota_check['required_quota'],
                'estimated_hours': estimated_hours,
                'debt': quota_check.get('debt', 0.0)
            }
        )
    
    return current_user
