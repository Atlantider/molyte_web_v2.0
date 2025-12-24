"""
计费和充值 API 端点
"""
from typing import List, Optional, Dict
from datetime import datetime
from fastapi import APIRouter, Depends, HTTPException, status, Request
from pydantic import BaseModel, Field
from sqlalchemy.orm import Session

from app.database import get_db
from app.dependencies import get_current_active_user, get_current_admin_user
from app.models.user import User
from app.models.billing import RechargeOrder, QuotaTransaction, PaymentStatus, PaymentMethod
from app.services.billing import BillingService
from app.utils.audit import create_audit_log
from app.core.logger import logger

router = APIRouter(prefix="/billing", tags=["billing"])


# ============ Schemas ============

class BalanceResponse(BaseModel):
    """用户余额响应"""
    balance: float
    frozen: float
    debt: float
    available: float
    price_per_hour: float
    
    class Config:
        from_attributes = True


class CreateOrderRequest(BaseModel):
    """创建充值订单请求"""
    amount: float = Field(..., ge=1, description="充值金额（元）")
    payment_method: str = Field(default="simulated", description="支付方式")


class OrderResponse(BaseModel):
    """订单响应"""
    id: int
    order_no: str
    amount: float
    cpu_hours: float
    price_per_hour: float
    payment_method: str
    payment_status: str
    created_at: datetime
    paid_at: Optional[datetime]
    
    class Config:
        from_attributes = True


class TransactionResponse(BaseModel):
    """交易记录响应"""
    id: int
    type: str
    amount: float
    balance_before: float
    balance_after: float
    description: Optional[str]
    created_at: datetime
    
    class Config:
        from_attributes = True


class SimulatePayRequest(BaseModel):
    """模拟支付请求"""
    order_id: int


class ConsumptionDetailResponse(BaseModel):
    """消费详情响应"""
    id: int
    user_id: int
    username: str
    cpu_hours: float
    amount: float  # 保持兼容性，等于cpu_hours
    money_amount: float  # 新增：金额（核时 × 单价）
    price_per_hour: float  # 新增：核时单价
    created_at: datetime

    class Config:
        from_attributes = True


# ============ 用户端 API ============

@router.get("/balance", response_model=BalanceResponse)
async def get_balance(
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """获取当前用户余额信息"""
    balance_info = BillingService.get_user_balance(db, current_user.id)
    if not balance_info:
        raise HTTPException(status_code=404, detail="用户不存在")

    # 获取用户的有效定价（优先使用用户自定义价格）
    balance_info["price_per_hour"] = BillingService.get_cpu_hour_price(db, current_user.id)
    return balance_info


@router.get("/can-submit")
async def check_can_submit(
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """检查用户是否可以提交任务"""
    can_submit, reason = BillingService.can_submit_job(db, current_user)
    return {
        "can_submit": can_submit,
        "reason": reason
    }


@router.post("/orders", response_model=OrderResponse)
async def create_order(
    request: CreateOrderRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """创建充值订单"""
    from app.services.payment import payment_service

    # 验证最低充值金额
    min_amount = float(BillingService.get_config(db, 'min_recharge_amount', '10'))
    if request.amount < min_amount:
        raise HTTPException(
            status_code=400,
            detail=f"最低充值金额为 ¥{min_amount}"
        )

    # 映射支付方式
    method_map = {
        "wechat": PaymentMethod.WECHAT,
        "alipay": PaymentMethod.ALIPAY,
        "simulated": PaymentMethod.SIMULATED,
    }
    payment_method = method_map.get(request.payment_method, PaymentMethod.SIMULATED)

    order = BillingService.create_recharge_order(
        db=db,
        user_id=current_user.id,
        amount=request.amount,
        payment_method=payment_method
    )

    # 调用支付服务创建支付
    success, payment_data = await payment_service.create_payment(order)

    response = OrderResponse(
        id=order.id,
        order_no=order.order_no,
        amount=order.amount,
        cpu_hours=order.cpu_hours,
        price_per_hour=order.price_per_hour,
        payment_method=order.payment_method,
        payment_status=order.payment_status,
        created_at=order.created_at,
        paid_at=order.paid_at
    )

    # 返回支付参数
    return {
        **response.model_dump(),
        "payment_data": payment_data if success else None,
    }


@router.get("/orders", response_model=List[OrderResponse])
async def get_orders(
    skip: int = 0,
    limit: int = 20,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """获取用户充值订单列表"""
    orders = BillingService.get_orders(db, current_user.id, skip, limit)
    return [OrderResponse(
        id=o.id,
        order_no=o.order_no,
        amount=o.amount,
        cpu_hours=o.cpu_hours,
        price_per_hour=o.price_per_hour,
        payment_method=o.payment_method,  # 已经是字符串
        payment_status=o.payment_status,  # 已经是字符串
        created_at=o.created_at,
        paid_at=o.paid_at
    ) for o in orders]


@router.get("/transactions", response_model=List[TransactionResponse])
async def get_transactions(
    skip: int = 0,
    limit: int = 20,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """获取用户配额变更记录"""
    transactions = BillingService.get_transactions(db, current_user.id, skip, limit)
    return [TransactionResponse(
        id=t.id,
        type=t.type,  # 已经是字符串
        amount=t.amount,
        balance_before=t.balance_before,
        balance_after=t.balance_after,
        description=t.description,
        created_at=t.created_at
    ) for t in transactions]


@router.post("/simulate-pay")
async def simulate_payment(
    request: SimulatePayRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """模拟支付（测试用）"""
    # 验证订单属于当前用户
    order = db.query(RechargeOrder).filter(
        RechargeOrder.id == request.order_id,
        RechargeOrder.user_id == current_user.id
    ).first()

    if not order:
        raise HTTPException(status_code=404, detail="订单不存在")

    if order.payment_status != PaymentStatus.PENDING:
        raise HTTPException(status_code=400, detail=f"订单状态异常: {order.payment_status.value}")

    success, message = BillingService.complete_payment(db, order.id)

    if not success:
        raise HTTPException(status_code=400, detail=message)

    # 获取更新后的余额
    balance_info = BillingService.get_user_balance(db, current_user.id)

    return {
        "success": True,
        "message": message,
        "balance": balance_info
    }


@router.get("/orders/{order_id}/status")
async def query_payment_status(
    order_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """查询支付状态"""
    from app.services.payment import payment_service

    order = db.query(RechargeOrder).filter(
        RechargeOrder.id == order_id,
        RechargeOrder.user_id == current_user.id
    ).first()

    if not order:
        raise HTTPException(status_code=404, detail="订单不存在")

    # 如果已支付，直接返回
    if order.payment_status == PaymentStatus.PAID:
        return {
            "paid": True,
            "order_no": order.order_no,
            "amount": order.amount,
            "cpu_hours": order.cpu_hours,
        }

    # 查询支付状态（仅对真实支付方式）
    if order.payment_method in ['wechat', 'alipay']:
        success, result = await payment_service.query_payment(order)

        if success and result.get("paid"):
            # 完成支付
            BillingService.complete_payment(db, order.id)
            return {
                "paid": True,
                "order_no": order.order_no,
                "amount": order.amount,
                "cpu_hours": order.cpu_hours,
            }

    return {
        "paid": False,
        "order_no": order.order_no,
        "payment_status": order.payment_status,
    }


@router.post("/callback/wechat")
async def wechat_payment_callback(
    request: Request,
    db: Session = Depends(get_db)
):
    """微信支付回调"""
    from app.services.payment import payment_service

    try:
        body = await request.json()
        success, order_no = payment_service.verify_callback("wechat", body)

        if success and order_no:
            order = db.query(RechargeOrder).filter(
                RechargeOrder.order_no == order_no
            ).first()

            if order and order.payment_status == PaymentStatus.PENDING:
                BillingService.complete_payment(db, order.id)

        return {"code": "SUCCESS", "message": "成功"}
    except Exception as e:
        logger.error(f"微信支付回调处理失败: {e}")
        return {"code": "FAIL", "message": str(e)}


@router.post("/callback/alipay")
async def alipay_payment_callback(
    request: Request,
    db: Session = Depends(get_db)
):
    """支付宝支付回调"""
    from app.services.payment import payment_service

    try:
        form_data = await request.form()
        data = dict(form_data)
        success, order_no = payment_service.verify_callback("alipay", data)

        if success and order_no:
            order = db.query(RechargeOrder).filter(
                RechargeOrder.order_no == order_no
            ).first()

            if order and order.payment_status == PaymentStatus.PENDING:
                BillingService.complete_payment(db, order.id)

        return "success"
    except Exception as e:
        logger.error(f"支付宝回调处理失败: {e}")
        return "fail"


# ============ 管理端 API ============

class PricingConfig(BaseModel):
    """定价配置"""
    cpu_hour_price: float = Field(..., ge=0.01, description="每核时单价（元）")
    min_recharge_amount: float = Field(..., ge=1, description="最低充值金额（元）")
    max_debt_cpu_hours: float = Field(..., ge=0, description="最大允许欠费机时")


class AdminAdjustRequest(BaseModel):
    """管理员调整余额请求"""
    user_id: int
    amount: float
    reason: str


class UserPricingConfig(BaseModel):
    """用户级别定价配置"""
    user_id: int
    custom_cpu_hour_price: Optional[float] = Field(None, ge=0.01, description="自定义核时单价（元），为 null 时使用全局定价")


class UserPricingResponse(BaseModel):
    """用户定价响应"""
    user_id: int
    username: str
    email: str
    custom_cpu_hour_price: Optional[float]
    global_price: float
    effective_price: float  # 实际生效的价格
    price_updated_at: Optional[str]
    price_updated_by: Optional[int]


@router.get("/admin/pricing", response_model=PricingConfig)
async def get_pricing_config(
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """获取定价配置（管理员）"""
    return PricingConfig(
        cpu_hour_price=float(BillingService.get_config(db, 'cpu_hour_price', '0.1')),
        min_recharge_amount=float(BillingService.get_config(db, 'min_recharge_amount', '10')),
        max_debt_cpu_hours=float(BillingService.get_config(db, 'max_debt_cpu_hours', '100'))
    )


@router.put("/admin/pricing")
async def update_pricing_config(
    config: PricingConfig,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """更新定价配置（管理员）"""
    BillingService.set_config(db, 'cpu_hour_price', str(config.cpu_hour_price),
                              '每核时单价（元）', admin.id)
    BillingService.set_config(db, 'min_recharge_amount', str(config.min_recharge_amount),
                              '最低充值金额（元）', admin.id)
    BillingService.set_config(db, 'max_debt_cpu_hours', str(config.max_debt_cpu_hours),
                              '最大允许欠费机时', admin.id)

    return {"success": True, "message": "定价配置已更新"}


@router.post("/admin/adjust")
async def admin_adjust_balance(
    request: AdminAdjustRequest,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """管理员手动调整用户余额"""
    success, message = BillingService.admin_adjust_balance(
        db=db,
        admin_id=admin.id,
        user_id=request.user_id,
        amount=request.amount,
        reason=request.reason
    )

    if not success:
        raise HTTPException(status_code=400, detail=message)

    return {"success": True, "message": message}


@router.post("/admin/grant")
async def admin_grant_cpu_hours(
    request: AdminAdjustRequest,
    http_request: Request,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """管理员赠送核时给用户（区别于调整余额）"""
    # 获取目标用户信息
    target_user = db.query(User).filter(User.id == request.user_id).first()
    if not target_user:
        raise HTTPException(status_code=404, detail="用户不存在")

    success, message = BillingService.admin_grant_cpu_hours(
        db=db,
        admin_id=admin.id,
        user_id=request.user_id,
        amount=request.amount,
        reason=request.reason
    )

    if not success:
        raise HTTPException(status_code=400, detail=message)

    # 创建审计日志
    create_audit_log(
        db=db,
        user=admin,
        action="grant_cpu_hours",
        resource_type="user",
        resource_id=request.user_id,
        details={
            "target_user": target_user.username,
            "amount": request.amount,
            "reason": request.reason,
            "admin": admin.username
        },
        request=http_request
    )

    return {"success": True, "message": message}


@router.get("/admin/transactions/{user_id}", response_model=List[TransactionResponse])
async def get_user_transactions_admin(
    user_id: int,
    skip: int = 0,
    limit: int = 50,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """获取用户的配额变更记录（管理员）"""
    user = db.query(User).filter(User.id == user_id).first()
    if not user:
        raise HTTPException(status_code=404, detail=f"用户 {user_id} 不存在")

    transactions = BillingService.get_transactions(db, user_id, skip, limit)
    return [TransactionResponse(
        id=t.id,
        type=t.type,
        amount=t.amount,
        balance_before=t.balance_before,
        balance_after=t.balance_after,
        description=t.description,
        created_at=t.created_at
    ) for t in transactions]


@router.get("/admin/user-pricing/{user_id}", response_model=UserPricingResponse)
async def get_user_pricing(
    user_id: int,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """获取用户的定价配置（管理员）"""
    user = db.query(User).filter(User.id == user_id).first()
    if not user:
        raise HTTPException(status_code=404, detail=f"用户 {user_id} 不存在")

    global_price = BillingService.get_cpu_hour_price(db)
    effective_price = BillingService.get_cpu_hour_price(db, user_id)

    return UserPricingResponse(
        user_id=user.id,
        username=user.username,
        email=user.email,
        custom_cpu_hour_price=user.custom_cpu_hour_price,
        global_price=global_price,
        effective_price=effective_price,
        price_updated_at=user.price_updated_at.isoformat() if user.price_updated_at else None,
        price_updated_by=user.price_updated_by
    )


@router.put("/admin/user-pricing/{user_id}")
async def set_user_pricing(
    user_id: int,
    config: UserPricingConfig,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """设置用户的自定义定价（管理员）"""
    success, message = BillingService.set_user_cpu_hour_price(
        db=db,
        user_id=user_id,
        price=config.custom_cpu_hour_price,
        admin_id=admin.id
    )

    if not success:
        raise HTTPException(status_code=400, detail=message)

    # 返回更新后的定价信息
    user = db.query(User).filter(User.id == user_id).first()
    global_price = BillingService.get_cpu_hour_price(db)
    effective_price = BillingService.get_cpu_hour_price(db, user_id)

    return {
        "success": True,
        "message": message,
        "pricing": UserPricingResponse(
            user_id=user.id,
            username=user.username,
            email=user.email,
            custom_cpu_hour_price=user.custom_cpu_hour_price,
            global_price=global_price,
            effective_price=effective_price,
            price_updated_at=user.price_updated_at.isoformat() if user.price_updated_at else None,
            price_updated_by=user.price_updated_by
        )
    }


@router.delete("/admin/user-pricing/{user_id}")
async def delete_user_pricing(
    user_id: int,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """删除用户的自定义定价，恢复为角色默认定价（管理员）"""
    success, message = BillingService.set_user_cpu_hour_price(
        db=db,
        user_id=user_id,
        price=None,  # 设置为 None 表示删除自定义定价
        admin_id=admin.id
    )

    if not success:
        raise HTTPException(status_code=400, detail=message)

    return {"success": True, "message": message}


class BatchUserPricingConfig(BaseModel):
    """批量用户定价配置"""
    user_ids: List[int]
    custom_cpu_hour_price: Optional[float] = None  # None 表示恢复全局定价


@router.put("/admin/user-pricing/batch")
async def batch_set_user_pricing(
    config: BatchUserPricingConfig,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """批量设置用户定价（管理员）"""
    results = []
    success_count = 0
    fail_count = 0

    for user_id in config.user_ids:
        try:
            success, message = BillingService.set_user_cpu_hour_price(
                db=db,
                user_id=user_id,
                price=config.custom_cpu_hour_price,
                admin_id=admin.id
            )
            if success:
                success_count += 1
                results.append({"user_id": user_id, "success": True, "message": message})
            else:
                fail_count += 1
                results.append({"user_id": user_id, "success": False, "message": message})
        except Exception as e:
            fail_count += 1
            results.append({"user_id": user_id, "success": False, "message": str(e)})

    return {
        "success": fail_count == 0,
        "message": f"成功更新 {success_count} 个用户，失败 {fail_count} 个",
        "success_count": success_count,
        "fail_count": fail_count,
        "results": results
    }


@router.get("/admin/orders", response_model=List[OrderResponse])
async def get_all_orders(
    skip: int = 0,
    limit: int = 50,
    status: Optional[str] = None,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """获取所有充值订单（管理员）"""
    query = db.query(RechargeOrder)

    if status:
        query = query.filter(RechargeOrder.payment_status == status)

    orders = query.order_by(RechargeOrder.created_at.desc()).offset(skip).limit(limit).all()

    return [OrderResponse(
        id=o.id,
        order_no=o.order_no,
        amount=o.amount,
        cpu_hours=o.cpu_hours,
        price_per_hour=o.price_per_hour,
        payment_method=o.payment_method,  # 已经是字符串
        payment_status=o.payment_status,  # 已经是字符串
        created_at=o.created_at,
        paid_at=o.paid_at
    ) for o in orders]


class RechargeRecordResponse(BaseModel):
    """充值记录响应"""
    id: int
    user_id: int
    username: str
    amount: float
    cpu_hours: float
    created_at: datetime

    class Config:
        from_attributes = True


@router.get("/admin/recharge-records", response_model=List[RechargeRecordResponse])
async def get_recharge_records(
    skip: int = 0,
    limit: int = 50,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """获取充值记录（管理员）"""
    orders = db.query(RechargeOrder).filter(
        RechargeOrder.payment_status == 'paid'
    ).order_by(RechargeOrder.created_at.desc()).offset(skip).limit(limit).all()

    result = []
    for order in orders:
        user = db.query(User).filter(User.id == order.user_id).first()
        if user:
            result.append(RechargeRecordResponse(
                id=order.id,
                user_id=order.user_id,
                username=user.username,
                amount=order.amount,
                cpu_hours=order.cpu_hours,
                created_at=order.created_at
            ))

    return result


@router.get("/admin/pricing-config")
async def get_admin_pricing_config(
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """获取定价配置（包括全局和角色定价）"""
    global_price = float(BillingService.get_config(db, 'cpu_hour_price', '0.1'))

    # 获取角色定价
    role_prices = {
        'ADMIN': float(BillingService.get_config(db, 'role_price_admin', '0.1')),
        'PREMIUM': float(BillingService.get_config(db, 'role_price_premium', '0.08')),
        'USER': float(BillingService.get_config(db, 'role_price_user', '0.1')),
        'GUEST': float(BillingService.get_config(db, 'role_price_guest', '0.15')),
    }

    # 获取自定义定价
    custom_prices = {}
    users = db.query(User).filter(User.custom_cpu_hour_price.isnot(None)).all()
    for user in users:
        custom_prices[str(user.id)] = user.custom_cpu_hour_price

    return {
        'global_price': global_price,
        'role_prices': role_prices,
        'custom_prices': custom_prices,
    }


class RolePricingUpdate(BaseModel):
    """角色定价更新"""
    role_prices: Dict[str, float] = Field(..., description="角色定价字典")


@router.put("/admin/role-pricing")
async def update_role_pricing(
    pricing: RolePricingUpdate,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """批量更新角色定价"""
    try:
        for role, price in pricing.role_prices.items():
            config_key = f'role_price_{role.lower()}'
            BillingService.set_config(
                db,
                config_key,
                str(price),
                f'{role}角色的核时单价（元/核时）',
                admin.id
            )

        return {
            "success": True,
            "message": "角色定价已更新"
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"更新失败: {str(e)}")


class UserPriceUpdate(BaseModel):
    """用户自定义价格更新"""
    user_id: int = Field(..., description="用户ID")
    price: Optional[float] = Field(None, ge=0.001, description="核时单价（元/核时），为null时删除自定义价格")


@router.put("/admin/user-price")
async def update_user_price(
    request: UserPriceUpdate,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """设置用户的自定义核时单价（管理员）"""
    try:
        success, message = BillingService.set_user_cpu_hour_price(
            db=db,
            user_id=request.user_id,
            price=request.price,
            admin_id=admin.id
        )

        if not success:
            raise HTTPException(status_code=400, detail=message)

        return {
            "success": True,
            "message": message
        }
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"更新失败: {str(e)}")


class ConsumptionRecord(BaseModel):
    """消费记录"""
    id: int
    user_id: int
    username: str
    amount: float
    cpu_hours: float
    created_at: datetime

    class Config:
        from_attributes = True


@router.get("/admin/consumption-records")
async def get_consumption_records(
    skip: int = 0,
    limit: int = 50,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """获取消费汇总记录（按用户聚合）"""
    # 查询所有消费交易
    transactions = db.query(QuotaTransaction).filter(
        QuotaTransaction.type == 'consume'
    ).order_by(QuotaTransaction.created_at.desc()).all()

    # 按用户聚合消费数据
    consumption_by_user = {}
    for t in transactions:
        if t.user_id not in consumption_by_user:
            consumption_by_user[t.user_id] = {
                'user_id': t.user_id,
                'total_consumption': 0.0,
                'total_cpu_hours': 0.0,
                'last_consumption_at': t.created_at,
            }
        consumption_by_user[t.user_id]['total_consumption'] += abs(t.amount)
        consumption_by_user[t.user_id]['total_cpu_hours'] += abs(t.amount)
        if t.created_at > consumption_by_user[t.user_id]['last_consumption_at']:
            consumption_by_user[t.user_id]['last_consumption_at'] = t.created_at

    # 转换为列表并应用分页
    result = list(consumption_by_user.values())
    result.sort(key=lambda x: x['last_consumption_at'], reverse=True)

    return result[skip:skip + limit]


@router.get("/admin/consumption-details", response_model=List[ConsumptionDetailResponse])
async def get_consumption_details(
    skip: int = 0,
    limit: int = 100,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """获取详细消费记录（管理员）"""
    transactions = db.query(QuotaTransaction).filter(
        QuotaTransaction.type == 'consume'
    ).order_by(QuotaTransaction.created_at.desc()).offset(skip).limit(limit).all()

    result = []
    for t in transactions:
        user = db.query(User).filter(User.id == t.user_id).first()
        if user:
            cpu_hours = abs(t.amount)  # 消费是负数，转为正数显示
            # 计算金额：核时 × 当时的价格（使用当前价格作为近似）
            price = BillingService.get_cpu_hour_price(db, t.user_id)
            money_amount = cpu_hours * price

            result.append(ConsumptionDetailResponse(
                id=t.id,
                user_id=t.user_id,
                username=user.username,
                cpu_hours=cpu_hours,
                amount=cpu_hours,  # 保持兼容性
                money_amount=round(money_amount, 2),  # 新增金额字段
                price_per_hour=price,  # 新增单价字段
                created_at=t.created_at
            ))

    return result

