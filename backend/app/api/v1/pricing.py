"""
定价管理API
Pricing Management API

提供任务类型定价和用户折扣的查询和管理接口
"""
from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy.orm import Session
from sqlalchemy import text
from typing import List, Optional
from pydantic import BaseModel, Field

from app.database import get_db
from app.dependencies import get_current_user, get_current_admin_user
from app.models.user import User
from app.services.pricing_service import PricingService
from app.core.logger import logger


router = APIRouter(prefix="/pricing", tags=["pricing"])


# ============================================================================
# Schemas
# ============================================================================

class TaskTypePriceSchema(BaseModel):
    """任务类型单价"""
    task_type: str = Field(..., description="任务类型")
    price_per_hour: float = Field(..., description="核时单价(元/核时)")
    
    class Config:
        from_attributes = True


class UserTypeDiscountSchema(BaseModel):
    """用户类型折扣"""
    user_type: str = Field(..., description="用户类型")
    discount_rate: float = Field(..., description="折扣率(0.9表示9折)")
    
    class Config:
        from_attributes = True


class UserTypePriceSchema(BaseModel):
    """用户类型核时单价"""
    user_type: str = Field(..., description="用户类型: GUEST/USER/PREMIUM/ADMIN")
    core_hour_price: float = Field(..., description="核时单价(元/核时)")
    
    class Config:
        from_attributes = True


class UserDiscountInfoSchema(BaseModel):
    """用户折扣信息"""
    discount_rate: float = Field(..., description="折扣率")
    discount_source: str = Field(..., description="折扣来源: custom/user_type/none")
    user_type: Optional[str] = Field(None, description="用户类型")


class PriceCalculationRequest(BaseModel):
    """价格计算请求"""
    task_type: str = Field(..., description="任务类型")
    cpu_hours: float = Field(..., description="核时数")


class PriceCalculationResponse(BaseModel):
    """价格计算响应"""
    task_type: str
    cpu_hours: float
    task_price: float = Field(..., description="任务类型单价")
    discount_rate: float = Field(..., description="用户折扣率")
    final_price: float = Field(..., description="最终价格")


class UpdateTaskTypePriceRequest(BaseModel):
    """更新任务类型单价请求"""
    price_per_hour: float = Field(..., gt=0, description="核时单价(元/核时)")


class UpdateUserTypePriceRequest(BaseModel):
    """更新用户类型核时单价请求"""
    core_hour_price: float = Field(..., ge=0, description="核时单价(元/核时)")


class SetUserCustomDiscountRequest(BaseModel):
    """设置用户自定义折扣请求"""
    discount_rate: Optional[float] = Field(None, gt=0, le=1, description="折扣率(None表示删除)")


class BillingConfigSchema(BaseModel):
    """计费配置"""
    pricing_mode: str = Field(..., description="计费模式: CORE_HOUR/TASK_TYPE")
    global_core_hour_price: float = Field(..., description="全局核时单价")


class UpdateBillingConfigRequest(BaseModel):
    """更新计费配置请求"""
    pricing_mode: Optional[str] = Field(None, description="计费模式: CORE_HOUR/TASK_TYPE")
    global_core_hour_price: Optional[float] = Field(None, description="全局核时单价")


# ============================================================================
# 用户端点
# ============================================================================

@router.get("/task-types", response_model=List[TaskTypePriceSchema])
async def get_task_type_prices(
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """获取所有任务类型单价"""
    try:
        result = db.execute(
            "SELECT task_type, price_per_hour FROM task_type_pricing WHERE is_active = TRUE ORDER BY task_type"
        ).fetchall()
        
        return [
            TaskTypePriceSchema(task_type=row[0], price_per_hour=row[1])
            for row in result
        ]
    except Exception as e:
        logger.error(f"Failed to get task type prices: {e}")
        # 返回默认值
        return [
            TaskTypePriceSchema(task_type=k, price_per_hour=v)
            for k, v in PricingService.DEFAULT_TASK_TYPE_PRICES.items()
        ]


@router.get("/my-discount", response_model=UserDiscountInfoSchema)
async def get_my_discount(
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """获取当前用户的折扣信息"""
    discount_rate = PricingService.get_user_discount_rate(db, current_user.id)
    
    # 判断折扣来源
    if current_user.custom_discount_rate is not None:
        source = "custom"
    elif discount_rate < 1.0:
        source = "user_type"
    else:
        source = "none"
    
    user_type_value = current_user.user_type.value if hasattr(current_user.user_type, 'value') else str(current_user.user_type)
    
    return UserDiscountInfoSchema(
        discount_rate=discount_rate,
        discount_source=source,
        user_type=user_type_value
    )


@router.post("/calculate", response_model=PriceCalculationResponse)
async def calculate_price(
    request: PriceCalculationRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """计算任务预估价格"""
    task_price = PricingService.get_task_type_price(db, request.task_type)
    discount_rate = PricingService.get_user_discount_rate(db, current_user.id)
    final_price = PricingService.calculate_final_price(
        db, request.task_type, request.cpu_hours, current_user.id
    )
    
    return PriceCalculationResponse(
        task_type=request.task_type,
        cpu_hours=request.cpu_hours,
        task_price=task_price,
        discount_rate=discount_rate,
        final_price=final_price
    )


@router.get("/config", response_model=BillingConfigSchema)
async def get_billing_config(
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """获取计费配置（用户端）"""
    try:
        result = db.execute(
            text("SELECT pricing_mode, global_core_hour_price FROM billing_config ORDER BY id DESC LIMIT 1")
        ).first()
        
        if result:
            return BillingConfigSchema(
                pricing_mode=result[0],
                global_core_hour_price=result[1]
            )
        else:
            # 返回默认配置
            return BillingConfigSchema(
                pricing_mode="CORE_HOUR",
                global_core_hour_price=0.1
            )
    except Exception as e:
        logger.error(f"Failed to get billing config: {e}")
        # 返回默认配置
        return BillingConfigSchema(
            pricing_mode="CORE_HOUR",
            global_core_hour_price=0.1
        )


# ============================================================================
# 管理员端点
# ============================================================================

@router.get("/admin/task-types", response_model=List[TaskTypePriceSchema])
async def admin_get_task_type_prices(
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_admin_user)
):
    """管理员获取所有任务类型单价(包括未激活的)"""
    try:
        result = db.execute(
            text("SELECT task_type, price_per_hour FROM task_type_pricing ORDER BY task_type")
        ).fetchall()
        
        return [
            TaskTypePriceSchema(task_type=row[0], price_per_hour=row[1])
            for row in result
        ]
    except Exception as e:
        logger.error(f"Failed to get task type prices: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.put("/admin/task-types/{task_type}")
async def admin_update_task_type_price(
    task_type: str,
    request: UpdateTaskTypePriceRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_admin_user)
):
    """管理员更新任务类型单价"""
    success, message = PricingService.admin_update_task_type_price(
        db, task_type, request.price_per_hour, current_user.id
    )
    
    if not success:
        raise HTTPException(status_code=400, detail=message)
    
    return {"success": True, "message": message}


@router.get("/admin/user-type-prices", response_model=List[UserTypePriceSchema])
async def admin_get_user_type_prices(
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_admin_user)
):
    """管理员获取所有用户类型核时单价"""
    try:
        result = db.execute(
            text("SELECT user_type, core_hour_price FROM user_type_prices ORDER BY user_type")
        ).fetchall()
        
        return [
            UserTypePriceSchema(user_type=row[0], core_hour_price=row[1])
            for row in result
        ]
    except Exception as e:
        logger.error(f"Failed to get user type prices: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.put("/admin/user-type-prices/{user_type}")
async def admin_update_user_type_price(
    user_type: str,
    request: UpdateUserTypePriceRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_admin_user)
):
    """管理员更新用户类型核时单价"""
    success, message = PricingService.admin_update_user_type_price(
        db, user_type, request.core_hour_price, current_user.id
    )
    
    if not success:
        raise HTTPException(status_code=400, detail=message)
    
    return {"success": True, "message": message}


@router.put("/admin/users/{user_id}/custom-discount")
async def admin_set_user_custom_discount(
    user_id: int,
    request: SetUserCustomDiscountRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_admin_user)
):
    """管理员设置特定用户自定义折扣"""
    success, message = PricingService.admin_set_user_custom_discount(
        db, user_id, request.discount_rate, current_user.id
    )
    
    if not success:
        raise HTTPException(status_code=400, detail=message)
    
    return {"success": True, "message": message}


@router.get("/admin/config", response_model=BillingConfigSchema)
async def admin_get_billing_config(
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_admin_user)
):
    """管理员获取计费配置"""
    try:
        result = db.execute(
            text("SELECT pricing_mode FROM billing_config ORDER BY id DESC LIMIT 1")
        ).first()
        
        if result:
            return BillingConfigSchema(
                pricing_mode=result[0]
            )
        else:
            # 返回默认配置
            return BillingConfigSchema(
                pricing_mode="CORE_HOUR"
            )
    except Exception as e:
        logger.error(f"Admin failed to get billing config: {e}")
        # 返回默认配置
        return BillingConfigSchema(
            pricing_mode="CORE_HOUR",
            global_core_hour_price=0.1
        )


@router.put("/admin/config")
async def admin_update_billing_config(
    config: UpdateBillingConfigRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_admin_user)
):
    """管理员更新计费配置"""
    try:
        # 检查billing_config表是否有数据
        result = db.execute(text("SELECT id FROM billing_config LIMIT 1")).first()
        
        update_fields = []
        if config.pricing_mode is not None:
            update_fields.append(f"pricing_mode = '{config.pricing_mode}'")
        if config.global_core_hour_price is not None:
            update_fields.append(f"global_core_hour_price = {config.global_core_hour_price}")
        
        if not update_fields:
            raise HTTPException(status_code=400, detail="没有要更新的字段")
        
        update_fields.append(f"updated_by = {current_user.id}")
        update_fields.append("updated_at = NOW()")
        
        if result:
            # 更新现有记录
            sql = text(f"UPDATE billing_config SET {', '.join(update_fields)} WHERE id = {result[0]}")
        else:
            # 插入新记录
            sql = text(f"""
                INSERT INTO billing_config (pricing_mode, updated_by)
                VALUES (
                    '{config.pricing_mode or 'CORE_HOUR'}',
                    {current_user.id}
                )
            """)
        
        db.execute(sql)
        db.commit()
        
        logger.info(f"Admin {current_user.id} updated billing config: {config}")
        
        return {
            "success": True,
            "message": "计费配置已更新"
        }
    except Exception as e:
        db.rollback()
        logger.error(f"Failed to update billing config: {e}")
        raise HTTPException(status_code=500, detail=f"更新失败: {str(e)}")
