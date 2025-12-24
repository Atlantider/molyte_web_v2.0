"""
充值套餐管理 API
"""
from typing import List, Optional
from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy.orm import Session
from pydantic import BaseModel, Field

from app.database import get_db
from app.dependencies import get_current_admin_user, get_optional_current_user
from app.models.user import User
from app.models.recharge_package import RechargePackage, PackageDiscount, UserTypeEnum

router = APIRouter(prefix="/recharge-packages", tags=["recharge-packages"])


# ============ Schemas ============

class RechargePackageResponse(BaseModel):
    """充值套餐响应"""
    id: int
    name: str
    description: Optional[str]
    user_type: str
    price: float
    cpu_hours: float
    display_order: int
    badge: Optional[str]
    color: str
    icon: Optional[str]
    is_active: bool

    class Config:
        from_attributes = True


class RechargePackageCreate(BaseModel):
    """创建充值套餐"""
    name: str = Field(..., min_length=1, max_length=100)
    description: Optional[str] = None
    user_type: str = Field(..., description="用户类型: STUDENT, RESEARCHER, COMPANY, TEAM")
    price: float = Field(..., gt=0)
    cpu_hours: float = Field(..., gt=0)
    display_order: int = 0
    badge: Optional[str] = None
    color: str = "#1890ff"
    icon: Optional[str] = None


class RechargePackageUpdate(BaseModel):
    """更新充值套餐"""
    name: Optional[str] = None
    description: Optional[str] = None
    price: Optional[float] = None
    cpu_hours: Optional[float] = None
    display_order: Optional[int] = None
    badge: Optional[str] = None
    color: Optional[str] = None
    icon: Optional[str] = None
    is_active: Optional[bool] = None


# ============ Public Endpoints ============

@router.get("/", response_model=List[RechargePackageResponse])
async def get_packages(
    user_type: Optional[str] = None,
    db: Session = Depends(get_db),
    current_user: Optional[User] = Depends(get_optional_current_user)
):
    """
    获取充值套餐列表

    如果提供了 user_type，则只返回该类型的套餐
    如果没有提供，则返回当前用户类型的套餐
    如果未登录，则返回所有套餐
    """
    query = db.query(RechargePackage).filter(RechargePackage.is_active == True)

    # 确定要查询的用户类型
    target_user_type = user_type
    if not target_user_type and current_user:
        target_user_type = current_user.user_type

    if target_user_type:
        query = query.filter(RechargePackage.user_type == target_user_type)

    # 按显示顺序排序
    packages = query.order_by(RechargePackage.display_order).all()
    return packages


@router.get("/{package_id}", response_model=RechargePackageResponse)
async def get_package(
    package_id: int,
    db: Session = Depends(get_db)
):
    """获取单个充值套餐详情"""
    package = db.query(RechargePackage).filter(
        RechargePackage.id == package_id,
        RechargePackage.is_active == True
    ).first()
    
    if not package:
        raise HTTPException(status_code=404, detail="Package not found")
    
    return package


# ============ Admin Endpoints ============

@router.post("/", response_model=RechargePackageResponse, status_code=status.HTTP_201_CREATED)
async def create_package(
    package_data: RechargePackageCreate,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """创建充值套餐（管理员）"""
    # 验证用户类型
    try:
        UserTypeEnum(package_data.user_type)
    except ValueError:
        raise HTTPException(
            status_code=400,
            detail=f"Invalid user_type. Must be one of: {', '.join([e.value for e in UserTypeEnum])}"
        )
    
    package = RechargePackage(**package_data.dict())
    db.add(package)
    db.commit()
    db.refresh(package)
    return package


@router.put("/{package_id}", response_model=RechargePackageResponse)
async def update_package(
    package_id: int,
    package_data: RechargePackageUpdate,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """更新充值套餐（管理员）"""
    package = db.query(RechargePackage).filter(RechargePackage.id == package_id).first()
    
    if not package:
        raise HTTPException(status_code=404, detail="Package not found")
    
    # 更新字段
    update_data = package_data.dict(exclude_unset=True)
    for field, value in update_data.items():
        setattr(package, field, value)
    
    db.commit()
    db.refresh(package)
    return package


@router.delete("/{package_id}", status_code=status.HTTP_204_NO_CONTENT)
async def delete_package(
    package_id: int,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """删除充值套餐（管理员）"""
    package = db.query(RechargePackage).filter(RechargePackage.id == package_id).first()
    
    if not package:
        raise HTTPException(status_code=404, detail="Package not found")
    
    db.delete(package)
    db.commit()


@router.get("/admin/all", response_model=List[RechargePackageResponse])
async def get_all_packages_admin(
    user_type: Optional[str] = None,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """获取所有充值套餐（包括已禁用的）（管理员）"""
    query = db.query(RechargePackage)
    
    if user_type:
        query = query.filter(RechargePackage.user_type == user_type)
    
    packages = query.order_by(RechargePackage.display_order).all()
    return packages

