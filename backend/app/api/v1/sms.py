"""
短信验证 API
"""
from fastapi import APIRouter, HTTPException, Depends
from pydantic import BaseModel, Field, validator
from sqlalchemy.orm import Session
import re

from app.database import get_db
from app.models.user import User
from app.services.sms import sms_service
from app.dependencies import get_current_user
from app.core.logger import logger

router = APIRouter(prefix="/sms", tags=["sms"])


class SendCodeRequest(BaseModel):
    """发送验证码请求"""
    phone: str = Field(..., description="手机号")
    purpose: str = Field(default="register", description="用途: register/bind/reset")
    
    @validator('phone')
    def validate_phone(cls, v):
        # 中国大陆手机号验证
        if not re.match(r'^1[3-9]\d{9}$', v):
            raise ValueError('请输入有效的手机号')
        return v


class VerifyCodeRequest(BaseModel):
    """验证验证码请求"""
    phone: str = Field(..., description="手机号")
    code: str = Field(..., min_length=4, max_length=6, description="验证码")
    purpose: str = Field(default="register", description="用途")


class BindPhoneRequest(BaseModel):
    """绑定手机号请求"""
    phone: str = Field(..., description="手机号")
    code: str = Field(..., description="验证码")


@router.post("/send-code")
async def send_verification_code(
    request: SendCodeRequest,
    db: Session = Depends(get_db)
):
    """发送短信验证码"""
    # 检查手机号是否已被注册（注册场景）
    if request.purpose == "register":
        existing_user = db.query(User).filter(User.phone == request.phone).first()
        if existing_user:
            raise HTTPException(status_code=400, detail="该手机号已被注册")
    
    # 发送验证码
    success, message = await sms_service.send_verification_code(
        phone=request.phone,
        purpose=request.purpose
    )
    
    if not success:
        raise HTTPException(status_code=400, detail=message)
    
    return {
        "success": True,
        "message": message,
        # 测试模式下返回验证码（生产环境应移除）
        "code": sms_service.get_code_for_testing(request.phone, request.purpose) if hasattr(sms_service, 'get_code_for_testing') else None
    }


@router.post("/verify-code")
async def verify_code(request: VerifyCodeRequest):
    """验证短信验证码"""
    success, message = sms_service.verify_code(
        phone=request.phone,
        code=request.code,
        purpose=request.purpose
    )
    
    if not success:
        raise HTTPException(status_code=400, detail=message)
    
    return {"success": True, "message": message}


@router.post("/bind-phone")
async def bind_phone(
    request: BindPhoneRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """绑定手机号到当前用户"""
    # 检查手机号是否已被其他用户绑定
    existing_user = db.query(User).filter(
        User.phone == request.phone,
        User.id != current_user.id
    ).first()
    if existing_user:
        raise HTTPException(status_code=400, detail="该手机号已被其他用户绑定")
    
    # 验证验证码
    success, message = sms_service.verify_code(
        phone=request.phone,
        code=request.code,
        purpose="bind"
    )
    
    if not success:
        raise HTTPException(status_code=400, detail=message)
    
    # 绑定手机号
    current_user.phone = request.phone
    current_user.phone_verified = True
    db.commit()
    
    logger.info(f"用户 {current_user.username} 绑定手机号 {request.phone}")
    
    return {
        "success": True,
        "message": "手机号绑定成功",
        "phone": request.phone
    }

