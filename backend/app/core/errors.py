"""
统一错误处理系统

提供标准化的错误码、错误类和错误响应格式
"""
from enum import Enum
from typing import Dict, Any, Optional
from fastapi import HTTPException
from fastapi.responses import JSONResponse


class ErrorCode(str, Enum):
    """标准错误码"""
    
    # ========== 客户端错误 (4xx) ==========
    
    # 输入验证错误
    INVALID_INPUT = "INVALID_INPUT"
    INVALID_SMILES = "INVALID_SMILES"
    INVALID_CONFIG = "INVALID_CONFIG"
    MISSING_REQUIRED_FIELD = "MISSING_REQUIRED_FIELD"
    
    # 资源错误
    NOT_FOUND = "NOT_FOUND"
    JOB_NOT_FOUND = "JOB_NOT_FOUND"
    USER_NOT_FOUND = "USER_NOT_FOUND"
    ELECTROLYTE_NOT_FOUND = "ELECTROLYTE_NOT_FOUND"
    
    # 认证授权错误
    UNAUTHORIZED = "UNAUTHORIZED"
    INVALID_CREDENTIALS = "INVALID_CREDENTIALS"
    TOKEN_EXPIRED = "TOKEN_EXPIRED"
    FORBIDDEN = "FORBIDDEN"
    INSUFFICIENT_PERMISSIONS = "INSUFFICIENT_PERMISSIONS"
    
    # 业务逻辑错误
    CONFLICT = "CONFLICT"
    DUPLICATE_RESOURCE = "DUPLICATE_RESOURCE"
    QUOTA_EXCEEDED = "QUOTA_EXCEEDED"
    INSUFFICIENT_BALANCE = "INSUFFICIENT_BALANCE"
    JOB_ALREADY_SUBMITTED = "JOB_ALREADY_SUBMITTED"
    JOB_CANNOT_BE_CANCELLED = "JOB_CANNOT_BE_CANCELLED"
    
    # 资源保护错误
    CONCURRENT_LIMIT_EXCEEDED = "CONCURRENT_LIMIT_EXCEEDED"
    DAILY_LIMIT_EXCEEDED = "DAILY_LIMIT_EXCEEDED"
    STORAGE_QUOTA_EXCEEDED = "STORAGE_QUOTA_EXCEEDED"
    
    # 状态错误
    INVALID_STATE = "INVALID_STATE"
    JOB_NOT_READY = "JOB_NOT_READY"
    
    # ========== 服务端错误 (5xx) ==========
    
    INTERNAL_ERROR = "INTERNAL_ERROR"
    DATABASE_ERROR = "DATABASE_ERROR"
    EXTERNAL_SERVICE_ERROR = "EXTERNAL_SERVICE_ERROR"
    SLURM_ERROR = "SLURM_ERROR"
    FILE_SYSTEM_ERROR = "FILE_SYSTEM_ERROR"
    CALCULATION_ERROR = "CALCULATION_ERROR"


class APIError(HTTPException):
    """统一的API错误类"""
    
    def __init__(
        self,
        code: ErrorCode,
        message: str,
        status_code: int = 400,
        details: Optional[Dict[str, Any]] = None
    ):
        """
        创建API错误
        
        Args:
            code: 错误码
            message: 用户友好的错误消息
            status_code: HTTP状态码
            details: 额外的错误详情
        """
        self.code = code
        self.message = message
        self.details = details or {}
        
        super().__init__(
            status_code=status_code,
            detail={
                "code": code.value,
                "message": message,
                "details": self.details
            }
        )


# ========== 便捷的错误创建函数 ===========

def invalid_input_error(message: str, field: str = None, details: Dict[str, Any] = None) -> APIError:
    """创建输入验证错误"""
    error_details = details or {}
    if field:
        error_details["field"] = field
    return APIError(ErrorCode.INVALID_INPUT, message, 400, error_details)


def already_exists_error(resource: str, details: Dict[str, Any] = None) -> APIError:
    """创建资源已存在错误"""
    message = f"{resource}已存在"
    return APIError(ErrorCode.DUPLICATE_RESOURCE, message, 400, details)


def not_found_error(resource: str, resource_id: Any = None) -> APIError:
    """创建资源未找到错误"""
    message = f"{resource}未找到"
    details = {"resource_id": resource_id} if resource_id else None
    return APIError(ErrorCode.NOT_FOUND, message, 404, details)


def unauthorized_error(message: str = "未授权，请先登录", details: Dict[str, Any] = None) -> APIError:
    """创建未授权错误"""
    return APIError(ErrorCode.UNAUTHORIZED, message, 401, details)


def forbidden_error(message: str = "权限不足") -> APIError:
    """创建权限不足错误"""
    return APIError(ErrorCode.FORBIDDEN, message, 403)


def quota_exceeded_error(quota_type: str, current: float, limit: float) -> APIError:
    """创建配额超限错误"""
    return APIError(
        ErrorCode.QUOTA_EXCEEDED,
        f"{quota_type}配额不足",
        400,
        {"current": current, "limit": limit}
    )


def internal_error(message: str = "服务器内部错误", details: Dict[str, Any] = None) -> APIError:
    """创建内部错误"""
    return APIError(ErrorCode.INTERNAL_ERROR, message, 500, details)


def database_error(operation: str, error: Exception) -> APIError:
    """创建数据库错误"""
    return APIError(
        ErrorCode.DATABASE_ERROR,
        f"数据库{operation}失败",
        500,
        {"error": str(error)}
    )


def slurm_error(message: str, details: Dict[str, Any] = None) -> APIError:
    """创建Slurm错误"""
    return APIError(ErrorCode.SLURM_ERROR, message, 500, details)


def concurrent_limit_error(current: int, limit: int) -> APIError:
    """创建并发限制超出错误"""
    return APIError(
        ErrorCode.CONCURRENT_LIMIT_EXCEEDED,
        f"并发任务数已达上限({current}/{limit})，请等待任务完成",
        400,
        {"current": current, "limit": limit}
    )


def daily_limit_error(current: int, limit: int) -> APIError:
    """创建每日限制超出错误"""
    return APIError(
        ErrorCode.DAILY_LIMIT_EXCEEDED,
        f"今日任务数已达上限({current}/{limit})，明日0点后可继续",
        400,
        {"current": current, "limit": limit}
    )


def insufficient_balance_error(required: float, available: float) -> APIError:
    """创建余额不足错误"""
    return APIError(
        ErrorCode.INSUFFICIENT_BALANCE,
        f"余额不足。需要{required:.2f}核时，可用{available:.2f}核时",
        400,
        {"required": required, "available": available}
    )


# ========== 错误处理装饰器 ==========

from functools import wraps
import logging

logger = logging.getLogger(__name__)


def handle_errors(func):
    """
    错误处理装饰器
    
    自动捕获并转换异常为APIError
    """
    @wraps(func)
    async def async_wrapper(*args, **kwargs):
        try:
            return await func(*args, **kwargs)
        except APIError:
            # APIError直接抛出
            raise
        except ValueError as e:
            # 值错误转换为输入验证错误
            raise invalid_input_error(str(e))
        except Exception as e:
            # 其他异常记录并转换为内部错误
            logger.exception(f"Unexpected error in {func.__name__}")
            raise internal_error(
                "操作失败，请稍后重试",
                {"function": func.__name__, "error": str(e)}
            )
    
    @wraps(func)
    def sync_wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except APIError:
            raise
        except ValueError as e:
            raise invalid_input_error(str(e))
        except Exception as e:
            logger.exception(f"Unexpected error in {func.__name__}")
            raise internal_error(
                "操作失败，请稍后重试",
                {"function": func.__name__, "error": str(e)}
            )
    
    # 根据函数类型返回合适的wrapper
    import inspect
    if inspect.iscoroutinefunction(func):
        return async_wrapper
    else:
        return sync_wrapper


# ========== 用户友好的错误消息映射 ==========

ERROR_MESSAGES = {
    ErrorCode.INVALID_INPUT: "输入参数不正确，请检查后重试",
    ErrorCode.INVALID_SMILES: "SMILES格式不正确，请检查分子结构",
    ErrorCode.NOT_FOUND: "请求的资源不存在",
    ErrorCode.JOB_NOT_FOUND: "任务不存在或已被删除",
    ErrorCode.UNAUTHORIZED: "请先登录",
    ErrorCode.FORBIDDEN: "您没有权限执行此操作",
    ErrorCode.QUOTA_EXCEEDED: "配额不足，请充值或联系管理员",
    ErrorCode.INSUFFICIENT_BALANCE: "余额不足，请充值",
    ErrorCode.JOB_ALREADY_SUBMITTED: "任务已提交，无法重复提交",
    ErrorCode.INTERNAL_ERROR: "服务器错误，请稍后重试",
    ErrorCode.DATABASE_ERROR: "数据库操作失败，请稍后重试",
    ErrorCode.SLURM_ERROR: "集群调度失败，请联系管理员",
}


def get_user_friendly_message(code: ErrorCode, default: str = None) -> str:
    """获取用户友好的错误消息"""
    return ERROR_MESSAGES.get(code, default or "操作失败")
