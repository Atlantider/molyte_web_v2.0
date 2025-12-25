"""
全局异常处理器

注册到FastAPI应用的全局异常处理
"""
from fastapi import Request, status
from fastapi.responses import JSONResponse
from fastapi.exceptions import RequestValidationError
from sqlalchemy.exc import SQLAlchemyError
import logging

from app.core.errors import APIError, ErrorCode, internal_error, database_error

logger = logging.getLogger(__name__)


async def api_error_handler(request: Request, exc: APIError) -> JSONResponse:
    """
    处理APIError异常
    
    返回标准化的错误响应格式
    """
    return JSONResponse(
        status_code=exc.status_code,
        content={
            "code": exc.code.value,
            "message": exc.message,
            "details": exc.details
        }
    )


async def validation_error_handler(request: Request, exc: RequestValidationError) -> JSONResponse:
    """
    处理Pydantic验证错误
    
    转换为统一的错误格式
    """
    errors = exc.errors()
    
    # 提取第一个错误作为主要消息
    first_error = errors[0] if errors else {}
    field = ".".join(str(loc) for loc in first_error.get("loc", []))
    msg = first_error.get("msg", "验证失败")
    
    return JSONResponse(
        status_code=status.HTTP_400_BAD_REQUEST,
        content={
            "code": ErrorCode.INVALID_INPUT.value,
            "message": f"参数 '{field}' {msg}",
            "details": {
                "validation_errors": errors
            }
        }
    )


async def database_error_handler(request: Request, exc: SQLAlchemyError) -> JSONResponse:
    """
    处理数据库错误
    """
    logger.exception("Database error occurred")
    
    error = database_error("操作", exc)
    
    return JSONResponse(
        status_code=500,
        content={
            "code": error.code.value,
            "message": error.message,
            "details": error.details
        }
    )


async def general_exception_handler(request: Request, exc: Exception) -> JSONResponse:
    """
    处理未捕获的异常
    
    防止敏感信息泄露
    """
    logger.exception(f"Unhandled exception: {exc}")
    
    error = internal_error()
    
    return JSONResponse(
        status_code=500,
        content={
            "code": error.code.value,
            "message": error.message,
            "details": {}  # 生产环境不返回详细错误
        }
    )


def register_exception_handlers(app):
    """
    注册所有异常处理器到FastAPI应用
    
    Args:
        app: FastAPI应用实例
    """
    app.add_exception_handler(APIError, api_error_handler)
    app.add_exception_handler(RequestValidationError, validation_error_handler)
    app.add_exception_handler(SQLAlchemyError, database_error_handler)
    app.add_exception_handler(Exception, general_exception_handler)
    
    logger.info("Exception handlers registered successfully")
