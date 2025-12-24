"""
FastAPI main application
"""
import logging
from fastapi import FastAPI, Request, status
from fastapi.middleware.cors import CORSMiddleware
from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse
from app.config import settings
from app.api.v1 import api_router
from app.core.logger import logger
from app.database import engine, Base
from app.init_db import init_db

# Disable all logging
logging.getLogger('sqlalchemy.engine').setLevel(logging.CRITICAL)
logging.getLogger('sqlalchemy.pool').setLevel(logging.CRITICAL)
logging.getLogger('sqlalchemy.dialects').setLevel(logging.CRITICAL)
logging.getLogger('uvicorn.access').setLevel(logging.CRITICAL)
logging.getLogger('uvicorn').setLevel(logging.CRITICAL)
logging.getLogger('app').setLevel(logging.CRITICAL)

# Create database tables and initialize database
Base.metadata.create_all(bind=engine)
init_db()

# Create FastAPI app
app = FastAPI(
    title=settings.APP_NAME,
    version=settings.APP_VERSION,
    description="Molyte Web - Electrolyte Molecular Dynamics Simulation Platform",
    docs_url="/docs",
    redoc_url="/redoc",
)

# Configure CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=settings.cors_origins_list,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


# Add validation error handler
@app.exception_handler(RequestValidationError)
async def validation_exception_handler(request: Request, exc: RequestValidationError):
    """Handle validation errors with detailed logging"""
    logger.error(f"Validation error for {request.method} {request.url}")
    try:
        body = await request.body()
        logger.error(f"Request body: {body}")
    except Exception:
        pass
    logger.error(f"Validation errors: {exc.errors()}")

    # 安全地处理 body，避免 FormData 序列化问题
    body_content = None
    if exc.body is not None:
        try:
            if hasattr(exc.body, '__dict__'):
                body_content = str(exc.body)
            else:
                body_content = exc.body
        except Exception:
            body_content = None

    return JSONResponse(
        status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
        content={
            "detail": exc.errors(),
            "body": body_content
        }
    )


@app.on_event("startup")
async def startup_event():
    """Application startup event"""
    pass


@app.on_event("shutdown")
async def shutdown_event():
    """Application shutdown event"""
    pass


@app.get("/")
async def root():
    """Root endpoint"""
    return {
        "name": settings.APP_NAME,
        "version": settings.APP_VERSION,
        "status": "running",
        "docs": "/docs"
    }


@app.get("/health")
async def health_check():
    """Health check endpoint"""
    return {"status": "healthy"}


# Include API router
app.include_router(api_router, prefix="/api/v1")


# Configure periodic tasks
@app.on_event("startup")
async def startup_event():
    """Configure periodic tasks on startup"""
    import logging
    logger = logging.getLogger(__name__)

    try:
        from app.celery_app import celery_app
        from celery.schedules import crontab

        # 每天凌晨 2 点检查核时过期
        celery_app.conf.beat_schedule = {
            'check-cpu-hours-expiration': {
                'task': 'app.tasks.compensation_scheduler.check_cpu_hours_expiration',
                'schedule': crontab(hour=2, minute=0),  # 每天 2:00 AM
            },
            'process-pending-compensations': {
                'task': 'app.tasks.compensation_scheduler.process_pending_compensations',
                'schedule': crontab(minute=0, hour='*/1'),  # 每小时执行一次
            },
            'cleanup-old-records': {
                'task': 'app.tasks.compensation_scheduler.cleanup_old_records',
                'schedule': crontab(hour=3, minute=0),  # 每天 3:00 AM
            },
            'aggregate-daily-stats': {
                'task': 'app.tasks.stats_aggregator.aggregate_daily_stats',
                'schedule': crontab(hour=1, minute=0),  # 每天 1:00 AM
            },
        }
        logger.info("Periodic tasks configured")
    except Exception as e:
        # Celery 导入失败，记录警告但继续启动
        logger.warning(f"⚠️  Celery 初始化失败: {e}. 后端将在没有异步任务的情况下运行。")


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(
        "app.main:app",
        host="0.0.0.0",
        port=8000,
        reload=settings.DEBUG
    )

