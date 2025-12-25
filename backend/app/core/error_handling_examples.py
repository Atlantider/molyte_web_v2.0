"""
使用新错误处理系统的API示例

展示如何在API endpoint中使用统一的错误处理
"""
from typing import List, Optional
from fastapi import APIRouter, Depends, Query
from sqlalchemy.orm import Session

from app.database import get_db
from app.core.errors import (
    APIError,
    ErrorCode,
    not_found_error,
    invalid_input_error,
    quota_exceeded_error,
    handle_errors,
)
from app.models.job import MDJob
from app.schemas.job import JobResponse, JobCreate

router = APIRouter()


# ========== 示例1: 基础错误处理 ==========

@router.get("/jobs/{job_id}", response_model=JobResponse)
@handle_errors  # 使用装饰器自动处理未捕获异常
async def get_job(
    job_id: int,
    db: Session = Depends(get_db)
):
    """获取任务详情"""
    job = db.query(MDJob).filter(MDJob.id == job_id).first()
    
    if not job:
        # 使用便捷函数创建错误
        raise not_found_error("任务", job_id)
    
    return job


# ========== 示例2: 带详情的错误 ==========

@router.post("/jobs", response_model=JobResponse, status_code=201)
async def create_job(
    job_data: JobCreate,
    user_id: int = Depends(get_current_user_id),
    db: Session = Depends(get_db)
):
    """创建新任务"""
    # 检查用户配额
    user = db.query(User).get(user_id)
    if user.cpu_hours < 10:
        # 创建带详情的配额错误
        raise quota_exceeded_error(
            quota_type="CPU核时",
            current=user.cpu_hours,
            limit=10
        )
    
    # 验证输入
    if not job_data.name or len(job_data.name) < 3:
        raise invalid_input_error(
            "任务名称至少需要3个字符",
            details={"field": "name", "min_length": 3}
        )
    
    # 创建任务
    job = MDJob(**job_data.dict(), user_id=user_id)
    db.add(job)
    db.commit()
    db.refresh(job)
    
    return job


# ========== 示例3: 自定义错误 ==========

@router.delete("/jobs/{job_id}")
async def cancel_job(
    job_id: int,
    db: Session = Depends(get_db)
):
    """取消任务"""
    job = db.query(MDJob).get(job_id)
    
    if not job:
        raise not_found_error("任务", job_id)
    
    # 检查任务状态
    if job.status not in ['CREATED', 'QUEUED']:
        # 创建自定义业务逻辑错误
        raise APIError(
            code=ErrorCode.JOB_CANNOT_BE_CANCELLED,
            message=f"任务状态为 {job.status}，无法取消",
            status_code=400,
            details={
                "current_status": job.status,
                "allowed_statuses": ['CREATED', 'QUEUED']
            }
        )
    
    job.status = 'CANCELLED'
    db.commit()
    
    return {"message": "任务已取消", "job_id": job_id}


# ========== 示例4: 数据库错误处理 ==========

from app.core.errors import database_error

@router.get("/jobs")
async def list_jobs(
    skip: int = Query(0, ge=0),
    limit: int = Query(10, ge=1, le=100),
    db: Session = Depends(get_db)
):
    """列出任务"""
    try:
        jobs = db.query(MDJob).offset(skip).limit(limit).all()
        total = db.query(MDJob).count()
        
        return {
            "items": jobs,
            "total": total,
            "skip": skip,
            "limit": limit
        }
    except SQLAlchemyError as e:
        # 数据库错误会被全局handler捕获
        # 这里手动处理是为了添加更多上下文
        raise database_error("查询", e)


# ========== 前端使用示例 ==========

"""
// 前端使用示例 (TypeScript)

import { extractAPIError, getUserFriendlyMessage, ErrorCode } from '@/utils/errorHandler';
import { message } from 'antd';

async function getJob(jobId: number) {
  try {
    const response = await axios.get(`/api/v1/jobs/${jobId}`);
    return response.data;
  } catch (error) {
    const apiError = extractAPIError(error);
    
    // 根据错误码处理
    if (apiError.code === ErrorCode.NOT_FOUND) {
      message.error('任务不存在');
      navigate('/jobs');
    } else if (apiError.code === ErrorCode.UNAUTHORIZED) {
      message.warning('请先登录');
      navigate('/login');
    } else {
      // 显示用户友好的消息
      message.error(getUserFriendlyMessage(apiError));
    }
    
    throw apiError;
  }
}

// 或使用React Hook
function JobDetail({ jobId }) {
  const { handleError } = useErrorHandler();
  
  const loadJob = async () => {
    try {
      const data = await getJob(jobId);
      setJob(data);
    } catch (error) {
      const { message, showNotification } = handleError(error);
      notification.error(showNotification());
    }
  };
}
"""
