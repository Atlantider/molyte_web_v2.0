"""
任务处理器基类

定义所有任务处理器的通用接口和方法
"""
import logging
from abc import ABC, abstractmethod
from typing import Dict, Any, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from worker.core.config import WorkerConfig
    from worker.core.client import APIClient
    from worker.coordinators.job_tracker import JobTracker


logger = logging.getLogger(__name__)


class BaseHandler(ABC):
    """任务处理器基类"""
    
    # 任务类型标识
    JOB_TYPE: str = "BASE"
    
    def __init__(
        self,
        config: 'WorkerConfig',
        client: 'APIClient',
        job_tracker: 'JobTracker'
    ):
        """
        初始化处理器
        
        Args:
            config: Worker 配置
            client: API 客户端
            job_tracker: 任务追踪器
        """
        self.config = config
        self.client = client
        self.job_tracker = job_tracker
        self.logger = logging.getLogger(self.__class__.__name__)
    
    @abstractmethod
    def process(self, job: Dict[str, Any]) -> bool:
        """
        处理任务
        
        Args:
            job: 任务信息
            
        Returns:
            是否成功开始处理
        """
        pass
    
    @abstractmethod
    def handle_completion(
        self,
        job_id: int,
        job_info: Dict[str, Any],
        slurm_status: str
    ):
        """
        处理任务完成
        
        Args:
            job_id: 任务 ID
            job_info: 任务信息
            slurm_status: Slurm 状态
        """
        pass
    
    def update_status(
        self,
        job_id: int,
        status: str,
        **kwargs
    ) -> bool:
        """
        更新任务状态
        
        Args:
            job_id: 任务 ID
            status: 新状态
            **kwargs: 其他参数
            
        Returns:
            是否成功
        """
        success, _ = self.client.update_job_status(
            job_id=job_id,
            job_type=self.JOB_TYPE,
            status=status,
            **kwargs
        )
        return success
    
    def handle_error(
        self,
        job_id: int,
        error: Exception,
        error_message: Optional[str] = None
    ):
        """
        处理错误
        
        Args:
            job_id: 任务 ID
            error: 异常对象
            error_message: 自定义错误消息
        """
        msg = error_message or str(error)
        self.logger.error(f"任务 {job_id} 错误: {msg}", exc_info=True)
        
        self.update_status(
            job_id=job_id,
            status='FAILED',
            error_message=msg
        )
        
        # 从追踪器移除
        self.job_tracker.remove_job(job_id)
    
    def mark_running(
        self,
        job_id: int,
        slurm_job_id: Optional[str] = None,
        work_dir: Optional[str] = None
    ):
        """
        标记任务为运行中
        
        Args:
            job_id: 任务 ID
            slurm_job_id: Slurm 任务 ID
            work_dir: 工作目录
        """
        self.update_status(
            job_id=job_id,
            status='RUNNING',
            slurm_job_id=slurm_job_id,
            work_dir=work_dir
        )
        
        # 添加到追踪器
        self.job_tracker.add_job(
            job_id=job_id,
            job_type=self.JOB_TYPE,
            slurm_job_id=slurm_job_id,
            work_dir=work_dir
        )
    
    def mark_completed(
        self,
        job_id: int,
        cpu_hours: Optional[float] = None,
        **kwargs
    ):
        """
        标记任务为完成
        
        Args:
            job_id: 任务 ID
            cpu_hours: CPU 核时
            **kwargs: 其他参数
        """
        self.update_status(
            job_id=job_id,
            status='COMPLETED',
            cpu_hours=cpu_hours,
            progress=100.0,
            **kwargs
        )
        
        # 从追踪器移除
        self.job_tracker.remove_job(job_id)
