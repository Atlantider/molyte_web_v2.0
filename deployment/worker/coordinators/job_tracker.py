"""
任务追踪模块

追踪运行中的任务
"""
import logging
from typing import Dict, Any, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from worker.core.config import WorkerConfig


logger = logging.getLogger(__name__)


class JobTracker:
    """任务追踪器"""
    
    def __init__(self, config: 'WorkerConfig'):
        """初始化任务追踪器"""
        self.config = config
        self.logger = logging.getLogger(self.__class__.__name__)
        
        # 运行中的任务 {job_id: job_info}
        self._running_jobs: Dict[int, Dict[str, Any]] = {}
    
    @property
    def running_jobs(self) -> Dict[int, Dict[str, Any]]:
        """获取运行中的任务"""
        return self._running_jobs
    
    @property
    def running_count(self) -> int:
        """获取运行中任务数量"""
        return len(self._running_jobs)
    
    def can_accept_new_job(self) -> bool:
        """检查是否可以接受新任务"""
        return self.running_count < self.config.max_concurrent_jobs
    
    def is_job_running(self, job_id: int) -> bool:
        """检查任务是否正在运行"""
        return job_id in self._running_jobs
    
    def add_job(
        self,
        job_id: int,
        job_type: str,
        slurm_job_id: Optional[str] = None,
        work_dir: Optional[str] = None,
        extra_info: Optional[Dict] = None
    ):
        """添加任务到追踪器"""
        import time
        
        self._running_jobs[job_id] = {
            'job_type': job_type,
            'slurm_job_id': slurm_job_id,
            'work_dir': work_dir,
            'start_time': time.time(),
            **(extra_info or {})
        }
        
        self.logger.debug(f"任务 {job_id} 已添加到追踪器")
    
    def remove_job(self, job_id: int):
        """从追踪器移除任务"""
        if job_id in self._running_jobs:
            del self._running_jobs[job_id]
            self.logger.debug(f"任务 {job_id} 已从追踪器移除")
    
    def get_job(self, job_id: int) -> Optional[Dict[str, Any]]:
        """获取任务信息"""
        return self._running_jobs.get(job_id)
    
    def update_job(self, job_id: int, **kwargs):
        """更新任务信息"""
        if job_id in self._running_jobs:
            self._running_jobs[job_id].update(kwargs)
