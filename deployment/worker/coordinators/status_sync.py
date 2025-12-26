"""
状态同步模块

负责与云端同步任务状态
"""
import logging
from typing import Dict, Any, TYPE_CHECKING

if TYPE_CHECKING:
    from worker.core.config import WorkerConfig
    from worker.core.client import APIClient


logger = logging.getLogger(__name__)


class StatusSynchronizer:
    """状态同步器"""
    
    def __init__(self, config: 'WorkerConfig', client: 'APIClient'):
        """初始化状态同步器"""
        self.config = config
        self.client = client
        self.logger = logging.getLogger(self.__class__.__name__)
    
    def sync_job_status(
        self,
        job_id: int,
        job_type: str,
        local_status: str
    ):
        """同步任务状态到云端"""
        self.client.update_job_status(
            job_id=job_id,
            job_type=job_type,
            status=local_status
        )
