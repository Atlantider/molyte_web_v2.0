"""
重试协调模块

协调失败任务的重试逻辑（与云端配合）
"""
import logging
from typing import Dict, Any, List, TYPE_CHECKING

if TYPE_CHECKING:
    from worker.core.config import WorkerConfig
    from worker.core.client import APIClient


logger = logging.getLogger(__name__)


class RetryCoordinator:
    """重试协调器
    
    工作原理：
    1. 云端检测到失败任务，决定是否重试
    2. 如果需要重试，云端重置任务状态为 SUBMITTED
    3. Polling Worker 在下一次轮询时获取这些任务
    4. Worker 使用调整后的参数重新处理
    """
    
    def __init__(self, config: 'WorkerConfig', client: 'APIClient'):
        """初始化重试协调器"""
        self.config = config
        self.client = client
        self.logger = logging.getLogger(self.__class__.__name__)
    
    def check_for_retry_jobs(self) -> List[Dict[str, Any]]:
        """检查是否有需要重试的任务"""
        # 这个功能主要由云端控制
        # Worker 只需要正常轮询 pending 任务
        # 云端会将需要重试的任务重置为 SUBMITTED
        return []
    
    def get_adjusted_parameters(
        self,
        job_type: str,
        original_config: Dict,
        retry_count: int,
        previous_error: str
    ) -> Dict[str, Any]:
        """
        获取调整后的参数（用于重试）
        
        Args:
            job_type: 任务类型
            original_config: 原始配置
            retry_count: 重试次数
            previous_error: 上次错误信息
            
        Returns:
            调整后的配置
        """
        adjusted = original_config.copy()
        
        if job_type == 'QC':
            # SCF 收敛问题
            if 'SCF' in previous_error or 'convergence' in previous_error.lower():
                if retry_count == 1:
                    # 第一次重试：放宽收敛标准
                    adjusted['scf_convergence'] = 'loose'
                    adjusted['max_scf_cycles'] = 200
                elif retry_count >= 2:
                    # 第二次重试：使用阻尼
                    adjusted['scf_algorithm'] = 'damping'
                    adjusted['damping_factor'] = 0.5
                    adjusted['max_scf_cycles'] = 300
            
            # 优化问题
            if 'optimization' in previous_error.lower():
                adjusted['opt_max_cycles'] = adjusted.get('opt_max_cycles', 100) + 50
        
        return adjusted
