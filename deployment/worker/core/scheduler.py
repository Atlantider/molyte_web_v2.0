"""
任务调度器模块

负责任务的调度和主循环管理
"""
import logging
import time
from typing import Dict, Any, List, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from worker.core.config import WorkerConfig
    from worker.core.client import APIClient
    from worker.handlers.base import BaseHandler
    from worker.coordinators.job_tracker import JobTracker


logger = logging.getLogger(__name__)


class TaskScheduler:
    """任务调度器"""
    
    def __init__(
        self,
        config: 'WorkerConfig',
        client: 'APIClient',
        handlers: Dict[str, 'BaseHandler'],
        job_tracker: 'JobTracker'
    ):
        """
        初始化调度器
        
        Args:
            config: Worker 配置
            client: API 客户端
            handlers: 任务处理器字典
            job_tracker: 任务追踪器
        """
        self.config = config
        self.client = client
        self.handlers = handlers
        self.job_tracker = job_tracker
        
        # 上次操作时间记录
        self._last_heartbeat = 0.0
        self._last_check_api_running = 0.0
        self._last_sync_anion_library = 0.0
        
        logger.info("任务调度器已初始化")
    
    def start_polling(self):
        """开始轮询主循环"""
        logger.info("开始轮询...")
        
        # 恢复运行中的任务
        self._recover_running_jobs()
        
        while True:
            try:
                self._poll_cycle()
                
                # 等待下一次轮询
                time.sleep(self.config.poll_interval)
                
            except KeyboardInterrupt:
                logger.info("收到中断信号，正在退出...")
                break
            except Exception as e:
                logger.error(f"主循环错误: {e}", exc_info=True)
                time.sleep(10)
    
    def _poll_cycle(self):
        """单次轮询周期"""
        current_time = time.time()
        
        # 1. 发送心跳
        if current_time - self._last_heartbeat > self.config.heartbeat_interval:
            self._send_heartbeat()
            self._last_heartbeat = current_time
        
        # 2. 检查运行中的任务
        self._check_running_jobs()
        
        # 3. 定期检查 API 中的 RUNNING 状态任务（容错机制）
        if current_time - self._last_check_api_running > 300:  # 5 分钟
            self._check_api_running_jobs()
            self._last_check_api_running = current_time
        
        # 4. 获取并处理新任务
        if self.job_tracker.can_accept_new_job():
            self._fetch_and_process_new_jobs()
        else:
            logger.debug(f"当前运行任务数已达最大并发数，跳过拉取新任务")
    
    def _recover_running_jobs(self):
        """恢复运行中的任务"""
        logger.info("正在恢复运行中的任务...")
        
        # 从 API 获取运行中的任务
        running_jobs = self.client.fetch_running_jobs()
        
        for job in running_jobs:
            job_id = job['id']
            job_type = job['type']
            
            # 添加到追踪器
            self.job_tracker.add_job(
                job_id=job_id,
                job_type=job_type,
                slurm_job_id=job.get('slurm_job_id'),
                work_dir=job.get('work_dir')
            )
            logger.info(f"恢复任务: {job_type} #{job_id}")
        
        logger.info(f"恢复了 {len(running_jobs)} 个运行中的任务")
    
    def _send_heartbeat(self):
        """发送心跳"""
        from worker.utils.slurm import SlurmManager
        
        # 获取分区信息
        slurm = SlurmManager()
        partitions = slurm.get_partitions()
        
        # 发送心跳
        self.client.send_heartbeat(
            running_jobs_count=self.job_tracker.running_count,
            partitions=partitions
        )
    
    def _check_running_jobs(self):
        """检查运行中的任务"""
        from worker.utils.slurm import SlurmManager
        
        slurm = SlurmManager()
        
        for job_id, job_info in list(self.job_tracker.running_jobs.items()):
            slurm_job_id = job_info.get('slurm_job_id')
            if not slurm_job_id:
                continue
            
            # 检查 Slurm 状态
            slurm_status = slurm.get_job_status(slurm_job_id)
            
            if slurm_status in ['COMPLETED', 'FAILED', 'CANCELLED', 'TIMEOUT']:
                # 任务完成，调用对应处理器处理完成逻辑
                job_type = job_info.get('job_type', 'MD')
                handler_key = job_type.lower()
                
                if handler_key in self.handlers:
                    self.handlers[handler_key].handle_completion(
                        job_id=job_id,
                        job_info=job_info,
                        slurm_status=slurm_status
                    )
                
                # 从追踪器移除
                self.job_tracker.remove_job(job_id)
    
    def _check_api_running_jobs(self):
        """检查 API 中的 RUNNING 状态任务（容错机制）"""
        # 这是为了处理 Worker 重启后漏掉的任务
        logger.debug("检查 API 中的 RUNNING 状态任务...")
        
        # 获取 API 中 RUNNING 状态的任务
        # 这个逻辑在各个 handler 中实现
        pass
    
    def _fetch_and_process_new_jobs(self):
        """获取并处理新任务"""
        
        # 按优先级获取各类型任务
        job_types = [
            ('MD', 'md'),
            ('QC', 'qc'),
            ('POSTPROCESS', 'postprocess'),
            ('BINDING', 'binding'),
            ('REDOX', 'redox'),
            ('REORG', 'reorg'),
            ('CLUSTER_ANALYSIS', 'cluster'),
            ('ANION_GENERATION', 'anion'),
        ]
        
        for api_type, handler_key in job_types:
            if not self.job_tracker.can_accept_new_job():
                break
            
            # 获取待处理任务
            jobs = self.client.fetch_pending_jobs(api_type, limit=1)
            
            for job in jobs:
                if not self.job_tracker.can_accept_new_job():
                    break
                
                # 交给对应处理器处理
                if handler_key in self.handlers:
                    try:
                        self.handlers[handler_key].process(job)
                    except Exception as e:
                        logger.error(f"处理 {api_type} 任务失败: {e}", exc_info=True)
    
    def stop(self):
        """停止调度器"""
        logger.info("正在停止调度器...")
        # 可以在这里添加清理逻辑
