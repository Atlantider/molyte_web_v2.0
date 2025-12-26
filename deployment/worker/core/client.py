"""
API 客户端模块

负责与腾讯云后端 API 通信
"""
import logging
import requests
from typing import Dict, Any, List, Optional, Tuple
from datetime import datetime

from worker.core.config import WorkerConfig


logger = logging.getLogger(__name__)


class APIClient:
    """腾讯云后端 API 客户端"""
    
    def __init__(self, config: WorkerConfig):
        """
        初始化 API 客户端
        
        Args:
            config: Worker 配置
        """
        self.config = config
        self.base_url = config.api_base_url
        self.headers = {
            'Authorization': f'Bearer {config.api_token}',
            'Content-Type': 'application/json'
        }
        self.timeout = 30
        logger.info(f"API 客户端已初始化: {self.base_url}")
    
    # ==================== 任务获取 ====================
    
    def fetch_pending_jobs(
        self, 
        job_type: str, 
        limit: int = 10,
        status_filter: Optional[str] = None
    ) -> List[Dict[str, Any]]:
        """
        获取待处理任务
        
        Args:
            job_type: 任务类型 (MD, QC, POSTPROCESS, etc.)
            limit: 最大返回数量
            status_filter: 状态过滤
            
        Returns:
            待处理任务列表
        """
        try:
            params = {
                'job_type': job_type,
                'limit': limit,
                'worker_name': self.config.worker_name
            }
            if status_filter:
                params['status_filter'] = status_filter
            
            response = requests.get(
                f"{self.base_url}/api/v1/workers/jobs/pending",
                headers=self.headers,
                params=params,
                timeout=self.timeout
            )
            
            if response.status_code == 200:
                jobs = response.json()
                if jobs:
                    logger.debug(f"获取到 {len(jobs)} 个 {job_type} 待处理任务")
                return jobs
            else:
                logger.warning(f"获取 {job_type} 任务失败: {response.status_code}")
                return []
                
        except Exception as e:
            logger.error(f"获取 {job_type} 任务异常: {e}")
            return []
    
    def fetch_running_jobs(self) -> List[Dict[str, Any]]:
        """获取运行中的任务（用于恢复）"""
        try:
            response = requests.get(
                f"{self.base_url}/api/v1/workers/jobs/running",
                headers=self.headers,
                timeout=self.timeout
            )
            
            if response.status_code == 200:
                return response.json()
            else:
                logger.warning(f"获取运行中任务失败: {response.status_code}")
                return []
                
        except Exception as e:
            logger.error(f"获取运行中任务异常: {e}")
            return []
    
    # ==================== 状态更新 ====================
    
    def update_job_status(
        self,
        job_id: int,
        job_type: str,
        status: str,
        slurm_job_id: Optional[str] = None,
        work_dir: Optional[str] = None,
        error_message: Optional[str] = None,
        progress: Optional[float] = None,
        cpu_hours: Optional[float] = None,
        result: Optional[Dict] = None,
        **kwargs
    ) -> Tuple[bool, str]:
        """
        更新任务状态
        
        Args:
            job_id: 任务 ID
            job_type: 任务类型
            status: 新状态
            slurm_job_id: Slurm 任务 ID
            work_dir: 工作目录
            error_message: 错误信息
            progress: 进度 (0-100)
            cpu_hours: CPU 核时
            result: 结果数据
            **kwargs: 其他参数
            
        Returns:
            (成功与否, 消息)
        """
        try:
            data = {
                'status': status,
                'job_type': job_type,
                'worker_name': self.config.worker_name
            }
            
            if slurm_job_id:
                data['slurm_job_id'] = slurm_job_id
            if work_dir:
                data['work_dir'] = work_dir
            if error_message is not None:
                data['error_message'] = error_message
            if progress is not None:
                data['progress'] = progress
            if cpu_hours is not None:
                data['cpu_hours'] = cpu_hours
            if result:
                data['result'] = result
            
            # 添加额外参数
            data.update(kwargs)
            
            response = requests.put(
                f"{self.base_url}/api/v1/workers/jobs/{job_id}/status",
                headers=self.headers,
                json=data,
                timeout=self.timeout
            )
            
            if response.status_code == 200:
                logger.info(f"任务 {job_id} 状态更新为 {status}")
                return True, "成功"
            else:
                msg = f"状态更新失败: {response.status_code} - {response.text}"
                logger.warning(msg)
                return False, msg
                
        except Exception as e:
            msg = f"状态更新异常: {e}"
            logger.error(msg)
            return False, msg
    
    # ==================== 心跳 ====================
    
    def send_heartbeat(
        self, 
        running_jobs_count: int,
        partitions: Optional[List[Dict]] = None
    ) -> bool:
        """
        发送心跳
        
        Args:
            running_jobs_count: 运行中任务数
            partitions: 分区信息
            
        Returns:
            是否成功
        """
        try:
            data = {
                'worker_name': self.config.worker_name,
                'status': 'running',
                'running_jobs': running_jobs_count,
                'timestamp': datetime.now().isoformat()
            }
            
            if partitions:
                data['partitions'] = partitions
            
            response = requests.post(
                f"{self.base_url}/api/v1/workers/heartbeat",
                headers=self.headers,
                json=data,
                timeout=self.timeout
            )
            
            if response.status_code == 200:
                logger.debug("心跳发送成功")
                return True
            else:
                logger.warning(f"心跳发送失败: {response.status_code}")
                return False
                
        except Exception as e:
            logger.error(f"心跳发送异常: {e}")
            return False
    
    # ==================== 结果提交 ====================
    
    def submit_qc_result(self, qc_job_id: int, result_data: Dict) -> bool:
        """提交 QC 结果"""
        try:
            response = requests.post(
                f"{self.base_url}/api/v1/workers/qc/{qc_job_id}/result",
                headers=self.headers,
                json=result_data,
                timeout=self.timeout
            )
            
            if response.status_code == 200:
                logger.info(f"QC 任务 {qc_job_id} 结果提交成功")
                return True
            else:
                logger.warning(f"QC 结果提交失败: {response.status_code}")
                return False
                
        except Exception as e:
            logger.error(f"QC 结果提交异常: {e}")
            return False
    
    def trigger_postprocess(self, md_job_id: int) -> bool:
        """触发后处理"""
        try:
            response = requests.post(
                f"{self.base_url}/api/v1/workers/md/{md_job_id}/postprocess",
                headers=self.headers,
                timeout=self.timeout
            )
            
            if response.status_code == 200:
                logger.info(f"MD 任务 {md_job_id} 后处理已触发")
                return True
            else:
                logger.warning(f"后处理触发失败: {response.status_code}")
                return False
                
        except Exception as e:
            logger.error(f"后处理触发异常: {e}")
            return False
    
    # ==================== Desolvation 相关 ====================
    
    def trigger_desolvation_phase(
        self, 
        job_id: int, 
        phase: int
    ) -> Tuple[bool, Dict]:
        """
        触发去溶剂化计算阶段
        
        Args:
            job_id: Postprocess 任务 ID
            phase: 阶段 (1 或 2)
            
        Returns:
            (成功与否, 结果)
        """
        try:
            response = requests.post(
                f"{self.base_url}/api/v1/workers/postprocess/{job_id}/execute",
                headers=self.headers,
                json={'phase': phase},
                timeout=60  # 较长超时
            )
            
            if response.status_code == 200:
                result = response.json()
                logger.info(f"Desolvation 任务 {job_id} 阶段 {phase} 执行成功")
                return True, result
            else:
                logger.warning(f"Desolvation 阶段执行失败: {response.status_code}")
                return False, {}
                
        except Exception as e:
            logger.error(f"Desolvation 阶段执行异常: {e}")
            return False, {}
    
    # ==================== 通用请求 ====================
    
    def get(self, endpoint: str, params: Optional[Dict] = None) -> Optional[Dict]:
        """通用 GET 请求"""
        try:
            response = requests.get(
                f"{self.base_url}{endpoint}",
                headers=self.headers,
                params=params,
                timeout=self.timeout
            )
            if response.status_code == 200:
                return response.json()
            return None
        except Exception as e:
            logger.error(f"GET {endpoint} 失败: {e}")
            return None
    
    def post(self, endpoint: str, data: Optional[Dict] = None) -> Optional[Dict]:
        """通用 POST 请求"""
        try:
            response = requests.post(
                f"{self.base_url}{endpoint}",
                headers=self.headers,
                json=data,
                timeout=self.timeout
            )
            if response.status_code == 200:
                return response.json()
            return None
        except Exception as e:
            logger.error(f"POST {endpoint} 失败: {e}")
            return None
