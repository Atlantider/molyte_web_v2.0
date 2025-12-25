"""
Binding 任务处理器

处理 Li-配体结合能计算
"""
import logging
from typing import Dict, Any

from worker.handlers.base import BaseHandler


logger = logging.getLogger(__name__)


class BindingHandler(BaseHandler):
    """Binding 任务处理器"""
    
    JOB_TYPE = "BINDING"
    
    def process(self, job: Dict[str, Any]) -> bool:
        """处理 Binding 任务"""
        job_id = job['id']
        self.logger.info(f"开始处理 Binding 任务 {job_id}")
        
        try:
            self.update_status(job_id, 'QUEUED')
            
            config = job.get('config', {})
            md_job_id = config.get('md_job_id')
            
            # 从现有 QC 结果计算 binding energy
            from worker.calculators.binding import BindingCalculator
            calculator = BindingCalculator(self.client)
            
            result = calculator.calculate_from_qc_results(
                md_job_id=md_job_id,
                config=config
            )
            
            if result.get('success'):
                self.update_status(
                    job_id, 'COMPLETED',
                    result=result,
                    qc_job_ids=result.get('qc_job_ids', [])
                )
                return True
            else:
                self.update_status(
                    job_id, 'FAILED',
                    error_message=result.get('error', 'Binding calculation failed')
                )
                return False
                
        except Exception as e:
            self.handle_error(job_id, e)
            return False
    
    def handle_completion(self, job_id: int, job_info: Dict, slurm_status: str):
        """Binding 任务通常不走 Slurm"""
        pass
