"""
Redox 任务处理器

处理氧化还原电位计算
"""
import logging
from typing import Dict, Any

from worker.handlers.base import BaseHandler


logger = logging.getLogger(__name__)


class RedoxHandler(BaseHandler):
    """Redox 任务处理器"""
    
    JOB_TYPE = "REDOX"
    
    def process(self, job: Dict[str, Any]) -> bool:
        """处理 Redox 任务"""
        job_id = job['id']
        self.logger.info(f"开始处理 Redox 任务 {job_id}")
        
        try:
            self.update_status(job_id, 'RUNNING')
            
            config = job.get('config', {})
            species_list = config.get('species_list', [])
            
            # 计算每个物种的氧化还原电位
            from worker.calculators.redox import RedoxCalculator
            calculator = RedoxCalculator(self.client)
            
            results = []
            failed = []
            
            for species in species_list:
                try:
                    result = calculator.calculate_redox_potential(
                        species_name=species.get('name'),
                        qc_job_id=species.get('qc_job_id'),
                        config=config
                    )
                    if result.get('success'):
                        results.append(result)
                    else:
                        failed.append({
                            'species': species.get('name'),
                            'error': result.get('error')
                        })
                except Exception as e:
                    failed.append({
                        'species': species.get('name'),
                        'error': str(e)
                    })
            
            # 判断是否为部分结果
            is_partial = len(failed) > 0 and len(results) > 0
            
            if results:
                self.update_status(
                    job_id, 'COMPLETED',
                    result={
                        'redox_potentials': results,
                        'is_partial': is_partial,
                        'failed_species': failed if failed else None
                    },
                    error_message=(
                        f"Partial: {len(failed)} species failed"
                        if is_partial else None
                    )
                )
                return True
            else:
                self.update_status(
                    job_id, 'FAILED',
                    error_message="All species calculations failed",
                    result={'failed_species': failed}
                )
                return False
                
        except Exception as e:
            self.handle_error(job_id, e)
            return False
    
    def handle_completion(self, job_id: int, job_info: Dict, slurm_status: str):
        """Redox 任务通常通过API处理"""
        pass
