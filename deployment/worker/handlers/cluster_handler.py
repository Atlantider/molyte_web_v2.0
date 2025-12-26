"""
Cluster Analysis 任务处理器

处理高级溶剂化结构分析任务
"""
import logging
from typing import Dict, Any

from worker.handlers.base import BaseHandler


logger = logging.getLogger(__name__)


class ClusterHandler(BaseHandler):
    """Cluster Analysis 任务处理器"""
    
    JOB_TYPE = "CLUSTER_ANALYSIS"
    
    def process(self, job: Dict[str, Any]) -> bool:
        """处理 Cluster Analysis 任务"""
        job_id = job['id']
        self.logger.info(f"开始处理 Cluster Analysis 任务 {job_id}")
        
        try:
            self.update_status(job_id, 'RUNNING')
            
            config = job.get('config', {})
            calc_types = config.get('calc_types', [])
            selected_structures = config.get('selected_structures', [])
            qc_config = config.get('qc_config', {})
            
            # 创建 QC 任务
            qc_job_ids = []
            
            for structure in selected_structures:
                for calc_type in calc_types:
                    # 创建 QC 任务
                    qc_result = self._create_qc_job(
                        structure=structure,
                        calc_type=calc_type,
                        qc_config=qc_config
                    )
                    
                    if qc_result.get('success'):
                        qc_job_ids.append(qc_result['qc_job_id'])
            
            if qc_job_ids:
                # 更新状态为等待 QC
                self.client.update_job_status(
                    job_id=job_id,
                    job_type='CLUSTER_ANALYSIS',
                    status='WAITING_QC',
                    qc_job_ids=qc_job_ids
                )
                
                self.logger.info(
                    f"Cluster Analysis {job_id} 创建了 {len(qc_job_ids)} 个 QC 任务"
                )
                return True
            else:
                self.update_status(
                    job_id, 'FAILED',
                    error_message="Failed to create any QC jobs"
                )
                return False
                
        except Exception as e:
            self.handle_error(job_id, e)
            return False
    
    def handle_completion(self, job_id: int, job_info: Dict, slurm_status: str):
        """Cluster Analysis 任务通过 QC 任务完成触发"""
        pass
    
    def _create_qc_job(
        self,
        structure: Dict,
        calc_type: str,
        qc_config: Dict
    ) -> Dict:
        """创建 QC 任务"""
        try:
            response = self.client.post(
                "/api/v1/qc/jobs",
                data={
                    'molecule_name': structure.get('name'),
                    'xyz_content': structure.get('xyz_content'),
                    'charge': structure.get('charge', 0),
                    'spin_multiplicity': structure.get('spin_multiplicity', 1),
                    **qc_config
                }
            )
            
            if response:
                return {
                    'success': True,
                    'qc_job_id': response.get('id')
                }
            return {'success': False, 'error': 'API returned no response'}
            
        except Exception as e:
            return {'success': False, 'error': str(e)}
