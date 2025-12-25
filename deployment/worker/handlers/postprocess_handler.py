"""
后处理任务处理器

处理去溶剂化能、Redox、Reorg等后处理任务
支持部分结果保存和容错
"""
import logging
from typing import Dict, Any, List, Optional, Tuple
from pathlib import Path

from worker.handlers.base import BaseHandler


logger = logging.getLogger(__name__)


class PostprocessHandler(BaseHandler):
    """后处理任务处理器"""
    
    JOB_TYPE = "POSTPROCESS"
    
    def process(self, job: Dict[str, Any]) -> bool:
        """
        处理后处理任务
        
        Args:
            job: 任务信息
            
        Returns:
            是否成功开始处理
        """
        job_id = job['id']
        config = job.get('config', {})
        job_type = config.get('job_type', 'UNKNOWN')
        
        self.logger.info(f"开始处理后处理任务 {job_id} (类型: {job_type})")
        
        try:
            # 更新状态为 QUEUED
            self.update_status(job_id, 'QUEUED')
            
            # 根据任务类型分发处理
            if job_type == 'DESOLVATION_ENERGY':
                return self._process_desolvation(job_id, config)
            else:
                raise ValueError(f"Unknown postprocess job type: {job_type}")
                
        except Exception as e:
            self.handle_error(job_id, e)
            return False
    
    def handle_completion(
        self,
        job_id: int,
        job_info: Dict[str, Any],
        slurm_status: str
    ):
        """处理任务完成（后处理任务通常通过API处理，不走Slurm）"""
        pass
    
    # ==================== 去溶剂化处理 ====================
    
    def _process_desolvation(self, job_id: int, config: Dict) -> bool:
        """
        处理去溶剂化能计算
        
        新增：
        1. 检查QC任务完成状态
        2. 支持部分结果计算和保存
        3. 失败任务自动跳过
        """
        self.logger.info(f"开始去溶剂化能计算任务 {job_id}")
        
        try:
            # 获取QC任务列表
            qc_job_ids = config.get('qc_job_ids', [])
            phase = config.get('phase', 1)
            
            if phase == 1:
                # Phase 1: 创建QC任务（通过API）
                return self._trigger_desolvation_phase1(job_id)
            
            elif phase == 2:
                # Phase 2: 检查QC结果并计算
                return self._process_desolvation_phase2(job_id, config)
            
            else:
                raise ValueError(f"Unknown phase: {phase}")
                
        except Exception as e:
            self.logger.error(f"去溶剂化任务 {job_id} 失败: {e}", exc_info=True)
            self.handle_error(job_id, e)
            return False
    
    def _trigger_desolvation_phase1(self, job_id: int) -> bool:
        """触发Phase 1: 创建QC任务"""
        success, result = self.client.trigger_desolvation_phase(job_id, phase=1)
        
        if success:
            self.logger.info(f"去溶剂化任务 {job_id} Phase 1 成功创建QC任务")
            return True
        else:
            self.logger.error(f"去溶剂化任务 {job_id} Phase 1 失败")
            return False
    
    def _process_desolvation_phase2(self, job_id: int, config: Dict) -> bool:
        """
        Phase 2: 检查QC结果并计算去溶剂化能
        
        支持部分结果：
        - 收集所有QC结果
        - 分类成功/失败
        - 只要有部分成功就计算部分结果
        """
        qc_job_ids = config.get('qc_job_ids', [])
        cluster_qc_id = config.get('cluster_qc_job_id')
        ligand_qc_jobs = config.get('ligand_qc_jobs', {})
        desolvation_mode = config.get('desolvation_mode', 'stepwise')
        
        # 1. 收集QC结果
        qc_results = self._collect_qc_results(qc_job_ids)
        
        # 2. 分类结果
        successful, failed = self._classify_qc_results(qc_results)
        
        self.logger.info(
            f"去溶剂化 {job_id}: {len(successful)} 成功, {len(failed)} 失败"
        )
        
        # 3. 检查Cluster是否成功（必需）
        if cluster_qc_id not in successful:
            self.logger.error(f"Cluster QC {cluster_qc_id} 失败，无法计算去溶剂化能")
            self.update_status(
                job_id, 'FAILED',
                error_message="Cluster QC calculation failed"
            )
            return False
        
        # 4. 根据模式处理
        if desolvation_mode == 'stepwise':
            return self._calculate_stepwise_desolvation(
                job_id, config, successful, failed
            )
        elif desolvation_mode == 'full':
            return self._calculate_full_desolvation(
                job_id, config, successful, failed
            )
        else:
            raise ValueError(f"Unknown desolvation mode: {desolvation_mode}")
    
    def _collect_qc_results(self, qc_job_ids: List[int]) -> Dict[int, Dict]:
        """收集QC任务结果"""
        results = {}
        
        for qc_id in qc_job_ids:
            result = self.client.get(f"/api/v1/qc/{qc_id}")
            if result:
                results[qc_id] = result
            else:
                results[qc_id] = {'status': 'NOT_FOUND', 'error': 'QC job not found'}
        
        return results
    
    def _classify_qc_results(
        self, qc_results: Dict[int, Dict]
    ) -> Tuple[List[int], List[Dict]]:
        """
        分类QC结果为成功和失败
        
        Returns:
            (成功的QC ID列表, 失败的详情列表)
        """
        successful = []
        failed = []
        
        for qc_id, result in qc_results.items():
            status = result.get('status', 'UNKNOWN')
            
            if status == 'COMPLETED':
                successful.append(qc_id)
            else:
                failed.append({
                    'qc_job_id': qc_id,
                    'status': status,
                    'error': result.get('error_message', 'Unknown error')
                })
        
        return successful, failed
    
    def _calculate_stepwise_desolvation(
        self,
        job_id: int,
        config: Dict,
        successful: List[int],
        failed: List[Dict]
    ) -> bool:
        """
        计算Stepwise模式去溶剂化能
        
        支持部分结果：跳过失败的配体，保存成功的部分
        """
        from worker.calculators.desolvation import DesolvationCalculator
        
        try:
            calculator = DesolvationCalculator(self.client)
            
            # 计算部分结果
            results = calculator.calculate_stepwise_partial(
                config=config,
                successful_qc_ids=successful,
                failed_details=failed
            )
            
            # 判断是否为部分结果
            is_partial = len(failed) > 0
            total_ligands = len(config.get('ligand_qc_jobs', {}))
            calculated_ligands = len([r for r in results.get('per_ligand_results', []) 
                                      if r.get('status') == 'calculated'])
            
            # 上传结果
            from worker.uploaders import PartialUploader
            uploader = PartialUploader(self.config)
            uploader.upload_desolvation_results(
                job_id=job_id,
                results=results,
                is_partial=is_partial
            )
            
            # 更新状态
            if calculated_ligands > 0:
                # 有成功的结果
                self.update_status(
                    job_id, 'COMPLETED',
                    result={
                        'is_partial': is_partial,
                        'total_ligands': total_ligands,
                        'calculated_ligands': calculated_ligands,
                        'failed_count': len(failed)
                    },
                    error_message=(
                        f"Partial: {len(failed)} ligands failed" 
                        if is_partial else None
                    )
                )
                self.logger.info(
                    f"去溶剂化 {job_id} 完成: {calculated_ligands}/{total_ligands} 计算成功"
                )
                return True
            else:
                # 全部失败
                self.update_status(
                    job_id, 'FAILED',
                    error_message=f"All {len(failed)} ligand calculations failed",
                    result={'failed_details': failed}
                )
                return False
                
        except Exception as e:
            self.logger.error(f"Stepwise计算失败: {e}", exc_info=True)
            self.handle_error(job_id, e)
            return False
    
    def _calculate_full_desolvation(
        self,
        job_id: int,
        config: Dict,
        successful: List[int],
        failed: List[Dict]
    ) -> bool:
        """
        计算Full模式去溶剂化能
        
        改进：支持部分结果
        - 如果有部分配体成功，计算部分去溶剂化能
        - 提供不完整结果的估算
        """
        from worker.calculators.desolvation import DesolvationCalculator
        
        try:
            calculator = DesolvationCalculator(self.client)
            
            # 检查中心离子是否成功
            center_ion_id = config.get('center_ion_job_id')
            if center_ion_id and center_ion_id not in successful:
                self.update_status(
                    job_id, 'FAILED',
                    error_message="Center ion QC calculation failed"
                )
                return False
            
            # 计算全模式结果（支持部分）
            results = calculator.calculate_full_partial(
                config=config,
                successful_qc_ids=successful,
                failed_details=failed
            )
            
            is_partial = results.get('is_partial', False)
            
            # 上传结果
            from worker.uploaders import PartialUploader
            uploader = PartialUploader(self.config)
            uploader.upload_desolvation_results(
                job_id=job_id,
                results=results,
                is_partial=is_partial
            )
            
            # 更新状态
            if results.get('total_desolvation_energy') is not None:
                self.update_status(
                    job_id, 'COMPLETED',
                    result={
                        'is_partial': is_partial,
                        'total_desolvation_energy': results['total_desolvation_energy'],
                        'calculated_ligands': results.get('calculated_ligand_count', 0),
                        'total_ligands': results.get('total_ligand_count', 0)
                    },
                    error_message=(
                        f"Partial: {len(failed)} ligands failed"
                        if is_partial else None
                    )
                )
                return True
            else:
                self.update_status(
                    job_id, 'FAILED',
                    error_message="Cannot calculate total desolvation energy",
                    result={'failed_details': failed}
                )
                return False
                
        except Exception as e:
            self.logger.error(f"Full模式计算失败: {e}", exc_info=True)
            self.handle_error(job_id, e)
            return False
