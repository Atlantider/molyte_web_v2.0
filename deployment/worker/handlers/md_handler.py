"""
MD 任务处理器

处理分子动力学模拟任务
"""
import logging
import time
from typing import Dict, Any, List, Optional
from pathlib import Path

from worker.handlers.base import BaseHandler


logger = logging.getLogger(__name__)


class MDHandler(BaseHandler):
    """MD 任务处理器"""
    
    JOB_TYPE = "MD"
    
    def process(self, job: Dict[str, Any]) -> bool:
        """
        处理 MD 任务
        
        Args:
            job: 任务信息
            
        Returns:
            是否成功开始处理
        """
        job_id = job['id']
        self.logger.info(f"开始处理 MD 任务 {job_id}")
        
        try:
            # 更新状态为 QUEUED
            self.update_status(job_id, 'QUEUED')
            
            # 获取任务配置
            job_data = job.get('config', {})
            charge_method = job_data.get('charge_method', 'ligpargen')
            
            # 初始化 MolyteWrapper
            wrapper = self._init_molyte_wrapper()
            
            # 检查是否需要 RESP 计算
            if charge_method == 'resp':
                solvents = job_data.get('solvents', [])
                solvents_needing_resp = wrapper.get_solvents_needing_resp(solvents)
                
                if solvents_needing_resp:
                    self.logger.info(
                        f"MD 任务 {job_id} 需要 RESP 计算: "
                        f"{[s['name'] for s in solvents_needing_resp]}"
                    )
                    return self._start_resp_calculations(
                        job_id, job, solvents_needing_resp, wrapper
                    )
            
            # 直接生成 LAMMPS 输入文件
            return self._continue_md_job(job_id, job, wrapper)
            
        except Exception as e:
            self.handle_error(job_id, e)
            return False
    
    def handle_completion(
        self,
        job_id: int,
        job_info: Dict[str, Any],
        slurm_status: str
    ):
        """
        处理 MD 任务完成
        
        Args:
            job_id: 任务 ID
            job_info: 任务信息
            slurm_status: Slurm 状态
        """
        work_dir = Path(job_info.get('work_dir', ''))
        slurm_job_id = job_info.get('slurm_job_id')
        
        self.logger.info(f"MD 任务 {job_id} Slurm 状态: {slurm_status}")
        
        if slurm_status == 'COMPLETED':
            try:
                # 检查 LAMMPS 输出
                if self._check_lammps_success(work_dir):
                    # 获取 CPU 核时
                    from worker.utils.slurm import SlurmManager
                    slurm = SlurmManager()
                    cpu_hours = slurm.get_job_cpu_hours(slurm_job_id)
                    
                    # 上传结果到 COS
                    from worker.uploaders import ResultUploader
                    uploader = ResultUploader(self.config)
                    uploader.upload_md_results(job_id, work_dir)
                    
                    # 触发后处理
                    self.client.trigger_postprocess(job_id)
                    
                    # 标记完成
                    self.mark_completed(job_id, cpu_hours=cpu_hours)
                else:
                    self.update_status(
                        job_id, 'FAILED',
                        error_message="LAMMPS simulation failed"
                    )
                    self.job_tracker.remove_job(job_id)
                    
            except Exception as e:
                self.handle_error(job_id, e)
        
        elif slurm_status in ['FAILED', 'TIMEOUT', 'CANCELLED']:
            self.update_status(
                job_id, 'FAILED',
                error_message=f"Slurm job {slurm_status}"
            )
            self.job_tracker.remove_job(job_id)
    
    # ==================== RESP 相关 ====================
    
    def _start_resp_calculations(
        self,
        job_id: int,
        job: Dict,
        solvents: List[Dict],
        wrapper
    ) -> bool:
        """启动 RESP 电荷计算"""
        try:
            # 创建 RESP 工作目录
            job_name = job.get('config', {}).get('name', f'MD-{job_id}')
            resp_base_dir = self.config.work_base_path / f"RESP_{job_name}"
            resp_base_dir.mkdir(parents=True, exist_ok=True)
            
            # 获取 Slurm 配置
            md_config = job.get('config', {})
            slurm_partition = md_config.get(
                'slurm_partition',
                self.config.default_partition
            )
            
            # 初始化 RESP 计算器
            from app.workers.resp_calculator import RESPCalculator
            resp_calculator = RESPCalculator(
                charge_save_path=self.config.work_base_path / 'charges',
                slurm_partition=slurm_partition
            )
            
            resp_jobs = []
            
            for solvent in solvents:
                name = solvent['name']
                smiles = solvent.get('smiles', '')
                
                if not smiles:
                    self.logger.warning(f"No SMILES for solvent {name}, skipping RESP")
                    continue
                
                # 为每个溶剂创建工作目录
                solvent_dir = resp_base_dir / name
                solvent_dir.mkdir(parents=True, exist_ok=True)
                
                # 启动 RESP 计算
                resp_result = resp_calculator.calculate(
                    molecule_name=name,
                    smiles=smiles,
                    work_dir=solvent_dir
                )
                
                resp_jobs.append({
                    'name': name,
                    'slurm_job_id': resp_result.get('slurm_job_id'),
                    'work_dir': str(solvent_dir)
                })
            
            # 保存 RESP 任务信息
            self.job_tracker.add_job(
                job_id=job_id,
                job_type='MD',
                extra_info={
                    'phase': 'resp',
                    'resp_jobs': resp_jobs,
                    'original_job': job
                }
            )
            
            return True
            
        except Exception as e:
            self.logger.error(f"启动 RESP 失败: {e}", exc_info=True)
            self.handle_error(job_id, e)
            return False
    
    def _continue_md_job(
        self,
        job_id: int,
        job: Dict,
        wrapper
    ) -> bool:
        """继续 MD 任务（RESP 完成后或不需要 RESP）"""
        try:
            job_data = job.get('config', {})
            job_name = job_data.get('name', f'MD-{job_id}')
            
            # 生成 LAMMPS 输入文件
            work_dir = wrapper.prepare_lammps_input(
                job_id=job_id,
                job_name=job_name,
                cations=job_data.get('cations', []),
                anions=job_data.get('anions', []),
                solvents=job_data.get('solvents', []),
                additives=job_data.get('additives', []),
                box_size=job_data.get('box_size', 50),
                temperature=job_data.get('temperature', 298.15),
                pressure=job_data.get('pressure', 1.0),
                slurm_config=job_data
            )
            
            # 提交到 Slurm
            from worker.utils.slurm import SlurmManager
            slurm = SlurmManager()
            slurm_result = slurm.submit_job(work_dir)
            
            if not slurm_result['success']:
                raise Exception(slurm_result.get('error', 'Failed to submit'))
            
            slurm_job_id = slurm_result['slurm_job_id']
            
            # 标记为运行中
            self.mark_running(
                job_id,
                slurm_job_id=slurm_job_id,
                work_dir=str(work_dir)
            )
            
            self.logger.info(f"MD 任务 {job_id} 已提交 (Slurm: {slurm_job_id})")
            return True
            
        except Exception as e:
            self.handle_error(job_id, e)
            return False
    
    # ==================== 辅助方法 ====================
    
    def _init_molyte_wrapper(self):
        """初始化 MolyteWrapper"""
        from app.workers.molyte_wrapper import MolyteWrapper
        
        local_config = self.config.get('local', {})
        
        return MolyteWrapper(
            work_base_path=Path(local_config.get('work_base_path', '/tmp')),
            initial_salts_path=Path(local_config.get('initial_salts_path', '')),
            ligpargen_path=Path(local_config.get('ligpargen_path', '')),
            packmol_path=Path(local_config.get('packmol_path', '')),
            ltemplify_path=Path(local_config.get('ltemplify_path', '')),
            moltemplate_path=Path(local_config.get('moltemplate_path', '')),
            charge_save_path=Path(local_config.get('charge_save_path', '')),
        )
    
    def _check_lammps_success(self, work_dir: Path) -> bool:
        """检查 LAMMPS 是否成功完成"""
        log_file = work_dir / 'log.lammps'
        
        if not log_file.exists():
            return False
        
        with open(log_file, 'r') as f:
            content = f.read()
        
        # 检查是否有 "Total wall time" 表示正常完成
        return 'Total wall time' in content
