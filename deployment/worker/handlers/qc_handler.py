"""
QC 任务处理器

处理 Gaussian/PySCF 量子化学计算任务
"""
import logging
import time
from typing import Dict, Any, Optional
from pathlib import Path

from worker.handlers.base import BaseHandler


logger = logging.getLogger(__name__)


class QCHandler(BaseHandler):
    """QC 任务处理器"""
    
    JOB_TYPE = "QC"
    
    def process(self, job: Dict[str, Any]) -> bool:
        """
        处理 QC 任务
        
        Args:
            job: 任务信息
            
        Returns:
            是否成功开始处理
        """
        job_id = job['id']
        self.logger.info(f"开始处理 QC 任务 {job_id}")
        
        try:
            # 检查是否已在处理中
            if self.job_tracker.is_job_running(job_id):
                self.logger.info(f"QC 任务 {job_id} 已在处理中，跳过")
                return False
            
            # 更新状态为 QUEUED
            self.update_status(job_id, 'QUEUED')
            
            # 获取任务配置
            config = job.get('config', {})
            molecule_name = config.get('molecule_name', f'QC_{job_id}')
            xyz_content = config.get('xyz_content', '')
            smiles = config.get('smiles', '')
            charge = config.get('charge', 0)
            spin_multiplicity = config.get('spin_multiplicity', 1)
            basis_set = config.get('basis_set', '6-31++g(d,p)')
            functional = config.get('functional', 'B3LYP')
            solvent_model = config.get('solvent_model', 'gas')
            solvent_name = config.get('solvent_name', '')
            slurm_partition = config.get('slurm_partition', 'hpc128c')
            slurm_cpus = config.get('slurm_cpus', 16)
            slurm_time = config.get('slurm_time', 7200)
            
            # 创建工作目录
            work_dir = self._create_work_dir(job_id, molecule_name, config)
            
            # 生成 Gaussian 输入文件
            gjf_path = self._generate_gaussian_input(
                work_dir, molecule_name,
                xyz_content=xyz_content,
                smiles=smiles,
                charge=charge,
                spin_multiplicity=spin_multiplicity,
                functional=functional,
                basis_set=basis_set,
                solvent_model=solvent_model,
                solvent_name=solvent_name,
                nprocs=slurm_cpus
            )
            
            # 生成 Slurm 作业脚本
            job_script = self._generate_job_script(
                work_dir, molecule_name,
                slurm_partition, slurm_cpus, slurm_time
            )
            
            # 提交到 Slurm
            from worker.utils.slurm import SlurmManager
            slurm = SlurmManager()
            slurm_result = slurm.submit_job(work_dir)
            
            if not slurm_result['success']:
                raise Exception(slurm_result.get('error', 'Failed to submit to Slurm'))
            
            slurm_job_id = slurm_result['slurm_job_id']
            
            # 标记为运行中
            self.mark_running(
                job_id,
                slurm_job_id=slurm_job_id,
                work_dir=str(work_dir)
            )
            
            self.logger.info(f"QC 任务 {job_id} 已提交到 Slurm (Job ID: {slurm_job_id})")
            return True
            
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
        处理 QC 任务完成
        
        Args:
            job_id: 任务 ID
            job_info: 任务信息
            slurm_status: Slurm 状态
        """
        work_dir = Path(job_info.get('work_dir', ''))
        slurm_job_id = job_info.get('slurm_job_id')
        
        self.logger.info(
            f"QC 任务 {job_id} Slurm 状态: {slurm_status}"
        )
        
        if slurm_status == 'COMPLETED':
            try:
                # 解析结果
                result = self._parse_gaussian_output(work_dir)
                
                if result.get('success'):
                    # 获取 CPU 核时
                    from worker.utils.slurm import SlurmManager
                    slurm = SlurmManager()
                    cpu_hours = slurm.get_job_cpu_hours(slurm_job_id)
                    
                    # 上传结果到 COS
                    from worker.uploaders import ResultUploader
                    uploader = ResultUploader(self.config)
                    uploader.upload_qc_results(job_id, work_dir, result)
                    
                    # 提交结果到后端
                    self.client.submit_qc_result(job_id, result)
                    
                    # 标记完成
                    self.mark_completed(job_id, cpu_hours=cpu_hours)
                    
                else:
                    # 计算失败
                    error_msg = result.get('error', 'Unknown Gaussian error')
                    self._handle_qc_failure(job_id, error_msg, work_dir)
                    
            except Exception as e:
                self._handle_qc_failure(job_id, str(e), work_dir)
        
        elif slurm_status in ['FAILED', 'TIMEOUT', 'CANCELLED']:
            self._handle_qc_failure(
                job_id,
                f"Slurm job {slurm_status}",
                work_dir
            )
    
    def _handle_qc_failure(
        self,
        job_id: int,
        error_message: str,
        work_dir: Path
    ):
        """
        处理 QC 失败（支持自动重试协调）
        
        云端会检查是否需要重试，如果需要会重置状态为 SUBMITTED
        """
        self.logger.error(f"QC 任务 {job_id} 失败: {error_message}")
        
        # 更新状态为失败（云端可能会触发重试）
        self.update_status(
            job_id, 'FAILED',
            error_message=error_message
        )
        
        # 从追踪器移除
        self.job_tracker.remove_job(job_id)
    
    # ==================== 辅助方法 ====================
    
    def _create_work_dir(
        self,
        job_id: int,
        molecule_name: str,
        config: Dict
    ) -> Path:
        """创建工作目录"""
        is_desolvation = bool(config.get('desolvation_job_type', ''))
        task_type = config.get('task_type', '')
        
        dir_name = f"QC-{job_id}-{molecule_name}"
        
        if is_desolvation and not task_type:
            base_path = self.config.work_base_path / 'cluster_qc'
        else:
            base_path = self.config.work_base_path / 'qc'
        
        work_dir = base_path / dir_name
        work_dir.mkdir(parents=True, exist_ok=True)
        
        return work_dir
    
    def _generate_gaussian_input(
        self,
        work_dir: Path,
        molecule_name: str,
        xyz_content: str = '',
        smiles: str = '',
        charge: int = 0,
        spin_multiplicity: int = 1,
        functional: str = 'B3LYP',
        basis_set: str = '6-31++g(d,p)',
        solvent_model: str = 'gas',
        solvent_name: str = '',
        nprocs: int = 16
    ) -> Path:
        """生成 Gaussian 输入文件"""
        from worker.utils.gaussian import GaussianUtils
        
        utils = GaussianUtils()
        gjf_path = work_dir / f"{self._sanitize_filename(molecule_name)}.gjf"
        
        utils.generate_input_file(
            output_path=gjf_path,
            molecule_name=molecule_name,
            xyz_content=xyz_content,
            smiles=smiles,
            charge=charge,
            spin_multiplicity=spin_multiplicity,
            functional=functional,
            basis_set=basis_set,
            solvent_model=solvent_model,
            solvent_name=solvent_name,
            nprocs=nprocs
        )
        
        return gjf_path
    
    def _generate_job_script(
        self,
        work_dir: Path,
        molecule_name: str,
        partition: str,
        cpus: int,
        time_limit: int
    ) -> Path:
        """生成 Slurm 作业脚本"""
        from worker.utils.slurm import SlurmManager
        
        slurm = SlurmManager()
        job_script = work_dir / "job.sh"
        
        slurm.generate_qc_job_script(
            output_path=job_script,
            job_name=f"QC-{self._sanitize_filename(molecule_name)}",
            partition=partition,
            cpus=cpus,
            time_limit=time_limit,
            work_dir=work_dir
        )
        
        return job_script
    
    def _parse_gaussian_output(self, work_dir: Path) -> Dict[str, Any]:
        """解析 Gaussian 输出"""
        from worker.utils.gaussian import GaussianUtils
        
        utils = GaussianUtils()
        return utils.parse_output(work_dir)
    
    def _sanitize_filename(self, name: str) -> str:
        """清理文件名中的特殊字符"""
        import re
        return re.sub(r'[^\w\-_.]', '_', name)
