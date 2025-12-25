"""
结果上传模块

负责上传 MD/QC 计算结果到 COS
"""
import logging
from typing import Dict, Any, TYPE_CHECKING
from pathlib import Path

if TYPE_CHECKING:
    from worker.core.config import WorkerConfig


logger = logging.getLogger(__name__)


class ResultUploader:
    """结果上传器"""
    
    def __init__(self, config: 'WorkerConfig'):
        """初始化结果上传器"""
        self.config = config
        self.logger = logging.getLogger(self.__class__.__name__)
        
        from worker.uploaders.cos_client import COSClient
        self.cos = COSClient(config)
    
    def upload_md_results(self, job_id: int, work_dir: Path):
        """上传 MD 结果"""
        cos_prefix = f"md_results/{job_id}/"
        
        # 上传关键文件
        files_to_upload = [
            'log.lammps',
            'trajectory.dcd',
            'final.data',
            '*.rdf',
            '*.msd'
        ]
        
        for pattern in files_to_upload:
            for file in work_dir.glob(pattern):
                key = f"{cos_prefix}{file.name}"
                self.cos.upload_file(file, key)
        
        self.logger.info(f"MD 结果已上传: {cos_prefix}")
    
    def upload_qc_results(
        self,
        job_id: int,
        work_dir: Path,
        result: Dict[str, Any]
    ):
        """上传 QC 结果"""
        cos_prefix = f"qc_results/{job_id}/"
        
        # 上传日志文件
        for log_file in work_dir.glob('*.log'):
            key = f"{cos_prefix}{log_file.name}"
            self.cos.upload_file(log_file, key)
        
        # 上传 checkpoint
        for chk_file in work_dir.glob('*.chk'):
            key = f"{cos_prefix}{chk_file.name}"
            self.cos.upload_file(chk_file, key)
        
        # 上传结果 JSON
        self.cos.upload_json(f"{cos_prefix}result.json", result)
        
        self.logger.info(f"QC 结果已上传: {cos_prefix}")
