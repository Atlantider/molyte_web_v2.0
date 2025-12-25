"""
QC后处理辅助函数 - 重试和恢复机制
"""
import logging
from pathlib import Path
from typing import Optional, Dict, Any
from datetime import datetime

logger = logging.getLogger(__name__)


def can_retry_postprocess(job) -> bool:
    """
    检查QC任务是否可以重试后处理
    
    Args:
        job: QCJob对象
        
    Returns:
        bool: 是否可以重试
    """
    # 检查重试次数
    if job.retry_count >= job.max_retries:
        logger.warning(f"Job {job.id} exceeded max retries ({job.max_retries})")
        return False
    
    # 检查工作目录是否存在
    if not job.work_dir or not Path(job.work_dir).exists():
        logger.warning(f"Job {job.id} work directory not found: {job.work_dir}")
        return False
    
    return True


def find_output_files(work_dir: Path, engine: str) -> Dict[str, Optional[Path]]:
    """
    在工作目录中查找输出文件
    
    Args:
        work_dir: 工作目录
        engine: QC引擎类型 (gaussian, pyscf)
        
    Returns:
        包含找到的文件路径的字典
    """
    files = {
        'log_file': None,
        'fchk_file': None,
        'results_file': None
    }
    
    if engine == 'gaussian':
        # 查找Gaussian输出文件
        # 优先查找 _out.log 文件
        log_files = list(work_dir.glob("*_out.log"))
        if not log_files:
            log_files = list(work_dir.glob("*.log"))
        
        if log_files:
            # 排除qc_out.log和qc_err.log（Slurm输出）
            log_files = [f for f in log_files if f.name not in ['qc_out.log', 'qc_err.log']]
            if log_files:
                files['log_file'] = log_files[0]
        
        # 查找fchk文件
        fchk_files = list(work_dir.glob("*.fchk"))
        if fchk_files:
            files['fchk_file'] = fchk_files[0]
    
    elif engine == 'pyscf':
        # 查找PySCF结果文件
        results_file = work_dir / "pyscf_results.json"
        if results_file.exists():
            files['results_file'] = results_file
        
        # 查找计算日志
        log_file = work_dir / "pyscf_calculation.log"
        if log_file.exists():
            files['log_file'] = log_file
    
    return files


def auto_resubmit_failed_job(job, db_session):
    """
    自动重新提交失败的QC任务
    
    Args:
        job: QCJob对象
        db_session: 数据库会话
        
    Returns:
        bool: 是否成功重新提交
    """
    try:
        from app.models.qc import QCJobStatus
        from app.tasks.qc_submission import submit_qc_job_task
        
        logger.info(f"Auto-resubmitting failed job {job.id}")
        
        # 重置状态
        job.status = QCJobStatus.CREATED
        job.retry_count += 1
        job.error_message = None
        job.slurm_job_id = None
        
        db_session.commit()
        
        # 异步提交任务
        submit_qc_job_task.delay(job.id)
        
        logger.info(f"Successfully resubmitted job {job.id} (retry {job.retry_count}/{job.max_retries})")
        return True
        
    except Exception as e:
        logger.error(f"Failed to auto-resubmit job {job.id}: {e}")
        return False


def extract_partial_results(work_dir: Path, engine: str) -> Optional[Dict[str, Any]]:
    """
    从部分完成的计算中提取可用结果
    
    对于未完全收敛但有能量输出的情况，尝试提取部分结果
    
    Args:
        work_dir: 工作目录
        engine: QC引擎
        
    Returns:
        部分结果字典或None
    """
    files = find_output_files(work_dir, engine)
    
    if engine == 'gaussian' and files['log_file']:
        try:
            from app.tasks.qc_postprocess import extract_gaussian_results
            results = extract_gaussian_results(str(files['log_file']))
            
            # 如果至少有能量值，认为可以使用
            if results.get('energy_au') is not None:
                logger.info(f"Extracted partial results from {files['log_file']}")
                return results
        except Exception as e:
            logger.warning(f"Failed to extract partial results: {e}")
    
    elif engine == 'pyscf' and files['results_file']:
        try:
            import json
            with open(files['results_file']) as f:
                results = json.load(f)
            
            if results.get('energy_au'):
                logger.info(f"Extracted PySCF results from {files['results_file']}")
                return results
        except Exception as e:
            logger.warning(f"Failed to extract PySCF results: {e}")
    
    return None
