"""
统一的路径配置管理 (Enhanced v2.0)

新特性:
- 按任务类型和状态组织 (active/completed/archived)
- 自包含任务目录 (input/work/output/logs)
- 生命周期管理支持
- 向后兼容旧API
"""
from pydantic_settings import BaseSettings
from pathlib import Path
from typing import Optional, Literal
import os
import logging

logger = logging.getLogger(__name__)

JobType = Literal["qc", "md", "anion", "cluster"]
JobStatus = Literal["active", "completed", "archived"]
JobSubdir = Literal["input", "work", "output", "logs"]


class PathSettings(BaseSettings):
    """统一的路径配置管理"""
    
    # 基础路径
    base_dir: str = "/public/home/xiaoji/molyte_web"
    
    # ========== 核心目录 ==========
    @property
    def data_dir(self) -> Path:
        """数据根目录"""
        return Path(self.base_dir) / "data"
    
    # ========== 新的任务目录结构 ==========
    
    def get_job_dir(
        self,
        job_type: JobType,
        job_id: str,
        status: JobStatus = "active"
    ) -> Path:
        """
        获取任务目录 (新API)
        
        Args:
            job_type: 任务类型 (qc, md, anion, cluster)
            job_id: 任务ID (如 QC-12345, MD-67890)
            status: 任务状态 (active, completed, archived)
            
        Returns:
            任务根目录路径
            
        Example:
            >>> paths.get_job_dir('qc', 'QC-12345', 'active')
            Path('/public/home/xiaoji/molyte_web/data/jobs/qc/active/QC-12345')
        """
        return self.data_dir / "jobs" / job_type / status / job_id
    
    def get_job_subdir(
        self,
        job_type: JobType,
        job_id: str,
        subdir: JobSubdir,
        status: JobStatus = "active"
    ) -> Path:
        """
        获取任务子目录 (新API)
        
        Args:
            job_type: 任务类型
            job_id: 任务ID
            subdir: 子目录名称 (input, work, output, logs)
            status: 任务状态
            
        Returns:
            任务子目录路径
            
        Example:
            >>> paths.get_job_subdir('qc', 'QC-12345', 'output')
            Path('/.../data/jobs/qc/active/QC-12345/output')
        """
        return self.get_job_dir(job_type, job_id, status) / subdir
    
    def create_job_structure(
        self,
        job_type: JobType,
        job_id: str,
        status: JobStatus = "active"
    ) -> dict[str, Path]:
        """
        创建完整的任务目录结构
        
        Returns:
            包含所有子目录的字典
        """
        job_dir = self.get_job_dir(job_type, job_id, status)
        
        subdirs = {
            'root': job_dir,
            'input': job_dir / 'input',
            'work': job_dir / 'work',
            'output': job_dir / 'output',
            'logs': job_dir / 'logs'
        }
        
        # 创建所有目录
        for path in subdirs.values():
            path.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"Created job structure for {job_type} {job_id}")
        return subdirs
    
    # ========== 库目录 (Library) ==========
    
    @property
    def library_dir(self) -> Path:
        """公共资源库目录"""
        return self.data_dir / "library"
    
    @property
    def force_fields_dir(self) -> Path:
        """力场库目录"""
        return self.library_dir / "force_fields"
    
    def get_force_field_dir(self, ff_type: Literal["cations", "anions"]) -> Path:
        """获取指定类型的力场目录"""
        return self.force_fields_dir / ff_type
    
    @property
    def charges_library_dir(self) -> Path:
        """RESP电荷库目录"""
        return self.library_dir / "charges"
    
    @property
    def molecules_library_dir(self) -> Path:
        """分子结构库目录"""
        return self.library_dir / "molecules"
    
    # ========== 用户目录 ==========
    
    @property
    def users_dir(self) -> Path:
        """用户数据根目录"""
        return self.data_dir / "users"
    
    def get_user_dir(self, user_id: int) -> Path:
        """获取用户目录"""
        return self.users_dir / f"user_{user_id:06d}"
    
    def get_user_uploads_dir(self, user_id: int) -> Path:
        """获取用户上传目录"""
        return self.get_user_dir(user_id) / "uploads"
    
    # ========== 临时目录 ==========
    
    @property
    def temp_dir(self) -> Path:
        """临时文件目录"""
        return self.data_dir / "temp"
    
    @property
    def temp_uploads_dir(self) -> Path:
        """临时上传目录"""
        return self.temp_dir / "uploads_pending"
    
    @property
    def temp_coord_gen_dir(self) -> Path:
        """坐标生成临时目录"""
        return self.temp_dir / "coord_gen"
    
    # ========== 备份目录 ==========
    
    @property
    def backup_dir(self) -> Path:
        """数据备份目录"""
        return self.data_dir / "backup"
    
    @property
    def db_backup_dir(self) -> Path:
        """数据库备份目录"""
        return self.backup_dir / "database"
    
    # ========== 日志目录 ==========
    
    @property
    def logs_dir(self) -> Path:
        """日志文件目录"""
        return Path(self.base_dir) / "logs"
    
    @property
    def api_logs_dir(self) -> Path:
        """API日志目录"""
        return self.logs_dir / "api"
    
    @property
    def worker_logs_dir(self) -> Path:
        """Worker日志目录"""
        return self.logs_dir / "workers"
    
    # ========== Slurm脚本 ==========
    
    @property
    def slurm_scripts_dir(self) -> Path:
        """Slurm提交脚本目录"""
        return Path(self.base_dir) / "slurm_scripts"
    
    # ========== 向后兼容的旧API (标记为废弃) ==========
    
    @property
    def qc_work_dir(self) -> Path:
        """
        QC计算工作目录 (废弃)
        
        推荐使用: get_job_dir('qc', job_id, 'active')
        """
        return self.data_dir / "jobs" / "qc" / "active"
    
    @property
    def qc_results_dir(self) -> Path:
        """
        QC计算结果目录 (废弃)
        
        推荐使用: get_job_subdir('qc', job_id, 'output')
        """
        return self.data_dir / "qc_results"
    
    @property
    def md_work_dir(self) -> Path:
        """
        MD计算工作目录 (废弃)
        
        推荐使用: get_job_dir('md', job_id, 'active')
        """
        return self.data_dir / "jobs" / "md" / "active"
    
    @property
    def md_results_dir(self) -> Path:
        """MD计算结果目录 (废弃)"""
        return self.data_dir / "md_results"
    
    @property
    def initial_salts_dir(self) -> Path:
        """
        初始盐结构目录 (重定向到library)
        
        推荐使用: force_fields_dir
        """
        return self.force_fields_dir
    
    @property
    def charges_dir(self) -> Path:
        """
        RESP电荷数据目录 (重定向到library)
        
        推荐使用: charges_library_dir
        """
        return self.charges_library_dir
    
    @property
    def cluster_work_dir(self) -> Path:
        """Cluster分析工作目录 (废弃)"""
        return self.data_dir / "jobs" / "cluster" / "active"
    
    @property
    def uploads_dir(self) -> Path:
        """用户上传文件目录 (废弃, 现在按用户隔离)"""
        return self.data_dir / "uploads"
    
    # ========== 工具路径 ==========
    
    @property
    def resp2_script_path(self) -> Path:
        """RESP2计算脚本路径"""
        default_path = "/public/home/xiaoji/molyte_v1/RESP/RESP2.sh"
        custom_path = os.getenv("MOLYTE_RESP2_SCRIPT", default_path)
        return Path(custom_path)
    
    # ========== 辅助方法 ==========
    
    def ensure_all_dirs(self):
        """确保所有必要的目录存在"""
        essential_dirs = [
            self.data_dir,
            self.library_dir,
            self.force_fields_dir,
            self.get_force_field_dir("cations"),
            self.get_force_field_dir("anions"),
            self.charges_library_dir,
            self.molecules_library_dir,
            self.users_dir,
            self.temp_dir,
            self.temp_uploads_dir,
            self.temp_coord_gen_dir,
            self.backup_dir,
            self.db_backup_dir,
            self.logs_dir,
            self.api_logs_dir,
            self.worker_logs_dir,
            self.slurm_scripts_dir,
        ]
        
        # 创建任务类型目录
        for job_type in ["qc", "md", "anion", "cluster"]:
            for status in ["active", "completed", "archived"]:
                job_base = self.data_dir / "jobs" / job_type / status
                essential_dirs.append(job_base)
        
        for dir_path in essential_dirs:
            try:
                dir_path.mkdir(parents=True, exist_ok=True)
            except Exception as e:
                logger.warning(f"Could not create directory {dir_path}: {e}")
    
    def migrate_old_structure(self):
        """
        迁移旧目录结构到新结构
        
        警告: 这是一个破坏性操作，应该谨慎使用
        建议先备份数据
        """
        logger.warning("migrate_old_structure() should be run carefully!")
        logger.warning("Please use the migration script instead")
        # 实际迁移逻辑应该在独立的迁移脚本中
        pass
    
    class Config:
        env_prefix = "MOLYTE_PATH_"
        case_sensitive = True


# 创建全局路径配置实例
paths = PathSettings()

