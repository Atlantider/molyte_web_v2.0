"""
Worker 配置管理模块

负责加载和管理 Worker 配置
"""
import logging
from pathlib import Path
from typing import Dict, Any, Optional
import yaml


logger = logging.getLogger(__name__)


class WorkerConfig:
    """Worker 配置管理类"""
    
    def __init__(self, config_path: str = "deployment/polling_worker_config.yaml"):
        """
        初始化配置
        
        Args:
            config_path: 配置文件路径
        """
        self._config = self._load_config(config_path)
        self._setup_logging()
        logger.info(f"配置已加载: {config_path}")
    
    def _load_config(self, config_path: str) -> Dict[str, Any]:
        """加载配置文件"""
        config_file = Path(config_path)
        if not config_file.exists():
            raise FileNotFoundError(f"配置文件不存在: {config_path}")
        
        with open(config_file, 'r', encoding='utf-8') as f:
            return yaml.safe_load(f)
    
    def _setup_logging(self):
        """设置日志"""
        log_file = self.worker_config.get('log_file', 'worker.log')
        log_level = getattr(logging, self.worker_config.get('log_level', 'INFO'))
        
        logging.basicConfig(
            level=log_level,
            format='%(asctime)s | %(levelname)s | %(name)s | %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler()
            ]
        )
    
    # ==================== 属性访问器 ====================
    
    @property
    def worker_config(self) -> Dict[str, Any]:
        """Worker 配置"""
        return self._config.get('worker', {})
    
    @property
    def worker_name(self) -> str:
        """Worker 名称"""
        return self.worker_config.get('name', 'default_worker')
    
    @property
    def max_concurrent_jobs(self) -> int:
        """最大并发任务数"""
        return self.worker_config.get('max_concurrent_jobs', 5)
    
    @property
    def heartbeat_interval(self) -> int:
        """心跳间隔（秒）"""
        return self.worker_config.get('heartbeat_interval', 60)
    
    @property
    def poll_interval(self) -> int:
        """轮询间隔（秒）"""
        return self._config.get('api', {}).get('poll_interval', 30)
    
    # ==================== API 配置 ====================
    
    @property
    def api_config(self) -> Dict[str, Any]:
        """API 配置"""
        return self._config.get('api', {})
    
    @property
    def api_base_url(self) -> str:
        """API 基础 URL"""
        return self.api_config.get('base_url', '')
    
    @property
    def api_token(self) -> str:
        """API Token"""
        return self.api_config.get('worker_token', '')
    
    # ==================== COS 配置 ====================
    
    @property
    def cos_config(self) -> Dict[str, Any]:
        """腾讯云 COS 配置"""
        return self._config.get('cos', {})
    
    @property
    def has_cos(self) -> bool:
        """是否配置了 COS"""
        return 'cos' in self._config
    
    @property
    def has_oss(self) -> bool:
        """是否配置了 OSS（阿里云）"""
        return 'oss' in self._config
    
    @property
    def oss_config(self) -> Dict[str, Any]:
        """阿里云 OSS 配置"""
        return self._config.get('oss', {})
    
    # ==================== 工作目录配置 ====================
    
    @property
    def work_base_path(self) -> Path:
        """工作基础路径"""
        path = self.worker_config.get('work_base_path', '/tmp/molyte_work')
        return Path(path)
    
    @property
    def molyte_path(self) -> Optional[Path]:
        """Molyte 项目路径"""
        path = self.worker_config.get('molyte_path')
        return Path(path) if path else None
    
    # ==================== 软件路径配置 ====================
    
    @property
    def gaussian_path(self) -> Optional[str]:
        """Gaussian 可执行文件路径"""
        return self._config.get('software', {}).get('gaussian_path')
    
    @property
    def pyscf_enabled(self) -> bool:
        """是否启用 PySCF"""
        return self._config.get('software', {}).get('pyscf_enabled', True)
    
    @property
    def multiwfn_path(self) -> Optional[str]:
        """Multiwfn 可执行文件路径"""
        return self._config.get('software', {}).get('multiwfn_path')
    
    # ==================== Slurm 配置 ====================
    
    @property
    def slurm_config(self) -> Dict[str, Any]:
        """Slurm 配置"""
        return self._config.get('slurm', {})
    
    @property
    def default_partition(self) -> str:
        """默认分区"""
        return self.slurm_config.get('default_partition', 'cpu')
    
    @property
    def supported_partitions(self) -> list:
        """支持的分区列表"""
        return self.slurm_config.get('supported_partitions', ['cpu'])
    
    # ==================== 原始配置访问 ====================
    
    def get(self, key: str, default: Any = None) -> Any:
        """获取配置项"""
        return self._config.get(key, default)
    
    def __getitem__(self, key: str) -> Any:
        """字典式访问"""
        return self._config[key]
    
    def __contains__(self, key: str) -> bool:
        """检查配置项是否存在"""
        return key in self._config
