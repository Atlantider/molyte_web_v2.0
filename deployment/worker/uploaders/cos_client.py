"""
COS 客户端模块

腾讯云对象存储客户端
"""
import logging
from typing import Dict, Any, Optional, TYPE_CHECKING
from pathlib import Path

if TYPE_CHECKING:
    from worker.core.config import WorkerConfig


logger = logging.getLogger(__name__)


class COSClient:
    """腾讯云 COS 客户端"""
    
    def __init__(self, config: 'WorkerConfig'):
        """初始化 COS 客户端"""
        self.config = config
        self.logger = logging.getLogger(self.__class__.__name__)
        self._client = None
        self._bucket = None
        self._init_client()
    
    def _init_client(self):
        """初始化 COS 客户端"""
        if self.config.has_cos:
            try:
                from qcloud_cos import CosConfig, CosS3Client
                
                cos_config = self.config.cos_config
                config = CosConfig(
                    Region=cos_config['region'],
                    SecretId=cos_config['secret_id'],
                    SecretKey=cos_config['secret_key'],
                    Scheme='https'
                )
                self._client = CosS3Client(config)
                self._bucket = cos_config['bucket']
                self.logger.info(f"COS 客户端已初始化 (Bucket: {self._bucket})")
            except Exception as e:
                self.logger.error(f"COS 客户端初始化失败: {e}")
        elif self.config.has_oss:
            self.logger.info("使用 OSS 模式（兼容）")
    
    def upload_file(
        self,
        local_path: Path,
        key: str,
        content_type: Optional[str] = None
    ) -> bool:
        """
        上传文件到 COS
        
        Args:
            local_path: 本地文件路径
            key: COS 对象键
            content_type: 内容类型
            
        Returns:
            是否成功
        """
        if not self._client:
            self.logger.warning("COS 客户端未初始化")
            return False
        
        try:
            self._client.upload_file(
                Bucket=self._bucket,
                Key=key,
                LocalFilePath=str(local_path)
            )
            self.logger.debug(f"上传文件成功: {key}")
            return True
        except Exception as e:
            self.logger.error(f"上传文件失败 {key}: {e}")
            return False
    
    def upload_json(self, key: str, data: Dict) -> bool:
        """上传 JSON 数据"""
        import json
        
        if not self._client:
            return False
        
        try:
            body = json.dumps(data, ensure_ascii=False, indent=2)
            self._client.put_object(
                Bucket=self._bucket,
                Key=key,
                Body=body.encode('utf-8'),
                ContentType='application/json'
            )
            return True
        except Exception as e:
            self.logger.error(f"上传 JSON 失败 {key}: {e}")
            return False
    
    def download_file(self, key: str, local_path: Path) -> bool:
        """下载文件"""
        if not self._client:
            return False
        
        try:
            self._client.download_file(
                Bucket=self._bucket,
                Key=key,
                DestFilePath=str(local_path)
            )
            return True
        except Exception as e:
            self.logger.error(f"下载文件失败 {key}: {e}")
            return False
