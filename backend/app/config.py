"""
Application configuration using Pydantic Settings
"""
from pydantic_settings import BaseSettings
from typing import List


class Settings(BaseSettings):
    """Application settings"""

    # Database
    DATABASE_URL: str = "postgresql://molyte_user:molyte2025@localhost:5432/molyte_web"

    # Security
    SECRET_KEY: str = "09d25e094faa6ca2556c818166b7a9563b93f7099f6f0f4caa6cf63b88e8d3e7"
    ALGORITHM: str = "HS256"
    ACCESS_TOKEN_EXPIRE_MINUTES: int = 1440  # 24小时 (原来是30分钟)

    # Application
    APP_NAME: str = "Molyte Web API"
    APP_VERSION: str = "1.0.0"
    DEBUG: bool = True

    # CORS - stored as comma-separated string in .env
    CORS_ORIGINS: str = "http://localhost:3000,http://localhost:5173,http://localhost:8080"

    # Logging
    LOG_LEVEL: str = "INFO"

    # Redis settings
    REDIS_HOST: str = "localhost"
    REDIS_PORT: int = 6379
    REDIS_PASSWORD: str = ""
    REDIS_DB: int = 0

    # Celery settings
    CELERY_BROKER_URL: str = "redis://localhost:6379/0"
    CELERY_RESULT_BACKEND: str = "redis://localhost:6379/1"

    # File paths - 使用统一路径配置
    # 保留这些属性以兼容现有代码,但从paths获取值
    @property
    def QC_WORK_DIR(self) -> str:
        from app.core.paths import paths
        return str(paths.qc_work_dir)
    
    @property
    def INITIAL_SALTS_DIR(self) -> str:
        from app.core.paths import paths
        return str(paths.initial_salts_dir)

    # SMS Settings (Aliyun)
    ALIYUN_SMS_ACCESS_KEY: str = ""
    ALIYUN_SMS_ACCESS_SECRET: str = ""
    ALIYUN_SMS_SIGN_NAME: str = "Molyte平台"
    ALIYUN_SMS_TEMPLATE_CODE: str = ""

    # SMS Settings (Tencent Cloud)
    TENCENT_SMS_APP_ID: str = ""
    TENCENT_SMS_APP_KEY: str = ""
    TENCENT_SMS_SIGN_NAME: str = "Molyte平台"
    TENCENT_SMS_TEMPLATE_ID: str = ""

    # Payment Settings (WeChat Pay)
    WECHAT_APP_ID: str = ""
    WECHAT_MCH_ID: str = ""  # 商户号
    WECHAT_API_KEY: str = ""  # API密钥
    WECHAT_API_V3_KEY: str = ""  # APIv3密钥
    WECHAT_SERIAL_NO: str = ""  # 证书序列号
    WECHAT_PRIVATE_KEY: str = ""  # 商户私钥路径
    WECHAT_NOTIFY_URL: str = "https://www.molyte.xyz/api/v1/billing/callback/wechat"

    # Payment Settings (Alipay)
    ALIPAY_APP_ID: str = ""
    ALIPAY_PRIVATE_KEY: str = ""  # 应用私钥
    ALIPAY_PUBLIC_KEY: str = ""  # 支付宝公钥
    ALIPAY_NOTIFY_URL: str = "https://www.molyte.xyz/api/v1/billing/callback/alipay"
    ALIPAY_RETURN_URL: str = "https://www.molyte.xyz/workspace/account-center?tab=recharge"

    @property
    def cors_origins_list(self) -> List[str]:
        """Convert CORS_ORIGINS string to list"""
        return [origin.strip() for origin in self.CORS_ORIGINS.split(",")]

    class Config:
        env_file = ".env"
        case_sensitive = True


# Create global settings instance
settings = Settings()

