"""
简化账户体系迁移脚本
目标：去掉组织管理，只保留主账号和子账号的简单关系
"""
import sys
import os
import logging
from pathlib import Path

# Add backend to path
sys.path.insert(0, str(Path(__file__).parent.parent))

import psycopg2
from app.config import settings

# 配置日志
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def migrate():
    """执行迁移"""
    try:
        # 解析数据库URL
        db_url = settings.DATABASE_URL
        if db_url.startswith("postgresql://"):
            db_url = db_url.replace("postgresql://", "")

        user_pass, host_db = db_url.split("@")
        user, password = user_pass.split(":")
        host_port, database = host_db.split("/")
        host, port = host_port.split(":")

        logger.info(f"连接到数据库: {host}:{port}/{database}")

        # 连接到数据库
        conn = psycopg2.connect(
            host=host,
            port=int(port),
            database=database,
            user=user,
            password=password
        )
        conn.autocommit = True
        cursor = conn.cursor()

        logger.info("开始执行账户体系简化迁移...")

        # 1. 为 sub_accounts 表添加分配配额字段
        logger.info("为 sub_accounts 表添加分配配额字段...")
        cursor.execute("""
            ALTER TABLE sub_accounts ADD COLUMN IF NOT EXISTS allocated_quota FLOAT DEFAULT 0.0 NOT NULL;
        """)
        cursor.execute("""
            ALTER TABLE sub_accounts ADD COLUMN IF NOT EXISTS allocated_used FLOAT DEFAULT 0.0 NOT NULL;
        """)

        # 2. 为 master_accounts 表添加个人配额字段
        logger.info("为 master_accounts 表添加个人配额字段...")
        cursor.execute("""
            ALTER TABLE master_accounts ADD COLUMN IF NOT EXISTS personal_quota FLOAT DEFAULT 0.0 NOT NULL;
        """)
        cursor.execute("""
            ALTER TABLE master_accounts ADD COLUMN IF NOT EXISTS personal_used FLOAT DEFAULT 0.0 NOT NULL;
        """)

        # 3. 创建索引
        logger.info("创建索引...")
        cursor.execute("""
            CREATE INDEX IF NOT EXISTS idx_sub_accounts_master_account_id ON sub_accounts(master_account_id);
        """)
        cursor.execute("""
            CREATE INDEX IF NOT EXISTS idx_sub_accounts_user_id ON sub_accounts(user_id);
        """)
        cursor.execute("""
            CREATE INDEX IF NOT EXISTS idx_master_accounts_user_id ON master_accounts(user_id);
        """)

        logger.info("✓ 迁移完成！")
        cursor.close()
        conn.close()
        return True

    except Exception as e:
        logger.error(f"✗ 迁移失败: {e}")
        return False

if __name__ == "__main__":
    success = migrate()
    sys.exit(0 if success else 1)

