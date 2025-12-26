"""
数据库初始化脚本
在应用启动时运行，确保所有必要的表和列都存在
"""
import logging
from sqlalchemy import text, inspect
from app.database import engine, Base

logger = logging.getLogger(__name__)

def init_db():
    """初始化数据库"""
    try:
        logger.info("开始初始化数据库...")
        
        # 1. 创建所有表（如果不存在）
        logger.info("创建表...")
        Base.metadata.create_all(bind=engine)
        
        # 2. 检查并添加缺失的列
        logger.info("检查并添加缺失的列...")
        _add_missing_columns()
        
        logger.info("✓ 数据库初始化完成")
        return True
        
    except Exception as e:
        logger.error(f"✗ 数据库初始化失败: {e}")
        return False

def _add_missing_columns():
    """添加缺失的列"""
    try:
        inspector = inspect(engine)

        # 检查 sub_accounts 表
        sub_accounts_columns = [col['name'] for col in inspector.get_columns('sub_accounts')]

        if 'allocated_quota' not in sub_accounts_columns:
            logger.info("添加 sub_accounts.allocated_quota 列...")
            try:
                with engine.connect() as conn:
                    conn.execute(text("""
                        ALTER TABLE sub_accounts
                        ADD COLUMN allocated_quota FLOAT DEFAULT 0.0 NOT NULL
                    """))
                    conn.commit()
            except Exception as e:
                logger.warning(f"无法添加 allocated_quota 列: {e}")

        if 'allocated_used' not in sub_accounts_columns:
            logger.info("添加 sub_accounts.allocated_used 列...")
            try:
                with engine.connect() as conn:
                    conn.execute(text("""
                        ALTER TABLE sub_accounts
                        ADD COLUMN allocated_used FLOAT DEFAULT 0.0 NOT NULL
                    """))
                    conn.commit()
            except Exception as e:
                logger.warning(f"无法添加 allocated_used 列: {e}")

        # 检查 master_accounts 表
        master_accounts_columns = [col['name'] for col in inspector.get_columns('master_accounts')]

        if 'personal_quota' not in master_accounts_columns:
            logger.info("添加 master_accounts.personal_quota 列...")
            try:
                with engine.connect() as conn:
                    conn.execute(text("""
                        ALTER TABLE master_accounts
                        ADD COLUMN personal_quota FLOAT DEFAULT 0.0 NOT NULL
                    """))
                    conn.commit()
            except Exception as e:
                logger.warning(f"无法添加 personal_quota 列: {e}")

        if 'personal_used' not in master_accounts_columns:
            logger.info("添加 master_accounts.personal_used 列...")
            try:
                with engine.connect() as conn:
                    conn.execute(text("""
                        ALTER TABLE master_accounts
                        ADD COLUMN personal_used FLOAT DEFAULT 0.0 NOT NULL
                    """))
                    conn.commit()
            except Exception as e:
                logger.warning(f"无法添加 personal_used 列: {e}")

        logger.info("✓ 列检查完成")
    except Exception as e:
        logger.warning(f"列检查过程中出现错误: {e}")

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    init_db()

