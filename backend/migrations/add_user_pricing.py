#!/usr/bin/env python3
"""
数据库迁移脚本：添加用户级别定价功能
"""
import sys
import os

# 添加项目根目录到 Python 路径
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from app.database import SessionLocal, engine
from sqlalchemy import text
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def run_migration():
    """执行迁移"""
    db = SessionLocal()
    
    try:
        # 1. 添加用户级别定价字段到 users 表
        logger.info("添加用户级别定价字段...")
        db.execute(text("""
            ALTER TABLE users ADD COLUMN IF NOT EXISTS custom_cpu_hour_price FLOAT;
        """))
        db.execute(text("""
            ALTER TABLE users ADD COLUMN IF NOT EXISTS price_updated_at TIMESTAMP WITH TIME ZONE;
        """))
        db.execute(text("""
            ALTER TABLE users ADD COLUMN IF NOT EXISTS price_updated_by INTEGER;
        """))
        
        # 2. 添加外键约束
        logger.info("添加外键约束...")
        try:
            db.execute(text("""
                ALTER TABLE users ADD CONSTRAINT fk_price_updated_by 
                  FOREIGN KEY (price_updated_by) REFERENCES users(id) ON DELETE SET NULL;
            """))
        except Exception as e:
            logger.warning(f"外键约束可能已存在: {e}")
        
        # 3. 创建索引
        logger.info("创建索引...")
        db.execute(text("""
            CREATE INDEX IF NOT EXISTS idx_users_custom_price ON users(custom_cpu_hour_price) 
              WHERE custom_cpu_hour_price IS NOT NULL;
        """))
        
        # 4. 添加系统配置
        logger.info("添加系统配置...")
        db.execute(text("""
            INSERT INTO system_configs (key, value, description) 
            VALUES ('user_pricing_enabled', 'true', '是否启用用户级别定价功能')
            ON CONFLICT (key) DO NOTHING;
        """))
        
        # 5. 创建审计日志表
        logger.info("创建审计日志表...")
        db.execute(text("""
            CREATE TABLE IF NOT EXISTS pricing_audit_log (
                id SERIAL PRIMARY KEY,
                admin_id INTEGER NOT NULL REFERENCES users(id) ON DELETE CASCADE,
                user_id INTEGER NOT NULL REFERENCES users(id) ON DELETE CASCADE,
                old_price FLOAT,
                new_price FLOAT,
                reason VARCHAR(255),
                created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
                
                CONSTRAINT check_prices CHECK (old_price IS NULL OR old_price >= 0),
                CONSTRAINT check_new_price CHECK (new_price IS NULL OR new_price >= 0)
            );
        """))
        
        # 6. 创建审计日志索引
        logger.info("创建审计日志索引...")
        db.execute(text("""
            CREATE INDEX IF NOT EXISTS idx_pricing_audit_user ON pricing_audit_log(user_id);
        """))
        db.execute(text("""
            CREATE INDEX IF NOT EXISTS idx_pricing_audit_admin ON pricing_audit_log(admin_id);
        """))
        db.execute(text("""
            CREATE INDEX IF NOT EXISTS idx_pricing_audit_created ON pricing_audit_log(created_at);
        """))
        
        db.commit()
        logger.info("✅ 迁移完成！")
        return True
        
    except Exception as e:
        logger.error(f"❌ 迁移失败: {e}")
        db.rollback()
        return False
    finally:
        db.close()

if __name__ == "__main__":
    success = run_migration()
    sys.exit(0 if success else 1)

