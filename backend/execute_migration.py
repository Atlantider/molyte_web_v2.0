#!/usr/bin/env python3
"""
Execute Phase 1 migration: Add account_type field to users table
"""
import sys
import os
import logging
from pathlib import Path

# Add backend to path
sys.path.insert(0, str(Path(__file__).parent))

import psycopg2
from psycopg2 import sql
from app.config import settings

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def execute_migration():
    """Execute the SQL migration using psycopg2"""
    try:
        # Parse database URL
        # Format: postgresql://user:password@host:port/database
        db_url = settings.DATABASE_URL
        if db_url.startswith("postgresql://"):
            db_url = db_url.replace("postgresql://", "")

        # Extract connection parameters
        user_pass, host_db = db_url.split("@")
        user, password = user_pass.split(":")
        host_port, database = host_db.split("/")
        host, port = host_port.split(":")

        logger.info(f"Connecting to database: {host}:{port}/{database}")

        # Connect to database
        conn = psycopg2.connect(
            host=host,
            port=int(port),
            database=database,
            user=user,
            password=password
        )
        conn.autocommit = True
        cursor = conn.cursor()

        logger.info("✓ Connected to database successfully")

        # Step 1: Create ENUM type
        logger.info("Step 1: Creating account_type_enum type...")
        try:
            cursor.execute("""
                CREATE TYPE account_type_enum AS ENUM (
                    'personal',
                    'organization_member',
                    'master_account',
                    'sub_account'
                );
            """)
            logger.info("✓ account_type_enum type created")
        except psycopg2.errors.DuplicateObject:
            logger.info("⚠ account_type_enum type already exists")

        # Step 2: Add account_type column
        logger.info("Step 2: Adding account_type column to users table...")
        try:
            cursor.execute("""
                ALTER TABLE users
                ADD COLUMN IF NOT EXISTS account_type account_type_enum DEFAULT 'personal' NOT NULL;
            """)
            logger.info("✓ account_type column added")
        except Exception as e:
            logger.warning(f"⚠ Column addition warning: {str(e)}")

        # Step 3: Create index
        logger.info("Step 3: Creating index on account_type...")
        try:
            cursor.execute("""
                CREATE INDEX IF NOT EXISTS idx_users_account_type ON users(account_type);
            """)
            logger.info("✓ Index created")
        except Exception as e:
            logger.warning(f"⚠ Index creation warning: {str(e)}")

        # Step 4: Verify migration
        logger.info("Step 4: Verifying migration...")
        cursor.execute("""
            SELECT COUNT(*) as total_users,
                   account_type,
                   COUNT(*) as count
            FROM users
            GROUP BY account_type;
        """)
        results = cursor.fetchall()
        logger.info("✓ User distribution by account_type:")
        for row in results:
            logger.info(f"  - {row[1]}: {row[2]} users")

        cursor.close()
        conn.close()

        logger.info("✓ Phase 1 migration completed successfully!")
        return True

    except Exception as e:
        logger.error(f"✗ Migration failed: {str(e)}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = execute_migration()
    sys.exit(0 if success else 1)

