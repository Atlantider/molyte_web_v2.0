#!/usr/bin/env python3
"""
Classify existing users into account types using direct SQL
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

import psycopg2
import logging
from app.config import settings

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def classify_users():
    """Classify all existing users into account types"""
    try:
        # Parse database URL
        db_url = settings.DATABASE_URL
        if db_url.startswith("postgresql://"):
            db_url = db_url.replace("postgresql://", "")
        
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
        
        # Step 1: Classify master account owners
        logger.info("Step 1: Classifying master account owners...")
        cursor.execute("""
            UPDATE users
            SET account_type = 'master_account'
            WHERE id IN (
                SELECT DISTINCT user_id FROM master_accounts
            )
            AND account_type = 'personal';
        """)
        master_count = cursor.rowcount
        logger.info(f"✓ Classified {master_count} master account owners")
        
        # Step 2: Classify sub account owners
        logger.info("Step 2: Classifying sub account owners...")
        cursor.execute("""
            UPDATE users
            SET account_type = 'sub_account'
            WHERE id IN (
                SELECT DISTINCT user_id FROM sub_accounts
            )
            AND account_type = 'personal';
        """)
        sub_count = cursor.rowcount
        logger.info(f"✓ Classified {sub_count} sub account owners")
        
        # Step 3: Classify organization members
        logger.info("Step 3: Classifying organization members...")
        cursor.execute("""
            UPDATE users
            SET account_type = 'organization_member'
            WHERE id IN (
                SELECT DISTINCT user_id FROM organization_members
            )
            AND account_type = 'personal';
        """)
        org_count = cursor.rowcount
        logger.info(f"✓ Classified {org_count} organization members")
        
        # Step 4: Verify classification
        logger.info("Step 4: Verifying classification...")
        cursor.execute("""
            SELECT account_type, COUNT(*) as count
            FROM users
            GROUP BY account_type
            ORDER BY account_type;
        """)
        results = cursor.fetchall()
        logger.info("✓ User distribution by account_type:")
        total = 0
        for account_type, count in results:
            logger.info(f"  - {account_type}: {count} users")
            total += count
        logger.info(f"  - Total: {total} users")
        
        cursor.close()
        conn.close()
        
        logger.info("✓ User classification completed successfully!")
        return True
        
    except Exception as e:
        logger.error(f"✗ Classification failed: {str(e)}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = classify_users()
    sys.exit(0 if success else 1)

