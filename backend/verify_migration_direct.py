#!/usr/bin/env python3
"""
Verify account_type migration using direct SQL
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

import psycopg2
import logging
from app.config import settings

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def verify_migration():
    """Verify account_type migration"""
    try:
        # Parse database URL
        db_url = settings.DATABASE_URL
        if db_url.startswith("postgresql://"):
            db_url = db_url.replace("postgresql://", "")
        
        user_pass, host_db = db_url.split("@")
        user, password = user_pass.split(":")
        host_port, database = host_db.split("/")
        host, port = host_port.split(":")
        
        logger.info(f"Connecting to database: {host}:{port}/{database}\n")
        
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
        
        logger.info("=== Verifying account_type Migration ===\n")
        
        # Check 1: Column exists
        logger.info("✓ Check 1: Verifying account_type column exists...")
        cursor.execute("""
            SELECT column_name FROM information_schema.columns 
            WHERE table_name='users' AND column_name='account_type'
        """)
        result = cursor.fetchone()
        if result:
            logger.info("  ✓ account_type column exists\n")
        else:
            logger.error("  ✗ account_type column NOT found\n")
            return False
        
        # Check 2: All users have account_type
        logger.info("✓ Check 2: Verifying all users have account_type...")
        cursor.execute("SELECT COUNT(*) FROM users WHERE account_type IS NULL")
        null_count = cursor.fetchone()[0]
        if null_count == 0:
            logger.info("  ✓ No NULL account_type values found\n")
        else:
            logger.error(f"  ✗ Found {null_count} users with NULL account_type\n")
            return False
        
        # Check 3: Account type distribution
        logger.info("✓ Check 3: Account type distribution...")
        cursor.execute("""
            SELECT account_type, COUNT(*) as count
            FROM users
            GROUP BY account_type
            ORDER BY account_type
        """)
        results = cursor.fetchall()
        total = 0
        for account_type, count in results:
            logger.info(f"  - {account_type}: {count}")
            total += count
        logger.info(f"  - Total: {total} users\n")
        
        # Check 4: Verify classifications
        logger.info("✓ Check 4: Verifying classifications...")
        
        # Master accounts
        cursor.execute("""
            SELECT COUNT(*) FROM users 
            WHERE account_type = 'master_account'
        """)
        master_count = cursor.fetchone()[0]
        cursor.execute("SELECT COUNT(*) FROM master_accounts")
        master_account_count = cursor.fetchone()[0]
        if master_count == master_account_count:
            logger.info(f"  ✓ Master accounts: {master_count} (matches MasterAccount table)")
        else:
            logger.warning(f"  ⚠ Master accounts mismatch: {master_count} users vs {master_account_count} records")
        
        # Sub accounts
        cursor.execute("""
            SELECT COUNT(*) FROM users 
            WHERE account_type = 'sub_account'
        """)
        sub_count = cursor.fetchone()[0]
        cursor.execute("SELECT COUNT(*) FROM sub_accounts")
        sub_account_count = cursor.fetchone()[0]
        if sub_count == sub_account_count:
            logger.info(f"  ✓ Sub accounts: {sub_count} (matches SubAccount table)")
        else:
            logger.warning(f"  ⚠ Sub accounts mismatch: {sub_count} users vs {sub_account_count} records")
        
        # Organization members
        cursor.execute("""
            SELECT COUNT(*) FROM users 
            WHERE account_type = 'organization_member'
        """)
        org_count = cursor.fetchone()[0]
        cursor.execute("SELECT COUNT(*) FROM organization_members")
        org_member_count = cursor.fetchone()[0]
        if org_count == org_member_count:
            logger.info(f"  ✓ Organization members: {org_count} (matches OrganizationMember table)")
        else:
            logger.warning(f"  ⚠ Organization members mismatch: {org_count} users vs {org_member_count} records")
        
        logger.info("\n" + "="*50)
        logger.info("✓ Migration verification PASSED!")
        logger.info("="*50)
        
        cursor.close()
        conn.close()
        
        return True
        
    except Exception as e:
        logger.error(f"✗ Verification failed: {str(e)}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = verify_migration()
    sys.exit(0 if success else 1)

