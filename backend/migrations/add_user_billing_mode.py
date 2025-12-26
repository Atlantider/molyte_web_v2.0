"""
Migration script: Add billing_mode and custom_task_prices to users table
Run with: python -m migrations.add_user_billing_mode
"""
from sqlalchemy import text
from app.database import engine

def run_migration():
    """Add billing_mode and custom_task_prices columns to users table."""
    with engine.begin() as conn:
        # Add billing_mode column
        conn.execute(text("""
            ALTER TABLE users 
            ADD COLUMN IF NOT EXISTS billing_mode VARCHAR(20) NOT NULL DEFAULT 'CORE_HOUR'
        """))
        print("✓ Added billing_mode column")
        
        # Add custom_task_prices column
        conn.execute(text("""
            ALTER TABLE users 
            ADD COLUMN IF NOT EXISTS custom_task_prices JSONB
        """))
        print("✓ Added custom_task_prices column")
        
        print("Migration completed successfully!")

if __name__ == "__main__":
    run_migration()
