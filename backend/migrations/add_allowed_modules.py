"""
Migration script to add allowed_modules field to users table
Run this script to add the new column to existing databases
"""

import sys
import os

# Add the backend directory to the path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/..')

from sqlalchemy import text
from app.database import engine

def migrate():
    """Add allowed_modules column to users table"""
    with engine.connect() as connection:
        try:
            # Check if column already exists
            result = connection.execute(
                text("SELECT column_name FROM information_schema.columns WHERE table_name='users' AND column_name='allowed_modules'")
            )
            if result.fetchone():
                print("✓ Column 'allowed_modules' already exists")
                return
            
            # Add the column
            connection.execute(
                text("""
                    ALTER TABLE users 
                    ADD COLUMN allowed_modules JSON DEFAULT '["electrolytes", "md", "analysis", "qc"]'
                """)
            )
            connection.commit()
            print("✓ Successfully added 'allowed_modules' column to users table")
            print("  Default value: ['electrolytes', 'md', 'analysis', 'qc']")
            
        except Exception as e:
            print(f"✗ Error adding column: {e}")
            sys.exit(1)

if __name__ == '__main__':
    print("Running migration: add_allowed_modules...")
    migrate()
    print("Migration completed!")

