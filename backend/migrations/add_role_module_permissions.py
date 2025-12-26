"""
Migration: Add role_module_permissions table
"""
import sys
import json
from pathlib import Path

# Add the backend directory to the path
sys.path.insert(0, str(Path(__file__).parent.parent))

from sqlalchemy import create_engine, Column, String, JSON, DateTime, text
from sqlalchemy.sql import func
from app.database import Base, engine
from app.models.role_permissions import RoleModulePermission, DEFAULT_ROLE_PERMISSIONS
from app.models.user import UserRole


def run_migration():
    """Run the migration to add role_module_permissions table"""

    # Create the table
    print("Creating role_module_permissions table...")
    Base.metadata.create_all(engine, tables=[RoleModulePermission.__table__])

    # Insert default role permissions
    print("Inserting default role permissions...")
    with engine.begin() as conn:
        # Check if table is empty
        result = conn.execute(text("SELECT COUNT(*) FROM role_module_permissions"))
        count = result.scalar()

        if count == 0:
            for role, config in DEFAULT_ROLE_PERMISSIONS.items():
                conn.execute(
                    text("""
                        INSERT INTO role_module_permissions (role, allowed_modules, description)
                        VALUES (:role, CAST(:allowed_modules AS jsonb), :description)
                    """),
                    {
                        "role": role.value,
                        "allowed_modules": json.dumps(config["allowed_modules"]),
                        "description": config["description"]
                    }
                )
            print(f"Inserted {len(DEFAULT_ROLE_PERMISSIONS)} default role permissions")
        else:
            print(f"Table already has {count} role permissions, skipping insertion")

    print("Migration completed successfully!")


if __name__ == "__main__":
    run_migration()

