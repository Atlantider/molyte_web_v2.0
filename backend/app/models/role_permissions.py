"""
Role-based module permissions model
"""
from sqlalchemy import Column, String, JSON, DateTime, Enum
from sqlalchemy.sql import func
from app.database import Base
from app.models.user import UserRole


class RoleModulePermission(Base):
    """
    Role-based module permission template
    Stores which modules are allowed for each user role
    """
    __tablename__ = "role_module_permissions"

    role = Column(String(50), primary_key=True, index=True)
    allowed_modules = Column(JSON, nullable=False, default=lambda: [])
    description = Column(String(500), nullable=True)
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    updated_at = Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now())

    def __repr__(self):
        return f"<RoleModulePermission(role={self.role}, modules={self.allowed_modules})>"


# Default role permissions
DEFAULT_ROLE_PERMISSIONS = {
    UserRole.ADMIN: {
        "allowed_modules": ["electrolytes", "md", "analysis", "qc", "ai-discovery", "anion-generation"],
        "description": "Administrator - full access to all modules"
    },
    UserRole.PREMIUM: {
        "allowed_modules": ["electrolytes", "md", "analysis", "qc", "ai-discovery", "anion-generation"],
        "description": "Premium user - access to all modules"
    },
    UserRole.USER: {
        "allowed_modules": ["electrolytes", "md", "ai-discovery"],
        "description": "Regular user - access to electrolytes, md, and ai-discovery modules"
    },
    UserRole.GUEST: {
        "allowed_modules": ["electrolytes", "ai-discovery"],
        "description": "Guest user - access to electrolytes and ai-discovery modules"
    }
}

