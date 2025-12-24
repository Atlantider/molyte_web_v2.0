"""
Module-level access control utilities
"""
from fastapi import HTTPException, status
from app.models.user import User, UserRole


def check_module_access(user: User, module_name: str) -> bool:
    """
    Check if a user has access to a specific module

    Args:
        user: The user object
        module_name: The module name to check (electrolytes, md, analysis, qc)

    Returns:
        bool: True if user has access, False otherwise
    """
    # Admin users always have access to all modules
    if user.role == UserRole.ADMIN:
        return True

    # If user has no allowed_modules, deny access (strict mode)
    # This ensures module permissions are enforced
    if not user.allowed_modules or len(user.allowed_modules) == 0:
        return False

    # Check if module is in allowed_modules
    return module_name in user.allowed_modules


def require_module_access(user: User, module_name: str) -> None:
    """
    Require that a user has access to a specific module
    Raises HTTPException if user doesn't have access
    
    Args:
        user: The user object
        module_name: The module name to check
    
    Raises:
        HTTPException: If user doesn't have access to the module
    """
    if not check_module_access(user, module_name):
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail=f"You don't have access to the {module_name} module"
        )


# Module name constants
MODULE_ELECTROLYTES = "electrolytes"
MODULE_MD = "md"
MODULE_ANALYSIS = "analysis"
MODULE_QC = "qc"
MODULE_AI_DISCOVERY = "ai-discovery"
MODULE_ANION_GENERATION = "anion-generation"

# All available modules
ALL_MODULES = [MODULE_ELECTROLYTES, MODULE_MD, MODULE_ANALYSIS, MODULE_QC, MODULE_AI_DISCOVERY, MODULE_ANION_GENERATION]

