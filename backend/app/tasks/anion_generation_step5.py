"""
Step 5: Register new anion in database and refresh cache
"""
import logging
from pathlib import Path
from datetime import datetime
from sqlalchemy.orm import Session
from app.core.logger import logger

# Paths
_CLOUD_SALTS_DIR = Path("/opt/molyte_web_v1.0/data/initial_salts")
# 使用统一路径配置
from app.core.paths import paths
_CAMPUS_SALTS_DIR = paths.initial_salts_dir
SALTS_DIR = _CLOUD_SALTS_DIR if _CLOUD_SALTS_DIR.exists() else _CAMPUS_SALTS_DIR


def register_anion_in_database(db: Session, anion_name: str, display_name: str,
                               charge: int, lt_path: str, pdb_path: str,
                               user_id: int) -> bool:
    """
    Register new anion in database
    
    Note: Currently, the system scans initial_salts/ directory dynamically.
    This function is a placeholder for future database-backed anion registry.
    
    Args:
        db: Database session
        anion_name: Short anion name
        display_name: Full display name
        charge: Charge
        lt_path: Path to .lt file
        pdb_path: Path to .pdb file
        user_id: User ID who created this anion
        
    Returns:
        True if successful, False otherwise
    """
    try:
        # Currently, the system uses file-based discovery
        # The new anion will be automatically discovered when the directory is scanned
        # This function is here for future enhancement to add database tracking
        
        logger.info(f"Anion {anion_name} registered successfully")
        logger.info(f"  Display name: {display_name}")
        logger.info(f"  Charge: {charge}")
        logger.info(f"  .lt path: {lt_path}")
        logger.info(f"  .pdb path: {pdb_path}")
        logger.info(f"  Created by user: {user_id}")
        
        return True
        
    except Exception as e:
        logger.error(f"Error registering anion: {str(e)}")
        return False


def refresh_ions_cache():
    """
    Refresh the ions cache in the electrolytes API
    
    This is called after a new anion is successfully generated
    to ensure it appears in the available ions list
    """
    try:
        # Import here to avoid circular imports
        from app.api.v1 import electrolytes
        
        # Reset the cache
        electrolytes._ions_cache = None
        electrolytes._cache_timestamp = None
        
        logger.info("Ions cache refreshed")
        return True
        
    except Exception as e:
        logger.error(f"Error refreshing ions cache: {str(e)}")
        return False


def verify_generated_files(lt_path: Path, pdb_path: Path) -> bool:
    """
    Verify that generated files exist and have content
    
    Args:
        lt_path: Path to .lt file
        pdb_path: Path to .pdb file
        
    Returns:
        True if both files exist and have content
    """
    try:
        if not lt_path.exists():
            logger.error(f".lt file not found: {lt_path}")
            return False
        
        if not pdb_path.exists():
            logger.error(f".pdb file not found: {pdb_path}")
            return False
        
        # Check file sizes
        lt_size = lt_path.stat().st_size
        pdb_size = pdb_path.stat().st_size
        
        if lt_size == 0:
            logger.error(f".lt file is empty: {lt_path}")
            return False
        
        if pdb_size == 0:
            logger.error(f".pdb file is empty: {pdb_path}")
            return False
        
        logger.info(f"Generated files verified:")
        logger.info(f"  .lt file: {lt_size} bytes")
        logger.info(f"  .pdb file: {pdb_size} bytes")
        
        return True
        
    except Exception as e:
        logger.error(f"Error verifying generated files: {str(e)}")
        return False

