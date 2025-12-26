"""
Force field and anion generation schemas
"""
from pydantic import BaseModel, Field
from typing import Optional, Dict, Any, List
from datetime import datetime


class AnionGenerationRequest(BaseModel):
    """Request to generate anion force field"""
    anion_name: str = Field(..., min_length=1, max_length=50, description="Anion short name (e.g., FSI)")
    display_name: Optional[str] = Field(None, max_length=255, description="Full display name")
    charge: int = Field(default=-1, description="Anion charge")
    identifier_type: str = Field(..., description="'smiles' or 'cas'")
    identifier_value: str = Field(..., description="SMILES string or CAS number")


class AnionGenerationResponse(BaseModel):
    """Response for anion generation submission"""
    job_id: str = Field(..., description="Unique job ID")
    status: str = Field(..., description="Job status")


class AnionGenerationStatusResponse(BaseModel):
    """Response for anion generation status query"""
    job_id: str
    status: str  # "pending", "running", "success", "failed"
    message: Optional[str] = None
    anion_key: Optional[str] = None  # anion_name when successful
    files: Optional[Dict[str, str]] = None  # {"lt_path": "...", "pdb_path": "..."}
    created_at: datetime
    updated_at: datetime
    started_at: Optional[datetime] = None
    finished_at: Optional[datetime] = None
    
    class Config:
        from_attributes = True


class AnionLibraryEntry(BaseModel):
    """Anion library entry"""
    anion_name: str
    display_name: str
    charge: int
    lt_path: str
    pdb_path: str
    source: str  # "manual" or "auto_generated_sob_gaussian"
    description: Optional[str] = None
    created_at: datetime
    
    class Config:
        from_attributes = True


class AnionLibraryListResponse(BaseModel):
    """List of available anions"""
    anions: List[AnionLibraryEntry]
    total: int

