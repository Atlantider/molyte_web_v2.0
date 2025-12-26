"""
Force field and anion generation models
"""
from sqlalchemy import Column, Integer, String, Text, DateTime, ForeignKey, Enum, Index, Boolean, Float
from sqlalchemy.orm import relationship
from sqlalchemy.sql import func
from sqlalchemy.dialects.postgresql import JSONB
from app.database import Base
import enum
import uuid


class AnionGenerationStatus(str, enum.Enum):
    """Anion generation job status"""
    PENDING = "pending"
    RUNNING = "running"
    QC_PENDING = "qc_pending"  # QC job created, waiting for polling_worker
    QC_RUNNING = "qc_running"  # QC job is running on campus cluster
    QC_COMPLETED = "qc_completed"  # QC job completed, processing results
    SUCCESS = "success"
    FAILED = "failed"


class AnionGenerationJob(Base):
    """Anion force field auto-generation job model"""
    __tablename__ = "anion_generation_jobs"

    id = Column(Integer, primary_key=True, index=True)
    job_id = Column(String(36), unique=True, nullable=False, index=True, default=lambda: str(uuid.uuid4()))
    user_id = Column(Integer, ForeignKey("users.id", ondelete="CASCADE"), nullable=False, index=True)
    
    # Anion information
    anion_name = Column(String(50), nullable=False, index=True)  # e.g., "FSI", "Cl"
    display_name = Column(String(255), nullable=True)  # e.g., "bis(fluorosulfonyl)imide"
    charge = Column(Integer, default=-1)
    
    # Input information
    identifier_type = Column(String(20), nullable=False)  # "smiles" or "cas"
    identifier_value = Column(String(500), nullable=False)  # SMILES string or CAS number
    
    # Job status
    status = Column(Enum(AnionGenerationStatus), default=AnionGenerationStatus.PENDING, nullable=False, index=True)
    message = Column(Text, nullable=True)  # Status message or error description
    
    # Generated files
    lt_path = Column(String(500), nullable=True)  # Path to .lt file
    pdb_path = Column(String(500), nullable=True)  # Path to .pdb file
    
    # Intermediate files
    work_dir = Column(String(500), nullable=True)  # Working directory
    gaussian_log = Column(String(500), nullable=True)  # Gaussian output log
    mol2_file = Column(String(500), nullable=True)  # Multiwfn mol2 output
    gromacs_top = Column(String(500), nullable=True)  # GROMACS topology file
    sob_output = Column(String(500), nullable=True)  # SOB tool output

    # Related QC job
    qc_job_id = Column(Integer, ForeignKey("qc_jobs.id", ondelete="SET NULL"), nullable=True)

    # Metadata
    config = Column(JSONB, default={})  # Additional configuration
    
    # Timestamps
    created_at = Column(DateTime(timezone=True), server_default=func.now(), nullable=False, index=True)
    updated_at = Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now(), nullable=False)
    started_at = Column(DateTime(timezone=True), nullable=True)
    finished_at = Column(DateTime(timezone=True), nullable=True)

    # CPU hours tracking
    cpu_hours = Column(Float, default=0.0, nullable=False)
    estimated_cpu_hours = Column(Float, default=0.0, nullable=False)

    # Soft delete
    is_deleted = Column(Boolean, default=False, index=True)
    deleted_at = Column(DateTime(timezone=True), nullable=True)
    deleted_by = Column(Integer, ForeignKey("users.id", ondelete="SET NULL"), nullable=True)
    delete_reason = Column(String(500), nullable=True)
    
    # Relationships
    user = relationship("User", back_populates="anion_generation_jobs", foreign_keys=[user_id])
    deleted_by_user = relationship("User", foreign_keys=[deleted_by])
    qc_job = relationship("QCJob", foreign_keys=[qc_job_id])
    
    # Indexes
    __table_args__ = (
        Index('idx_anion_gen_user_id', 'user_id'),
        Index('idx_anion_gen_status', 'status'),
        Index('idx_anion_gen_anion_name', 'anion_name'),
        Index('idx_anion_gen_created_at', 'created_at'),
    )
    
    def __repr__(self):
        return f"<AnionGenerationJob(job_id={self.job_id}, anion_name={self.anion_name}, status={self.status})>"


class AnionLibrary(Base):
    """Registered anion library entry"""
    __tablename__ = "anion_library"

    id = Column(Integer, primary_key=True, index=True)
    anion_name = Column(String(50), unique=True, nullable=False, index=True)  # e.g., "FSI"
    display_name = Column(String(255), nullable=False)  # e.g., "bis(fluorosulfonyl)imide"
    charge = Column(Integer, default=-1)
    
    # File paths
    lt_path = Column(String(500), nullable=False)  # Path to .lt file
    pdb_path = Column(String(500), nullable=False)  # Path to .pdb file
    
    # Source information
    source = Column(String(100), default="manual")  # "manual" or "auto_generated_sob_gaussian"
    generation_job_id = Column(Integer, ForeignKey("anion_generation_jobs.id", ondelete="SET NULL"), nullable=True)
    
    # Metadata
    description = Column(Text, nullable=True)
    config = Column(JSONB, default={})
    
    # Timestamps
    created_at = Column(DateTime(timezone=True), server_default=func.now(), nullable=False)
    created_by = Column(Integer, ForeignKey("users.id", ondelete="SET NULL"), nullable=True)
    updated_at = Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now(), nullable=False)
    
    # Soft delete
    is_deleted = Column(Boolean, default=False, index=True)
    deleted_at = Column(DateTime(timezone=True), nullable=True)
    
    # Relationships
    generation_job = relationship("AnionGenerationJob", foreign_keys=[generation_job_id])
    created_by_user = relationship("User", foreign_keys=[created_by])
    
    # Indexes
    __table_args__ = (
        Index('idx_anion_lib_name', 'anion_name'),
        Index('idx_anion_lib_created_at', 'created_at'),
    )
    
    def __repr__(self):
        return f"<AnionLibrary(anion_name={self.anion_name}, source={self.source})>"

