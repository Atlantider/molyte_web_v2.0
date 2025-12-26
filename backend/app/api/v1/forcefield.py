"""
Force field and anion generation API endpoints
"""
from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy.orm import Session
from app.database import get_db
from app.dependencies import get_current_user
from app.models import User, AnionGenerationJob, AnionLibrary, AnionGenerationStatus
from app.utils.permissions import require_module_access, MODULE_ANION_GENERATION
from app.schemas.forcefield import (
    AnionGenerationRequest,
    AnionGenerationResponse,
    AnionGenerationStatusResponse,
    AnionLibraryListResponse,
    AnionLibraryEntry
)
from sqlalchemy import desc
from datetime import datetime
import logging

router = APIRouter(tags=["forcefield"])
logger = logging.getLogger(__name__)


@router.post("/anions/auto-generate", response_model=AnionGenerationResponse, status_code=status.HTTP_201_CREATED)
async def auto_generate_anion(
    request: AnionGenerationRequest,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Submit a job to auto-generate anion force field

    The job will be picked up by polling_worker and processed asynchronously.

    Parameters:
    - anion_name: Short name for the anion (e.g., "FSI", "Cl")
    - display_name: Full display name (optional)
    - charge: Anion charge (default: -1)
    - identifier_type: "smiles" or "cas"
    - identifier_value: SMILES string or CAS number
    """
    # Check module access
    require_module_access(current_user, MODULE_ANION_GENERATION)

    try:
        # Check if anion already exists
        existing = db.query(AnionLibrary).filter(
            AnionLibrary.anion_name == request.anion_name,
            AnionLibrary.is_deleted == False
        ).first()

        if existing:
            raise HTTPException(
                status_code=status.HTTP_409_CONFLICT,
                detail=f"Anion '{request.anion_name}' already exists in library"
            )

        # Create generation job
        # Status is PENDING, polling_worker will pick it up
        job = AnionGenerationJob(
            user_id=current_user.id,
            anion_name=request.anion_name,
            display_name=request.display_name or request.anion_name,
            charge=request.charge,
            identifier_type=request.identifier_type,
            identifier_value=request.identifier_value,
            status=AnionGenerationStatus.PENDING,
            message="Job submitted, waiting for polling_worker to process"
        )

        db.add(job)
        db.commit()
        db.refresh(job)

        logger.info(f"Anion generation job {job.job_id} created by user {current_user.id}, "
                   f"waiting for polling_worker to pick up")

        return AnionGenerationResponse(
            job_id=job.job_id,
            status=job.status.value
        )

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error submitting anion generation job: {e}", exc_info=True)
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to submit anion generation job"
        )


@router.get("/anions/auto-generate/{job_id}", response_model=AnionGenerationStatusResponse)
async def get_anion_generation_status(
    job_id: str,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Query the status of an anion generation job
    """
    # Check module access
    require_module_access(current_user, MODULE_ANION_GENERATION)

    try:
        job = db.query(AnionGenerationJob).filter(
            AnionGenerationJob.job_id == job_id,
            AnionGenerationJob.user_id == current_user.id
        ).first()
        
        if not job:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail="Job not found"
            )
        
        response = AnionGenerationStatusResponse(
            job_id=job.job_id,
            status=job.status.value,
            message=job.message,
            anion_key=job.anion_name if job.status == AnionGenerationStatus.SUCCESS else None,
            files={
                "lt_path": job.lt_path,
                "pdb_path": job.pdb_path
            } if job.status == AnionGenerationStatus.SUCCESS else None,
            created_at=job.created_at,
            updated_at=job.updated_at,
            started_at=job.started_at,
            finished_at=job.finished_at
        )
        
        return response
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error querying anion generation status: {e}", exc_info=True)
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to query job status"
        )


@router.get("/anions/jobs")
async def get_anion_generation_jobs(
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Get list of anion generation jobs for current user
    """
    # Check module access
    require_module_access(current_user, MODULE_ANION_GENERATION)

    try:
        jobs = db.query(AnionGenerationJob).filter(
            AnionGenerationJob.user_id == current_user.id,
            AnionGenerationJob.is_deleted == False
        ).order_by(desc(AnionGenerationJob.created_at)).all()

        jobs_data = []
        for job in jobs:
            jobs_data.append({
                "id": job.id,
                "job_id": job.job_id,
                "anion_name": job.anion_name,
                "display_name": job.display_name,
                "charge": job.charge,
                "identifier_type": job.identifier_type,
                "identifier_value": job.identifier_value,
                "status": job.status.value,
                "message": job.message,
                "lt_path": job.lt_path,
                "pdb_path": job.pdb_path,
                "cpu_hours_used": job.cpu_hours or 0,  # 添加核时字段
                "estimated_cpu_hours": job.estimated_cpu_hours or 0,
                "created_at": job.created_at.isoformat() if job.created_at else None,
                "started_at": job.started_at.isoformat() if job.started_at else None,
                "finished_at": job.finished_at.isoformat() if job.finished_at else None,
            })

        return {"jobs": jobs_data}

    except Exception as e:
        logger.error(f"Error fetching anion generation jobs: {e}", exc_info=True)
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to fetch anion generation jobs"
        )


@router.delete("/anions/jobs/{job_id}")
async def delete_anion_generation_job(
    job_id: int,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Delete an anion generation job (soft delete)
    """
    # Check module access
    require_module_access(current_user, MODULE_ANION_GENERATION)

    try:
        job = db.query(AnionGenerationJob).filter(
            AnionGenerationJob.id == job_id,
            AnionGenerationJob.user_id == current_user.id
        ).first()

        if not job:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail="Job not found"
            )

        # Soft delete
        job.is_deleted = True
        job.deleted_at = datetime.utcnow()
        job.deleted_by = current_user.id
        db.commit()

        logger.info(f"Anion generation job {job_id} deleted by user {current_user.id}")

        return {"message": "Job deleted successfully"}

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error deleting anion generation job: {e}", exc_info=True)
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to delete job"
        )


@router.get("/anions/library", response_model=AnionLibraryListResponse)
async def get_anion_library(
    db: Session = Depends(get_db)
):
    """
    Get list of available anions in the library
    """
    try:
        anions = db.query(AnionLibrary).filter(
            AnionLibrary.is_deleted == False
        ).order_by(AnionLibrary.anion_name).all()

        entries = [AnionLibraryEntry.from_orm(anion) for anion in anions]

        return AnionLibraryListResponse(
            anions=entries,
            total=len(entries)
        )

    except Exception as e:
        logger.error(f"Error fetching anion library: {e}", exc_info=True)
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to fetch anion library"
        )


@router.delete("/anions/{anion_name}")
async def delete_anion_from_library(
    anion_name: str,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Delete anion from library (soft delete)

    This endpoint is called by polling_worker when syncing filesystem with database.
    Only worker users can call this endpoint.
    """
    try:
        # Find the anion
        anion = db.query(AnionLibrary).filter(
            AnionLibrary.anion_name == anion_name,
            AnionLibrary.is_deleted == False
        ).first()

        if not anion:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"Anion '{anion_name}' not found in library"
            )

        # Soft delete
        anion.is_deleted = True
        anion.deleted_at = datetime.utcnow()

        db.commit()

        logger.info(f"Anion {anion_name} soft deleted by user {current_user.username}")

        return {
            "status": "ok",
            "message": f"Anion '{anion_name}' deleted successfully"
        }

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error deleting anion {anion_name}: {e}", exc_info=True)
        db.rollback()
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to delete anion"
        )


@router.post("/anions/manual-register")
async def manual_register_anion(
    anion_data: dict,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Manually register an anion in the library

    This endpoint is called by polling_worker when syncing filesystem with database.
    It registers anions that exist in initial_salts folder but not in database.
    """
    try:
        anion_name = anion_data.get('anion_name')
        display_name = anion_data.get('display_name', anion_name)
        charge = anion_data.get('charge', -1)
        lt_path = anion_data.get('lt_path')
        pdb_path = anion_data.get('pdb_path')
        source = anion_data.get('source', 'manual')

        if not anion_name or not lt_path or not pdb_path:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail="Missing required fields: anion_name, lt_path, pdb_path"
            )

        # Check if already exists
        existing = db.query(AnionLibrary).filter(
            AnionLibrary.anion_name == anion_name
        ).first()

        if existing:
            if existing.is_deleted:
                # Restore soft-deleted entry
                existing.is_deleted = False
                existing.deleted_at = None
                existing.lt_path = lt_path
                existing.pdb_path = pdb_path
                existing.updated_at = datetime.utcnow()
                db.commit()

                logger.info(f"Anion {anion_name} restored from soft delete")

                return {
                    "status": "ok",
                    "message": f"Anion '{anion_name}' restored successfully"
                }
            else:
                raise HTTPException(
                    status_code=status.HTTP_409_CONFLICT,
                    detail=f"Anion '{anion_name}' already exists in library"
                )

        # Create new entry
        anion = AnionLibrary(
            anion_name=anion_name,
            display_name=display_name,
            charge=charge,
            lt_path=lt_path,
            pdb_path=pdb_path,
            source=source,
            created_by=current_user.id
        )

        db.add(anion)
        db.commit()
        db.refresh(anion)

        logger.info(f"Anion {anion_name} manually registered by user {current_user.username}")

        return {
            "status": "ok",
            "anion_name": anion_name,
            "message": f"Anion '{anion_name}' registered successfully"
        }

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error registering anion: {e}", exc_info=True)
        db.rollback()
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to register anion"
        )

