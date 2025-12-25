"""
Helper functions for job CPU hours estimation
"""
from typing import Any, Dict


def estimate_job_cpu_hours(job_data: Any) -> float:
    """
    Estimate CPU hours required for a job
    
    Args:
        job_data: Job creation data (MDJobCreate, QCJobCreate, etc.)
        
    Returns:
        Estimated CPU hours
    """
    # Simple estimation based on job type and parameters
    # This is a placeholder - you should adjust based on actual resource usage
    
    # For MD jobs
    if hasattr(job_data, 'nsteps_npt'):
        # Estimate based on simulation steps
        total_steps = (job_data.nsteps_npt or 500000) + (job_data.nsteps_nvt or 500000)
        cores = (job_data.slurm_ntasks or 8) * (job_data.slurm_cpus_per_task or 8)
        
        # Rough estimate: 1M steps on 64 cores takes ~2 hours
        hours_per_million_steps = 2.0
        estimated_hours = (total_steps / 1000000) * hours_per_million_steps * (64 / cores)
        
        # Add QC overhead if enabled
        if hasattr(job_data, 'qc_options') and job_data.qc_options and job_data.qc_options.enabled:
            estimated_hours += 5.0  # Rough QC overhead
        
        return max(1.0, estimated_hours)
    
    # For QC jobs
    if hasattr(job_data, 'basis_set'):
        # Estimate based on basis set and functional
        # Larger basis sets take more time
        basis_set = job_data.basis_set or "6-31g(d,p)"
        if "++" in basis_set or "aug" in basis_set:
            base_hours = 2.0  # Diffuse functions are expensive
        elif "cc-pv" in basis_set.lower():
            base_hours = 3.0  # Correlation consistent basis sets
        else:
            base_hours = 1.0
        
        return base_hours
    
    # Default fallback
    return 1.0
