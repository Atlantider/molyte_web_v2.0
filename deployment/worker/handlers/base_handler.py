from abc import ABC, abstractmethod
from typing import Dict, Any

class BaseHandler(ABC):
    """Abstract base class for job handlers"""
    
    def __init__(self, worker):
        """
        Initialize handler.
        
        Args:
            worker: The parent worker instance (provides logger, config, API client)
        """
        self.worker = worker
        self.logger = worker.logger
        self.config = worker.config
        
    @abstractmethod
    def handle_job(self, job: Dict) -> None:
        """
        Process a job.
        
        Args:
            job: Job dictionary
        """
        pass
    
    def update_status(self, job_id, status, job_type, **kwargs):
        """Helper to update job status via worker"""
        self.worker._update_job_status(job_id, status, job_type, **kwargs)
