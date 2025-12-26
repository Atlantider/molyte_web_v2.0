from pathlib import Path
from typing import Dict
import time
from .base_handler import BaseHandler
from ..engines.gaussian import GaussianEngine
from ..engines.pyscf_engine import PySCFEngine

class QCHandler(BaseHandler):
    """Handler for QC tasks"""
    
    def __init__(self, worker):
        super().__init__(worker)
        self.engines = {
            'gaussian': GaussianEngine(self.config, self.logger),
            'pyscf': PySCFEngine(self.config, self.logger)
        }
    
    def handle_job(self, job: Dict) -> None:
        job_id = job['id']
        self.logger.info(f"Start handling QC Job {job_id}")
        
        try:
            # Check dup
            if job_id in self.worker.running_jobs:
                self.logger.info(f"QC Job {job_id} already running, skipping")
                return

            self.update_status(job_id, 'QUEUED', 'qc')
            
            # Select engine
            config = job.get('config', {})
            engine_name = config.get('engine', 'gaussian').lower() # Default to Gaussian
            
            if engine_name not in self.engines:
                self.logger.warning(f"Unknown engine {engine_name}, falling back to Gaussian")
                engine_name = 'gaussian'
                
            engine = self.engines[engine_name]
            self.logger.info(f"Using engine: {engine_name}")

            # Prepare work dir
            qc_work_base = Path(self.config['local']['qc_work_base_path'])
            molecule_name = config.get('molecule_name', f'QC_{job_id}')
            dir_name = f"QC-{job_id}-{molecule_name}"
            work_dir = qc_work_base / dir_name
            work_dir.mkdir(parents=True, exist_ok=True)
            
            # Generate input
            input_path = engine.generate_input(job, work_dir)
            self.logger.info(f"Generated input: {input_path}")
            
            # Generate script
            job_script = engine.generate_job_script(job, input_path, work_dir)
            
            # Submit
            submit_result = engine.submit_job(job_script, work_dir)
            
            if not submit_result['success']:
                raise Exception(f"Submission failed: {submit_result.get('error')}")
                
            slurm_job_id = submit_result['slurm_job_id']
            
            # Update status
            self.update_status(
                job_id, 'RUNNING', 'qc',
                slurm_job_id=slurm_job_id,
                work_dir=str(work_dir)
            )
            
            # Track
            self.worker.running_jobs[job_id] = {
                'type': 'qc',
                'slurm_job_id': slurm_job_id,
                'work_dir': work_dir,
                'start_time': time.time()
            }
            
            self.logger.info(f"QC Job {job_id} submitted (Slurm ID: {slurm_job_id})")
            
        except Exception as e:
            self.logger.error(f"QC Job {job_id} failed: {e}", exc_info=True)
            self.update_status(job_id, 'FAILED', 'qc', error_message=str(e))
