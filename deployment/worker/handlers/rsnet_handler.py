from pathlib import Path
from typing import Dict, List
import json
import time
import requests
from .base_handler import BaseHandler

class RSNetHandler(BaseHandler):
    """Handler for Reaction Network (RS-Net) tasks"""

    def handle_job(self, job: Dict) -> None:
        job_id = job['id']
        self.logger.info(f"Processing Reaction Network Job {job_id}")
        
        try:
            if job_id in self.worker.running_jobs:
                self.logger.info(f"Job {job_id} already running, skipping")
                return

            self.update_status(job_id, 'QUEUED', 'reaction_network')
            
            config = job.get('config', {})
            job_name = config.get('job_name', f'RN_{job_id}')
            initial_smiles = config.get('initial_smiles', [])
            
            if not initial_smiles:
                raise ValueError("Initial SMILES list cannot be empty")
            
            # config params
            temperature = config.get('temperature', 300.0)
            electrode_type = config.get('electrode_type', 'anode')
            voltage = config.get('voltage', 0.1)
            max_generations = config.get('max_generations', 3)
            max_species = config.get('max_species', 50)
            energy_cutoff = config.get('energy_cutoff', 80.0)
            
            # Slurm
            slurm_partition = config.get('slurm_partition', 'cpu')
            slurm_cpus = config.get('slurm_cpus', 16)
            slurm_time = config.get('slurm_time', 7200)
            
            # Work dir
            rn_work_base = Path(self.config['local'].get('rn_work_base_path', 
                                                          self.config['local']['work_base_path']))
            work_dir = rn_work_base / f"RN-{job_id}-{job_name}"
            work_dir.mkdir(parents=True, exist_ok=True)
            
            # Generate Python Driver Script
            script_path = work_dir / "run_rsnet.py"
            self._generate_rsnet_script(
                script_path,
                initial_smiles=initial_smiles,
                temperature=temperature,
                electrode_type=electrode_type,
                voltage=voltage,
                max_generations=max_generations,
                max_species=max_species,
                energy_cutoff=energy_cutoff,
                output_dir=str(work_dir)
            )
            
            # Generate Slurm Script
            job_script = work_dir / "job.sh"
            self._generate_rsnet_slurm_script(
                job_script,
                work_dir=work_dir,
                cpus=slurm_cpus,
                time_limit=slurm_time,
                partition=slurm_partition
            )
            
            # Submit
            slurm_res = self._submit_to_slurm(work_dir)
            if not slurm_res['success']:
                raise Exception(f"Slurm submission failed: {slurm_res.get('error')}")
                
            slurm_job_id = slurm_res['slurm_job_id']
            
            self.update_status(
                job_id, 'RUNNING', 'reaction_network',
                slurm_job_id=slurm_job_id,
                work_dir=str(work_dir)
            )
            
            self.worker.running_jobs[job_id] = {
                'type': 'reaction_network',
                'slurm_job_id': slurm_job_id,
                'work_dir': str(work_dir),
                'start_time': time.time()
            }
            
            self.logger.info(f"RN Job {job_id} submitted (Slurm: {slurm_job_id})")

        except Exception as e:
            self.logger.error(f"RN Job {job_id} failed: {e}", exc_info=True)
            self.update_status(job_id, 'FAILED', 'reaction_network', error_message=str(e))

    def _generate_rsnet_script(self, script_path: Path, **kwargs):
        """Generate RSNet driver script"""
        smiles_str = ',\n        '.join([f'"{s}"' for s in kwargs['initial_smiles']])
        content = f'''#!/usr/bin/env python3
import sys
import json
import os
from pathlib import Path

# Add backend to Python path
backend_path = Path(__file__).parent.parent.parent.parent / "backend"
sys.path.insert(0, str(backend_path))

try:
    # Import RSNet service from backend
    from app.services.rsnet import RSNetService
    
    # Initialize RSNet service and generate network
    service = RSNetService()
    result = service.generate_network(
        smiles_list=[
            {smiles_str}
        ],
        temperature={kwargs['temperature']},
        electrode_type="{kwargs['electrode_type']}",
        voltage={kwargs['voltage']},
        max_generations={kwargs['max_generations']},
        max_species={kwargs['max_species']},
        energy_cutoff={kwargs['energy_cutoff']},
        visualize=True,
        save_results=True,
        output_dir="{kwargs['output_dir']}"
    )
    
    output_file = Path("{kwargs['output_dir']}") / "network_result.json"
    with open(output_file, "w", encoding='utf-8') as f:
        # Simplify serialization
        json.dump(result, f, default=lambda x: str(x), indent=2)
        
except Exception as e:
    print(f"RSNet Execution Error: {{e}}")
    sys.exit(1)
'''
        with open(script_path, 'w', encoding='utf-8') as f:
            f.write(content)
        # os.chmod(script_path, 0o755) # Windows/WSL might not matter but good practice


    def _generate_rsnet_slurm_script(self, script_path: Path, work_dir: Path, cpus, time_limit, partition):
        content = f'''#!/bin/bash
#SBATCH --job-name=RSNet
#SBATCH --output=rsnet_out.log
#SBATCH --error=rsnet_err.log
#SBATCH --partition={partition}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={cpus}
#SBATCH --time={time_limit}
#SBATCH --mem=32G

cd "{work_dir}"

# Activate conda
source /opt/miniconda3/etc/profile.d/conda.sh || source /public/software/anaconda3/etc/profile.d/conda.sh
conda activate molyte

# PYTHONPATH will be set by the Python script itself
export OMP_NUM_THREADS={cpus}

python run_rsnet.py
'''
        with open(script_path, 'w', encoding='utf-8') as f:
            f.write(content)

    def _submit_to_slurm(self, work_dir: Path) -> Dict:
        # Re-implement or use BaseHandler/Worker helper if available
        # But BaseHandler doesn't have it, Engine does. 
        # Actually worker.runner has _check_slurm_status but not submit.
        # We need to use subprocess directly or move submit logic to BaseHandler.
        # For now, duplicate simple submit logic.
        import subprocess
        import re
        try:
            res = subprocess.run(['sbatch', 'job.sh'], cwd=str(work_dir), capture_output=True, text=True)
            if res.returncode != 0:
                return {'success': False, 'error': res.stderr}
            match = re.search(r'Submitted batch job (\d+)', res.stdout)
            if match:
                 return {'success': True, 'slurm_job_id': match.group(1)}
            return {'success': False, 'error': 'No job ID parsed'}
        except Exception as e:
            return {'success': False, 'error': str(e)}
