from pathlib import Path
from typing import Dict, List, Tuple
import time
import shutil
import subprocess
from .base_handler import BaseHandler

class MDHandler(BaseHandler):
    """Handler for MD and RESP tasks"""

    def handle_job(self, job: Dict) -> None:
        job_id = job['id']
        self.logger.info(f"Start processing MD Job {job_id}")
        
        try:
            self.update_status(job_id, 'QUEUED', 'md')
            
            # Late imports to avoid circular dependency issues if any
            from app.workers.molyte_wrapper import MolyteWrapper
            from app.workers.resp_calculator import RESPCalculator
            
            wrapper = MolyteWrapper(
                work_base_path=Path(self.config['local']['work_base_path']),
                initial_salts_path=Path(self.config['local']['initial_salts_path']),
                ligpargen_path=Path(self.config['local']['ligpargen_path']),
                packmol_path=Path(self.config['local']['packmol_path']),
                ltemplify_path=Path(self.config['local']['ltemplify_path']),
                moltemplate_path=Path(self.config['local']['moltemplate_path']),
                charge_save_path=Path(self.config['local']['charge_save_path']),
            )
            
            job_data = job['config']
            charge_method = job_data.get("charge_method", "ligpargen")
            
            if charge_method == "resp":
                solvents = job_data.get("solvents", [])
                solvents_needing_resp = wrapper.get_solvents_needing_resp(solvents)
                
                if solvents_needing_resp:
                    self.logger.info(f"MD Job {job_id} requires RESP: {[s['name'] for s in solvents_needing_resp]}")
                    self._start_resp_calculations(job_id, job, solvents_needing_resp, wrapper)
                    return

            self._continue_md_job(job_id, job, wrapper)

        except Exception as e:
            self.logger.error(f"MD Job {job_id} failed: {e}", exc_info=True)
            self.update_status(job_id, 'FAILED', 'md', error_message=str(e))

    def _start_resp_calculations(self, job_id, job, solvents, wrapper):
        from app.workers.resp_calculator import RESPCalculator
        
        try:
            resp_base_dir = Path(self.config['local']['work_base_path']) / f"RESP_{job['config'].get('name', job_id)}"
            resp_base_dir.mkdir(parents=True, exist_ok=True)
            
            md_config = job.get('config', {})
            slurm_partition = md_config.get('slurm_partition', self.config.get('slurm', {}).get('partition', 'cpu'))
            
            resp_calculator = RESPCalculator(
                charge_save_path=Path(self.config['local']['charge_save_path']),
                slurm_partition=slurm_partition
            )
            
            resp_jobs = []
            for solvent in solvents:
                name = solvent['name']
                smiles = solvent.get('smiles', '')
                if not smiles: continue
                
                solvent_dir = resp_base_dir / name
                solvent_dir.mkdir(parents=True, exist_ok=True)
                
                # Run LigParGen
                ligpargen_path = self.config['local']['ligpargen_path']
                cmd = f"{ligpargen_path}/ligpargen -s '{smiles}' -n {name} -r MOL -c 0 -o 0 -cgen CM1A"
                res = subprocess.run(cmd, shell=True, cwd=str(solvent_dir), capture_output=True, text=True)
                if res.returncode != 0:
                    self.logger.error(f"LigParGen failed for {name}: {res.stderr}")
                    continue
                    
                pdb_file = f"{name}.charmm.pdb"
                script_path = resp_calculator.generate_resp_slurm_script(
                    work_dir=solvent_dir,
                    pdb_file=pdb_file,
                    molecule_name=name,
                    charge=0,
                    spin_multiplicity=1,
                    solvent="water",
                    cpus=16,
                    time_limit_hours=24
                )
                
                success, slurm_id, err = resp_calculator.submit_resp_job(solvent_dir, script_path)
                if success:
                    resp_jobs.append({
                        'molecule_name': name,
                        'slurm_job_id': slurm_id,
                        'work_dir': str(solvent_dir),
                        'status': 'RUNNING',
                        'cpu_hours': 0.0
                    })
            
            if resp_jobs:
                self.worker.running_jobs[job_id] = {
                    'type': 'md_waiting_resp',
                    'job': job,
                    'wrapper': wrapper,
                    'resp_jobs': resp_jobs,
                    'resp_base_dir': str(resp_base_dir),
                    'start_time': time.time()
                }
                self.logger.info(f"MD Job {job_id} waiting for {len(resp_jobs)} RESP jobs")
            else:
                self.logger.warning(f"No RESP jobs started for {job_id}, falling back")
                job['config']['charge_method'] = 'ligpargen'
                self._continue_md_job(job_id, job, wrapper)
                
        except Exception as e:
            self.logger.error(f"Failed to start RESP: {e}")
            self.update_status(job_id, 'FAILED', 'md', error_message=f"RESP failed: {e}")

    def _continue_md_job(self, job_id, job, wrapper):
        try:
            job_data = job['config']
            result = wrapper.generate_lammps_input(job_data)
            if not result['success']:
                raise Exception(result.get('error', 'Unknown error'))
            
            work_dir = result['work_dir']
            slurm_res = wrapper.submit_to_slurm(work_dir)
            if not slurm_res['success']:
                 raise Exception(slurm_res.get('error'))
            
            slurm_job_id = slurm_res['slurm_job_id']
            
            resp_cpu = 0.0
            if job_id in self.worker.running_jobs:
                 resp_cpu = self.worker.running_jobs[job_id].get('total_resp_cpu_hours', 0.0)
            
            self.update_status(job_id, 'RUNNING', 'md', slurm_job_id=slurm_job_id, work_dir=str(work_dir))
            
            self.worker.running_jobs[job_id] = {
                'type': 'md',
                'slurm_job_id': slurm_job_id,
                'work_dir': work_dir,
                'start_time': time.time(),
                'resp_cpu_hours': resp_cpu
            }
        except Exception as e:
            self.logger.error(f"MD Continuation failed: {e}")
            self.update_status(job_id, 'FAILED', 'md', error_message=str(e))
    
    def check_resp_jobs(self, job_id, job_info) -> Tuple[bool, bool]:
        from app.workers.resp_calculator import RESPCalculator
        resp_calc = RESPCalculator(charge_save_path=Path(self.config['local']['charge_save_path']))
        
        all_completed = True
        any_failed = False
        total_cpu = 0.0
        
        for rjob in job_info['resp_jobs']:
            if rjob['status'] in ['COMPLETED', 'FAILED']:
                total_cpu += rjob.get('cpu_hours', 0.0)
                if rjob['status'] == 'FAILED': any_failed = True
                continue
            
            status = resp_calc.check_job_status(rjob['slurm_job_id'])
            if status in ['PENDING', 'RUNNING']:
                all_completed = False
            elif status == 'COMPLETED':
                rjob['status'] = 'COMPLETED'
                rjob['cpu_hours'] = resp_calc.get_job_cpu_hours(rjob['slurm_job_id'])
                total_cpu += rjob['cpu_hours']
            else:
                rjob['status'] = 'FAILED'
                any_failed = True
        
        job_info['total_resp_cpu_hours'] = total_cpu
        return all_completed, any_failed
