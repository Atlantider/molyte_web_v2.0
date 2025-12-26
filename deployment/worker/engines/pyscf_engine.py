import json
from pathlib import Path
from typing import Dict, Any, Optional
import subprocess
from .base_engine import BaseQCEngine
import sys

# Import project utils for path resolution
sys.path.append(str(Path(__file__).parents[3] / 'backend'))
from app.utils.qc_parameter_converter import ParameterConverter

class PySCFEngine(BaseQCEngine):
    """PySCF execution engine"""

    def generate_input(self, job: Dict, work_dir: Path) -> Path:
        """Generate PySCF Python script"""
        config = job.get('config', {})
        molecule_name = config.get('molecule_name', f'QC_{job["id"]}')
        
        # Convert parameters
        pyscf_config = ParameterConverter.prepare_for_engine(config, 'pyscf')
        
        smiles = pyscf_config.get('smiles', '')
        charge = pyscf_config.get('charge', 0)
        spin = pyscf_config.get('spin', 0)
        basis = pyscf_config.get('basis_set', 'sto-3g')
        functional = pyscf_config.get('functional', 'b3lyp')
        solvent_model = pyscf_config.get('solvent_model')
        solvent_name = pyscf_config.get('solvent_name', 'water')
        solvent_eps = pyscf_config.get('solvent_dielectric')
        
        # Geometric optimization?
        job_type = config.get('job_type', 'SP').upper() # SP, OPT, FREQ
        
        xyz_content = config.get('xyz_content')
        if config.get('initial_xyz'):
             xyz_content = config.get('initial_xyz')

        safe_name = self._sanitize_filename(molecule_name)
        py_script_path = work_dir / f"{safe_name}.py"
        
        script_content = self._build_pyscf_script(
            mol_name=safe_name,
            smiles=smiles,
            xyz=xyz_content,
            charge=charge,
            spin=spin,
            basis=basis,
            functional=functional,
            job_type=job_type,
            solvent_model=solvent_model,
            solvent_eps=solvent_eps,
            nprocs=int(config.get('slurm_cpus', 4))
        )
        
        with open(py_script_path, 'w') as f:
            f.write(script_content)
            
        return py_script_path

    def _build_pyscf_script(self, mol_name, smiles, xyz, charge, spin, basis, functional, job_type, solvent_model, solvent_eps, nprocs):
        """Construct the Python script content"""
        
        script = f"""
import numpy as np
import json
from pyscf import gto, scf, dft, solvent
from pyscf.geomopt import geometric_solver

# Redirect stdout to capture results cleanly
import sys

def run_calculation():
    results = {{}}
    try:
        mol = gto.Mole()
        mol.verbose = 4
        # mol.output = '{mol_name}.log'
        mol.atom = '''{xyz}''' if '''{xyz}''' else '{smiles}'
        # If no xyz/smiles provided, this will fail, handle in outer logic? 
        # Assuming one is present.
        
        mol.charge = {charge}
        mol.spin = {spin} 
        mol.basis = '{basis}'
        mol.max_memory = 4000 
        mol.build()

        # Mean Field
        if '{functional}' == 'hf':
            mf = scf.RHF(mol) if mol.spin == 0 else scf.UHF(mol)
        else:
            mf = dft.RKS(mol) if mol.spin == 0 else dft.UKS(mol)
            mf.xc = '{functional}'

        # Solvent
        # PySCF solvent support is evolving. 
        # solvent.DDCOSMO is standard
        if '{solvent_model}' in ['pcm', 'smd', 'custom'] and {solvent_eps} is not None:
             # Using DDCOSMO as approximation/implementation for implicit solvent
             # This is a simplification. Real SMD needs pyscf-properties or plug-ins
             # Just setting eps for DDCOSMO for now as 'implicit' generic
             mf = mf.ddCOSMO()
             mf.with_solvent.eps = {solvent_eps}

        # Run 
        if '{job_type}' == 'OPT':
            mol_eq = geometric_solver.optimize(mf, maxsteps=100)
            results['energy'] = mol_eq.e_tot
            results['xyz'] = mol_eq.atom_coords() # Needs formatting
            # Recalculate energy at optimized geometry to be sure?
            # geometric_solver returns the optimized mole object
            
        elif '{job_type}' == 'FREQ':
            mf.kernel()
            hessian = mf.Hessian().kernel()
            results['energy'] = mf.e_tot
            # Freq analysis logic... complex to dump to JSON for now
            
        else: # SP
            mf.kernel()
            results['energy'] = mf.e_tot

        results['converged'] = mf.converged
        
        # Save results
        with open('results.json', 'w') as f:
            json.dump(results, f)
            
    except Exception as e:
        with open('error.log', 'w') as f:
            f.write(str(e))
        sys.exit(1)

if __name__ == "__main__":
    run_calculation()
"""
        return script

    def generate_job_script(self, job: Dict, input_path: Path, work_dir: Path) -> Path:
        """Generate Slurm script for PySCF"""
        config = job.get('config', {})
        safe_name = self._sanitize_filename(input_path.stem)
        slurm_partition = config.get('slurm_partition', 'cpu')
        slurm_cpus = config.get('slurm_cpus', 4)
        slurm_time = config.get('slurm_time', 7200)

        script_content = f"""#!/bin/bash
#SBATCH --job-name={safe_name}
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --partition={slurm_partition}
#SBATCH --ntasks={slurm_cpus}
#SBATCH --time={slurm_time // 60}:{(slurm_time % 60):02d}

# Activate environment
source activate molyte

# Run PySCF
export OMP_NUM_THREADS={slurm_cpus}
python {input_path.name} > {safe_name}.log
"""
        job_script = work_dir / "job.sh"
        with open(job_script, "w", encoding="utf-8") as f:
            f.write(script_content)
        return job_script

    def submit_job(self, job_script: Path, work_dir: Path) -> Dict[str, Any]:
        """Submit via sbatch"""
        try:
            result = subprocess.run(
                ["sbatch", str(job_script)],
                cwd=str(work_dir),
                capture_output=True,
                text=True,
                timeout=30
            )
            if result.returncode == 0:
                slurm_job_id = result.stdout.strip().split()[-1]
                return {"success": True, "slurm_job_id": slurm_job_id}
            else:
                return {"success": False, "error": result.stderr.strip()}
        except Exception as e:
            return {"success": False, "error": str(e)}

    def parse_result(self, work_dir: Path) -> Dict[str, Any]:
        """Parse results.json"""
        result_file = work_dir / "results.json"
        if not result_file.exists():
             return None
        
        try:
            with open(result_file, 'r') as f:
                data = json.load(f)
            return data
        except Exception:
            return None
