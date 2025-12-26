"""
Anion force field auto-generation task
"""
import os
import subprocess
import logging
import time
from datetime import datetime
from pathlib import Path
from celery import shared_task
from sqlalchemy.orm import Session
from app.database import SessionLocal
from app.models import AnionGenerationJob, AnionLibrary, AnionGenerationStatus
from app.core.config import settings

logger = logging.getLogger(__name__)

# Try to import RDKit
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    logger.warning("RDKit not available - SMILES parsing will fail")


def process_anion_generation_job(job_id: int):
    """
    Backend preparation for anion generation
    Called by polling_worker via API

    This function only does the minimal preparation work on Tencent Cloud:
    1. Update job status to QC_PENDING
    2. Set up work directory path
    3. Let polling_worker handle all RDKit/computation work

    Args:
        job_id: Database ID of the AnionGenerationJob
    """
    db = SessionLocal()
    try:
        job = db.query(AnionGenerationJob).filter(AnionGenerationJob.id == job_id).first()
        if not job:
            logger.error(f"Job {job_id} not found")
            return False

        # Update status to QC_PENDING directly
        # All actual work (RDKit, QC job creation) will be done by polling_worker
        job.status = AnionGenerationStatus.QC_PENDING
        job.started_at = datetime.utcnow()
        job.message = "Ready for polling_worker to process (parse input, create QC job, run calculations)"

        # Set up work directory path (but don't create        # 使用统一路径配置
        from app.core.paths import paths
        work_base_path = paths.qc_work_dir / f"anion_{job.job_id}"
        job.work_dir = str(work_base_path) # Convert Path object to string for database storage

        db.commit()

        logger.info(f"Anion generation job {job.job_id} (anion: {job.anion_name}) marked as QC_PENDING for polling_worker")
        return True

    except Exception as e:
        logger.error(f"Error in anion generation preparation {job_id}: {e}", exc_info=True)
        job.status = AnionGenerationStatus.FAILED
        job.finished_at = datetime.utcnow()
        job.message = f"Error in preparation: {str(e)}"
        db.commit()
        return False

    finally:
        db.close()


def _parse_input_and_generate_3d(job: AnionGenerationJob) -> dict:
    """
    Step 1: Parse user input and generate 3D structure

    For salt inputs (CAS numbers), automatically extracts the anion portion.
    For example: CAS 678966-16-0 (LiDFOP) -> extracts DFOP anion

    Returns: dict with coordinates data
    """
    if not RDKIT_AVAILABLE:
        raise RuntimeError("RDKit is not available. Please install it: pip install rdkit")

    try:
        if job.identifier_type == "smiles":
            smiles = job.identifier_value
        elif job.identifier_type == "cas":
            # Query PubChem API to get SMILES
            full_smiles = _resolve_cas_to_smiles(job.identifier_value)
            logger.info(f"CAS {job.identifier_value} resolved to SMILES: {full_smiles}")

            # Extract anion from salt (e.g., LiDFOP -> DFOP-)
            smiles = _extract_anion_from_salt(full_smiles, expected_charge=job.charge)
            logger.info(f"Extracted anion SMILES: {smiles}")
        else:
            raise ValueError(f"Unknown identifier type: {job.identifier_type}")

        # 使用渐进式坐标生成
        try:
            from app.utils.coordinate_generator import generate_3d_coordinates
            
            logger.info(f"Generating 3D coordinates for anion: {smiles}")
            
            coord_result = generate_3d_coordinates(
                smiles=smiles,
                molecule_name=job.anion_name,
                charge=job.charge,
                multiplicity=1,  # 大多数阴离子是闭壳层
                enable_xtb=True
            )
            
            if coord_result.source == 'random':
                logger.warning(
                    f"⚠ Anion {job.anion_name} 使用随机坐标生成\n"
                    f"   质量: {coord_result.quality}\n"
                    f"   最小距离: {coord_result.min_distance:.2f} Å\n"
                    f"   建议: 人工检查结构合理性"
                )
            else:
                logger.info(f"✓ Anion坐标生成成功: {coord_result.source} (quality: {coord_result.quality})")
            
            # 从XYZ内容解析坐标
            xyz_lines = coord_result.xyz_content.strip().split('\n')
            coords = []
            mol = Chem.MolFromSmiles(smiles)
            mol = Chem.AddHs(mol)
            
            for i, line in enumerate(xyz_lines[2:]):  # 跳过前两行
                parts = line.split()
                if len(parts) >= 4:
                    coords.append({
                        'idx': i,
                        'symbol': parts[0],
                        'x': float(parts[1]),
                        'y': float(parts[2]),
                        'z': float(parts[3])
                    })

            return {
                'mol': mol,
                'smiles': smiles,
                'coords': coords,
                'coordinate_source': coord_result.source,
                'coordinate_quality': coord_result.quality
            }
            
        except Exception as e:
            raise RuntimeError(f"坐标生成失败: {str(e)}")
            'coords': coords
        }

    except Exception as e:
        raise RuntimeError(f"SMILES_PARSE_FAILED: {str(e)}")


def _resolve_cas_to_smiles(cas_number: str) -> str:
    """
    Resolve CAS number to SMILES using PubChem API
    Requires network access
    """
    try:
        import requests

        # Query PubChem for CAS number
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{cas_number}/JSON"
        response = requests.get(url, timeout=10)

        if response.status_code != 200:
            raise ValueError(f"CAS number not found: {cas_number}")

        data = response.json()
        smiles = data['PC_Compounds'][0]['props'][5]['ival']['sval']
        return smiles

    except Exception as e:
        raise RuntimeError(f"CAS_RESOLVE_FAILED: {str(e)}")


def _extract_anion_from_salt(full_smiles: str, expected_charge: int = -1) -> str:
    """
    Extract anion from salt SMILES (e.g., LiDFOP -> DFOP-)

    Strategies:
    1. If SMILES contains '.', split and find fragment with expected charge
    2. If covalent structure, identify and remove metal ions (Li, Na, K, Mg, Ca)
    3. Verify final charge matches expected

    Args:
        full_smiles: SMILES of salt (e.g., "[Li+].O=P(F)(F)([O-])..." or "O=P(F)(F)([O-])...")
        expected_charge: Expected charge of anion (default -1)

    Returns:
        SMILES of anion only
    """
    from rdkit import Chem

    logger.info(f"Extracting anion from salt SMILES: {full_smiles}")

    # Strategy 1: Split by '.' and find anion fragment
    if '.' in full_smiles:
        logger.info("SMILES contains '.', attempting to split fragments")
        fragments = full_smiles.split('.')

        for frag in fragments:
            try:
                frag_mol = Chem.MolFromSmiles(frag)
                if frag_mol is None:
                    continue

                # Calculate charge
                charge = sum(atom.GetFormalCharge() for atom in frag_mol.GetAtoms())

                logger.info(f"Fragment: {frag}, charge: {charge}")

                # Check if this is the anion
                if charge == expected_charge:
                    logger.info(f"Found anion fragment: {frag}")
                    return frag
            except Exception as e:
                logger.warning(f"Failed to parse fragment {frag}: {e}")
                continue

    # Strategy 2: Remove metal ions from covalent structure
    logger.info("Attempting to remove metal ions from structure")
    mol = Chem.MolFromSmiles(full_smiles)
    if mol is None:
        raise ValueError(f"Failed to parse SMILES: {full_smiles}")

    # List of metal ions to remove
    metal_symbols = ['Li', 'Na', 'K', 'Mg', 'Ca', 'Al', 'Zn']

    # Find metal atoms
    editable_mol = Chem.RWMol(mol)
    atoms_to_remove = []

    for atom in editable_mol.GetAtoms():
        if atom.GetSymbol() in metal_symbols:
            logger.info(f"Found metal ion: {atom.GetSymbol()} at index {atom.GetIdx()}")
            atoms_to_remove.append(atom.GetIdx())

    # Remove metal atoms (in reverse order to maintain indices)
    for idx in sorted(atoms_to_remove, reverse=True):
        editable_mol.RemoveAtom(idx)

    anion_mol = editable_mol.GetMol()
    anion_smiles = Chem.MolToSmiles(anion_mol)

    # Verify charge
    charge = sum(atom.GetFormalCharge() for atom in anion_mol.GetAtoms())
    logger.info(f"Extracted anion SMILES: {anion_smiles}, charge: {charge}")

    if charge != expected_charge:
        logger.warning(
            f"Extracted anion charge ({charge}) != expected ({expected_charge}). "
            f"This may be normal if the salt structure is complex."
        )

    return anion_smiles


def _coords_to_xyz(coords: list) -> str:
    """
    Convert coordinates list to XYZ format string

    Args:
        coords: List of dicts with 'symbol', 'x', 'y', 'z'

    Returns:
        XYZ format string
    """
    num_atoms = len(coords)
    xyz_lines = [str(num_atoms), "Anion structure"]

    for coord in coords:
        symbol = coord['symbol']
        x = coord['x']
        y = coord['y']
        z = coord['z']
        xyz_lines.append(f"{symbol:2s}  {x:12.8f}  {y:12.8f}  {z:12.8f}")

    return '\n'.join(xyz_lines)


def _create_gaussian_qc_job(job: AnionGenerationJob, coords_data: dict) -> int:
    """
    Step 2: Create QC job for Gaussian optimization

    Creates a QC job in the database with SUBMITTED status.
    The polling_worker will pick it up and process it.

    Returns: QC job ID
    """
    from app.models import QCJob, QCJobStatus
    from app.database import SessionLocal

    try:
        db = SessionLocal()

        # Create QC job
        qc_job = QCJob(
            user_id=job.user_id,
            molecule_name=f"anion_{job.anion_name}",
            smiles=coords_data['smiles'],
            charge=job.charge,
            spin_multiplicity=1,  # Most anions are closed-shell
            functional="B3LYP",
            basis_set="6-31++G(d,p)",  # Diffuse functions for anions!
            task_type="opt",
            solvent_model="gas",
            molecule_type="anion",
            config={
                "slurm_partition": "hpc128c",
                "slurm_cpus": 16,
                "slurm_time": 7200,  # 2 hours
                "custom_keywords": "Pop=MK IOp(6/33=2,6/42=6)",  # MK charges
                "initial_xyz": _coords_to_xyz(coords_data['coords'])  # Provide initial geometry
            },
            status=QCJobStatus.SUBMITTED  # Set to SUBMITTED so polling_worker picks it up
        )

        db.add(qc_job)
        db.commit()
        db.refresh(qc_job)

        qc_job_id = qc_job.id
        logger.info(f"Created QC job {qc_job_id} for anion {job.anion_name}")

        db.close()
        return qc_job_id

    except Exception as e:
        logger.error(f"Failed to create QC job for anion {job.anion_name}: {e}")
        raise RuntimeError(f"QC_JOB_CREATION_FAILED: {str(e)}")


def _wait_for_qc_completion(job: AnionGenerationJob, qc_job_id: int) -> str:
    """
    Step 3: Wait for QC job to complete

    Polls the QC job status until it completes or fails.
    Returns the path to the Gaussian log file.
    """
    from app.models import QCJob, QCJobStatus
    from app.database import SessionLocal

    try:
        max_polls = 1800  # 30 minutes
        poll_interval = 10  # 10 seconds
        poll_count = 0

        while poll_count < max_polls:
            time.sleep(poll_interval)
            poll_count += 1

            # Query QC job status
            db = SessionLocal()
            qc_job = db.query(QCJob).filter(QCJob.id == qc_job_id).first()
            db.close()

            if not qc_job:
                raise RuntimeError(f"QC job {qc_job_id} not found")

            if qc_job.status == QCJobStatus.COMPLETED:
                logger.info(f"QC job {qc_job_id} completed successfully")
                # Return path to Gaussian log file
                work_dir = Path(qc_job.work_dir) if qc_job.work_dir else Path(settings.QC_WORK_BASE_PATH) / f"QC-{qc_job_id}"
                gaussian_log = work_dir / f"{qc_job.molecule_name}.log"
                return str(gaussian_log)

            elif qc_job.status == QCJobStatus.FAILED:
                raise RuntimeError(f"QC job {qc_job_id} failed: {qc_job.error_message}")

            elif qc_job.status == QCJobStatus.CANCELLED:
                raise RuntimeError(f"QC job {qc_job_id} was cancelled")

            # Log progress every 60 seconds
            if poll_count % 6 == 0:
                logger.info(f"Waiting for QC job {qc_job_id}... ({poll_count * poll_interval}s elapsed)")

        raise RuntimeError(f"QC job {qc_job_id} timeout (30 minutes)")

    except Exception as e:
        if "QC_JOB" in str(e):
            raise
        raise RuntimeError(f"QC_JOB_WAIT_FAILED: {str(e)}")





def _write_gaussian_input(gjf_path: Path, job: AnionGenerationJob, coords_data: dict):
    """Write Gaussian input file for anion calculation

    Important: Use diffuse functions (6-31+G(d) or 6-31++G(d,p)) for anions
    to properly describe the electron distribution
    """
    coords = coords_data['coords']

    # Determine spin multiplicity
    spin_multiplicity = 1  # Default for closed-shell
    # Could add logic here to detect odd electrons

    # For anions, use basis set with diffuse functions
    # 6-31++G(d,p): double diffuse (on heavy atoms and H) + polarization
    # This is critical for accurate charge distribution in anions
    basis_set = "6-31++G(d,p)"

    gjf_content = f"""%chk=anion_{job.anion_name}.chk
%mem=4GB
%nprocshared=8
#p B3LYP/{basis_set} Opt TightSCF Pop=MK IOp(6/33=2,6/42=6)

{job.display_name} anion auto-generated by Molyte

{job.charge} {spin_multiplicity}
"""

    # Add coordinates
    for coord in coords:
        gjf_content += f"{coord['symbol']:2s}  {coord['x']:12.8f}  {coord['y']:12.8f}  {coord['z']:12.8f}\n"

    gjf_content += "\n"

    with open(gjf_path, 'w') as f:
        f.write(gjf_content)


def _run_multiwfn(job: AnionGenerationJob, gaussian_log: str) -> str:
    """
    Step 3: Run Multiwfn to generate mol2 file from Gaussian output
    Returns: path to mol2 file
    """
    try:
        work_dir = Path(job.work_dir)
        mol2_file = work_dir / f"anion_{job.anion_name}.mol2"

        # Create Multiwfn input script
        # This script extracts charges and generates mol2 file
        mwfn_script = work_dir / f"multiwfn_{job.anion_name}.txt"

        script_content = f"""
{gaussian_log}
7
1
{mol2_file}
q
"""

        with open(mwfn_script, 'w') as f:
            f.write(script_content)

        # Run Multiwfn
        result = subprocess.run(
            ["Multiwfn", str(mwfn_script)],
            cwd=str(work_dir),
            capture_output=True,
            text=True,
            timeout=600
        )

        if result.returncode != 0:
            logger.error(f"Multiwfn failed: {result.stderr}")
            raise RuntimeError(f"MULTIWFN_FAILED: {result.stderr[:500]}")

        if not mol2_file.exists():
            raise RuntimeError("MULTIWFN_FAILED: mol2 file not created")

        logger.info(f"Multiwfn generated mol2 file: {mol2_file}")
        return str(mol2_file)

    except Exception as e:
        if "MULTIWFN_FAILED" in str(e):
            raise
        raise RuntimeError(f"MULTIWFN_FAILED: {str(e)}")


def _run_sobtop(job: AnionGenerationJob, mol2_file: str) -> str:
    """
    Step 4: Run Sobtop to generate GROMACS topology from mol2 file
    Returns: path to GROMACS topology file
    """
    try:
        work_dir = Path(job.work_dir)
        gromacs_top = work_dir / f"anion_{job.anion_name}.top"

        # Prepare sobtop input script
        # Sobtop is interactive, so we need to provide input via stdin
        sobtop_input = f"""{mol2_file}
{gromacs_top}
"""

        # Run Sobtop wrapper script
        # The wrapper handles the sobtop.ini file location
        result = subprocess.run(
            ["/public/software/sobtop_wrapper.sh"],
            input=sobtop_input,
            cwd=str(work_dir),
            capture_output=True,
            text=True,
            timeout=600
        )

        if result.returncode != 0:
            logger.error(f"Sobtop failed: {result.stderr}")
            raise RuntimeError(f"SOBTOP_FAILED: {result.stderr[:500]}")

        if not gromacs_top.exists():
            raise RuntimeError("SOBTOP_FAILED: GROMACS topology file not created")

        logger.info(f"Sobtop generated GROMACS topology: {gromacs_top}")
        return str(gromacs_top)

    except FileNotFoundError:
        raise RuntimeError("SOBTOP_FAILED: Sobtop tool not found. Please install it.")
    except Exception as e:
        if "SOBTOP_FAILED" in str(e):
            raise
        raise RuntimeError(f"SOBTOP_FAILED: {str(e)}")


def _convert_gromacs_to_lammps(job: AnionGenerationJob, gromacs_top: str, coords_data: dict) -> tuple:
    """
    Step 5: Convert GROMACS topology to .lt and .pdb files
    Returns: (lt_path, pdb_path)
    """
    try:
        # Create anion directory
        anion_dir = Path(settings.INITIAL_SALTS_DIR) / job.anion_name
        anion_dir.mkdir(parents=True, exist_ok=True)

        lt_path = anion_dir / f"{job.anion_name}.lt"
        pdb_path = anion_dir / f"{job.anion_name}.pdb"

        # Parse GROMACS topology and convert to LAMMPS format
        _generate_lt_from_gromacs(lt_path, job, gromacs_top)

        # Generate .pdb file from coordinates
        _generate_pdb_file(pdb_path, job, coords_data)

        logger.info(f"Generated force field files: {lt_path}, {pdb_path}")
        return (str(lt_path), str(pdb_path))

    except Exception as e:
        raise RuntimeError(f"OUTPUT_FILE_INVALID: {str(e)}")


def _build_lt_content(job, atoms_data, bonds_data, angles_data, dihedrals_data,
                      atom_types_params, bond_types_params, angle_types_params, dihedral_types_params) -> str:
    """
    Build complete LAMMPS .lt file content
    """
    # Get unique atom types and elements
    unique_atom_types = {}
    unique_elements = set()
    for atom_id, atom_info in atoms_data.items():
        atom_type = atom_info['type']
        element = atom_info['element']
        unique_atom_types[atom_type] = atom_info
        unique_elements.add(element)

    # Build masses section - 使用原子类型而不是元素符号
    masses_section = "  write_once(\"Data Masses\") {\n"
    for atom_type, atom_info in sorted(unique_atom_types.items()):
        element = atom_info['element']
        mass = _get_element_mass(element)
        masses_section += f"    @atom:{atom_type} {mass}\n"
    masses_section += "  }\n"

    # Build pair coefficients section
    pair_coeff_section = "  write_once(\"In Settings\") {\n"
    pair_coeff_section += "    # Pair coefficients (LJ parameters)\n"
    for atom_type, params in atom_types_params.items():
        sigma = params['sigma']
        epsilon = params['epsilon']
        pair_coeff_section += f"    pair_coeff @atom:{atom_type} @atom:{atom_type} {epsilon:.6f} {sigma:.6f}\n"
    pair_coeff_section += "  }\n"

    # Build bond coefficients section
    bond_coeff_section = "  write_once(\"In Settings\") {\n"
    bond_coeff_section += "    # Bond coefficients\n"
    for bond_key, params in bond_types_params.items():
        kr = params['kr']
        r0 = params['r0']
        bond_coeff_section += f"    bond_coeff @bond:{bond_key} {kr:.4f} {r0:.4f}\n"
    bond_coeff_section += "  }\n"

    # Build angle coefficients section
    angle_coeff_section = "  write_once(\"In Settings\") {\n"
    angle_coeff_section += "    # Angle coefficients\n"
    for angle_key, params in angle_types_params.items():
        ktheta = params['ktheta']
        theta0 = params['theta0']
        angle_coeff_section += f"    angle_coeff @angle:{angle_key} {ktheta:.4f} {theta0:.2f}\n"
    angle_coeff_section += "  }\n"

    # Build dihedral coefficients section
    dihedral_coeff_section = ""
    if dihedral_types_params:
        dihedral_coeff_section = "  write_once(\"In Settings\") {\n"
        dihedral_coeff_section += "    # Dihedral coefficients\n"
        for dihedral_key, params in dihedral_types_params.items():
            for i, param_set in enumerate(params['params']):
                phi = param_set['phi']
                kphi = param_set['kphi']
                n = param_set['n']
                dihedral_coeff_section += f"    dihedral_coeff @dihedral:{dihedral_key}_{i} {kphi:.4f} {n} {phi:.2f}\n"
        dihedral_coeff_section += "  }\n"

    # Build atoms section
    atoms_section = "  write(\"Data Atoms\") {\n"
    for atom_id, atom_info in sorted(atoms_data.items(), key=lambda x: int(x[0])):
        atom_type = atom_info['type']
        charge = atom_info['charge']
        # Note: coordinates will be added from PDB file by moltemplate
        atoms_section += f"    $atom:{atom_id} $mol @atom:{atom_type} {charge:.10f}\n"
    atoms_section += "  }\n"

    # Build bonds section
    bonds_section = "  write(\"Data Bonds\") {\n"
    for idx, bond in enumerate(bonds_data, 1):
        bond_key = f"{bond['i']}-{bond['j']}"
        bonds_section += f"    $bond:id{idx} @bond:{bond_key} $atom:{bond['i']} $atom:{bond['j']}\n"
    bonds_section += "  }\n"

    # Build angles section
    angles_section = ""
    if angles_data:
        angles_section = "  write(\"Data Angles\") {\n"
        for idx, angle in enumerate(angles_data, 1):
            angle_key = f"{angle['i']}-{angle['j']}-{angle['k']}"
            angles_section += f"    $angle:id{idx} @angle:{angle_key} $atom:{angle['i']} $atom:{angle['j']} $atom:{angle['k']}\n"
        angles_section += "  }\n"

    # Build dihedrals section
    dihedrals_section = ""
    if dihedrals_data:
        dihedrals_section = "  write(\"Data Dihedrals\") {\n"
        for idx, dihedral in enumerate(dihedrals_data, 1):
            dihedral_key = f"{dihedral['i']}-{dihedral['j']}-{dihedral['k']}-{dihedral['l']}"
            dihedrals_section += f"    $dihedral:id{idx} @dihedral:{dihedral_key}_0 $atom:{dihedral['i']} $atom:{dihedral['j']} $atom:{dihedral['k']} $atom:{dihedral['l']}\n"
        dihedrals_section += "  }\n"

    # Build charges section
    charges_section = "  write_once(\"In Charges\") {\n"
    for atom_id, atom_info in sorted(atoms_data.items(), key=lambda x: int(x[0])):
        atom_type = atom_info['type']
        charge = atom_info['charge']
        charges_section += f"    set type @atom:{atom_type} charge {charge:.10f}\n"
    charges_section += "  }\n"

    # Build list_salt section
    atom_types_list = " ".join(sorted(unique_elements))
    list_salt_section = f"""  write_once("In List_salt") {{
    group        {job.anion_name}           type   {' '.join(f'@atom:{e}' for e in sorted(unique_elements))}
    variable     {job.anion_name}_list      index  "{atom_types_list}"
  }}
"""

    # Assemble complete file
    lt_content = f"""{job.anion_name} {{

### LAMMPS force field for {job.display_name}
### Auto-generated from GROMACS topology by Molyte
### Source: Gaussian {job.identifier_value}

{masses_section}

{pair_coeff_section}

{bond_coeff_section}

{angle_coeff_section}

{dihedral_coeff_section}

{atoms_section}

{bonds_section}

{angles_section}

{dihedrals_section}

{charges_section}

{list_salt_section}

}} # end of "{job.anion_name}" type definition
"""

    return lt_content


def _get_element_mass(element: str) -> float:
    """Get atomic mass for element"""
    masses = {
        'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999,
        'F': 18.998, 'S': 32.060, 'P': 30.974, 'Cl': 35.453,
        'Br': 79.904, 'I': 126.904, 'B': 10.811, 'Si': 28.086
    }
    return masses.get(element, 12.0)


def _generate_lt_from_gromacs(lt_path: Path, job: AnionGenerationJob, gromacs_top: str):
    """
    Generate LAMMPS .lt file from GROMACS topology
    Parses GROMACS .top file and converts to LAMMPS format
    """
    try:
        # Parse GROMACS topology file
        with open(gromacs_top, 'r') as f:
            gromacs_content = f.read()

        # Extract all data from GROMACS format
        atoms_data = _parse_gromacs_atoms_detailed(gromacs_content)
        bonds_data = _parse_gromacs_bonds_detailed(gromacs_content)
        angles_data = _parse_gromacs_angles_detailed(gromacs_content)
        dihedrals_data = _parse_gromacs_dihedrals_detailed(gromacs_content)
        atom_types_params = _parse_gromacs_atom_types(gromacs_content)
        bond_types_params = _parse_gromacs_bond_types(gromacs_content)
        angle_types_params = _parse_gromacs_angle_types(gromacs_content)
        dihedral_types_params = _parse_gromacs_dihedral_types(gromacs_content)

        # Generate LAMMPS .lt file
        lt_content = _build_lt_content(
            job, atoms_data, bonds_data, angles_data, dihedrals_data,
            atom_types_params, bond_types_params, angle_types_params, dihedral_types_params
        )

        with open(lt_path, 'w') as f:
            f.write(lt_content)

        logger.info(f"Generated LAMMPS .lt file: {lt_path}")

    except Exception as e:
        logger.error(f"Error generating .lt file: {e}", exc_info=True)
        raise RuntimeError(f"LT_GENERATION_FAILED: {str(e)}")


def _parse_gromacs_atoms_detailed(gromacs_content: str) -> dict:
    """
    Parse detailed atom information from GROMACS topology
    Returns: dict with atom_id -> {type, charge, mass, element}
    """
    atoms = {}
    in_atoms = False

    for line in gromacs_content.split('\n'):
        if '[ atoms ]' in line:
            in_atoms = True
            continue
        if in_atoms and line.startswith('['):
            break
        if in_atoms and line.strip() and not line.startswith(';'):
            parts = line.split()
            if len(parts) >= 7:
                atom_id = parts[0]
                atom_type = parts[1]
                res_num = parts[2]
                res_name = parts[3]
                atom_name = parts[4]
                charge = float(parts[6])
                mass = float(parts[7]) if len(parts) > 7 else 0.0

                # Infer element from atom type or name
                element = _infer_element(atom_type, atom_name)

                atoms[atom_id] = {
                    'type': atom_type,
                    'name': atom_name,
                    'charge': charge,
                    'mass': mass,
                    'element': element,
                    'res_name': res_name
                }

    return atoms


def _infer_element(atom_type: str, atom_name: str) -> str:
    """Infer element symbol from atom type or name"""
    # Try to extract element from atom type (e.g., 'C1' -> 'C', 'O2' -> 'O')
    for char in atom_type:
        if char.isalpha():
            element = char
            # Check if next char is also alpha (for two-letter elements)
            idx = atom_type.index(char)
            if idx + 1 < len(atom_type) and atom_type[idx + 1].isalpha() and atom_type[idx + 1].islower():
                element = atom_type[idx:idx+2]
            return element

    # Fallback to atom name
    for char in atom_name:
        if char.isalpha():
            return char

    return 'X'  # Unknown element


def _parse_gromacs_bonds_detailed(gromacs_content: str) -> list:
    """
    Parse detailed bond information from GROMACS topology
    Returns: list of {i, j, type, func_type}
    """
    bonds = []
    in_bonds = False

    for line in gromacs_content.split('\n'):
        if '[ bonds ]' in line:
            in_bonds = True
            continue
        if in_bonds and line.startswith('['):
            break
        if in_bonds and line.strip() and not line.startswith(';'):
            parts = line.split()
            if len(parts) >= 3:
                bonds.append({
                    'i': parts[0],
                    'j': parts[1],
                    'func_type': parts[2] if len(parts) > 2 else '1',
                    'type': f"bond_{parts[0]}_{parts[1]}"
                })

    return bonds


def _parse_gromacs_angles_detailed(gromacs_content: str) -> list:
    """
    Parse detailed angle information from GROMACS topology
    Returns: list of {i, j, k, func_type, type}
    """
    angles = []
    in_angles = False

    for line in gromacs_content.split('\n'):
        if '[ angles ]' in line:
            in_angles = True
            continue
        if in_angles and line.startswith('['):
            break
        if in_angles and line.strip() and not line.startswith(';'):
            parts = line.split()
            if len(parts) >= 4:
                angles.append({
                    'i': parts[0],
                    'j': parts[1],
                    'k': parts[2],
                    'func_type': parts[3] if len(parts) > 3 else '1',
                    'type': f"angle_{parts[0]}_{parts[1]}_{parts[2]}"
                })

    return angles


def _parse_gromacs_dihedrals_detailed(gromacs_content: str) -> list:
    """
    Parse detailed dihedral information from GROMACS topology
    Returns: list of {i, j, k, l, func_type, type}
    """
    dihedrals = []
    in_dihedrals = False

    for line in gromacs_content.split('\n'):
        if '[ dihedrals ]' in line:
            in_dihedrals = True
            continue
        if in_dihedrals and line.startswith('['):
            break
        if in_dihedrals and line.strip() and not line.startswith(';'):
            parts = line.split()
            if len(parts) >= 5:
                dihedrals.append({
                    'i': parts[0],
                    'j': parts[1],
                    'k': parts[2],
                    'l': parts[3],
                    'func_type': parts[4] if len(parts) > 4 else '1',
                    'type': f"dihedral_{parts[0]}_{parts[1]}_{parts[2]}_{parts[3]}"
                })

    return dihedrals


def _parse_gromacs_atom_types(gromacs_content: str) -> dict:
    """Parse atom type parameters from GROMACS topology"""
    atom_types = {}
    in_atomtypes = False

    for line in gromacs_content.split('\n'):
        if '[ atomtypes ]' in line:
            in_atomtypes = True
            continue
        if in_atomtypes and line.startswith('['):
            break
        if in_atomtypes and line.strip() and not line.startswith(';'):
            parts = line.split()
            if len(parts) >= 6:
                atom_type = parts[0]
                mass = float(parts[2])
                charge = float(parts[3])
                sigma = float(parts[5])  # in nm
                epsilon = float(parts[6]) if len(parts) > 6 else 0.0  # in kJ/mol

                # Convert to LAMMPS units (Angstrom and kcal/mol)
                sigma_angstrom = sigma * 10.0
                epsilon_kcal = epsilon / 4.184

                atom_types[atom_type] = {
                    'mass': mass,
                    'charge': charge,
                    'sigma': sigma_angstrom,
                    'epsilon': epsilon_kcal
                }

    return atom_types


def _parse_gromacs_bond_types(gromacs_content: str) -> dict:
    """Parse bond type parameters from GROMACS topology"""
    bond_types = {}
    in_bondtypes = False

    for line in gromacs_content.split('\n'):
        if '[ bondtypes ]' in line:
            in_bondtypes = True
            continue
        if in_bondtypes and line.startswith('['):
            break
        if in_bondtypes and line.strip() and not line.startswith(';'):
            parts = line.split()
            if len(parts) >= 5:
                atom1 = parts[0]
                atom2 = parts[1]
                func_type = parts[2]
                r0 = float(parts[3])  # in nm
                kr = float(parts[4])  # in kJ/(mol*nm^2)

                # Convert to LAMMPS units
                r0_angstrom = r0 * 10.0
                kr_kcal = kr / 4.184 * 100.0  # Convert to kcal/(mol*Angstrom^2)

                key = f"{atom1}-{atom2}"
                bond_types[key] = {
                    'r0': r0_angstrom,
                    'kr': kr_kcal,
                    'func_type': func_type
                }

    return bond_types


def _parse_gromacs_angle_types(gromacs_content: str) -> dict:
    """Parse angle type parameters from GROMACS topology"""
    angle_types = {}
    in_angletypes = False

    for line in gromacs_content.split('\n'):
        if '[ angletypes ]' in line:
            in_angletypes = True
            continue
        if in_angletypes and line.startswith('['):
            break
        if in_angletypes and line.strip() and not line.startswith(';'):
            parts = line.split()
            if len(parts) >= 5:
                atom1 = parts[0]
                atom2 = parts[1]
                atom3 = parts[2]
                func_type = parts[3]
                theta0 = float(parts[4])  # in degrees
                ktheta = float(parts[5]) if len(parts) > 5 else 0.0  # in kJ/(mol*rad^2)

                # Convert to LAMMPS units
                ktheta_kcal = ktheta / 4.184

                key = f"{atom1}-{atom2}-{atom3}"
                angle_types[key] = {
                    'theta0': theta0,
                    'ktheta': ktheta_kcal,
                    'func_type': func_type
                }

    return angle_types


def _parse_gromacs_dihedral_types(gromacs_content: str) -> dict:
    """Parse dihedral type parameters from GROMACS topology"""
    dihedral_types = {}
    in_dihedraltypes = False

    for line in gromacs_content.split('\n'):
        if '[ dihedraltypes ]' in line:
            in_dihedraltypes = True
            continue
        if in_dihedraltypes and line.startswith('['):
            break
        if in_dihedraltypes and line.strip() and not line.startswith(';'):
            parts = line.split()
            if len(parts) >= 6:
                atom1 = parts[0]
                atom2 = parts[1]
                atom3 = parts[2]
                atom4 = parts[3]
                func_type = parts[4]

                # Parse dihedral parameters (can be multiple sets)
                params = []
                for i in range(5, len(parts), 3):
                    if i + 2 < len(parts):
                        phi = float(parts[i])
                        kphi = float(parts[i+1])
                        n = int(parts[i+2])
                        params.append({'phi': phi, 'kphi': kphi / 4.184, 'n': n})

                key = f"{atom1}-{atom2}-{atom3}-{atom4}"
                dihedral_types[key] = {
                    'params': params,
                    'func_type': func_type
                }

    return dihedral_types


def _generate_pdb_file(pdb_path: Path, job: AnionGenerationJob, coords_data: dict):
    """Generate PDB file

    PDB 格式规范（参考 FSI.pdb）：
    - 列 1-6: 记录名称（ATOM）
    - 列 7-11: 原子序号
    - 列 13-16: 原子名称
    - 列 18-20: 残基名称
    - 列 23-26: 残基序号（整数，右对齐）
    - 列 31-38: X 坐标
    - 列 39-46: Y 坐标
    - 列 47-54: Z 坐标

    注意：不包含链标识符，以确保与 Packmol 兼容
    """
    coords = coords_data['coords']

    pdb_content = "REMARK   Auto-generated by Molyte\n"
    pdb_content += f"REMARK   Anion: {job.display_name}\n"

    for i, coord in enumerate(coords, 1):
        element = coord['symbol']
        x = coord['x']
        y = coord['y']
        z = coord['z']
        # 参考 FSI.pdb 的格式：ATOM      1  S1  MOL     1      -1.492  -0.011   0.142  1.00  0.00           S
        pdb_content += f"ATOM  {i:5d}  {element:2s}  MOL     {1:1d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}  \n"

    pdb_content += "TER\n"
    pdb_content += "END\n"

    with open(pdb_path, 'w') as f:
        f.write(pdb_content)


def _register_anion_in_library(job: AnionGenerationJob, db: Session):
    """
    Step 5: Register anion in library database
    """
    try:
        # Check if already exists
        existing = db.query(AnionLibrary).filter(
            AnionLibrary.anion_name == job.anion_name
        ).first()
        
        if existing:
            # Update existing entry
            existing.lt_path = job.lt_path
            existing.pdb_path = job.pdb_path
            existing.source = "auto_generated_sobtop_gaussian"
            existing.generation_job_id = job.id
        else:
            # Create new entry
            anion = AnionLibrary(
                anion_name=job.anion_name,
                display_name=job.display_name,
                charge=job.charge,
                lt_path=job.lt_path,
                pdb_path=job.pdb_path,
                source="auto_generated_sobtop_gaussian",
                generation_job_id=job.id,
                created_by=job.user_id
            )
            db.add(anion)
        
        db.commit()
        
    except Exception as e:
        raise RuntimeError(f"DB_SAVE_FAILED: {str(e)}")


def _coords_to_xyz(coords: list) -> str:
    """
    Convert coordinates to XYZ format string

    Args:
        coords: List of [atom_symbol, x, y, z] tuples

    Returns:
        XYZ format string
    """
    lines = [str(len(coords)), "Generated by RDKit"]
    for atom_symbol, x, y, z in coords:
        lines.append(f"{atom_symbol:2s}  {x:12.6f}  {y:12.6f}  {z:12.6f}")
    return "\n".join(lines)

