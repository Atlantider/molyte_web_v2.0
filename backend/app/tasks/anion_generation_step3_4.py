"""
Step 3-4: sob tool invocation and .lt/.pdb file generation
"""
import logging
import subprocess
import re
from pathlib import Path
from typing import Dict, Any, Optional, List, Tuple
from app.core.logger import logger

# Paths
_CLOUD_SALTS_DIR = Path("/opt/molyte_web_v1.0/data/initial_salts")
_CAMPUS_SALTS_DIR = Path("/public/home/xiaoji/molyte_web/data/initial_salts")
SALTS_DIR = _CLOUD_SALTS_DIR if _CLOUD_SALTS_DIR.exists() else _CAMPUS_SALTS_DIR


def _run_sob(log_path: Path, work_dir: Path, anion_name: str) -> Optional[Dict[str, Any]]:
    """
    Run sob tool to generate force field parameters from Gaussian output
    
    Args:
        log_path: Path to Gaussian .log file
        work_dir: Working directory
        anion_name: Anion name
        
    Returns:
        Dict with sob output info, or None if failed
    """
    try:
        output_prefix = f"anion_{anion_name}"
        
        logger.info(f"Running sob for {anion_name}...")
        result = subprocess.run(
            ["sob", "-i", str(log_path), "-o", output_prefix],
            cwd=str(work_dir),
            capture_output=True,
            text=True,
            timeout=600  # 10 minute timeout
        )
        
        if result.returncode != 0:
            logger.error(f"sob failed with return code {result.returncode}")
            logger.error(f"stderr: {result.stderr}")
            return None
        
        logger.info(f"sob completed successfully")
        logger.info(f"sob output:\n{result.stdout}")
        
        # Look for generated files
        # sob typically generates files like anion_<name>.itp, anion_<name>.gro, etc.
        generated_files = list(work_dir.glob(f"{output_prefix}.*"))
        logger.info(f"Generated files: {[f.name for f in generated_files]}")
        
        return {
            "output_prefix": output_prefix,
            "generated_files": [str(f) for f in generated_files],
            "stdout": result.stdout,
            "stderr": result.stderr,
        }
        
    except subprocess.TimeoutExpired:
        logger.error(f"sob calculation timed out for {anion_name}")
        return None
    except FileNotFoundError:
        logger.error("sob command not found - sob tool may not be installed")
        return None
    except Exception as e:
        logger.error(f"Error running sob: {str(e)}")
        return None


def _parse_sob_output(sob_output: Dict[str, Any], work_dir: Path, 
                     anion_name: str) -> Optional[Dict[str, Any]]:
    """
    Parse sob output to extract force field parameters
    
    Args:
        sob_output: Output from _run_sob
        work_dir: Working directory
        anion_name: Anion name
        
    Returns:
        Dict with parsed parameters, or None if failed
    """
    try:
        # Look for .itp file (GROMACS format) which contains force field parameters
        itp_files = list(work_dir.glob(f"anion_{anion_name}.itp"))
        
        if not itp_files:
            logger.warning(f"No .itp file found for {anion_name}")
            # Try to parse from stdout
            return {
                "atoms": [],
                "bonds": [],
                "angles": [],
                "dihedrals": [],
                "pairs": [],
            }
        
        itp_path = itp_files[0]
        logger.info(f"Parsing {itp_path}")
        
        with open(itp_path, 'r') as f:
            content = f.read()
        
        # Parse atoms section
        atoms = []
        atom_section = re.search(r'\[\s*atoms\s*\](.*?)(?:\[|$)', content, re.DOTALL)
        if atom_section:
            for line in atom_section.group(1).strip().split('\n'):
                line = line.strip()
                if line and not line.startswith(';'):
                    parts = line.split()
                    if len(parts) >= 5:
                        atoms.append({
                            "nr": parts[0],
                            "type": parts[1],
                            "resnr": parts[2],
                            "resname": parts[3],
                            "atom": parts[4],
                            "charge": float(parts[6]) if len(parts) > 6 else 0.0,
                        })
        
        # Parse bonds section
        bonds = []
        bond_section = re.search(r'\[\s*bonds\s*\](.*?)(?:\[|$)', content, re.DOTALL)
        if bond_section:
            for line in bond_section.group(1).strip().split('\n'):
                line = line.strip()
                if line and not line.startswith(';'):
                    parts = line.split()
                    if len(parts) >= 3:
                        bonds.append({
                            "ai": parts[0],
                            "aj": parts[1],
                            "funct": parts[2],
                        })
        
        logger.info(f"Parsed {len(atoms)} atoms and {len(bonds)} bonds")
        
        return {
            "atoms": atoms,
            "bonds": bonds,
            "itp_path": str(itp_path),
        }
        
    except Exception as e:
        logger.error(f"Error parsing sob output: {str(e)}")
        return None


def _get_template_anion() -> Optional[Path]:
    """
    Get a template anion .lt file to use as reference
    
    Returns:
        Path to template .lt file
    """
    try:
        # Look for any .lt file in initial_salts
        lt_files = list(SALTS_DIR.glob("*/*.lt"))
        if lt_files:
            return lt_files[0]
        
        # Fallback: look in direct subdirectories
        lt_files = list(SALTS_DIR.glob("*.lt"))
        if lt_files:
            return lt_files[0]
        
        logger.warning("No template .lt file found")
        return None
        
    except Exception as e:
        logger.error(f"Error finding template anion: {str(e)}")
        return None


def _generate_lt_file(anion_name: str, display_name: str, charge: int,
                     parsed_params: Dict[str, Any], work_dir: Path,
                     output_dir: Path) -> Optional[Path]:
    """
    Generate .lt file following project conventions
    
    Args:
        anion_name: Short anion name (e.g., FSI)
        display_name: Full display name
        charge: Charge
        parsed_params: Parsed parameters from sob
        work_dir: Working directory
        output_dir: Output directory for .lt file
        
    Returns:
        Path to generated .lt file
    """
    try:
        # Create output directory
        output_dir.mkdir(parents=True, exist_ok=True)
        lt_path = output_dir / f"{anion_name}.lt"
        
        # Get template for reference
        template_path = _get_template_anion()
        template_content = ""
        if template_path:
            with open(template_path, 'r') as f:
                template_content = f.read()
        
        # Generate .lt file
        lines = []
        lines.append(f"# {display_name} ({anion_name})")
        lines.append(f"# Auto-generated by Molyte")
        lines.append("")
        lines.append(f"{anion_name} {{")
        lines.append("")
        lines.append('  write_once("Data Masses") {')
        
        # Add atom masses (simplified - would need actual atomic masses)
        atom_symbols = set()
        for atom in parsed_params.get("atoms", []):
            atom_symbols.add(atom.get("type", "C"))
        
        # Common atomic masses
        masses = {
            "H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999,
            "F": 18.998, "S": 32.060, "Cl": 35.453, "P": 30.974,
        }
        
        for i, symbol in enumerate(sorted(atom_symbols), 1):
            mass = masses.get(symbol, 12.0)
            lines.append(f'    @atom:{symbol}{i}  {mass}')
        
        lines.append("  }")
        lines.append("")
        lines.append('  write_once("In Settings") {')
        lines.append("    # Force field parameters would go here")
        lines.append("    # (parsed from sob output)")
        lines.append("  }")
        lines.append("")
        lines.append('  write("Data Atoms") {')
        
        # Add atoms (simplified)
        for i, atom in enumerate(parsed_params.get("atoms", []), 1):
            atom_type = atom.get("type", "C")
            charge_val = atom.get("charge", 0.0)
            lines.append(f'    ${"{"}atom:{i}{"}"}  ${"{"}mol{"}"}  @atom:{atom_type}  {charge_val}  0  0  0')
        
        lines.append("  }")
        lines.append("")
        lines.append('  write("Data Bonds") {')
        lines.append("    # Bonds would go here")
        lines.append("  }")
        lines.append("")
        lines.append('  write_once("In Charges") {')
        lines.append("    # Charge assignments")
        lines.append("  }")
        lines.append("")
        lines.append('  write_once("In List_salt") {')
        lines.append(f'    group  {anion_name}  type  @atom:*')
        lines.append(f'    variable  {anion_name}_list  index  "*"')
        lines.append("  }")
        lines.append("")
        lines.append("}")
        
        # Write to file
        with open(lt_path, 'w') as f:
            f.write('\n'.join(lines))
        
        logger.info(f"Generated .lt file: {lt_path}")
        return lt_path
        
    except Exception as e:
        logger.error(f"Error generating .lt file: {str(e)}")
        return None


def _generate_pdb_file(anion_name: str, struct_data: Dict[str, Any],
                      output_dir: Path) -> Optional[Path]:
    """
    Generate .pdb file from structure data

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

    Args:
        anion_name: Short anion name
        struct_data: Structure data from RDKit
        output_dir: Output directory

    Returns:
        Path to generated .pdb file
    """
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
        pdb_path = output_dir / f"{anion_name}.pdb"

        lines = []
        lines.append("REMARK   Auto-generated by Molyte")

        # Add atoms
        atoms = struct_data.get("atoms", [])
        coords = struct_data.get("coords", [])

        for i, (atom, coord) in enumerate(zip(atoms, coords), 1):
            symbol = atom.get("symbol", "C")
            x, y, z = coord
            # PDB format: ATOM record (参考 FSI.pdb)
            # ATOM      1  S1  MOL     1      -1.492  -0.011   0.142  1.00  0.00           S
            line = f"ATOM  {i:5d}  {symbol:2s}  {anion_name:3s}     {1:1d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {symbol:>2s}  "
            lines.append(line)

        # Add connectivity
        lines.append("TER")
        lines.append("END")

        # Write to file
        with open(pdb_path, 'w') as f:
            f.write('\n'.join(lines))

        logger.info(f"Generated .pdb file: {pdb_path}")
        return pdb_path
        
    except Exception as e:
        logger.error(f"Error generating .pdb file: {str(e)}")
        return None

