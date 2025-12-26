"""
Structure analysis utilities for RSNet.

This module provides functions for analyzing molecular structures
and identifying reactive sites and functional groups.
"""

from typing import List, Dict, Set, Tuple, Any
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Descriptors, Crippen
from ..core.molecule import Molecule


def detect_heteroatoms(mol: Molecule) -> List[int]:
    """
    Detect heteroatoms (non-carbon, non-hydrogen atoms) in the molecule.
    
    Args:
        mol: Molecule object
        
    Returns:
        List of atom indices that are heteroatoms
    """
    heteroatoms = []
    rdkit_mol = mol.rdkit_mol
    
    for atom in rdkit_mol.GetAtoms():
        if atom.GetAtomicNum() not in [1, 6]:  # Not H or C
            heteroatoms.append(atom.GetIdx())
    
    return heteroatoms


def detect_polar_bonds(mol: Molecule) -> List[Tuple[int, int]]:
    """
    Detect polar bonds (bonds involving heteroatoms) in the molecule.
    
    Args:
        mol: Molecule object
        
    Returns:
        List of tuples (atom1_idx, atom2_idx) representing polar bonds
    """
    polar_bonds = []
    rdkit_mol = mol.rdkit_mol
    heteroatoms = set(detect_heteroatoms(mol))
    
    for bond in rdkit_mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        
        # Bond is polar if it involves at least one heteroatom
        if begin_idx in heteroatoms or end_idx in heteroatoms:
            polar_bonds.append((begin_idx, end_idx))
    
    return polar_bonds


def detect_acidic_hydrogens(mol: Molecule) -> List[int]:
    """
    Detect potentially acidic hydrogen atoms (attached to heteroatoms).
    
    Args:
        mol: Molecule object
        
    Returns:
        List of hydrogen atom indices that are potentially acidic
    """
    acidic_hydrogens = []
    rdkit_mol = mol.rdkit_mol
    heteroatoms = set(detect_heteroatoms(mol))
    
    for atom in rdkit_mol.GetAtoms():
        if atom.GetAtomicNum() == 1:  # Hydrogen
            # Check if hydrogen is bonded to a heteroatom
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() in heteroatoms:
                    acidic_hydrogens.append(atom.GetIdx())
                    break
    
    return acidic_hydrogens


def detect_hydrogen_bond_acceptors(mol: Molecule) -> List[int]:
    """
    Detect potential hydrogen bond acceptor atoms (O, N with lone pairs).
    
    Args:
        mol: Molecule object
        
    Returns:
        List of atom indices that can accept hydrogen bonds
    """
    acceptors = []
    rdkit_mol = mol.rdkit_mol
    
    for atom in rdkit_mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        
        # Oxygen or nitrogen atoms are potential acceptors
        if atomic_num in [7, 8]:  # N, O
            # Check if atom has lone pairs (simplified check)
            valence = atom.GetTotalValence()
            if (atomic_num == 8 and valence <= 2) or (atomic_num == 7 and valence <= 3):
                acceptors.append(atom.GetIdx())
    
    return acceptors


def detect_weak_bonds(mol: Molecule, energy_threshold: float = 60.0) -> List[Tuple[int, int]]:
    """
    Detect potentially weak bonds that could break easily.
    
    Args:
        mol: Molecule object
        energy_threshold: Energy threshold for considering a bond weak (kcal/mol)
        
    Returns:
        List of tuples (atom1_idx, atom2_idx) representing weak bonds
    """
    weak_bonds = []
    rdkit_mol = mol.rdkit_mol
    
    for bond in rdkit_mol.GetBonds():
        begin_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()
        
        # Heuristic rules for weak bonds
        is_weak = False
        
        # Single bonds involving heteroatoms are often weaker
        if bond.GetBondType() == Chem.BondType.SINGLE:
            if (begin_atom.GetAtomicNum() != 6 or end_atom.GetAtomicNum() != 6):
                is_weak = True
        
        # Bonds in small rings are strained
        if bond.IsInRingSize(3) or bond.IsInRingSize(4):
            is_weak = True
        
        # Long bonds (C-C single bonds in extended chains)
        if (bond.GetBondType() == Chem.BondType.SINGLE and 
            begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6):
            # Check if both atoms have multiple neighbors (branching)
            if len(begin_atom.GetNeighbors()) > 2 and len(end_atom.GetNeighbors()) > 2:
                is_weak = True
        
        if is_weak:
            weak_bonds.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
    
    return weak_bonds


def detect_small_rings(mol: Molecule, max_size: int = 5) -> List[List[int]]:
    """
    Detect small rings in the molecule.
    
    Args:
        mol: Molecule object
        max_size: Maximum ring size to consider as "small"
        
    Returns:
        List of lists, each containing atom indices forming a small ring
    """
    rdkit_mol = mol.rdkit_mol
    ring_info = rdkit_mol.GetRingInfo()
    
    small_rings = []
    for ring in ring_info.AtomRings():
        if len(ring) <= max_size:
            small_rings.append(list(ring))
    
    return small_rings


def detect_pi_systems(mol: Molecule) -> List[List[int]]:
    """
    Detect π-electron systems (double bonds, aromatic rings).
    
    Args:
        mol: Molecule object
        
    Returns:
        List of lists, each containing atom indices participating in π-systems
    """
    rdkit_mol = mol.rdkit_mol
    pi_systems = []
    
    # Find aromatic rings
    ring_info = rdkit_mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        # Check if ring is aromatic
        is_aromatic = True
        for atom_idx in ring:
            atom = rdkit_mol.GetAtomWithIdx(atom_idx)
            if not atom.GetIsAromatic():
                is_aromatic = False
                break
        
        if is_aromatic:
            pi_systems.append(list(ring))
    
    # Find isolated double bonds
    for bond in rdkit_mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE and not bond.GetIsAromatic():
            pi_systems.append([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()])
    
    return pi_systems


def get_structure_tags(mol: Molecule) -> Dict[str, Any]:
    """
    Get comprehensive structure tags for a molecule.
    
    Args:
        mol: Molecule object
        
    Returns:
        Dictionary containing various structural features
    """
    tags = {
        'heteroatoms': detect_heteroatoms(mol),
        'polar_bonds': detect_polar_bonds(mol),
        'acidic_hydrogens': detect_acidic_hydrogens(mol),
        'h_bond_acceptors': detect_hydrogen_bond_acceptors(mol),
        'weak_bonds': detect_weak_bonds(mol),
        'small_rings': detect_small_rings(mol),
        'pi_systems': detect_pi_systems(mol),
        'molecular_weight': mol.molecular_weight,
        'num_rotatable_bonds': rdMolDescriptors.CalcNumRotatableBonds(mol.rdkit_mol),
        'num_hbd': rdMolDescriptors.CalcNumHBD(mol.rdkit_mol),  # H-bond donors
        'num_hba': rdMolDescriptors.CalcNumHBA(mol.rdkit_mol),  # H-bond acceptors
        'logp': Crippen.MolLogP(mol.rdkit_mol),
        'tpsa': rdMolDescriptors.CalcTPSA(mol.rdkit_mol),  # Topological polar surface area
    }
    
    # Add boolean flags for easy checking
    tags['has_heteroatoms'] = len(tags['heteroatoms']) > 0
    tags['has_polar_bonds'] = len(tags['polar_bonds']) > 0
    tags['has_acidic_hydrogens'] = len(tags['acidic_hydrogens']) > 0
    tags['has_small_rings'] = len(tags['small_rings']) > 0
    tags['has_pi_systems'] = len(tags['pi_systems']) > 0
    tags['has_weak_bonds'] = len(tags['weak_bonds']) > 0
    
    return tags


def find_hydrogen_transfer_sites(mol: Molecule) -> List[Dict[str, Any]]:
    """
    Find potential hydrogen transfer sites in a molecule.
    
    Args:
        mol: Molecule object
        
    Returns:
        List of dictionaries describing potential H-transfer sites
    """
    sites = []
    tags = get_structure_tags(mol)
    
    # Find donor-acceptor pairs for hydrogen transfer
    donors = tags['acidic_hydrogens']
    acceptors = tags['h_bond_acceptors']
    
    rdkit_mol = mol.rdkit_mol
    
    for donor_h_idx in donors:
        donor_h = rdkit_mol.GetAtomWithIdx(donor_h_idx)
        # Find the heavy atom bonded to this hydrogen
        donor_heavy = None
        for neighbor in donor_h.GetNeighbors():
            if neighbor.GetAtomicNum() != 1:
                donor_heavy = neighbor
                break
        
        if donor_heavy is None:
            continue
        
        for acceptor_idx in acceptors:
            if acceptor_idx == donor_heavy.GetIdx():
                continue  # Skip self-transfer

            acceptor = rdkit_mol.GetAtomWithIdx(acceptor_idx)
            
            # Calculate distance (if 3D coordinates available)
            distance = None
            coords = mol.get_coordinates()
            if coords is not None:
                donor_pos = coords[donor_heavy.GetIdx()]
                acceptor_pos = coords[acceptor_idx]
                distance = np.linalg.norm(donor_pos - acceptor_pos)
            
            sites.append({
                'type': 'hydrogen_transfer',
                'donor_hydrogen': donor_h_idx,
                'donor_heavy': donor_heavy.GetIdx(),
                'acceptor': acceptor_idx,
                'distance': distance,
                'donor_element': donor_heavy.GetSymbol(),
                'acceptor_element': acceptor.GetSymbol(),
            })
    
    return sites
