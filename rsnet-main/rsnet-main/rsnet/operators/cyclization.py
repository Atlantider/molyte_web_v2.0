"""
Cyclization operator for RSNet.

This module implements cyclization reactions where linear or branched molecules
form cyclic structures through intramolecular bond formation.
"""

import logging
from typing import List, Dict, Set, Optional, Tuple, Any
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, rdchem

from ..core.molecule import Molecule
from ..core.environment import Environment
from ..core.reaction import Reaction, ReactionTemplate
from .base import BaseOperator
# from ..utils.structure_analysis import StructureAnalyzer  # Not needed for basic functionality


logger = logging.getLogger(__name__)


class CyclizationTemplate(ReactionTemplate):
    """
    Template for cyclization reactions.
    
    Represents the formation of a ring by connecting two atoms
    in a linear or branched molecule.
    """
    
    def __init__(
        self,
        atom1_idx: int,
        atom2_idx: int,
        bond_type: str = 'single',
        ring_size: int = 5,
        **kwargs
    ):
        """
        Initialize cyclization template.
        
        Args:
            atom1_idx: Index of first atom to connect
            atom2_idx: Index of second atom to connect
            bond_type: Type of bond to form ('single', 'double', 'triple')
            ring_size: Size of ring being formed
        """
        super().__init__(**kwargs)
        self.atom1_idx = atom1_idx
        self.atom2_idx = atom2_idx
        self.bond_type = bond_type
        self.ring_size = ring_size
        
        # Bond type mapping
        self.bond_type_map = {
            'single': Chem.BondType.SINGLE,
            'double': Chem.BondType.DOUBLE,
            'triple': Chem.BondType.TRIPLE
        }
    
    def apply(self, molecules: List[Molecule]) -> List[Reaction]:
        """
        Apply cyclization template to molecules.
        
        Args:
            molecules: List of molecules to apply template to
            
        Returns:
            List of cyclization reactions
        """
        if len(molecules) != 1:
            return []
        
        molecule = molecules[0]
        
        try:
            # Create cyclized product
            cyclized_mol = self._form_cycle(molecule)
            
            if cyclized_mol is None:
                return []
            
            # Create reaction
            reaction = Reaction(
                reactants=[molecule],
                products=[cyclized_mol],
                name=f"Cyclization: {molecule.name} -> cyclic product",
                template=self
            )
            
            return [reaction]
            
        except Exception as e:
            logger.warning(f"Cyclization template failed: {e}")
            return []
    
    def _form_cycle(self, molecule: Molecule) -> Optional[Molecule]:
        """
        Form a cycle by connecting two atoms.
        
        Args:
            molecule: Molecule to cyclize
            
        Returns:
            Cyclized molecule or None if failed
        """
        try:
            # Get RDKit molecule with explicit hydrogens
            rdkit_mol = Chem.AddHs(molecule.rdkit_mol)
            
            # Check if atoms exist and are valid for cyclization
            if (self.atom1_idx >= rdkit_mol.GetNumAtoms() or 
                self.atom2_idx >= rdkit_mol.GetNumAtoms()):
                return None
            
            atom1 = rdkit_mol.GetAtomWithIdx(self.atom1_idx)
            atom2 = rdkit_mol.GetAtomWithIdx(self.atom2_idx)
            
            # Check if atoms are already bonded
            if rdkit_mol.GetBondBetweenAtoms(self.atom1_idx, self.atom2_idx):
                return None
            
            # Check valence constraints
            if not self._can_form_bond(atom1, atom2):
                return None
            
            # Create editable molecule
            editable = Chem.EditableMol(rdkit_mol)
            
            # Add the cyclization bond
            bond_type = self.bond_type_map.get(self.bond_type, Chem.BondType.SINGLE)
            editable.AddBond(self.atom1_idx, self.atom2_idx, bond_type)
            
            # Get the cyclized molecule
            cyclized_rdkit = editable.GetMol()
            
            # Sanitize and remove explicit hydrogens
            try:
                Chem.SanitizeMol(cyclized_rdkit)
                cyclized_rdkit = Chem.RemoveHs(cyclized_rdkit)
            except:
                return None
            
            # Create new Molecule object
            cyclized_smiles = Chem.MolToSmiles(cyclized_rdkit)
            cyclized_mol = Molecule(
                smiles=cyclized_smiles,
                name=f"{molecule.name}_cyclized",
                rdkit_mol=cyclized_rdkit
            )
            
            return cyclized_mol
            
        except Exception as e:
            logger.debug(f"Cyclization failed: {e}")
            return None
    
    def _can_form_bond(self, atom1: rdchem.Atom, atom2: rdchem.Atom) -> bool:
        """
        Check if two atoms can form a bond based on valence.
        
        Args:
            atom1: First atom
            atom2: Second atom
            
        Returns:
            True if bond can be formed
        """
        try:
            # Get current valence and maximum valence
            val1 = atom1.GetTotalValence()
            val2 = atom2.GetTotalValence()
            
            max_val1 = Chem.GetPeriodicTable().GetValenceList(atom1.GetAtomicNum())[0]
            max_val2 = Chem.GetPeriodicTable().GetValenceList(atom2.GetAtomicNum())[0]
            
            # Check if atoms can accommodate additional bond
            bond_order = 1 if self.bond_type == 'single' else (2 if self.bond_type == 'double' else 3)
            
            return (val1 + bond_order <= max_val1 and val2 + bond_order <= max_val2)
            
        except:
            return False


class CyclizationOperator(BaseOperator):
    """
    Operator for cyclization reactions.
    
    Identifies molecules that can undergo intramolecular cyclization
    and generates the corresponding cyclic products.
    """
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """
        Initialize cyclization operator.
        
        Args:
            config: Configuration parameters
        """
        super().__init__(
            name="Cyclization",
            description="Generates cyclization reactions",
            config=config
        )

        # Configuration parameters
        self.min_chain_length = self.config.get('min_chain_length', 4)
        self.max_chain_length = self.config.get('max_chain_length', 12)
        self.preferred_ring_sizes = self.config.get('preferred_ring_sizes', [5, 6, 7])
        self.allow_strained_rings = self.config.get('allow_strained_rings', False)
        self.max_templates_per_molecule = self.config.get('max_templates_per_molecule', 5)
    
    def is_applicable(self, molecules: List[Molecule]) -> bool:
        """
        Check if cyclization operator is applicable.
        
        Args:
            molecules: List of molecules to check
            
        Returns:
            True if operator can be applied
        """
        if len(molecules) != 1:
            return False
        
        molecule = molecules[0]
        
        # Check if molecule is already cyclic
        if self._has_rings(molecule):
            return False
        
        # Check chain length
        chain_length = self._get_longest_chain_length(molecule)
        if chain_length < self.min_chain_length or chain_length > self.max_chain_length:
            return False
        
        # Check if cyclization sites exist
        cyclization_sites = self._find_cyclization_sites(molecule)
        
        return len(cyclization_sites) > 0

    def can_apply(self, molecules: List[Molecule], environment: Environment) -> bool:
        """
        Check if cyclization operator can be applied.

        Args:
            molecules: List of molecules to check
            environment: Reaction environment

        Returns:
            True if operator can be applied
        """
        if not self.is_applicable(molecules):
            return False

        # Get driving forces from environment
        drives = environment.get_active_drives()

        # Cyclization reactions are favored by:
        # - Thermal conditions
        # - Ring strain relief
        # - π-system formation

        required_drives = [
            drives.get('thermal', False),
            drives.get('ring_strain', False),
            drives.get('pi_system_reaction', False)
        ]

        return any(required_drives)

    def apply(self, molecules: List[Molecule], environment: Environment) -> List[Reaction]:
        """
        Apply cyclization operator to generate reactions.
        
        Args:
            molecules: List of molecules to apply operator to
            environment: Reaction environment
            
        Returns:
            List of cyclization reactions
        """
        if not self.is_applicable(molecules):
            return []
        
        molecule = molecules[0]
        reactions = []
        
        # Find cyclization sites
        cyclization_sites = self._find_cyclization_sites(molecule)
        
        # Limit number of templates
        cyclization_sites = cyclization_sites[:self.max_templates_per_molecule]
        
        # Generate reactions for each site
        for site in cyclization_sites:
            template = CyclizationTemplate(
                atom1_idx=site['atom1_idx'],
                atom2_idx=site['atom2_idx'],
                bond_type=site['bond_type'],
                ring_size=site['ring_size']
            )
            
            template_reactions = template.apply([molecule])
            reactions.extend(template_reactions)
        
        logger.debug(f"Cyclization operator generated {len(reactions)} reactions for {molecule.name}")
        
        return reactions
    
    def _has_rings(self, molecule: Molecule) -> bool:
        """Check if molecule already has rings."""
        return rdMolDescriptors.CalcNumRings(molecule.rdkit_mol) > 0
    
    def _get_longest_chain_length(self, molecule: Molecule) -> int:
        """Get the length of the longest carbon chain."""
        try:
            # Simple approximation: count heavy atoms
            return molecule.rdkit_mol.GetNumHeavyAtoms()
        except:
            return 0
    
    def _find_cyclization_sites(self, molecule: Molecule) -> List[Dict[str, Any]]:
        """
        Find potential cyclization sites in a molecule.
        
        Args:
            molecule: Molecule to analyze
            
        Returns:
            List of cyclization site dictionaries
        """
        sites = []
        rdkit_mol = molecule.rdkit_mol
        
        try:
            # Get all heavy atoms
            heavy_atoms = [atom for atom in rdkit_mol.GetAtoms() if atom.GetAtomicNum() > 1]
            
            # Find pairs of atoms that could form rings
            for i, atom1 in enumerate(heavy_atoms):
                for j, atom2 in enumerate(heavy_atoms[i+1:], i+1):
                    
                    # Skip if atoms are already bonded
                    if rdkit_mol.GetBondBetweenAtoms(atom1.GetIdx(), atom2.GetIdx()):
                        continue
                    
                    # Calculate potential ring size
                    ring_size = self._estimate_ring_size(rdkit_mol, atom1.GetIdx(), atom2.GetIdx())
                    
                    if ring_size is None:
                        continue
                    
                    # Check if ring size is preferred
                    if ring_size not in self.preferred_ring_sizes:
                        if not self.allow_strained_rings or ring_size < 3 or ring_size > 8:
                            continue
                    
                    # Check if atoms can form bond
                    if not self._can_cyclize(atom1, atom2):
                        continue
                    
                    # Determine bond type (default to single)
                    bond_type = self._determine_bond_type(atom1, atom2)
                    
                    sites.append({
                        'atom1_idx': atom1.GetIdx(),
                        'atom2_idx': atom2.GetIdx(),
                        'ring_size': ring_size,
                        'bond_type': bond_type,
                        'strain_energy': self._estimate_strain_energy(ring_size)
                    })
            
            # Sort by strain energy (prefer low strain)
            sites.sort(key=lambda x: x['strain_energy'])
            
        except Exception as e:
            logger.warning(f"Error finding cyclization sites: {e}")
        
        return sites
    
    def _estimate_ring_size(self, mol: Chem.Mol, atom1_idx: int, atom2_idx: int) -> Optional[int]:
        """
        Estimate the ring size that would be formed by connecting two atoms.
        
        Args:
            mol: RDKit molecule
            atom1_idx: Index of first atom
            atom2_idx: Index of second atom
            
        Returns:
            Estimated ring size or None if no reasonable path exists
        """
        try:
            # Find shortest path between atoms
            path = Chem.GetShortestPath(mol, atom1_idx, atom2_idx)
            
            if path is None or len(path) < 3:
                return None
            
            # Ring size would be the path length
            ring_size = len(path)
            
            # Only consider reasonable ring sizes
            if 3 <= ring_size <= 12:
                return ring_size
            
            return None
            
        except:
            return None
    
    def _can_cyclize(self, atom1: rdchem.Atom, atom2: rdchem.Atom) -> bool:
        """
        Check if two atoms can participate in cyclization.
        
        Args:
            atom1: First atom
            atom2: Second atom
            
        Returns:
            True if cyclization is possible
        """
        try:
            # Check if atoms are carbon (most common cyclization)
            if atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 6:
                return True
            
            # Allow cyclization with heteroatoms (N, O, S)
            heteroatoms = {6, 7, 8, 16}  # C, N, O, S
            if (atom1.GetAtomicNum() in heteroatoms and 
                atom2.GetAtomicNum() in heteroatoms):
                return True
            
            return False
            
        except:
            return False
    
    def _determine_bond_type(self, atom1: rdchem.Atom, atom2: rdchem.Atom) -> str:
        """
        Determine the type of bond to form in cyclization.
        
        Args:
            atom1: First atom
            atom2: Second atom
            
        Returns:
            Bond type string
        """
        # For now, default to single bonds
        # Could be enhanced with more sophisticated logic
        return 'single'
    
    def _estimate_strain_energy(self, ring_size: int) -> float:
        """
        Estimate ring strain energy based on ring size.
        
        Args:
            ring_size: Size of the ring
            
        Returns:
            Estimated strain energy in kcal/mol
        """
        # Approximate strain energies for different ring sizes
        strain_energies = {
            3: 27.5,  # Cyclopropane
            4: 26.3,  # Cyclobutane
            5: 6.2,   # Cyclopentane
            6: 0.0,   # Cyclohexane (reference)
            7: 6.2,   # Cycloheptane
            8: 9.7,   # Cyclooctane
            9: 12.8,  # Cyclononane
            10: 12.4, # Cyclodecane
            11: 11.3, # Cycloundecane
            12: 4.6   # Cyclododecane
        }
        
        return strain_energies.get(ring_size, 20.0)  # Default high strain for unusual sizes
