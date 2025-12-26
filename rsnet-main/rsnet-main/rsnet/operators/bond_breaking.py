"""
Bond breaking operator for RSNet.

This module implements bond breaking reactions, which are fundamental
in many chemical processes including thermal decomposition and radical reactions.
"""

from typing import List, Dict, Any, Optional, Tuple
from rdkit import Chem
from rdkit.Chem import AllChem
import copy

from .base import BaseOperator
from ..core.molecule import Molecule
from ..core.reaction import Reaction, ReactionTemplate
from ..core.environment import Environment
# from ..utils.structure_analysis import get_structure_tags, detect_weak_bonds


class BondBreakingTemplate(ReactionTemplate):
    """
    Template for bond breaking reactions.
    """
    
    def __init__(
        self,
        bond_idx: int,
        atom1_idx: int,
        atom2_idx: int,
        bond_type: str = "single",
        name: Optional[str] = None
    ):
        """
        Initialize bond breaking template.
        
        Args:
            bond_idx: Index of bond to break
            atom1_idx: Index of first atom in bond
            atom2_idx: Index of second atom in bond
            bond_type: Type of bond (single, double, etc.)
            name: Optional template name
        """
        self.bond_idx = bond_idx
        self.atom1_idx = atom1_idx
        self.atom2_idx = atom2_idx
        self.bond_type = bond_type
        
        if name is None:
            name = f"bond_break_{atom1_idx}_{atom2_idx}"
        
        super().__init__(
            name=name,
            description=f"Break {bond_type} bond between atoms {atom1_idx} and {atom2_idx}",
            pattern=None
        )
    
    def apply(self, molecules: List[Molecule]) -> List[Reaction]:
        """
        Apply bond breaking template to generate reactions.
        
        Args:
            molecules: List of molecules (should contain exactly one)
            
        Returns:
            List of generated reactions
        """
        if len(molecules) != 1:
            return []
        
        reactant = molecules[0]
        
        try:
            # Create products by breaking the bond
            products = self._break_bond(reactant)
            if not products:
                return []
            
            # Create reaction
            reaction = Reaction(
                reactants=[reactant],
                products=products,
                name=f"Bond breaking: {reactant.name} -> fragments"
            )
            
            return [reaction]
            
        except Exception as e:
            print(f"Warning: Bond breaking failed: {e}")
            return []
    
    def _break_bond(self, mol: Molecule) -> List[Molecule]:
        """
        Break the specified bond and return fragment molecules.
        
        Args:
            mol: Input molecule
            
        Returns:
            List of fragment molecules, or empty list if failed
        """
        try:
            # Create a copy of the RDKit molecule
            rdkit_mol = Chem.Mol(mol.rdkit_mol)
            
            # Remove the bond
            editable_mol = Chem.EditableMol(rdkit_mol)
            editable_mol.RemoveBond(self.atom1_idx, self.atom2_idx)
            
            # Get the modified molecule
            broken_mol = editable_mol.GetMol()
            
            # Get fragments
            fragments = Chem.GetMolFrags(broken_mol, asMols=True, sanitizeFrags=True)
            
            if len(fragments) < 2:
                return []  # Bond breaking should create at least 2 fragments
            
            # Create Molecule objects for fragments
            products = []
            for i, frag in enumerate(fragments):
                frag_name = f"{mol.name}_frag_{i+1}"
                try:
                    product = Molecule(frag, name=frag_name)
                    products.append(product)
                except Exception as e:
                    print(f"Warning: Could not create fragment {i+1}: {e}")
                    continue
            
            return products
            
        except Exception as e:
            print(f"Error in bond breaking: {e}")
            return []


class BondBreakingOperator(BaseOperator):
    """
    Operator for bond breaking reactions.
    
    This operator identifies weak bonds that can break under given conditions
    and generates corresponding fragmentation reactions.
    """
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """
        Initialize bond breaking operator.
        
        Args:
            config: Configuration parameters
        """
        super().__init__(
            name="BondBreaking",
            description="Generates bond breaking reactions",
            config=config
        )
        
        # Configuration parameters
        self.min_fragment_size = self.config.get('min_fragment_size', 2)  # Minimum heavy atoms
        self.max_bonds_to_break = self.config.get('max_bonds_to_break', 5)
        self.prefer_weak_bonds = self.config.get('prefer_weak_bonds', True)
        self.energy_threshold = self.config.get('energy_threshold', 80.0)  # kcal/mol
    
    def is_applicable(self, molecules: List[Molecule]) -> bool:
        """
        Check if bond breaking is applicable.
        
        Args:
            molecules: Input molecules
            
        Returns:
            True if applicable
        """
        if not super().is_applicable(molecules):
            return False
        
        # Only handle single molecules
        if len(molecules) != 1:
            return False
        
        mol = molecules[0]
        
        # Need at least 2 heavy atoms to break bonds
        if mol.num_heavy_atoms < 2:
            return False
        
        # Check if molecule has breakable bonds
        rdkit_mol = mol.rdkit_mol
        num_bonds = rdkit_mol.GetNumBonds()
        
        return num_bonds > 0

    def can_apply(self, molecules: List[Molecule], environment: Environment) -> bool:
        """
        Check if bond breaking operator can be applied.

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

        # Bond breaking reactions are favored by:
        # - High temperature
        # - Thermal conditions
        # - Radical environments

        required_drives = [
            drives.get('thermal', False),
            drives.get('high_temperature', False),
            drives.get('radical_environment', False)
        ]

        return any(required_drives)

    def find_reaction_sites(self, molecules: List[Molecule]) -> List[Dict[str, Any]]:
        """
        Find potential bond breaking reaction sites.

        Args:
            molecules: List of molecules to analyze

        Returns:
            List of reaction site information
        """
        if len(molecules) != 1:
            return []

        sites = []
        molecule = molecules[0]
        mol = molecule.rdkit_mol

        # Find breakable bonds
        for bond in mol.GetBonds():
            bond_type = bond.GetBondType()
            sites.append({
                'molecule': molecule,
                'bond_idx': bond.GetIdx(),
                'atom1_idx': bond.GetBeginAtomIdx(),
                'atom2_idx': bond.GetEndAtomIdx(),
                'bond_type': str(bond_type),
                'site_type': 'breakable_bond'
            })

        return sites

    def apply(
        self, 
        molecules: List[Molecule], 
        environment: Environment
    ) -> List[Reaction]:
        """
        Apply bond breaking operator to generate reactions.
        
        Args:
            molecules: Input molecules
            environment: Reaction environment
            
        Returns:
            List of generated reactions
        """
        if not self.is_applicable(molecules):
            return []
        
        mol = molecules[0]
        reactions = []
        
        # Get breakable bonds
        breakable_bonds = self._find_breakable_bonds(mol, environment)
        
        # Limit number of bonds to consider
        if len(breakable_bonds) > self.max_bonds_to_break:
            # Sort by priority (weak bonds first if preferred)
            if self.prefer_weak_bonds:
                from ..utils.structure_analysis import get_structure_tags
                tags = get_structure_tags(mol)
                weak_bond_pairs = set(tags.get('weak_bonds', []))
                breakable_bonds.sort(key=lambda x: (x[0], x[1]) not in weak_bond_pairs)
            
            breakable_bonds = breakable_bonds[:self.max_bonds_to_break]
        
        # Create templates for each breakable bond
        for bond_info in breakable_bonds:
            bond_idx, atom1_idx, atom2_idx, bond_type = bond_info
            
            template = BondBreakingTemplate(
                bond_idx=bond_idx,
                atom1_idx=atom1_idx,
                atom2_idx=atom2_idx,
                bond_type=bond_type
            )
            
            # Apply template
            bond_reactions = template.apply(molecules)
            
            # Filter products by minimum fragment size
            filtered_reactions = []
            for rxn in bond_reactions:
                valid = True
                for product in rxn.products:
                    if product.num_heavy_atoms < self.min_fragment_size:
                        valid = False
                        break
                if valid:
                    filtered_reactions.append(rxn)
            
            reactions.extend(filtered_reactions)
        
        # Filter reactions by energy threshold
        filtered_reactions = self.filter_reactions(reactions)
        
        return filtered_reactions
    
    def _find_breakable_bonds(
        self, 
        mol: Molecule, 
        environment: Environment
    ) -> List[Tuple[int, int, int, str]]:
        """
        Find bonds that can be broken under the given conditions.
        
        Args:
            mol: Input molecule
            environment: Reaction environment
            
        Returns:
            List of tuples (bond_idx, atom1_idx, atom2_idx, bond_type)
        """
        breakable_bonds = []
        rdkit_mol = mol.rdkit_mol
        
        # Get structure tags for weak bond identification
        from ..utils.structure_analysis import get_structure_tags, detect_weak_bonds
        try:
            tags = get_structure_tags(mol)
            weak_bond_pairs = set(tags.get('weak_bonds', []))
        except Exception:
            # Fallback if structure analysis fails
            tags = {}
            weak_bond_pairs = set()
        
        for bond in rdkit_mol.GetBonds():
            bond_idx = bond.GetIdx()
            atom1_idx = bond.GetBeginAtomIdx()
            atom2_idx = bond.GetEndAtomIdx()
            bond_type = str(bond.GetBondType()).lower()
            
            # Check if bond is breakable
            is_breakable = False
            
            # Always consider weak bonds
            if (atom1_idx, atom2_idx) in weak_bond_pairs or (atom2_idx, atom1_idx) in weak_bond_pairs:
                is_breakable = True
            
            # Consider single bonds at high temperature
            elif bond.GetBondType() == Chem.BondType.SINGLE and environment.temperature > 400:
                is_breakable = True
            
            # Consider bonds in small rings (ring strain)
            elif bond.IsInRingSize(3) or bond.IsInRingSize(4):
                is_breakable = True
            
            # Consider bonds involving heteroatoms
            elif (bond.GetBeginAtom().GetAtomicNum() not in [1, 6] or 
                  bond.GetEndAtom().GetAtomicNum() not in [1, 6]):
                is_breakable = True
            
            if is_breakable:
                breakable_bonds.append((bond_idx, atom1_idx, atom2_idx, bond_type))
        
        return breakable_bonds
    
    def get_templates(self) -> List[ReactionTemplate]:
        """
        Get reaction templates for this operator.
        
        Note: This operator generates templates dynamically based on input molecules.
        
        Returns:
            Empty list (templates are generated dynamically)
        """
        return []
    
    def get_required_drives(self) -> List[str]:
        """
        Get the driving forces required for this operator.
        
        Returns:
            List of required driving force names
        """
        return ['thermal', 'radical_environment']
    
    def check_structure_conditions(self, tags: Dict[str, Any]) -> bool:
        """
        Check if molecular structure satisfies conditions for this operator.
        
        Args:
            tags: Structure tags from get_structure_tags()
            
        Returns:
            True if structure conditions are met
        """
        # Need at least one bond to break
        return tags.get('has_weak_bonds', False) or len(tags.get('weak_bonds', [])) > 0
