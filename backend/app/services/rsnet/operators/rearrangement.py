"""
Rearrangement reaction operator for RSNet.

This module implements molecular rearrangement reactions including
sigmatropic rearrangements, ring expansions/contractions, and
skeletal rearrangements.
"""

from typing import List, Dict, Any, Optional, Tuple
from rdkit import Chem
from rdkit.Chem import AllChem
import copy

from .base import BaseOperator
from ..core.molecule import Molecule
from ..core.reaction import Reaction, ReactionTemplate
from ..core.environment import Environment


class RearrangementTemplate(ReactionTemplate):
    """Template for rearrangement reactions."""
    
    def __init__(
        self,
        reactant_pattern: str,
        product_pattern: str,
        rearrangement_type: str,
        name: str = "Rearrangement"
    ):
        """
        Initialize rearrangement template.
        
        Args:
            reactant_pattern: SMARTS pattern for reactant
            product_pattern: SMARTS pattern for product
            rearrangement_type: Type of rearrangement
            name: Template name
        """
        super().__init__(
            name=name,
            description=f"{rearrangement_type} rearrangement",
            pattern=f"{reactant_pattern}>>{product_pattern}"
        )
        self.reactant_pattern = reactant_pattern
        self.product_pattern = product_pattern
        self.rearrangement_type = rearrangement_type
    
    def apply(self, molecules: List[Molecule]) -> List[Reaction]:
        """
        Apply rearrangement template to generate reactions.
        
        Args:
            molecules: List of molecules (should contain exactly one)
            
        Returns:
            List of generated reactions
        """
        if len(molecules) != 1:
            return []
        
        reactant = molecules[0]
        
        try:
            product = self._perform_rearrangement(reactant)
            if product is not None:
                reaction = Reaction(
                    reactants=[reactant],
                    products=[product],
                    name=f"{self.name}: {reactant.name} rearrangement"
                )
                return [reaction]
        except Exception as e:
            print(f"Warning: Rearrangement failed: {e}")
        
        return []
    
    def _perform_rearrangement(self, reactant: Molecule) -> Optional[Molecule]:
        """Perform the rearrangement reaction."""
        mol = reactant.rdkit_mol
        
        # Check if molecule matches the reactant pattern
        pattern = Chem.MolFromSmarts(self.reactant_pattern)
        if pattern is None:
            return None
        
        matches = mol.GetSubstructMatches(pattern)
        if not matches:
            return None
        
        # Perform rearrangement based on type
        if self.rearrangement_type == "ring_expansion":
            return self._ring_expansion(mol, matches[0])
        elif self.rearrangement_type == "ring_contraction":
            return self._ring_contraction(mol, matches[0])
        elif self.rearrangement_type == "sigmatropic":
            return self._sigmatropic_rearrangement(mol, matches[0])
        elif self.rearrangement_type == "skeletal":
            return self._skeletal_rearrangement(mol, matches[0])
        else:
            return None
    
    def _ring_expansion(self, mol: Chem.Mol, match: Tuple[int, ...]) -> Optional[Molecule]:
        """Perform ring expansion rearrangement."""
        # Simplified ring expansion - insert a carbon into a ring
        editable_mol = Chem.EditableMol(mol)
        
        # Find a ring bond to break and insert carbon
        ring_info = mol.GetRingInfo()
        for ring in ring_info.AtomRings():
            if len(ring) < 7:  # Only expand small rings
                # Break a bond in the ring and insert carbon
                atom1_idx = ring[0]
                atom2_idx = ring[1]
                
                # Remove the bond
                editable_mol.RemoveBond(atom1_idx, atom2_idx)
                
                # Add new carbon
                new_carbon_idx = editable_mol.AddAtom(Chem.Atom(6))  # Carbon
                
                # Connect new carbon to both atoms
                editable_mol.AddBond(atom1_idx, new_carbon_idx, Chem.BondType.SINGLE)
                editable_mol.AddBond(new_carbon_idx, atom2_idx, Chem.BondType.SINGLE)
                
                break
        
        try:
            product_mol = editable_mol.GetMol()
            Chem.SanitizeMol(product_mol)
            return Molecule(product_mol, name="ring_expanded")
        except:
            return None
    
    def _ring_contraction(self, mol: Chem.Mol, match: Tuple[int, ...]) -> Optional[Molecule]:
        """Perform ring contraction rearrangement."""
        # Simplified ring contraction - remove an atom from a ring
        editable_mol = Chem.EditableMol(mol)
        
        ring_info = mol.GetRingInfo()
        for ring in ring_info.AtomRings():
            if len(ring) > 5:  # Only contract larger rings
                # Remove an atom from the ring
                atom_to_remove = ring[0]
                
                # Get neighbors of the atom to remove
                atom = mol.GetAtomWithIdx(atom_to_remove)
                neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
                
                if len(neighbors) == 2:
                    # Connect the neighbors directly
                    editable_mol.AddBond(neighbors[0], neighbors[1], Chem.BondType.SINGLE)
                
                # Remove the atom
                editable_mol.RemoveAtom(atom_to_remove)
                break
        
        try:
            product_mol = editable_mol.GetMol()
            Chem.SanitizeMol(product_mol)
            return Molecule(product_mol, name="ring_contracted")
        except:
            return None
    
    def _sigmatropic_rearrangement(self, mol: Chem.Mol, match: Tuple[int, ...]) -> Optional[Molecule]:
        """Perform sigmatropic rearrangement (simplified)."""
        # This is a very simplified implementation
        # Real sigmatropic rearrangements are much more complex
        
        editable_mol = Chem.EditableMol(mol)
        
        # Find a suitable bond to migrate
        for bond in mol.GetBonds():
            if bond.GetBondType() == Chem.BondType.SINGLE:
                atom1_idx = bond.GetBeginAtomIdx()
                atom2_idx = bond.GetEndAtomIdx()
                
                # Remove the bond
                editable_mol.RemoveBond(atom1_idx, atom2_idx)
                
                # Try to form a new bond elsewhere (simplified)
                atom1 = mol.GetAtomWithIdx(atom1_idx)
                for neighbor in atom1.GetNeighbors():
                    if neighbor.GetIdx() != atom2_idx:
                        neighbor_neighbors = [n.GetIdx() for n in neighbor.GetNeighbors()]
                        for nn_idx in neighbor_neighbors:
                            if nn_idx != atom1_idx and nn_idx != atom2_idx:
                                # Form new bond
                                editable_mol.AddBond(atom2_idx, nn_idx, Chem.BondType.SINGLE)
                                break
                        break
                break
        
        try:
            product_mol = editable_mol.GetMol()
            Chem.SanitizeMol(product_mol)
            return Molecule(product_mol, name="sigmatropic_product")
        except:
            return None
    
    def _skeletal_rearrangement(self, mol: Chem.Mol, match: Tuple[int, ...]) -> Optional[Molecule]:
        """Perform skeletal rearrangement."""
        # Simplified skeletal rearrangement - rearrange carbon skeleton
        editable_mol = Chem.EditableMol(mol)
        
        # Find carbon-carbon bonds to rearrange
        carbon_bonds = []
        for bond in mol.GetBonds():
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            if atom1.GetSymbol() == 'C' and atom2.GetSymbol() == 'C':
                carbon_bonds.append(bond)
        
        if len(carbon_bonds) >= 2:
            # Break one bond and form another
            bond1 = carbon_bonds[0]
            bond2 = carbon_bonds[1]
            
            # Remove first bond
            editable_mol.RemoveBond(bond1.GetBeginAtomIdx(), bond1.GetEndAtomIdx())
            
            # Form new bond between different carbons
            editable_mol.AddBond(
                bond1.GetBeginAtomIdx(), 
                bond2.GetBeginAtomIdx(), 
                Chem.BondType.SINGLE
            )
        
        try:
            product_mol = editable_mol.GetMol()
            Chem.SanitizeMol(product_mol)
            return Molecule(product_mol, name="skeletal_rearranged")
        except:
            return None


class RearrangementOperator(BaseOperator):
    """
    Operator for rearrangement reactions.
    
    This operator identifies potential rearrangement reactions including
    ring expansions/contractions, sigmatropic rearrangements, and
    skeletal rearrangements.
    """
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """
        Initialize rearrangement operator.
        
        Args:
            config: Configuration parameters
        """
        super().__init__(
            name="Rearrangement",
            description="Generates molecular rearrangement reactions",
            config=config
        )
        
        # Configuration parameters
        self.min_ring_size = self.config.get('min_ring_size', 3)
        self.max_ring_size = self.config.get('max_ring_size', 8)
        self.allow_skeletal = self.config.get('allow_skeletal', True)
        
        # Define rearrangement templates
        self._setup_templates()
    
    def _setup_templates(self):
        """Set up rearrangement reaction templates."""
        
        # Ring expansion (small ring to larger ring)
        ring_expansion = RearrangementTemplate(
            reactant_pattern="[C,c]1[C,c][C,c][C,c]1",  # 4-membered ring
            product_pattern="[C,c]1[C,c][C,c][C,c][C,c]1",  # 5-membered ring
            rearrangement_type="ring_expansion",
            name="Ring Expansion"
        )
        self.add_template(ring_expansion)
        
        # Ring contraction (larger ring to smaller ring)
        ring_contraction = RearrangementTemplate(
            reactant_pattern="[C,c]1[C,c][C,c][C,c][C,c][C,c]1",  # 6-membered ring
            product_pattern="[C,c]1[C,c][C,c][C,c][C,c]1",  # 5-membered ring
            rearrangement_type="ring_contraction",
            name="Ring Contraction"
        )
        self.add_template(ring_contraction)
        
        # Sigmatropic rearrangement (simplified)
        sigmatropic = RearrangementTemplate(
            reactant_pattern="[C]=[C][C][C]=[C]",  # Conjugated system
            product_pattern="[C][C]=[C][C][C]",  # Rearranged system
            rearrangement_type="sigmatropic",
            name="Sigmatropic Rearrangement"
        )
        self.add_template(sigmatropic)
        
        if self.allow_skeletal:
            # Skeletal rearrangement
            skeletal = RearrangementTemplate(
                reactant_pattern="[C][C][C][C]",  # Carbon chain
                product_pattern="[C]([C])[C][C]",  # Branched chain
                rearrangement_type="skeletal",
                name="Skeletal Rearrangement"
            )
            self.add_template(skeletal)
    
    def can_apply(self, molecules: List[Molecule], environment: Environment) -> bool:
        """
        Check if rearrangement operator can be applied.
        
        Args:
            molecules: List of molecules to check
            environment: Reaction environment
            
        Returns:
            True if operator can be applied
        """
        if len(molecules) != 1:
            return False
        
        # Get driving forces from environment
        drives = environment.get_active_drives()
        
        # Rearrangement reactions are favored by:
        # - High temperature (thermal rearrangements)
        # - Radical environments
        # - Ring strain relief
        
        required_drives = [
            drives.get('thermal', False),
            drives.get('high_temperature', False),
            drives.get('radical_environment', False)
        ]
        
        return any(required_drives)
    
    def find_reaction_sites(self, molecules: List[Molecule]) -> List[Dict[str, Any]]:
        """
        Find potential rearrangement reaction sites.
        
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
        
        # Look for rings that can expand or contract
        ring_info = mol.GetRingInfo()
        for ring in ring_info.AtomRings():
            ring_size = len(ring)
            
            if ring_size < 6:  # Small rings can expand
                sites.append({
                    'molecule': molecule,
                    'ring_atoms': ring,
                    'ring_size': ring_size,
                    'rearrangement_type': 'ring_expansion',
                    'driving_force': 'ring_strain_relief'
                })
            
            if ring_size > 5:  # Larger rings can contract
                sites.append({
                    'molecule': molecule,
                    'ring_atoms': ring,
                    'ring_size': ring_size,
                    'rearrangement_type': 'ring_contraction',
                    'driving_force': 'entropy_increase'
                })
        
        # Look for conjugated systems (potential sigmatropic sites)
        for atom in mol.GetAtoms():
            if atom.GetIsAromatic() or atom.GetHybridization() == Chem.HybridizationType.SP2:
                sites.append({
                    'molecule': molecule,
                    'atom_idx': atom.GetIdx(),
                    'rearrangement_type': 'sigmatropic',
                    'driving_force': 'conjugation_stabilization'
                })
        
        return sites
    
    def apply(self, molecules: List[Molecule], environment: Environment) -> List[Reaction]:
        """
        Apply rearrangement operator to generate reactions.
        
        Args:
            molecules: List of molecules (should contain exactly one)
            environment: Reaction environment
            
        Returns:
            List of generated reactions
        """
        if not self.can_apply(molecules, environment):
            return []
        
        reactions = []
        
        # Apply each template
        for template in self.get_templates():
            template_reactions = template.apply(molecules)
            reactions.extend(template_reactions)
        
        return reactions
