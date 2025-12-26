"""
Hydrogen transfer operator for RSNet.

This module implements hydrogen transfer reactions, which are fundamental
in many chemical processes including acid-base reactions and proton transfers.
"""

from typing import List, Dict, Any, Optional
from rdkit import Chem
from rdkit.Chem import AllChem
import copy

from .base import BaseOperator
from ..core.molecule import Molecule
from ..core.reaction import Reaction, ReactionTemplate
from ..core.environment import Environment
# from ..utils.structure_analysis import find_hydrogen_transfer_sites, get_structure_tags


class HydrogenTransferTemplate(ReactionTemplate):
    """
    Template for hydrogen transfer reactions.
    """
    
    def __init__(
        self,
        donor_heavy_idx: int,
        donor_hydrogen_idx: int,
        acceptor_idx: int,
        name: Optional[str] = None
    ):
        """
        Initialize hydrogen transfer template.
        
        Args:
            donor_heavy_idx: Index of heavy atom donating hydrogen
            donor_hydrogen_idx: Index of hydrogen being transferred
            acceptor_idx: Index of acceptor atom
            name: Optional template name
        """
        self.donor_heavy_idx = donor_heavy_idx
        self.donor_hydrogen_idx = donor_hydrogen_idx
        self.acceptor_idx = acceptor_idx
        
        if name is None:
            name = f"H_transfer_{donor_heavy_idx}_{acceptor_idx}"
        
        super().__init__(
            name=name,
            description=f"Hydrogen transfer from atom {donor_heavy_idx} to atom {acceptor_idx}",
            pattern=None  # We'll use programmatic approach instead of SMARTS
        )
    
    def apply(self, molecules: List[Molecule]) -> List[Reaction]:
        """
        Apply hydrogen transfer template to generate reactions.
        
        Args:
            molecules: List of molecules (should contain exactly one for intramolecular transfer)
            
        Returns:
            List of generated reactions
        """
        if len(molecules) != 1:
            return []  # For now, only handle intramolecular transfers
        
        reactant = molecules[0]
        
        try:
            # Create product by transferring hydrogen
            product = self._transfer_hydrogen(reactant)
            if product is None:
                return []
            
            # Create reaction
            reaction = Reaction(
                reactants=[reactant],
                products=[product],
                name=f"H-transfer: {reactant.name} -> {product.name}"
            )
            
            return [reaction]
            
        except Exception as e:
            print(f"Warning: Hydrogen transfer failed: {e}")
            return []
    
    def _transfer_hydrogen(self, mol: Molecule) -> Optional[Molecule]:
        """
        Perform the actual hydrogen transfer on the molecule.

        This is a simplified implementation that uses SMARTS-based transformations
        for common hydrogen transfer patterns.

        Args:
            mol: Input molecule

        Returns:
            Product molecule with transferred hydrogen, or None if failed
        """
        try:
            # For now, we'll implement a very simple case: intramolecular proton transfer
            # This is a placeholder - a full implementation would need more sophisticated
            # reaction templates and proper handling of electron movement

            rdkit_mol = mol.rdkit_mol
            donor_atom = rdkit_mol.GetAtomWithIdx(self.donor_heavy_idx)
            acceptor_atom = rdkit_mol.GetAtomWithIdx(self.acceptor_idx)

            # Check if this is a reasonable transfer
            # For now, only allow O-H → O transfers (like in carboxylic acids)
            if (donor_atom.GetSymbol() == 'O' and acceptor_atom.GetSymbol() == 'O'):
                # This would be a tautomerization - very complex to implement correctly
                # For now, return None to indicate this needs more sophisticated handling
                return None

            # For other cases, also return None for now
            # A proper implementation would need:
            # 1. Reaction templates for specific transfer types
            # 2. Proper handling of formal charges
            # 3. Consideration of resonance structures
            # 4. Energy calculations to determine feasibility

            return None

        except Exception as e:
            print(f"Error in hydrogen transfer: {e}")
            return None


class HydrogenTransferOperator(BaseOperator):
    """
    Operator for hydrogen transfer reactions.
    
    This operator identifies potential hydrogen transfer sites and generates
    corresponding reactions.
    """
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """
        Initialize hydrogen transfer operator.
        
        Args:
            config: Configuration parameters
        """
        super().__init__(
            name="HydrogenTransfer",
            description="Generates hydrogen transfer reactions",
            config=config
        )
        
        # Configuration parameters
        self.max_distance = self.config.get('max_distance', 4.0)  # Angstroms
        self.min_distance = self.config.get('min_distance', 1.5)  # Angstroms
        self.require_3d = self.config.get('require_3d', False)

    def can_apply(self, molecules: List[Molecule], environment: Environment) -> bool:
        """
        Check if hydrogen transfer operator can be applied.

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

        # Hydrogen transfer reactions are favored by:
        # - Thermal conditions
        # - Radical environments
        # - High temperature

        required_drives = [
            drives.get('thermal', False),
            drives.get('radical_environment', False),
            drives.get('high_temperature', False)
        ]

        return any(required_drives)

    def find_reaction_sites(self, molecules: List[Molecule]) -> List[Dict[str, Any]]:
        """
        Find potential hydrogen transfer reaction sites.

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

        # Find atoms with transferable hydrogens
        for atom in mol.GetAtoms():
            if atom.GetTotalNumHs() > 0:
                sites.append({
                    'molecule': molecule,
                    'atom_idx': atom.GetIdx(),
                    'site_type': 'hydrogen_donor',
                    'num_hydrogens': atom.GetTotalNumHs()
                })

        return sites

    def is_applicable(self, molecules: List[Molecule]) -> bool:
        """
        Check if hydrogen transfer is applicable.
        
        Args:
            molecules: Input molecules
            
        Returns:
            True if applicable
        """
        if not super().is_applicable(molecules):
            return False
        
        # For now, only handle single molecules (intramolecular transfer)
        if len(molecules) != 1:
            return False
        
        mol = molecules[0]
        from ..utils.structure_analysis import get_structure_tags
        tags = get_structure_tags(mol)
        
        # Need both acidic hydrogens and acceptors
        return tags['has_acidic_hydrogens'] and tags['has_heteroatoms']
    
    def apply(
        self, 
        molecules: List[Molecule], 
        environment: Environment
    ) -> List[Reaction]:
        """
        Apply hydrogen transfer operator to generate reactions.
        
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
        
        # Find hydrogen transfer sites
        from ..utils.structure_analysis import find_hydrogen_transfer_sites
        sites = find_hydrogen_transfer_sites(mol)
        
        for site in sites:
            # Check distance constraint if 3D coordinates available
            if site['distance'] is not None:
                if (site['distance'] < self.min_distance or 
                    site['distance'] > self.max_distance):
                    continue
            elif self.require_3d:
                continue  # Skip if 3D required but not available
            
            # Create template for this site
            template = HydrogenTransferTemplate(
                donor_heavy_idx=site['donor_heavy'],
                donor_hydrogen_idx=site['donor_hydrogen'],
                acceptor_idx=site['acceptor']
            )
            
            # Apply template
            site_reactions = template.apply(molecules)
            reactions.extend(site_reactions)
        
        # Filter reactions
        filtered_reactions = self.filter_reactions(reactions)
        
        return filtered_reactions
    
    def get_templates(self) -> List[ReactionTemplate]:
        """
        Get reaction templates for this operator.
        
        Note: This operator generates templates dynamically based on input molecules,
        so this method returns an empty list.
        
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
        return ['thermal', 'polar_environment']
    
    def check_structure_conditions(self, tags: Dict[str, Any]) -> bool:
        """
        Check if molecular structure satisfies conditions for this operator.
        
        Args:
            tags: Structure tags from get_structure_tags()
            
        Returns:
            True if structure conditions are met
        """
        return tags['has_acidic_hydrogens'] and tags['has_heteroatoms']
