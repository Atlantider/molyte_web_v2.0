"""
Addition reaction operator for RSNet.

This module implements addition reactions where two molecules combine
to form a single product, including nucleophilic and electrophilic additions.
"""

from typing import List, Dict, Any, Optional, Tuple
from rdkit import Chem
from rdkit.Chem import AllChem
import copy

from .base import BaseOperator
from ..core.molecule import Molecule
from ..core.reaction import Reaction, ReactionTemplate
from ..core.environment import Environment


class AdditionTemplate(ReactionTemplate):
    """Template for addition reactions."""
    
    def __init__(
        self,
        nucleophile_pattern: str,
        electrophile_pattern: str,
        product_pattern: str,
        name: str = "Addition"
    ):
        """
        Initialize addition template.
        
        Args:
            nucleophile_pattern: SMARTS pattern for nucleophile
            electrophile_pattern: SMARTS pattern for electrophile
            product_pattern: SMARTS pattern for expected product
            name: Template name
        """
        super().__init__(
            name=name,
            description=f"Addition of nucleophile to electrophile",
            pattern=f"{nucleophile_pattern}.{electrophile_pattern}>>{product_pattern}"
        )
        self.nucleophile_pattern = nucleophile_pattern
        self.electrophile_pattern = electrophile_pattern
        self.product_pattern = product_pattern
    
    def apply(self, molecules: List[Molecule]) -> List[Reaction]:
        """
        Apply addition template to generate reactions.
        
        Args:
            molecules: List of molecules (should contain exactly two)
            
        Returns:
            List of generated reactions
        """
        if len(molecules) != 2:
            return []
        
        mol1, mol2 = molecules
        reactions = []
        
        # Try both orientations (mol1 as nucleophile, mol2 as electrophile and vice versa)
        for nucleophile, electrophile in [(mol1, mol2), (mol2, mol1)]:
            try:
                product = self._perform_addition(nucleophile, electrophile)
                if product is not None:
                    reaction = Reaction(
                        reactants=[nucleophile, electrophile],
                        products=[product],
                        name=f"{self.name}: {nucleophile.name} + {electrophile.name}"
                    )
                    reactions.append(reaction)
            except Exception as e:
                print(f"Warning: Addition failed: {e}")
                continue
        
        return reactions
    
    def _perform_addition(self, nucleophile: Molecule, electrophile: Molecule) -> Optional[Molecule]:
        """Perform the addition reaction."""
        # This is a simplified implementation
        # In practice, you would use more sophisticated reaction templates
        
        # Check if molecules match the patterns
        nuc_mol = nucleophile.rdkit_mol
        elec_mol = electrophile.rdkit_mol
        
        nuc_pattern = Chem.MolFromSmarts(self.nucleophile_pattern)
        elec_pattern = Chem.MolFromSmarts(self.electrophile_pattern)
        
        if nuc_pattern is None or elec_pattern is None:
            return None
        
        nuc_matches = nuc_mol.GetSubstructMatches(nuc_pattern)
        elec_matches = elec_mol.GetSubstructMatches(elec_pattern)
        
        if not nuc_matches or not elec_matches:
            return None
        
        # Simple combination - just combine the molecules
        # In practice, you would perform the actual bond formation
        combined_mol = Chem.CombineMols(nuc_mol, elec_mol)
        
        # Add a bond between reactive sites (simplified)
        editable_mol = Chem.EditableMol(combined_mol)
        
        # Get the first matching sites
        nuc_site = nuc_matches[0][0]  # First atom of first match
        elec_site = elec_matches[0][0] + nuc_mol.GetNumAtoms()  # Offset by nucleophile size
        
        # Add single bond
        editable_mol.AddBond(nuc_site, elec_site, Chem.BondType.SINGLE)
        
        product_mol = editable_mol.GetMol()
        
        if product_mol is None:
            return None
        
        # Sanitize the molecule
        try:
            Chem.SanitizeMol(product_mol)
            return Molecule(product_mol, name=f"addition_product")
        except:
            return None


class AdditionOperator(BaseOperator):
    """
    Operator for addition reactions.
    
    This operator identifies potential addition reactions between nucleophiles
    and electrophiles, including Michael additions, nucleophilic additions to
    carbonyls, and other addition mechanisms.
    """
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """
        Initialize addition operator.
        
        Args:
            config: Configuration parameters
        """
        super().__init__(
            name="Addition",
            description="Generates addition reactions between nucleophiles and electrophiles",
            config=config
        )
        
        # Define common addition templates
        self._setup_templates()
    
    def _setup_templates(self):
        """Set up common addition reaction templates."""
        
        # Nucleophilic addition to carbonyl
        carbonyl_addition = AdditionTemplate(
            nucleophile_pattern="[N,O,S,C-]",  # Nucleophiles with lone pairs or negative charge
            electrophile_pattern="[CX3]=[OX1]",  # Carbonyl carbon
            product_pattern="[N,O,S,C][CX4][OX2]",  # Addition product
            name="Carbonyl Addition"
        )
        self.add_template(carbonyl_addition)
        
        # Michael addition (1,4-addition to α,β-unsaturated carbonyls)
        michael_addition = AdditionTemplate(
            nucleophile_pattern="[C,N,O,S-]",  # Michael donors
            electrophile_pattern="[CX3]=[CX3][CX3]=[OX1]",  # α,β-unsaturated carbonyl
            product_pattern="[C,N,O,S][CX4][CX3]=[CX3][OX1]",  # 1,4-addition product
            name="Michael Addition"
        )
        self.add_template(michael_addition)
        
        # Nucleophilic addition to alkenes (with electron-withdrawing groups)
        alkene_addition = AdditionTemplate(
            nucleophile_pattern="[N,O,S-]",  # Nucleophiles
            electrophile_pattern="[CX3]=[CX3][C,N,O,S]",  # Activated alkene
            product_pattern="[N,O,S][CX4][CX4][C,N,O,S]",  # Addition product
            name="Alkene Addition"
        )
        self.add_template(alkene_addition)
        
        # Cycloaddition (simplified Diels-Alder type)
        cycloaddition = AdditionTemplate(
            nucleophile_pattern="[CX3]=[CX3][CX3]=[CX3]",  # Diene
            electrophile_pattern="[CX3]=[CX3]",  # Dienophile
            product_pattern="[CX4]1[CX4][CX4][CX4][CX4][CX4]1",  # Cyclohexene ring
            name="Cycloaddition"
        )
        self.add_template(cycloaddition)
    
    def can_apply(self, molecules: List[Molecule], environment: Environment) -> bool:
        """
        Check if addition operator can be applied.
        
        Args:
            molecules: List of molecules to check
            environment: Reaction environment
            
        Returns:
            True if operator can be applied
        """
        if len(molecules) < 1:
            return False
        
        # Get driving forces from environment
        drives = environment.get_active_drives()
        
        # Addition reactions are favored by:
        # - Nucleophilic/electrophilic environments
        # - Moderate temperatures
        # - Solution phase reactions
        
        required_drives = [
            drives.get('solution_phase', False),
            drives.get('thermal', False) or drives.get('electrochemical', False)
        ]
        
        return any(required_drives)
    
    def find_reaction_sites(self, molecules: List[Molecule]) -> List[Dict[str, Any]]:
        """
        Find potential addition reaction sites.
        
        Args:
            molecules: List of molecules to analyze
            
        Returns:
            List of reaction site information
        """
        if len(molecules) != 2:
            return []
        
        sites = []
        mol1, mol2 = molecules
        
        # Look for nucleophilic sites in mol1 and electrophilic sites in mol2
        nuc_sites = self._find_nucleophilic_sites(mol1)
        elec_sites = self._find_electrophilic_sites(mol2)
        
        for nuc_site in nuc_sites:
            for elec_site in elec_sites:
                sites.append({
                    'nucleophile_mol': mol1,
                    'electrophile_mol': mol2,
                    'nucleophile_site': nuc_site,
                    'electrophile_site': elec_site,
                    'reaction_type': 'addition'
                })
        
        # Also check the reverse orientation
        nuc_sites = self._find_nucleophilic_sites(mol2)
        elec_sites = self._find_electrophilic_sites(mol1)
        
        for nuc_site in nuc_sites:
            for elec_site in elec_sites:
                sites.append({
                    'nucleophile_mol': mol2,
                    'electrophile_mol': mol1,
                    'nucleophile_site': nuc_site,
                    'electrophile_site': elec_site,
                    'reaction_type': 'addition'
                })
        
        return sites
    
    def _find_nucleophilic_sites(self, molecule: Molecule) -> List[Dict[str, Any]]:
        """Find nucleophilic sites in a molecule."""
        sites = []
        mol = molecule.rdkit_mol
        
        for atom in mol.GetAtoms():
            is_nucleophilic = False
            reasons = []
            
            # Atoms with lone pairs
            if atom.GetSymbol() in ['N', 'O', 'S', 'P']:
                is_nucleophilic = True
                reasons.append('lone_pairs')
            
            # Negatively charged atoms
            if atom.GetFormalCharge() < 0:
                is_nucleophilic = True
                reasons.append('negative_charge')
            
            # Carbanions (carbon with negative charge or high electron density)
            if atom.GetSymbol() == 'C' and atom.GetFormalCharge() <= 0:
                # Check for electron-donating neighbors
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() in ['Li', 'Na', 'K', 'Mg']:
                        is_nucleophilic = True
                        reasons.append('organometallic')
            
            if is_nucleophilic:
                sites.append({
                    'atom_idx': atom.GetIdx(),
                    'atom_symbol': atom.GetSymbol(),
                    'reasons': reasons
                })
        
        return sites
    
    def _find_electrophilic_sites(self, molecule: Molecule) -> List[Dict[str, Any]]:
        """Find electrophilic sites in a molecule."""
        sites = []
        mol = molecule.rdkit_mol
        
        for atom in mol.GetAtoms():
            is_electrophilic = False
            reasons = []
            
            # Positively charged atoms
            if atom.GetFormalCharge() > 0:
                is_electrophilic = True
                reasons.append('positive_charge')
            
            # Carbonyl carbons
            if atom.GetSymbol() == 'C':
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'O':
                        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                        if bond.GetBondType() == Chem.BondType.DOUBLE:
                            is_electrophilic = True
                            reasons.append('carbonyl_carbon')
            
            # Alkene carbons with electron-withdrawing groups
            if atom.GetSymbol() == 'C' and atom.GetHybridization() == Chem.HybridizationType.SP2:
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() in ['O', 'N'] and neighbor.GetFormalCharge() >= 0:
                        is_electrophilic = True
                        reasons.append('activated_alkene')
            
            if is_electrophilic:
                sites.append({
                    'atom_idx': atom.GetIdx(),
                    'atom_symbol': atom.GetSymbol(),
                    'reasons': reasons
                })
        
        return sites
    
    def apply(self, molecules: List[Molecule], environment: Environment) -> List[Reaction]:
        """
        Apply addition operator to generate reactions.
        
        Args:
            molecules: List of molecules (should contain exactly two)
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
