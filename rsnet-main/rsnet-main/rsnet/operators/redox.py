"""
Redox (oxidation-reduction) reaction operator for RSNet.

This module implements oxidation and reduction reactions, including
electron transfer, radical formation, and electrochemical processes.
"""

from typing import List, Dict, Any, Optional, Tuple
from rdkit import Chem
from rdkit.Chem import AllChem
import copy

from .base import BaseOperator
from ..core.molecule import Molecule
from ..core.reaction import Reaction, ReactionTemplate
from ..core.environment import Environment


class RedoxTemplate(ReactionTemplate):
    """Template for redox reactions."""
    
    def __init__(
        self,
        reactant_pattern: str,
        product_pattern: str,
        redox_type: str,
        electron_change: int,
        name: str = "Redox"
    ):
        """
        Initialize redox template.
        
        Args:
            reactant_pattern: SMARTS pattern for reactant
            product_pattern: SMARTS pattern for product
            redox_type: Type of redox reaction ("oxidation" or "reduction")
            electron_change: Number of electrons transferred (positive for oxidation)
            name: Template name
        """
        super().__init__(
            name=name,
            description=f"{redox_type} reaction with {abs(electron_change)} electron(s)",
            pattern=f"{reactant_pattern}>>{product_pattern}"
        )
        self.reactant_pattern = reactant_pattern
        self.product_pattern = product_pattern
        self.redox_type = redox_type
        self.electron_change = electron_change
    
    def apply(self, molecules: List[Molecule]) -> List[Reaction]:
        """
        Apply redox template to generate reactions.
        
        Args:
            molecules: List of molecules (should contain exactly one)
            
        Returns:
            List of generated reactions
        """
        if len(molecules) != 1:
            return []
        
        reactant = molecules[0]
        
        try:
            product = self._perform_redox(reactant)
            if product is not None:
                reaction = Reaction(
                    reactants=[reactant],
                    products=[product],
                    name=f"{self.name}: {reactant.name} {self.redox_type}"
                )
                return [reaction]
        except Exception as e:
            print(f"Warning: Redox reaction failed: {e}")
        
        return []
    
    def _perform_redox(self, reactant: Molecule) -> Optional[Molecule]:
        """Perform the redox reaction."""
        mol = reactant.rdkit_mol
        
        # Check if molecule matches the reactant pattern
        pattern = Chem.MolFromSmarts(self.reactant_pattern)
        if pattern is None:
            return None
        
        matches = mol.GetSubstructMatches(pattern)
        if not matches:
            return None
        
        # Perform redox reaction
        if self.redox_type == "oxidation":
            return self._perform_oxidation(mol, matches[0])
        elif self.redox_type == "reduction":
            return self._perform_reduction(mol, matches[0])
        else:
            return None
    
    def _perform_oxidation(self, mol: Chem.Mol, match: Tuple[int, ...]) -> Optional[Molecule]:
        """Perform oxidation reaction."""
        editable_mol = Chem.EditableMol(mol)
        
        # Simple oxidation: increase formal charge or form double bonds
        target_atom_idx = match[0]
        target_atom = mol.GetAtomWithIdx(target_atom_idx)
        
        # Create new atom with increased formal charge
        new_atom = Chem.Atom(target_atom.GetAtomicNum())
        new_atom.SetFormalCharge(target_atom.GetFormalCharge() + self.electron_change)
        new_atom.SetNumRadicalElectrons(target_atom.GetNumRadicalElectrons())
        
        # Replace the atom
        editable_mol.ReplaceAtom(target_atom_idx, new_atom)
        
        try:
            product_mol = editable_mol.GetMol()
            Chem.SanitizeMol(product_mol)
            return Molecule(product_mol, name="oxidized")
        except:
            # If sanitization fails, try alternative oxidation
            return self._alternative_oxidation(mol, match)
    
    def _perform_reduction(self, mol: Chem.Mol, match: Tuple[int, ...]) -> Optional[Molecule]:
        """Perform reduction reaction."""
        editable_mol = Chem.EditableMol(mol)
        
        # Simple reduction: decrease formal charge or break double bonds
        target_atom_idx = match[0]
        target_atom = mol.GetAtomWithIdx(target_atom_idx)
        
        # Create new atom with decreased formal charge
        new_atom = Chem.Atom(target_atom.GetAtomicNum())
        new_atom.SetFormalCharge(target_atom.GetFormalCharge() - abs(self.electron_change))
        new_atom.SetNumRadicalElectrons(target_atom.GetNumRadicalElectrons())
        
        # Replace the atom
        editable_mol.ReplaceAtom(target_atom_idx, new_atom)
        
        try:
            product_mol = editable_mol.GetMol()
            Chem.SanitizeMol(product_mol)
            return Molecule(product_mol, name="reduced")
        except:
            # If sanitization fails, try alternative reduction
            return self._alternative_reduction(mol, match)
    
    def _alternative_oxidation(self, mol: Chem.Mol, match: Tuple[int, ...]) -> Optional[Molecule]:
        """Alternative oxidation mechanism (bond formation)."""
        editable_mol = Chem.EditableMol(mol)
        
        # Try to form a double bond (oxidation often involves bond formation)
        target_atom_idx = match[0]
        target_atom = mol.GetAtomWithIdx(target_atom_idx)
        
        for neighbor in target_atom.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(target_atom_idx, neighbor.GetIdx())
            if bond.GetBondType() == Chem.BondType.SINGLE:
                # Convert single bond to double bond
                editable_mol.RemoveBond(target_atom_idx, neighbor.GetIdx())
                editable_mol.AddBond(target_atom_idx, neighbor.GetIdx(), Chem.BondType.DOUBLE)
                break
        
        try:
            product_mol = editable_mol.GetMol()
            Chem.SanitizeMol(product_mol)
            return Molecule(product_mol, name="oxidized_alt")
        except:
            return None
    
    def _alternative_reduction(self, mol: Chem.Mol, match: Tuple[int, ...]) -> Optional[Molecule]:
        """Alternative reduction mechanism (bond breaking)."""
        editable_mol = Chem.EditableMol(mol)
        
        # Try to break a double bond (reduction often involves bond breaking)
        target_atom_idx = match[0]
        target_atom = mol.GetAtomWithIdx(target_atom_idx)
        
        for neighbor in target_atom.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(target_atom_idx, neighbor.GetIdx())
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                # Convert double bond to single bond
                editable_mol.RemoveBond(target_atom_idx, neighbor.GetIdx())
                editable_mol.AddBond(target_atom_idx, neighbor.GetIdx(), Chem.BondType.SINGLE)
                break
        
        try:
            product_mol = editable_mol.GetMol()
            Chem.SanitizeMol(product_mol)
            return Molecule(product_mol, name="reduced_alt")
        except:
            return None


class RedoxOperator(BaseOperator):
    """
    Operator for redox (oxidation-reduction) reactions.
    
    This operator identifies potential redox reactions including
    electron transfer, radical formation, and electrochemical processes.
    """
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """
        Initialize redox operator.
        
        Args:
            config: Configuration parameters
        """
        super().__init__(
            name="Redox",
            description="Generates oxidation and reduction reactions",
            config=config
        )
        
        # Configuration parameters
        self.max_electron_transfer = self.config.get('max_electron_transfer', 2)
        self.allow_radical_formation = self.config.get('allow_radical_formation', True)
        
        # Define redox templates
        self._setup_templates()
    
    def _setup_templates(self):
        """Set up redox reaction templates."""
        
        # Oxidation of alcohols to aldehydes/ketones
        alcohol_oxidation = RedoxTemplate(
            reactant_pattern="[CX4][OH]",  # Primary/secondary alcohol
            product_pattern="[CX3]=[OX1]",  # Aldehyde/ketone
            redox_type="oxidation",
            electron_change=2,
            name="Alcohol Oxidation"
        )
        self.add_template(alcohol_oxidation)
        
        # Reduction of carbonyls to alcohols
        carbonyl_reduction = RedoxTemplate(
            reactant_pattern="[CX3]=[OX1]",  # Aldehyde/ketone
            product_pattern="[CX4][OH]",  # Alcohol
            redox_type="reduction",
            electron_change=-2,
            name="Carbonyl Reduction"
        )
        self.add_template(carbonyl_reduction)
        
        # Oxidation of alkenes to epoxides
        alkene_oxidation = RedoxTemplate(
            reactant_pattern="[CX3]=[CX3]",  # Alkene
            product_pattern="[CX4]1[OX2][CX4]1",  # Epoxide
            redox_type="oxidation",
            electron_change=2,
            name="Alkene Epoxidation"
        )
        self.add_template(alkene_oxidation)
        
        # Reduction of nitro groups to amines
        nitro_reduction = RedoxTemplate(
            reactant_pattern="[NX3+](=O)[O-]",  # Nitro group
            product_pattern="[NX3]",  # Amine
            redox_type="reduction",
            electron_change=-6,
            name="Nitro Reduction"
        )
        self.add_template(nitro_reduction)
        
        if self.allow_radical_formation:
            # Single electron oxidation (radical formation)
            radical_oxidation = RedoxTemplate(
                reactant_pattern="[C,N,O,S]",  # Various atoms
                product_pattern="[C,N,O,S+]",  # Radical cation
                redox_type="oxidation",
                electron_change=1,
                name="Radical Formation"
            )
            self.add_template(radical_oxidation)
    
    def can_apply(self, molecules: List[Molecule], environment: Environment) -> bool:
        """
        Check if redox operator can be applied.
        
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
        
        # Redox reactions are favored by:
        # - Oxidizing or reducing environments
        # - Electrochemical conditions
        # - High temperature (for some redox reactions)
        
        required_drives = [
            drives.get('oxidation', False),
            drives.get('reduction', False),
            drives.get('electrochemical', False),
            drives.get('radical_environment', False)
        ]
        
        return any(required_drives)
    
    def find_reaction_sites(self, molecules: List[Molecule]) -> List[Dict[str, Any]]:
        """
        Find potential redox reaction sites.
        
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
        
        # Look for oxidizable sites
        oxidizable_sites = self._find_oxidizable_sites(mol)
        for site in oxidizable_sites:
            sites.append({
                'molecule': molecule,
                'atom_idx': site['atom_idx'],
                'redox_type': 'oxidation',
                'site_type': site['site_type'],
                'electron_change': site.get('electron_change', 1)
            })
        
        # Look for reducible sites
        reducible_sites = self._find_reducible_sites(mol)
        for site in reducible_sites:
            sites.append({
                'molecule': molecule,
                'atom_idx': site['atom_idx'],
                'redox_type': 'reduction',
                'site_type': site['site_type'],
                'electron_change': site.get('electron_change', -1)
            })
        
        return sites
    
    def _find_oxidizable_sites(self, mol: Chem.Mol) -> List[Dict[str, Any]]:
        """Find sites that can be oxidized."""
        sites = []

        for atom in mol.GetAtoms():
            site_type = None
            electron_change = 1

            # Alcohols (can be oxidized to aldehydes/ketones)
            if atom.GetSymbol() == 'C':
                # Check for C-OH pattern (alcohol carbon)
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'O' and neighbor.GetTotalNumHs() == 1:
                        # This carbon is connected to OH
                        site_type = 'alcohol'
                        electron_change = 2
                        break
            
            # Alkenes (can be oxidized to epoxides)
            elif atom.GetSymbol() == 'C' and atom.GetHybridization() == Chem.HybridizationType.SP2:
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'C' and neighbor.GetHybridization() == Chem.HybridizationType.SP2:
                        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                        if bond.GetBondType() == Chem.BondType.DOUBLE:
                            site_type = 'alkene'
                            electron_change = 2
                            break
            
            # Amines (can be oxidized)
            elif atom.GetSymbol() == 'N' and atom.GetFormalCharge() == 0:
                site_type = 'amine'
                electron_change = 1
            
            # Sulfides (can be oxidized to sulfoxides/sulfones)
            elif atom.GetSymbol() == 'S' and atom.GetFormalCharge() == 0:
                site_type = 'sulfide'
                electron_change = 2
            
            if site_type:
                sites.append({
                    'atom_idx': atom.GetIdx(),
                    'site_type': site_type,
                    'electron_change': electron_change
                })
        
        return sites
    
    def _find_reducible_sites(self, mol: Chem.Mol) -> List[Dict[str, Any]]:
        """Find sites that can be reduced."""
        sites = []
        
        for atom in mol.GetAtoms():
            site_type = None
            electron_change = -1
            
            # Carbonyls (can be reduced to alcohols)
            if atom.GetSymbol() == 'C':
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'O':
                        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                        if bond.GetBondType() == Chem.BondType.DOUBLE:
                            site_type = 'carbonyl'
                            electron_change = -2
                            break
            
            # Nitro groups (can be reduced to amines)
            elif atom.GetSymbol() == 'N' and atom.GetFormalCharge() > 0:
                oxygen_neighbors = 0
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'O':
                        oxygen_neighbors += 1
                
                if oxygen_neighbors >= 2:
                    site_type = 'nitro'
                    electron_change = -6
            
            # Imines (can be reduced to amines)
            elif atom.GetSymbol() == 'N' and atom.GetHybridization() == Chem.HybridizationType.SP2:
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'C':
                        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                        if bond.GetBondType() == Chem.BondType.DOUBLE:
                            site_type = 'imine'
                            electron_change = -2
                            break
            
            # Positive charges (can be reduced)
            elif atom.GetFormalCharge() > 0:
                site_type = 'cation'
                electron_change = -atom.GetFormalCharge()
            
            if site_type:
                sites.append({
                    'atom_idx': atom.GetIdx(),
                    'site_type': site_type,
                    'electron_change': electron_change
                })
        
        return sites
    
    def apply(self, molecules: List[Molecule], environment: Environment) -> List[Reaction]:
        """
        Apply redox operator to generate reactions.
        
        Args:
            molecules: List of molecules (should contain exactly one)
            environment: Reaction environment
            
        Returns:
            List of generated reactions
        """
        if not self.can_apply(molecules, environment):
            return []
        
        reactions = []
        drives = environment.get_active_drives()
        
        # Apply templates based on environment
        for template in self.get_templates():
            # Only apply oxidation templates in oxidizing environments
            if template.redox_type == "oxidation" and not drives.get('oxidation', False):
                continue
            
            # Only apply reduction templates in reducing environments
            if template.redox_type == "reduction" and not drives.get('reduction', False):
                continue
            
            template_reactions = template.apply(molecules)
            reactions.extend(template_reactions)
        
        return reactions
