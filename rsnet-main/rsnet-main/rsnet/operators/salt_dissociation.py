"""
Salt Dissociation Operator
==========================

Handles the heterolytic dissociation of salts and complex anions.
Specifically targets:
1. Complex anion dissociation: PF6- -> PF5 + F-
2. Contact Ion Pair dissociation: LiF -> Li+ + F- (reverse of precipitation)
"""

from typing import List, Dict, Any, Optional
from rdkit import Chem

from ..core.molecule import Molecule
from ..core.environment import Environment
from ..core.reaction import Reaction
from .base import BaseOperator

class SaltDissociationOperator(BaseOperator):
    """
    Operator for salt dissociation/decomposition.
    """
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        super().__init__(
            name="SaltDissociation",
            description="Salt dissociation (PF6- -> PF5 + F-)",
            config=config
        )
        self.arity = 1
        
    def can_apply(self, molecules: List[Molecule], environment: Environment) -> bool:
        """Check if any molecule is a complex anion capable of dissociation."""
        for mol in molecules:
            # Check for anions with dissociable ligands
            if self._is_complex_anion(mol):
                return True
        return False
        
    def apply(self, molecules: List[Molecule], environment: Environment) -> List[Reaction]:
        reactions = []
        for mol in molecules:
            if self._is_complex_anion(mol):
                # Try to dissociate ligands
                new_rxns = self._dissociate_anion(mol)
                reactions.extend(new_rxns)
        return reactions

    def _is_complex_anion(self, mol: Molecule) -> bool:
        """
        Check if molecule is a complex anion (e.g. PF6-, BF4-, AsF6-).
        Criteria:
        1. Net charge < 0
        2. Contains a central atom with hypervariance or high coordination
        3. Connected to electronegative ligands (F, Cl, O, etc.)
        """
        rdkit_mol = mol.rdkit_mol
        if Chem.GetFormalCharge(rdkit_mol) >= 0:
            return False
            
        # Find the negative center
        for atom in rdkit_mol.GetAtoms():
            if atom.GetFormalCharge() < 0:
                # Check neighbors for potential leaving groups
                # E.g. P-F bond
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() in ['F', 'Cl', 'Br', 'I', 'O']:
                        return True
        return False

    def _dissociate_anion(self, mol: Molecule) -> List[Reaction]:
        """
        Generalized dissociation: [M-X]- -> M + X-
        Strategy:
        1. Identify negative center M-
        2. Identify ligand X attached to M-
        3. Break M-X bond heterolytically: electrons go to X, making X-
        """
        reactions = []
        rdkit_mol = mol.rdkit_mol
        
        # Identify bonds to break
        bonds_to_break = []
        
        for atom in rdkit_mol.GetAtoms():
            if atom.GetFormalCharge() < 0: # Negative center (e.g., P in PF6-)
                center_idx = atom.GetIdx()
                
                for neighbor in atom.GetNeighbors():
                    # Check if neighbor is a valid leaving group (halogens, etc.)
                    # and not a carbon backbone (usually)
                    sym = neighbor.GetSymbol()
                    if sym in ['F', 'Cl', 'Br', 'I', 'O', 'S']:
                        ligand_idx = neighbor.GetIdx()
                        bonds_to_break.append((center_idx, ligand_idx, sym))

        # Generate reactions for each unique dissociation path
        # Limit to 1 per distinct ligand type to avoid explosion? 
        # For PF6, all F are equivalent (topologically), RDKit might handle symmetry or we get 6 identical reactions.
        # We can unique-ify by SMILES check later.
        
        processed_smiles = set()
        
        for center_idx, ligand_idx, ligand_sym in bonds_to_break:
            try:
                rw_mol = Chem.RWMol(rdkit_mol)
                rw_mol.RemoveBond(center_idx, ligand_idx)
                
                # Charge Re-distribution for Heterolysis
                # [M-]-X -> M(0) + [X-]
                # Current: M is -1, X is 0. 
                # After break: M becomes 0 (loses electron density to X?), X becomes -1.
                
                # Logic: Heterolytic bond cleavage where the pair goes to the ligand.
                # Since M was -1 and linked to neutral X. M effectively had "extra" electron.
                # If bond breaks and electrons go to X:
                # M (-1) -> loses pairing electron? No.
                # Let's count electrons.
                # PF6-: P is -1 (formal). 6 bonds. Valence 5 + 1 = 6.
                # Dissociate F-: P becomes PF5 (Valence 5, 0 charge). F becomes F- (Valence 7+1=8, -1 charge).
                # So:
                # Center: -1 -> 0
                # Ligand: 0 -> -1
                
                center_atom = rw_mol.GetAtomWithIdx(center_idx)
                ligand_atom = rw_mol.GetAtomWithIdx(ligand_idx)
                
                new_center_charge = center_atom.GetFormalCharge() + 1
                new_ligand_charge = ligand_atom.GetFormalCharge() - 1
                
                center_atom.SetFormalCharge(new_center_charge)
                ligand_atom.SetFormalCharge(new_ligand_charge)
                
                # Sanitize
                try:
                    Chem.SanitizeMol(rw_mol)
                except:
                   # Sometimes partial sanitization is needed for transition states
                   continue
                   
                frags = Chem.GetMolFrags(rw_mol, asMols=True, sanitizeFrags=True)
                if len(frags) != 2:
                    continue
                    
                # Store Product
                products = []
                prod_smiles = []
                for f in frags:
                    s = Chem.MolToSmiles(f)
                    prod_smiles.append(s)
                    # Naming
                    name = s # Placeholder
                    if "P" in s and "F" in s: name = "PF5"
                    if s == "[F-]": name = "F-"
                    if s == "[Li+]": name = "Li+"
                    products.append(Molecule.from_smiles(s, name=name))
                
                # Unique Check
                products_str = sorted(prod_smiles)
                products_hash = "_".join(products_str)
                if products_hash in processed_smiles:
                    continue
                processed_smiles.add(products_hash)
                
                rxn = Reaction(
                    reactants=[mol],
                    products=products,
                    operator_name="SaltDissociation"
                )
                rxn.reaction_energy = 10.0 # Endo
                reactions.append(rxn)
                
            except Exception as e:
                pass
                
        return reactions
