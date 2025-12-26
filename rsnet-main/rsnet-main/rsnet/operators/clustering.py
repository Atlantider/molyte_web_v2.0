"""
Clustering Operator - Multi-Molecular Aggregation Strategy
===========================================================

This module implements the "Anchor-Based Aggregation" strategy for handling
multi-molecular interactions (e.g., Li+ + 4EC -> Li(EC)4) without O(N^3) complexity.
"""

from typing import List, Dict, Any, Optional, Tuple
from rdkit import Chem

from ..core.molecule import Molecule
from ..core.environment import Environment
from ..core.reaction import Reaction
from .base import BaseOperator
# from ..utils.structure_analysis import get_structure_tags

class ClusteringOperator(BaseOperator):
    """
    Clustering Operator for multi-molecular aggregation.
    
    Arity: N (Variable)
    Strategy: Anchor-Based
    - Identify "Anchors" (Metal ions, Nucleation centers)
    - Identify "Ligands" (Solvent molecules, Monomers)
    - Form stable clusters based on coordination number or solvation limits
    """
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        super().__init__(
            name="Clustering",
            description="Multi-molecular clustering (Solvation/Coordination)",
            config=config
        )
        self.arity = 0  # Special flag for N-body
        self.reaction_type = "clustering"
        self.max_coordination = self.config.get('max_coordination', 4)
        
        # Anchor definitions
        self.anchors = ['Li', 'Na', 'K', 'Mg', 'Ca', 'Al', 'Zn', 'Ni', 'Mn', 'Co', 'Fe']
        
    def can_apply(self, molecules: List[Molecule], environment: Environment) -> bool:
        """
        Check if clustering is applicable.
        For clustering, we usually check the whole pool, but here we can check if 
        at least one molecule is an anchor.
        """
        for mol in molecules:
            if self._is_anchor(mol):
                return True
        return False

    def apply(self, molecules: List[Molecule], environment: Environment) -> List[Reaction]:
        """
        Apply clustering logic.
        
        Args:
            molecules: List of available molecules (Context: The entire SpeciesPool usually)
            environment: Reaction environment
            
        Returns:
            List of reactions (e.g., Li+ + 4EC -> Li(EC)4)
        """
        reactions = []
        
        # 1. Identify Anchors and Ligands
        anchors = []
        ligands = []
        
        for mol in molecules:
            if self._is_anchor(mol):
                anchors.append(mol)
            elif self._is_ligand(mol):
                ligands.append(mol)
                
        if not anchors or not ligands:
            return []
            
        # 2. Form Clusters
        # Simple strategy: Greedy allocation
        # For each anchor, try to find enough ligands
        
        used_ligands = set()
        
        for anchor in anchors:
            # Determine coordination number
            target_cn = self._get_target_cn(anchor)
            
            # Find available ligands (ignore used ones for strict stoichiometry, 
            # or allow reuse for "representative" reaction generation)
            # In network generation, we want to generate the *possibility* of reaction.
            
            # Form pairs first to sort by affinity? 
            # Simplified: Just grab the most abundant/compatible ligands
            
            # We take up to `target_cn` ligands
            current_ligands = []
            for l in ligands:
                if len(current_ligands) >= target_cn:
                    break
                current_ligands.append(l)
            
            if len(current_ligands) > 0:
                # Generate Reaction: Anchor + n*Ligand -> Complex
                # Check reaction uniqueness is handled by Evolver
                
                # Check intermediate steps?
                # Li+ + EC -> Li(EC)+
                # Li(EC)+ + EC -> Li(EC)2+...
                # We can generate the "Full Shell" directly or sequential
                
                # Option 1: Full Shell (Direct)
                if len(current_ligands) == target_cn:
                    complex_mol = self._create_complex(anchor, current_ligands)
                    if complex_mol:
                        rxn = Reaction(
                            reactants=[anchor] + current_ligands,
                            products=[complex_mol],
                            operator_name="solvation_shell_formation"
                        )
                        rxn.reaction_energy = -15.0 * len(current_ligands) # Mock energy
                        reactions.append(rxn)

                # Option 2: Stepwise (1 by 1) - Often better for detailed network
                # But user asked for "Multi-molecular", preventing N^3.
                # So we support the "Direct N-body" assembly.
        
        return reactions

    def _is_anchor(self, mol: Molecule) -> bool:
        """Check if molecule is a coordination center."""
        # Check atomic symbols
        rdkit_mol = mol.rdkit_mol
        if rdkit_mol.GetNumAtoms() == 1:
            sym = rdkit_mol.GetAtomWithIdx(0).GetSymbol()
            if sym in self.anchors:
                return True
                
        # Also check for existing complexes that aren't saturated
        # (Identifying "unsaturated complex" is harder, skipping for now)
        return False

    def _is_ligand(self, mol: Molecule) -> bool:
        """Check if molecule can act as a ligand."""
        # Has heteroatoms?
        from ..utils.structure_analysis import get_structure_tags
        tags = get_structure_tags(mol)
        has_hetero = any(atom.GetSymbol() in ['O', 'N', 'S', 'P', 'F'] for atom in mol.rdkit_mol.GetAtoms())
        
        # Check formal charge: Ligands should not be positively charged (avoids Cation-Cation clustering)
        charge = Chem.GetFormalCharge(mol.rdkit_mol)
        if charge > 0:
            return False
            
        return has_hetero and not self._is_anchor(mol)

    def _get_target_cn(self, anchor: Molecule) -> int:
        """Get target coordination number."""
        sym = anchor.rdkit_mol.GetAtomWithIdx(0).GetSymbol()
        if sym == 'Li': return 4
        if sym in ['Mg', 'Ca', 'Al']: return 6
        return 4
        
    def _create_complex(self, anchor: Molecule, ligands: List[Molecule]) -> Optional[Molecule]:
        """Create the complex molecule object (merged)."""
        try:
            # Combine all
            combo = anchor.rdkit_mol
            for lig in ligands:
                combo = Chem.CombineMols(combo, lig.rdkit_mol)
            
            # Note: We don't draw explicit bonds to metal usually in force fields,
            # but for connectivity in RDKit we might want to?
            # For now, just a disconnected salt/complex object is fine, 
            # Or we add dummy bonds.
            
            # Let's add dummy zero-order bonds for visualization if possible,
            # or just leave them disconnected but in one Mol object.
            # Disconnected is safer for valence checks.
            
            # Name
            ligand_name = ligands[0].name.split('_')[0]
            name = f"{anchor.name}({ligand_name}){len(ligands)}"
            
            smiles = Chem.MolToSmiles(combo)
            return Molecule.from_smiles(smiles, name=name)
        except:
            return None
