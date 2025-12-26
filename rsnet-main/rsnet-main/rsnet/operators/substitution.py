"""
Nucleophilic substitution operator for RSNet.

This module implements nucleophilic substitution reactions (Sn2) and 
transesterification (addition-elimination) relevant for electrolyte decomposition.
"""

from typing import List, Dict, Any, Optional
from rdkit import Chem
from rdkit.Chem import AllChem

from .base import BaseOperator
from ..core.molecule import Molecule
from ..core.reaction import Reaction, ReactionTemplate
from ..core.environment import Environment

class SubstitutionTemplate(ReactionTemplate):
    """
    Template for substitution reactions using SMARTS transformations.
    """
    
    def __init__(
        self,
        name: str,
        reaction_smarts: str,
        description: str,
        forward_barrier: float = 25.0
    ):
        """
        Initialize substitution template.
        
        Args:
            name: Template name
            reaction_smarts: SMARTS string defining the reaction (reactants>>products)
            description: Description of the reaction
            forward_barrier: Estimated activation energy (kcal/mol)
        """
        super().__init__(
            name=name,
            description=description,
            pattern=reaction_smarts
        )
        self.reaction_smarts = reaction_smarts
        self.forward_barrier = forward_barrier
        self._rxn_obj = AllChem.ReactionFromSmarts(reaction_smarts)

    def apply(self, molecules: List[Molecule]) -> List[Reaction]:
        """
        Apply the SMARTS transformation to valid reactant combinations.
        
        Args:
            molecules: List of potential reactant molecules (expecting 2 for substitution)
            
        Returns:
            List of generated Reaction objects
        """
        if len(molecules) != 2:
            return []
            
        reactions = []
        
        # Try both orderings: [A, B] and [B, A] since SMARTS is order-dependent
        # SMARTS usually defines: nucleophile.substrate >> product1.product2
        
        combinations = [
            (molecules[0], molecules[1]),
            (molecules[1], molecules[0])
        ]
        
        unique_products = set()

        for r1, r2 in combinations:
            # Prepare RDKit reactants
            rd_reactants = (r1.rdkit_mol, r2.rdkit_mol)
            
            try:
                # Run the reaction
                product_sets = self._rxn_obj.RunReactants(rd_reactants)
                
                for product_tuple in product_sets:
                    generated_mols = []
                    valid_product_set = True
                    product_smiles_list = []
                    
                    for prod in product_tuple:
                        try:
                            # Sanitize and create Molecule object
                            Chem.SanitizeMol(prod)
                            # Convert to simplified smiles for uniqueness check
                            smi = Chem.MolToSmiles(prod, isomericSmiles=True)
                            product_smiles_list.append(smi)
                            
                            mol_name = f"product_{smi[:10]}" # Temp name
                            generated_mols.append(Molecule(prod, name=mol_name))
                        except Exception:
                            valid_product_set = False
                            break
                    
                    if not valid_product_set:
                        continue
                        
                    # Tuple of sorted SMILES strings to identify unique product sets
                    prod_sig = tuple(sorted(product_smiles_list))
                    if prod_sig in unique_products:
                        continue
                    unique_products.add(prod_sig)
                    
                    # Create Reaction Object
                    reaction_name = f"{self.name}: {r1.name} + {r2.name}"
                    rxn = Reaction(
                        reactants=[r1, r2],
                        products=generated_mols,
                        name=reaction_name,
                        operator_name="Substitution",
                        activation_energy=self.forward_barrier
                    )
                    reactions.append(rxn)
                    
            except Exception as e:
                # Silent fail for individual reaction attempts (common in library generation)
                continue
                
        return reactions

class SubstitutionOperator(BaseOperator):
    """
    Operator for nucleophilic substitution reactions.
    
    Includes:
    - Sn2 Alkylation (Attack on Carbonate/Alkyl carbon)
    - Transesterification (Attack on Carbonyl carbon)
    """
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        super().__init__(
            name="Substitution",
            description="Generates nucleophilic substitution and transesterification reactions",
            config=config
        )
        self.arity = 2
        self.reaction_type = "substitution"
        self._setup_templates()
        
    def _setup_templates(self):
        """Define SMARTS templates for substitutions."""
        
        # 1. Transesterification (Alkoxide + Carbonate/Ester -> New Carbonate/Ester + Alkoxide)
        # Nucleophile: [O-] (Alkoxide)
        # Leaving Group: [O-] (Alkoxide)
        # Explicitly set charge to 0 for new ester oxygen, and -1 for leaving oxygen.
        transesterification = SubstitutionTemplate(
            name="Transesterification",
            # [O-:1] attacks [C:2], releasing [O:4]
            reaction_smarts="[O&-1:1].[C:2](=[O:3])[O:4][C:5]>>[O&+0:1][C:2](=[O:3]).[O&-1:4][C:5]",
            description="Nucleophilic acyl substitution at carbonyl carbon",
            forward_barrier=15.0
        )
        self.add_template(transesterification)
        
        # 2. Sn2 Alkylation (Alkoxide + Alkyl Carbonate -> Ether + Carbonate Anion)
        # [O-:1] attacks [C:2], releasing [O:3]
        sn2_alkylation = SubstitutionTemplate(
            name="Sn2_Alkylation",
            # Attack on Alkyl Carbon [C:2] attached to Oxygen [O:3]
            reaction_smarts="[O&-1:1].[C;X4:2][O:3][C:4](=[O:5])>>[O&+0:1][C:2].[O&-1:3][C:4](=[O:5])",
            description="Sn2 attack on alkyl group of ester/carbonate",
            forward_barrier=30.0
        )
        self.add_template(sn2_alkylation)
        
        # 3. Hydrolysis (Water + Ester -> Acid + Alcohol)
        # [O:1] (water) attacks, becomes protonated ester initially -> separate step?
        # Simplified: H2O + RCOOR' -> RCOOH + R'OH
        # [O&H2:1].[C:2](=[O:3])[O:4][C:5]>>[O:1][C:2](=[O:3]).[O:4][C:5]
        # This one is trickier with proton transfer. 
        # Let's use Hydroxide [OH-] for basic hydrolysis which is common in batteries.
        basic_hydrolysis = SubstitutionTemplate(
            name="Basic_Hydrolysis",
            reaction_smarts="[O&-1:1].[C:2](=[O:3])[O:4][C:5]>>[O&+0:1][C:2](=[O:3]).[O&-1:4][C:5]",
            description="Hydrolysis of ester/carbonate by hydroxide",
            forward_barrier=20.0
        )
        self.add_template(basic_hydrolysis)

    def can_apply(self, molecules: List[Molecule], environment: Environment) -> bool:
        """Check if substitution/transesterification is favorable."""
        if len(molecules) != 2:
            return False
            
        # Substitution generally requires:
        # 1. Thermal driving force (for barrier crossing)
        # 2. Solution phase (for ion mobility)
        # 3. Presence of nucleophiles/electrophiles (handled by templates)
        
        drives = environment.get_active_drives()
        
        return drives.get('thermal', False) or drives.get('solution_phase', False)
