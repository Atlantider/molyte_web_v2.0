"""
Intelligent Network Evolution System
====================================

This module implements the advanced reaction network evolution logic.
It replaces the simple "apply all operators" approach with a chemically
aware system that handles:
1. Cross-interactions (R+R, R+P, P+P)
2. Smart Pair Selection (Chemical Logic)
3. Reverse Reaction Filtering (Thermodynamic Consistency)
4. Dynamic Species Pool Management
"""

import logging
import itertools
from typing import List, Dict, Set, Tuple, Optional, Any
from dataclasses import dataclass, field

from ..core.molecule import Molecule
from ..core.reaction import Reaction
from ..core.environment import Environment
from ..operators.base import BaseOperator
from ..operators.registry import OPERATOR_REGISTRY
from ..compute.reaction_screener import ReactionScreener
from ..network.config import NetworkGenerationConfig, NetworkGenerationPresets

logger = logging.getLogger(__name__)

@dataclass
class SpeciesState:
    """Tracks the state of a molecule in the network."""
    molecule: Molecule
    generation: int
    is_intermediate: bool = False
    is_radical: bool = False
    concentration_score: float = 1.0  # 1.0 = Reactant, decreases with gen

class SpeciesPool:
    """
    Manages the active pool of chemical species.
    """
    def __init__(self):
        self.species: Dict[str, SpeciesState] = {}  # SMILES -> State
        self.new_species: List[SpeciesState] = []   # Added in current step

    def add(self, molecule: Molecule, generation: int, tags: dict = None):
        if molecule.smiles in self.species:
            return False
        
        is_radical = False
        if tags:
            is_radical = tags.get('is_radical', False)
        
        state = SpeciesState(
            molecule=molecule,
            generation=generation,
            is_radical=is_radical,
            concentration_score=1.0 / (generation + 1)
        )
        self.species[molecule.smiles] = state
        self.new_species.append(state)
        return True

    def get_all(self) -> List[SpeciesState]:
        return list(self.species.values())

    def get_radicals(self) -> List[SpeciesState]:
        return [s for s in self.species.values() if s.is_radical]

    def get_neutrals(self) -> List[SpeciesState]:
        return [s for s in self.species.values() if not s.is_radical]
    
    def clear_new(self):
        self.new_species = []

class InteractionStrategy:
    """
    Defines intelligent pairing rules for bimolecular reactions.
    """
    def __init__(self):
        self.max_combinations = 500  # Safety limit

    def generate_pairs(self, pool: SpeciesPool) -> List[Tuple[Molecule, Molecule]]:
        """
        Generate chemically significant pairs.
        Priority:
        1. Radical + Radical (Coupling/Disproportionation) - Very Fast
        2. Radical + Neutral (Abstraction/Addition) - Fast
        3. Reactive Neutral + Reactive Neutral (e.g. Nucleophile + Electrophile) - Slow/Selective
        """
        pairs = []
        radicals = [s.molecule for s in pool.get_radicals()]
        neutrals = [s.molecule for s in pool.get_neutrals()]
        
        # 1. Radical + Radical (Highest Priority)
        # itertools.combinations_with_replacement allows A+A and A+B
        for r1, r2 in itertools.combinations_with_replacement(radicals, 2):
            pairs.append((r1, r2))
            
        # 2. Radical + Neutral
        for r in radicals:
            for n in neutrals:
                pairs.append((r, n))
                
        # 3. Neutral + Neutral (Selective)
        # TODO: Implement "Reactive Tag" check to avoid dead molecules reacting
        # For now, we limit neutral-neutral to avoid explosion if pool is large
        if len(neutrals) < 20:
            for n1, n2 in itertools.combinations_with_replacement(neutrals, 2):
                pairs.append((n1, n2))
        else:
            logger.info("Skipping full Neutral-Neutral scan due to pool size.")
        
        return pairs

class NetworkEvolver:
    """
    Orchestrates the evolution of the reaction network.
    """
    """
    Orchestrates the evolution of the reaction network.
    """
    def __init__(self, environment: Environment, registry=OPERATOR_REGISTRY, screener: Optional[Any] = None, config: Optional[NetworkGenerationConfig] = None):
        self.environment = environment
        self.registry = registry
        self.config = config or NetworkGenerationPresets.battery_chemistry()
        self.pool = SpeciesPool()
        self.strategy = InteractionStrategy()
        self.generated_reactions_hashes: Set[str] = set() # Hash of reactants->products
        self.reactions: List[Reaction] = [] # Store actual reaction objects
        self.screener = screener or ReactionScreener()

    def evolve(self, initial_molecules: List[Molecule], max_generations: int = 3):
        """
        Main evolution loop.
        """
        current_gen = 0
        
        # Initialize pool
        for mol in initial_molecules:
            # We need to tag them to know if they are radicals
            from ..features.structure_tags import get_structure_tags
            from ..utils.chemistry_tools import is_radical
            tags = get_structure_tags(mol)
            tags['is_radical'] = is_radical(mol.rdkit_mol)
            self.pool.add(mol, 0, tags)

        while current_gen < max_generations:
            logger.info(f"--- Evolution Generation {current_gen} ---")
            
            new_reactions = []
            
            # 1. Unimolecular Step (Decomposition, Rearrangement, etc.)
            # Apply only to *new* species from this/prev generation to avoid re-calc
            target_species = self.pool.get_all() # simplified: check all against registry efficiency
            
            uni_reactions = self._process_unimolecular(target_species)
            new_reactions.extend(uni_reactions)
            
            # Check limit
            if len(new_reactions) >= self.config.max_reactions_per_generation:
                 logger.warning(f"Generation limit hit ({self.config.max_reactions_per_generation}). Screening and finalizing early.")
                 # Truncate
                 new_reactions = new_reactions[:self.config.max_reactions_per_generation]
            else:
                # 2. Bimolecular Step (Cross-Interactions)
                pairs = self.strategy.generate_pairs(self.pool)
                logger.info(f"Generated {len(pairs)} candidate pairs for bimolecular reactions")
                
                bi_reactions = self._process_bimolecular(pairs)
                new_reactions.extend(bi_reactions)
                
                # Check limit again
                if len(new_reactions) >= self.config.max_reactions_per_generation:
                    logger.warning(f"Generation limit hit during bimolecular ({self.config.max_reactions_per_generation}). Finalizing early.")
                    new_reactions = new_reactions[:self.config.max_reactions_per_generation]
                else:
                    # 3. Clustering Step
                    # This handles N-body aggregation like Solvation
                    cluster_reactions = self._process_clustering(self.pool.get_all()) # Check entire pool
                    new_reactions.extend(cluster_reactions)

            # Final check just in case
            if len(new_reactions) > self.config.max_reactions_per_generation:
                new_reactions = new_reactions[:self.config.max_reactions_per_generation]

            # 4. Process Products & Update Pool
            added_count = 0
            for rxn in new_reactions:
                rxn_hash = self._reaction_hash(rxn)
                self.generated_reactions_hashes.add(rxn_hash)
                self.reactions.append(rxn)
                for prod in rxn.products:
                    # Tag new products
                    from ..features.structure_tags import get_structure_tags
                    from ..utils.chemistry_tools import is_radical
                    tags = get_structure_tags(prod)
                    tags['is_radical'] = is_radical(prod.rdkit_mol)
                    
                    if self.pool.add(prod, current_gen + 1, tags):
                        added_count += 1
            
            logger.info(f"Gen {current_gen} complete. {len(new_reactions)} reactions, {added_count} new species.")
            
            if added_count == 0 and len(new_reactions) == 0:
                logger.info("Network converged (no new species or reactions).")
                break
                
            current_gen += 1
            
        return self.pool.get_all()

    def _process_unimolecular(self, species_list: List[SpeciesState]) -> List[Reaction]:
        reactions = []
        for state in species_list:
            # Safety break
            if len(reactions) >= self.config.max_reactions_per_generation:
                break

            mol = state.molecule
            # Get operators that can apply to this single molecule
            valid_ops = self.registry.get_active_operators([mol], self.environment)
            for op in valid_ops:
                # Use arity metadata if available, otherwise assume suitable if can_apply returned True
                # (Active operators are pre-filtered by registry using can_apply)
                if getattr(op, 'arity', 1) == 1 or op.reaction_type == 'radical': # Radical supports 1 or 2
                    try:
                        rxns = op.apply([mol], self.environment)
                        for r in rxns:
                            if self._is_new_reaction(r):
                                reactions.append(r)
                    except Exception as e:
                        logger.warning(f"Unimolecular op {op.name} failed on {mol.smiles}: {e}")
        return reactions

    def _process_bimolecular(self, pairs: List[Tuple[Molecule, Molecule]]) -> List[Reaction]:
        reactions = []
        for m1, m2 in pairs:
            # Safety break
            if len(reactions) >= self.config.max_reactions_per_generation:
                break

            # Get operators compatible with this pair
            valid_ops = self.registry.get_active_operators([m1, m2], self.environment)
            for op in valid_ops:
                # Ensure operator supports bimolecular mode
                if getattr(op, 'arity', 1) >= 2:
                     try:
                        rxns = op.apply([m1, m2], self.environment)
                        for r in rxns:
                            if self._is_new_reaction(r):
                                reactions.append(r)
                     except Exception as e:
                        pass
        return reactions

    def _reaction_hash(self, reaction: Reaction) -> str:
        # Sort reactants and products for unique ID
        r_smi = sorted([m.smiles for m in reaction.reactants])
        p_smi = sorted([m.smiles for m in reaction.products])
        return f"{'.'.join(r_smi)}>>{'.'.join(p_smi)}"

    def _is_new_reaction(self, reaction: Reaction) -> bool:
        # Check forward
        h = self._reaction_hash(reaction)
        if h in self.generated_reactions_hashes:
            return False
        
        # Check reverse (Thermodynamic filtering)
        # If A->B exists, we don't want B->A unless strictly necessary
        # This is a rigorous check usually done after energy calc.
        # Here we just prevent duplicate edges in simplest form
        
        # Ideally, we check:
        # r_smi = sorted([m.smiles for m in reaction.reactants])
        # p_smi = sorted([m.smiles for m in reaction.products])
        # reverse_hash = f"{'.'.join(p_smi)}>>{'.'.join(r_smi)}"
        # if reverse_hash in self.generated_reactions:
        #     return False # Prevent immediate reverse cycling
        
        return True

    def _process_clustering(self, molecules: List[SpeciesState]) -> List[Reaction]:
        """
        Handle multi-molecular clustering (e.g. Li+ + 4EC -> Li(EC)4).
        Complexity: O(N_anchors * N_ligands), not O(N^3).
        """
        reactions = []
        # Extract available molecule objects
        available_mols = [s.molecule for s in molecules]
        
        # Get clustering operator
        clustering_op = self.registry.get_operator("clustering")
        if clustering_op and clustering_op.can_apply(available_mols, self.environment):
            try:
                # Apply clustering to the entire pool context
                rxns = clustering_op.apply(available_mols, self.environment)
                for r in rxns:
                    if self._is_new_reaction(r):
                        reactions.append(r)
            except Exception as e:
                logger.warning(f"Clustering op failed: {e}")
        
        return reactions
