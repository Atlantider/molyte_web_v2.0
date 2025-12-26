"""
Reaction network generator for RSNet.

This module implements the core network generation algorithm that:
1. Takes initial molecules and environment
2. Applies operators to generate reactions
3. Screens reactions using xTB calculations
4. Grows the network iteratively
5. Applies pruning strategies
"""

import time
import logging
from typing import List, Dict, Set, Optional, Tuple, Any
from collections import deque
import networkx as nx

from ..core.molecule import Molecule
from ..core.environment import Environment
from ..core.reaction import Reaction
from ..operators.base import BaseOperator
# from ..operators.hydrogen_transfer import HydrogenTransferOperator  <-- Removed to avoid circular import
# from ..operators.bond_breaking import BondBreakingOperator      <-- Removed to avoid circular import
from ..operators.registry import OPERATOR_REGISTRY, OperatorRegistry
from ..features.driving_forces import get_driving_forces
from ..features.structure_tags import get_structure_tags
from ..compute.reaction_screener import ReactionScreener
from .config import NetworkGenerationConfig, NetworkGenerationStrategy, OperatorSelectionStrategy


logger = logging.getLogger(__name__)


class ReactionNetwork:
    """
    Represents a reaction network as a directed graph.
    
    Nodes: Molecules
    Edges: Reactions
    """
    
    def __init__(self):
        """Initialize an empty reaction network."""
        self.graph = nx.DiGraph()
        self.molecules = {}  # SMILES -> Molecule mapping
        self.reactions = {}  # reaction_id -> Reaction mapping
        self.generation_info = {}  # molecule_id -> generation number
        
    def add_molecule(self, molecule: Molecule, generation: int = 0) -> None:
        """Add a molecule to the network."""
        if molecule.smiles not in self.molecules:
            self.molecules[molecule.smiles] = molecule
            self.graph.add_node(molecule.smiles, molecule=molecule)
            self.generation_info[molecule.smiles] = generation
    
    def add_reaction(self, reaction: Reaction) -> bool:
        """
        Add a reaction to the network.
        
        Returns:
            True if reaction was added, False if it already exists
        """
        # Create reaction signature for deduplication
        reactant_smiles = tuple(sorted([r.smiles for r in reaction.reactants]))
        product_smiles = tuple(sorted([p.smiles for p in reaction.products]))
        reaction_signature = (reactant_smiles, product_smiles)
        
        # Check if reaction already exists
        for existing_reaction in self.reactions.values():
            existing_reactants = tuple(sorted([r.smiles for r in existing_reaction.reactants]))
            existing_products = tuple(sorted([p.smiles for p in existing_reaction.products]))
            if (existing_reactants, existing_products) == reaction_signature:
                return False
        
        # Add reaction
        self.reactions[reaction.id] = reaction
        
        # Add molecules if not present
        for mol in reaction.reactants + reaction.products:
            if mol.smiles not in self.molecules:
                self.add_molecule(mol)
        
        # Add edges
        for reactant in reaction.reactants:
            for product in reaction.products:
                self.graph.add_edge(
                    reactant.smiles, 
                    product.smiles,
                    reaction=reaction,
                    weight=reaction.reaction_energy or 0.0
                )
        
        return True
    
    def get_molecules_by_generation(self, generation: int) -> List[Molecule]:
        """Get all molecules from a specific generation."""
        return [
            mol for smiles, mol in self.molecules.items()
            if self.generation_info.get(smiles, 0) == generation
        ]
    
    def get_statistics(self) -> Dict[str, Any]:
        """Get network statistics."""
        return {
            'num_molecules': len(self.molecules),
            'num_reactions': len(self.reactions),
            'num_edges': self.graph.number_of_edges(),
            'max_generation': max(self.generation_info.values()) if self.generation_info else 0,
            'molecules_by_generation': {
                gen: len([1 for g in self.generation_info.values() if g == gen])
                for gen in set(self.generation_info.values())
            }
        }


class NetworkGenerator:
    """
    Main network generation engine.
    
    Implements the growth-based reaction network generation algorithm:
    1. Start with seed molecules
    2. Apply operators to generate reactions
    3. Screen reactions for feasibility
    4. Add feasible reactions to network
    5. Repeat with new molecules until convergence
    """
    
    def __init__(
        self,
        operators: Optional[List[BaseOperator]] = None,
        operator_registry: Optional[OperatorRegistry] = None,
        screener: Optional[ReactionScreener] = None,
        config: Optional[NetworkGenerationConfig] = None,
        legacy_config: Optional[Dict[str, Any]] = None
    ):
        """
        Initialize the network generator.

        Args:
            operators: List of reaction operators to use (deprecated, use operator_registry)
            operator_registry: Operator registry for intelligent operator selection
            screener: Reaction screener for filtering
            config: Network generation configuration object
            legacy_config: Legacy configuration dictionary (deprecated)
        """
        # Handle configuration
        if config is not None:
            self.config = config
        elif legacy_config is not None:
            # Convert legacy config to new format
            self.config = NetworkGenerationConfig()
            for key, value in legacy_config.items():
                if hasattr(self.config, key):
                    setattr(self.config, key, value)
        else:
            self.config = NetworkGenerationConfig()

        # Use operator registry for intelligent selection
        if operator_registry is not None:
            self.operator_registry = operator_registry
            self.use_intelligent_selection = True
            self.operators = None  # Will be selected dynamically
        elif operators is not None:
            # Legacy mode: use provided operators
            self.operators = operators
            self.operator_registry = None
            self.use_intelligent_selection = False
        else:
            # Default: use global registry
            self.operator_registry = OPERATOR_REGISTRY
            self.use_intelligent_selection = True
            self.operators = None

        # Override intelligent selection based on config
        if self.config.generation_strategy == NetworkGenerationStrategy.EXHAUSTIVE:
            self.use_intelligent_selection = False
            if self.operators is None:
                # Get all operators from registry
                self.operators = self.operator_registry.get_all_operators()
        
        # Default screener
        if screener is None:
            screener_config = {
                'optimize_geometries': False,
                'max_workers': self.config.max_screening_workers if self.config.parallel_screening else 1
            }
            screener = ReactionScreener(**screener_config)
        self.screener = screener

        # Statistics
        self.generation_stats = []
        
    def generate_network(
        self,
        seed_molecules: List[Molecule],
        environment: Environment,
        max_time: Optional[float] = None
    ) -> ReactionNetwork:
        """
        Generate a reaction network from seed molecules.
        
        Args:
            seed_molecules: Initial molecules to start from
            environment: Reaction environment
            max_time: Maximum time in seconds (None for no limit)
            
        Returns:
            Generated reaction network
        """
        logger.info(f"Starting network generation with {len(seed_molecules)} seed molecules")
        start_time = time.time()
        
        # Initialize network
        network = ReactionNetwork()
        
        # Add seed molecules
        for mol in seed_molecules:
            network.add_molecule(mol, generation=0)
        
        # Generation queue - molecules to process in next generation
        current_generation = 0
        molecules_to_process = deque(seed_molecules)
        
        while (current_generation < self.config.max_generations and
               len(network.molecules) < self.config.max_species and
               molecules_to_process):
            
            # Check time limit
            if max_time and (time.time() - start_time) > max_time:
                logger.warning(f"Time limit reached after {time.time() - start_time:.1f}s")
                break
            
            logger.info(f"Processing generation {current_generation} with {len(molecules_to_process)} molecules")
            
            # Process current generation
            generation_start = time.time()
            new_reactions = self._process_generation(
                list(molecules_to_process), 
                environment, 
                network,
                current_generation
            )
            generation_time = time.time() - generation_start
            
            # Statistics
            gen_stats = {
                'generation': current_generation,
                'input_molecules': len(molecules_to_process),
                'new_reactions': len(new_reactions),
                'total_molecules': len(network.molecules),
                'total_reactions': len(network.reactions),
                'time': generation_time
            }
            self.generation_stats.append(gen_stats)
            
            logger.info(f"Generation {current_generation}: {len(new_reactions)} new reactions, "
                       f"{len(network.molecules)} total molecules, {generation_time:.2f}s")
            
            # Check convergence
            if len(new_reactions) < self.config.min_new_reactions_per_generation:
                logger.info(f"Convergence reached: only {len(new_reactions)} new reactions")
                break
            
            # Prepare next generation
            molecules_to_process.clear()
            next_generation = current_generation + 1
            
            # Add new molecules from this generation's reactions
            for reaction in new_reactions:
                for product in reaction.products:
                    if product.smiles not in network.molecules:
                        network.add_molecule(product, generation=next_generation)
                        molecules_to_process.append(product)
            
            current_generation = next_generation
        
        total_time = time.time() - start_time
        logger.info(f"Network generation completed in {total_time:.2f}s")
        logger.info(f"Final network: {len(network.molecules)} molecules, {len(network.reactions)} reactions")
        
        return network
    
    def _process_generation(
        self,
        molecules: List[Molecule],
        environment: Environment,
        network: ReactionNetwork,
        generation: int
    ) -> List[Reaction]:
        """
        Process a single generation of molecules.

        Args:
            molecules: Molecules to process
            environment: Reaction environment
            network: Current network state
            generation: Current generation number

        Returns:
            List of new feasible reactions
        """
        all_reactions = []

        if self.use_intelligent_selection:
            # Use intelligent operator selection
            all_reactions = self._generate_reactions_intelligent(molecules, environment, generation)
        else:
            # Legacy mode: apply all operators
            all_reactions = self._generate_reactions_legacy(molecules, environment)

        logger.debug(f"Generated {len(all_reactions)} raw reactions")
        
        # Screen reactions for feasibility
        if all_reactions:
            # Enhanced screening with driving force evaluation
            screening_results = self._screen_reactions_enhanced(all_reactions, environment, generation)

            # Filter feasible reactions
            feasible_reactions = []
            for result in screening_results:
                if (result['success'] and
                    result['feasible'] and
                    abs(result['reaction_energy']) <= self.config.energy_cutoff):

                    reaction = result['reaction']
                    reaction.reaction_energy = result['reaction_energy']
                    reaction.activation_energy = result['activation_energy']
                    reaction.driving_force_score = result.get('driving_force_score', 0.0)
                    feasible_reactions.append(reaction)

            logger.debug(f"Screened to {len(feasible_reactions)} feasible reactions")

            # Sort reactions by combined score (energy + driving forces)
            if self.use_intelligent_selection:
                feasible_reactions = self._rank_reactions_by_relevance(feasible_reactions, environment)

            # Add feasible reactions to network
            new_reactions = []
            for reaction in feasible_reactions:
                if network.add_reaction(reaction):
                    new_reactions.append(reaction)

            return new_reactions
        
        return []

    def _generate_reactions_intelligent(
        self,
        molecules: List[Molecule],
        environment: Environment,
        generation: int
    ) -> List[Reaction]:
        """
        Generate reactions using intelligent operator selection.

        Args:
            molecules: Molecules to process
            environment: Reaction environment
            generation: Current generation number

        Returns:
            List of generated reactions
        """
        all_reactions = []

        # Group molecules by similar characteristics for efficient processing
        molecule_groups = self._group_molecules_by_features(molecules)

        for group_molecules in molecule_groups:
            # Get driving forces for this group
            driving_forces = get_driving_forces(group_molecules, environment)
            strong_forces = {k: v for k, v in driving_forces.items() if v > self.config.driving_force_threshold}

            if not strong_forces:
                logger.debug(f"No strong driving forces for {len(group_molecules)} molecules, skipping intelligent selection")
                continue

            # Select operators based on strategy
            if self.config.operator_selection_strategy == OperatorSelectionStrategy.ADAPTIVE:
                operators = self._select_operators_adaptive(group_molecules, environment, strong_forces)
            elif self.config.operator_selection_strategy == OperatorSelectionStrategy.TOP_N:
                operators = self.operator_registry.get_recommended_operators(
                    group_molecules, environment, max_operators=self.config.max_operators_per_generation
                )
            elif self.config.operator_selection_strategy == OperatorSelectionStrategy.THRESHOLD:
                all_operators = self.operator_registry.get_active_operators(group_molecules, environment)
                operators = [op for op in all_operators
                           if self.operator_registry.get_operator_priority(op, group_molecules, environment)
                           >= self.config.operator_priority_threshold]
            else:  # ALL
                operators = self.operator_registry.get_active_operators(group_molecules, environment)

            logger.debug(f"Selected {len(operators)} operators for {len(group_molecules)} molecules: "
                        f"{[op.name for op in operators]}")

            # Apply selected operators
            for operator in operators:
                for molecule in group_molecules:
                     # Safety check
                    if len(all_reactions) >= self.config.max_reactions_per_generation:
                        logger.warning(f"Reaction limit reached ({self.config.max_reactions_per_generation}). Stopping generation early.")
                        return all_reactions

                    if operator.can_apply([molecule], environment):
                        try:
                            reactions = operator.apply([molecule], environment)
                            all_reactions.extend(reactions)
                        except Exception as e:
                            logger.warning(f"Operator {operator.name} failed on {molecule.smiles}: {e}")

        return all_reactions

    def _generate_reactions_legacy(
        self,
        molecules: List[Molecule],
        environment: Environment
    ) -> List[Reaction]:
        """
        Generate reactions using legacy operator application.

        Args:
            molecules: Molecules to process
            environment: Reaction environment

        Returns:
            List of generated reactions
        """
        all_reactions = []

        # Apply each operator to each molecule
        for operator in self.operators:
            for molecule in molecules:
                # Safety check
                if len(all_reactions) >= self.config.max_reactions_per_generation:
                    logger.warning(f"Reaction limit reached ({self.config.max_reactions_per_generation}). Stopping generation early.")
                    return all_reactions

                if hasattr(operator, 'can_apply'):
                    can_apply = operator.can_apply([molecule], environment)
                else:
                    can_apply = operator.is_applicable([molecule])

                if can_apply:
                    try:
                        reactions = operator.apply([molecule], environment)
                        all_reactions.extend(reactions)
                    except Exception as e:
                        logger.warning(f"Operator {operator.name} failed on {molecule.smiles}: {e}")

        return all_reactions

    def _group_molecules_by_features(self, molecules: List[Molecule]) -> List[List[Molecule]]:
        """
        Group molecules by similar structural features for efficient processing.

        Args:
            molecules: List of molecules to group

        Returns:
            List of molecule groups
        """
        if not self.config.use_structure_based_filtering or len(molecules) <= 5:
            # For small sets or when disabled, treat all as one group
            return [molecules]

        # Group by basic structural features
        groups = {}
        for mol in molecules:
            try:
                tags = get_structure_tags([mol])
                # Create a simple grouping key based on key features
                key = (
                    tags.get('has_aromatic_rings', False),
                    tags.get('has_small_rings', False),
                    len(tags.get('functional_groups', [])),
                    tags.get('molecular_weight', 0) // 50  # Group by weight ranges
                )
                if key not in groups:
                    groups[key] = []
                groups[key].append(mol)
            except Exception as e:
                logger.warning(f"Failed to get structure tags for {mol.smiles}: {e}")
                # Put in default group
                if 'default' not in groups:
                    groups['default'] = []
                groups['default'].append(mol)

        return list(groups.values())

    def _select_operators_adaptive(
        self,
        molecules: List[Molecule],
        environment: Environment,
        strong_forces: Dict[str, float]
    ) -> List[BaseOperator]:
        """
        Adaptively select operators based on driving forces and molecular features.

        Args:
            molecules: Molecules to process
            environment: Reaction environment
            strong_forces: Strong driving forces (> threshold)

        Returns:
            List of selected operators
        """
        # Get all active operators
        active_operators = self.operator_registry.get_active_operators(molecules, environment)

        if not active_operators:
            return []

        # Score operators based on driving forces and molecular compatibility
        operator_scores = {}
        for operator in active_operators:
            score = self._calculate_operator_score(operator, molecules, environment, strong_forces)
            operator_scores[operator] = score

        # Sort by score and select top operators
        sorted_operators = sorted(operator_scores.items(), key=lambda x: x[1], reverse=True)
        selected_operators = [op for op, score in sorted_operators[:self.config.max_operators_per_generation]]

        return selected_operators

    def _calculate_operator_score(
        self,
        operator: BaseOperator,
        molecules: List[Molecule],
        environment: Environment,
        strong_forces: Dict[str, float]
    ) -> float:
        """
        Calculate a score for an operator based on its relevance to current conditions.

        Args:
            operator: Operator to score
            molecules: Current molecules
            environment: Reaction environment
            strong_forces: Strong driving forces

        Returns:
            Operator score (higher is better)
        """
        base_score = 1.0

        # Bonus for operators that match strong driving forces
        force_bonus = 0.0
        if 'thermal' in strong_forces and operator.name.lower() in ['bond_breaking', 'rearrangement']:
            force_bonus += strong_forces['thermal'] * 0.5
        if 'electrochemical' in strong_forces and operator.name.lower() in ['redox', 'addition']:
            force_bonus += strong_forces['electrochemical'] * 0.5
        if 'oxidation' in strong_forces and operator.name.lower() == 'redox':
            force_bonus += strong_forces['oxidation'] * 0.7

        # Bonus for operators that can apply to many molecules
        applicability_bonus = 0.0
        applicable_count = 0
        for mol in molecules:
            if operator.can_apply([mol], environment):
                applicable_count += 1
        applicability_bonus = (applicable_count / len(molecules)) * 0.3

        return base_score + force_bonus + applicability_bonus

    def _screen_reactions_enhanced(
        self,
        reactions: List[Reaction],
        environment: Environment,
        generation: int
    ) -> List[Dict[str, Any]]:
        """
        Enhanced reaction screening with driving force evaluation.

        Args:
            reactions: Reactions to screen
            environment: Reaction environment
            generation: Current generation number

        Returns:
            List of screening results with driving force scores
        """
        # Standard xTB screening
        screening_results = self.screener.screen_reactions(reactions, environment)

        # Add driving force evaluation
        for result in screening_results:
            if result['success']:
                reaction = result['reaction']

                # Evaluate driving forces for reactants and products
                reactant_forces = get_driving_forces(reaction.reactants, environment)
                product_forces = get_driving_forces(reaction.products, environment)

                # Calculate driving force score
                driving_force_score = self._calculate_driving_force_score(
                    reactant_forces, product_forces, environment
                )
                result['driving_force_score'] = driving_force_score

                # Adjust feasibility based on driving forces
                if driving_force_score > 0.7:
                    # Strong driving forces can make marginally unfeasible reactions feasible
                    if not result['feasible'] and abs(result.get('reaction_energy', 0)) < self.config.energy_cutoff * 1.2:
                        result['feasible'] = True
                        logger.debug(f"Reaction made feasible by strong driving forces: {driving_force_score:.2f}")

        return screening_results

    def _calculate_driving_force_score(
        self,
        reactant_forces: Dict[str, float],
        product_forces: Dict[str, float],
        environment: Environment
    ) -> float:
        """
        Calculate overall driving force score for a reaction.

        Args:
            reactant_forces: Driving forces for reactants
            product_forces: Driving forces for products
            environment: Reaction environment

        Returns:
            Combined driving force score (0.0-1.0)
        """
        # Weight different driving forces
        force_weights = {
            'thermal': 0.3,
            'electrochemical': 0.4,
            'oxidation': 0.3,
            'reduction': 0.3,
            'li_coordination': 0.2,
            'surface_reaction': 0.2,
            'ring_strain': 0.4,
            'polar_bond_cleavage': 0.3,
            'weak_bond_breaking': 0.5
        }

        total_score = 0.0
        total_weight = 0.0

        # Consider forces that favor the reaction
        for force_name, weight in force_weights.items():
            reactant_force = reactant_forces.get(force_name, 0.0)
            product_force = product_forces.get(force_name, 0.0)

            # Driving force is the difference (products should be more stable or reactive)
            force_contribution = max(0, reactant_force - product_force * 0.5)  # Favor reactant instability

            total_score += force_contribution * weight
            total_weight += weight

        # Normalize score
        if total_weight > 0:
            return min(1.0, total_score / total_weight)
        return 0.0

    def _rank_reactions_by_relevance(
        self,
        reactions: List[Reaction],
        environment: Environment
    ) -> List[Reaction]:
        """
        Rank reactions by combined energy and driving force relevance.

        Args:
            reactions: Reactions to rank
            environment: Reaction environment

        Returns:
            Sorted list of reactions (best first)
        """
        def reaction_score(reaction):
            # Energy component (lower energy is better)
            energy_score = max(0, 1.0 - abs(reaction.reaction_energy or 0) / self.config.energy_cutoff)

            # Driving force component
            driving_force_score = getattr(reaction, 'driving_force_score', 0.0)

            # Combined score (weighted average)
            return 0.6 * energy_score + 0.4 * driving_force_score

        return sorted(reactions, key=reaction_score, reverse=True)

    def get_generation_statistics(self) -> List[Dict[str, Any]]:
        """Get statistics for each generation."""
        return self.generation_stats.copy()

    def reset_statistics(self) -> None:
        """Reset generation statistics."""
        self.generation_stats.clear()
