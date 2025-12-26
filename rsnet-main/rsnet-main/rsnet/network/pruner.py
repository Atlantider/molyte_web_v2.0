"""
Network pruning strategies for RSNet.

This module implements various pruning strategies to:
1. Remove high-energy reactions
2. Filter out unstable intermediates
3. Remove dead-end pathways
4. Simplify complex networks
"""

import logging
from typing import List, Dict, Set, Optional, Tuple, Any, Callable
import networkx as nx
from collections import defaultdict

from ..core.molecule import Molecule
from ..core.reaction import Reaction
from .generator import ReactionNetwork


logger = logging.getLogger(__name__)


class NetworkPruner:
    """
    Pruning strategies for reaction networks.
    
    Implements various pruning methods to simplify and clean networks:
    - Energy-based pruning
    - Connectivity-based pruning
    - Chemical stability pruning
    - Pathway relevance pruning
    """
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """
        Initialize the pruner.
        
        Args:
            config: Pruning configuration parameters
        """
        self.config = config or {}
        
        # Default pruning thresholds
        self.energy_threshold = self.config.get('energy_threshold', 50.0)  # kcal/mol
        self.barrier_threshold = self.config.get('barrier_threshold', 80.0)  # kcal/mol
        self.min_pathway_length = self.config.get('min_pathway_length', 2)
        self.max_pathway_length = self.config.get('max_pathway_length', 20)
        self.min_degree = self.config.get('min_degree', 1)
        
    def prune_network(
        self,
        network: ReactionNetwork,
        strategies: Optional[List[str]] = None,
        preserve_seeds: bool = True
    ) -> ReactionNetwork:
        """
        Apply multiple pruning strategies to a network.
        
        Args:
            network: Network to prune
            strategies: List of pruning strategies to apply
            preserve_seeds: Whether to preserve seed molecules
            
        Returns:
            Pruned network
        """
        if strategies is None:
            strategies = [
                'energy_filter',
                'barrier_filter',
                'dead_end_removal',
                'isolated_node_removal',
                'low_degree_removal'
            ]
        
        logger.info(f"Pruning network with strategies: {strategies}")
        original_stats = network.get_statistics()
        
        # Create a copy of the network
        pruned_network = self._copy_network(network)
        
        # Apply each pruning strategy
        for strategy in strategies:
            if hasattr(self, f'_{strategy}'):
                logger.debug(f"Applying {strategy}")
                pruned_network = getattr(self, f'_{strategy}')(pruned_network, preserve_seeds)
            else:
                logger.warning(f"Unknown pruning strategy: {strategy}")
        
        # Log pruning results
        final_stats = pruned_network.get_statistics()
        logger.info(f"Pruning complete: {original_stats['num_molecules']} -> {final_stats['num_molecules']} molecules, "
                   f"{original_stats['num_reactions']} -> {final_stats['num_reactions']} reactions")
        
        return pruned_network
    
    def _copy_network(self, network: ReactionNetwork) -> ReactionNetwork:
        """Create a deep copy of a network."""
        new_network = ReactionNetwork()
        
        # Copy molecules
        for smiles, molecule in network.molecules.items():
            generation = network.generation_info.get(smiles, 0)
            new_network.add_molecule(molecule, generation)
        
        # Copy reactions
        for reaction in network.reactions.values():
            new_network.add_reaction(reaction)
        
        return new_network
    
    def _energy_filter(self, network: ReactionNetwork, preserve_seeds: bool = True) -> ReactionNetwork:
        """
        Remove reactions with high reaction energies.
        
        Args:
            network: Network to filter
            preserve_seeds: Whether to preserve seed molecules
            
        Returns:
            Filtered network
        """
        reactions_to_remove = []
        
        for reaction in network.reactions.values():
            if (reaction.reaction_energy is not None and 
                abs(reaction.reaction_energy) > self.energy_threshold):
                reactions_to_remove.append(reaction.id)
        
        logger.debug(f"Removing {len(reactions_to_remove)} high-energy reactions")
        
        # Remove reactions
        for reaction_id in reactions_to_remove:
            self._remove_reaction(network, reaction_id)
        
        # Clean up orphaned molecules
        if not preserve_seeds:
            self._remove_orphaned_molecules(network)
        
        return network
    
    def _barrier_filter(self, network: ReactionNetwork, preserve_seeds: bool = True) -> ReactionNetwork:
        """
        Remove reactions with high activation barriers.
        
        Args:
            network: Network to filter
            preserve_seeds: Whether to preserve seed molecules
            
        Returns:
            Filtered network
        """
        reactions_to_remove = []
        
        for reaction in network.reactions.values():
            if (reaction.activation_energy is not None and 
                reaction.activation_energy > self.barrier_threshold):
                reactions_to_remove.append(reaction.id)
        
        logger.debug(f"Removing {len(reactions_to_remove)} high-barrier reactions")
        
        # Remove reactions
        for reaction_id in reactions_to_remove:
            self._remove_reaction(network, reaction_id)
        
        # Clean up orphaned molecules
        if not preserve_seeds:
            self._remove_orphaned_molecules(network)
        
        return network
    
    def _dead_end_removal(self, network: ReactionNetwork, preserve_seeds: bool = True) -> ReactionNetwork:
        """
        Remove dead-end pathways that don't lead anywhere.
        
        Args:
            network: Network to prune
            preserve_seeds: Whether to preserve seed molecules
            
        Returns:
            Pruned network
        """
        # Find molecules with no outgoing reactions (sinks)
        sinks = [
            smiles for smiles in network.graph.nodes()
            if network.graph.out_degree(smiles) == 0
        ]
        
        # Find molecules that only lead to dead ends
        molecules_to_remove = set()
        
        for smiles in network.graph.nodes():
            generation = network.generation_info.get(smiles, 0)
            
            # Skip seed molecules if preserving them
            if preserve_seeds and generation == 0:
                continue
            
            # Check if this molecule only leads to sinks or other dead ends
            if self._leads_only_to_dead_ends(network, smiles, sinks):
                molecules_to_remove.add(smiles)
        
        logger.debug(f"Removing {len(molecules_to_remove)} dead-end molecules")
        
        # Remove molecules and associated reactions
        for smiles in molecules_to_remove:
            self._remove_molecule(network, smiles)
        
        return network
    
    def _isolated_node_removal(self, network: ReactionNetwork, preserve_seeds: bool = True) -> ReactionNetwork:
        """
        Remove isolated molecules with no connections.
        
        Args:
            network: Network to prune
            preserve_seeds: Whether to preserve seed molecules
            
        Returns:
            Pruned network
        """
        isolated_molecules = []
        
        for smiles in network.graph.nodes():
            generation = network.generation_info.get(smiles, 0)
            
            # Skip seed molecules if preserving them
            if preserve_seeds and generation == 0:
                continue
            
            # Check if molecule is isolated
            if network.graph.degree(smiles) == 0:
                isolated_molecules.append(smiles)
        
        logger.debug(f"Removing {len(isolated_molecules)} isolated molecules")
        
        # Remove isolated molecules
        for smiles in isolated_molecules:
            self._remove_molecule(network, smiles)
        
        return network
    
    def _low_degree_removal(self, network: ReactionNetwork, preserve_seeds: bool = True) -> ReactionNetwork:
        """
        Remove molecules with very low connectivity.
        
        Args:
            network: Network to prune
            preserve_seeds: Whether to preserve seed molecules
            
        Returns:
            Pruned network
        """
        low_degree_molecules = []
        
        for smiles in network.graph.nodes():
            generation = network.generation_info.get(smiles, 0)
            
            # Skip seed molecules if preserving them
            if preserve_seeds and generation == 0:
                continue
            
            # Check if molecule has low degree
            if network.graph.degree(smiles) < self.min_degree:
                low_degree_molecules.append(smiles)
        
        logger.debug(f"Removing {len(low_degree_molecules)} low-degree molecules")
        
        # Remove low-degree molecules
        for smiles in low_degree_molecules:
            self._remove_molecule(network, smiles)
        
        return network
    
    def _leads_only_to_dead_ends(
        self, 
        network: ReactionNetwork, 
        smiles: str, 
        sinks: List[str],
        visited: Optional[Set[str]] = None
    ) -> bool:
        """
        Check if a molecule only leads to dead ends.
        
        Args:
            network: Network to check
            smiles: Molecule SMILES to check
            sinks: List of sink molecules
            visited: Set of visited molecules (for cycle detection)
            
        Returns:
            True if molecule only leads to dead ends
        """
        if visited is None:
            visited = set()
        
        if smiles in visited:
            return True  # Cycle detected, consider as dead end
        
        if smiles in sinks:
            return True  # This is a sink
        
        visited.add(smiles)
        
        # Check all successors
        successors = list(network.graph.successors(smiles))
        if not successors:
            return True  # No successors, this is a dead end
        
        # All successors must lead to dead ends
        for successor in successors:
            if not self._leads_only_to_dead_ends(network, successor, sinks, visited.copy()):
                return False
        
        return True
    
    def _remove_reaction(self, network: ReactionNetwork, reaction_id: str) -> None:
        """Remove a reaction from the network."""
        if reaction_id in network.reactions:
            reaction = network.reactions[reaction_id]
            
            # Remove edges from graph
            for reactant in reaction.reactants:
                for product in reaction.products:
                    if network.graph.has_edge(reactant.smiles, product.smiles):
                        edge_data = network.graph.get_edge_data(reactant.smiles, product.smiles)
                        if edge_data and edge_data.get('reaction') == reaction:
                            network.graph.remove_edge(reactant.smiles, product.smiles)
            
            # Remove from reactions dict
            del network.reactions[reaction_id]
    
    def _remove_molecule(self, network: ReactionNetwork, smiles: str) -> None:
        """Remove a molecule and all associated reactions from the network."""
        if smiles in network.molecules:
            # Find all reactions involving this molecule
            reactions_to_remove = []
            
            for reaction in network.reactions.values():
                reactant_smiles = [r.smiles for r in reaction.reactants]
                product_smiles = [p.smiles for p in reaction.products]
                
                if smiles in reactant_smiles or smiles in product_smiles:
                    reactions_to_remove.append(reaction.id)
            
            # Remove reactions
            for reaction_id in reactions_to_remove:
                self._remove_reaction(network, reaction_id)
            
            # Remove molecule
            if network.graph.has_node(smiles):
                network.graph.remove_node(smiles)
            
            if smiles in network.molecules:
                del network.molecules[smiles]
            
            if smiles in network.generation_info:
                del network.generation_info[smiles]
    
    def _remove_orphaned_molecules(self, network: ReactionNetwork) -> None:
        """Remove molecules that are no longer connected to any reactions."""
        orphaned_molecules = []
        
        for smiles in network.molecules:
            if network.graph.degree(smiles) == 0:
                generation = network.generation_info.get(smiles, 0)
                if generation > 0:  # Don't remove seed molecules
                    orphaned_molecules.append(smiles)
        
        logger.debug(f"Removing {len(orphaned_molecules)} orphaned molecules")
        
        for smiles in orphaned_molecules:
            self._remove_molecule(network, smiles)
    
    def apply_custom_filter(
        self,
        network: ReactionNetwork,
        filter_func: Callable[[Reaction], bool],
        preserve_seeds: bool = True
    ) -> ReactionNetwork:
        """
        Apply a custom filter function to reactions.
        
        Args:
            network: Network to filter
            filter_func: Function that returns True for reactions to keep
            preserve_seeds: Whether to preserve seed molecules
            
        Returns:
            Filtered network
        """
        reactions_to_remove = []
        
        for reaction in network.reactions.values():
            if not filter_func(reaction):
                reactions_to_remove.append(reaction.id)
        
        logger.debug(f"Custom filter removing {len(reactions_to_remove)} reactions")
        
        # Remove reactions
        for reaction_id in reactions_to_remove:
            self._remove_reaction(network, reaction_id)
        
        # Clean up orphaned molecules
        if not preserve_seeds:
            self._remove_orphaned_molecules(network)
        
        return network
