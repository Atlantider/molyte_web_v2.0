"""
Network analysis tools for RSNet.

This module provides tools for analyzing reaction networks:
1. Dominant pathway identification
2. Network topology analysis
3. Key intermediate identification
4. Thermodynamic analysis
"""

import logging
from typing import List, Dict, Set, Optional, Tuple, Any
import networkx as nx
import numpy as np
from collections import defaultdict, Counter

from ..core.molecule import Molecule
from ..core.reaction import Reaction
from .generator import ReactionNetwork


logger = logging.getLogger(__name__)


class NetworkAnalyzer:
    """
    Analyzer for reaction networks.
    
    Provides methods for:
    - Finding dominant reaction pathways
    - Analyzing network topology
    - Identifying key intermediates
    - Computing thermodynamic properties
    """
    
    def __init__(self, network: ReactionNetwork):
        """
        Initialize the analyzer.
        
        Args:
            network: Reaction network to analyze
        """
        self.network = network
        self.graph = network.graph
        
    def find_dominant_pathways(
        self,
        source_molecules: List[str],
        target_molecules: Optional[List[str]] = None,
        max_pathways: int = 10,
        energy_weighted: bool = True
    ) -> List[Dict[str, Any]]:
        """
        Find dominant reaction pathways in the network.
        
        Args:
            source_molecules: Starting molecules (SMILES)
            target_molecules: Target molecules (SMILES). If None, find paths to all sinks
            max_pathways: Maximum number of pathways to return
            energy_weighted: Whether to use energy-weighted path finding
            
        Returns:
            List of pathway dictionaries with path info and energetics
        """
        pathways = []
        
        # If no targets specified, find all sink nodes (no outgoing edges)
        if target_molecules is None:
            target_molecules = [
                node for node in self.graph.nodes()
                if self.graph.out_degree(node) == 0
            ]
        
        logger.info(f"Finding pathways from {len(source_molecules)} sources to {len(target_molecules)} targets")
        
        # Find paths between each source-target pair
        for source in source_molecules:
            if source not in self.graph.nodes():
                logger.warning(f"Source molecule {source} not found in network")
                continue
                
            for target in target_molecules:
                if target not in self.graph.nodes():
                    continue
                
                try:
                    if energy_weighted:
                        # Use energy-weighted shortest paths
                        paths = list(nx.shortest_simple_paths(
                            self.graph, source, target, weight='weight'
                        ))
                    else:
                        # Use unweighted shortest paths
                        paths = list(nx.shortest_simple_paths(
                            self.graph, source, target
                        ))
                    
                    # Analyze each path
                    for i, path in enumerate(paths[:max_pathways]):
                        pathway_info = self._analyze_pathway(path)
                        pathway_info.update({
                            'source': source,
                            'target': target,
                            'rank': i + 1,
                            'path_length': len(path) - 1
                        })
                        pathways.append(pathway_info)
                        
                except nx.NetworkXNoPath:
                    logger.debug(f"No path found from {source} to {target}")
                    continue
                except Exception as e:
                    logger.warning(f"Error finding paths from {source} to {target}: {e}")
                    continue
        
        # Sort pathways by total energy
        pathways.sort(key=lambda x: x.get('total_energy', float('inf')))
        
        return pathways[:max_pathways]
    
    def _analyze_pathway(self, path: List[str]) -> Dict[str, Any]:
        """
        Analyze a single reaction pathway.
        
        Args:
            path: List of molecule SMILES in the pathway
            
        Returns:
            Dictionary with pathway analysis
        """
        reactions = []
        total_energy = 0.0
        max_barrier = 0.0
        
        # Analyze each step in the pathway
        for i in range(len(path) - 1):
            source_smiles = path[i]
            target_smiles = path[i + 1]
            
            # Find the reaction between these molecules
            edge_data = self.graph.get_edge_data(source_smiles, target_smiles)
            if edge_data and 'reaction' in edge_data:
                reaction = edge_data['reaction']
                reactions.append(reaction)
                
                if reaction.reaction_energy is not None:
                    total_energy += reaction.reaction_energy
                
                if reaction.activation_energy is not None:
                    max_barrier = max(max_barrier, reaction.activation_energy)
        
        return {
            'path': path,
            'reactions': reactions,
            'total_energy': total_energy,
            'max_activation_energy': max_barrier,
            'num_steps': len(reactions),
            'molecules': [self.network.molecules[smiles] for smiles in path]
        }
    
    def identify_key_intermediates(self, top_n: int = 10) -> List[Dict[str, Any]]:
        """
        Identify key intermediate molecules in the network.
        
        Uses centrality measures to identify important intermediates.
        
        Args:
            top_n: Number of top intermediates to return
            
        Returns:
            List of intermediate info dictionaries
        """
        # Calculate centrality measures
        betweenness = nx.betweenness_centrality(self.graph)
        closeness = nx.closeness_centrality(self.graph)
        degree = dict(self.graph.degree())
        
        # Combine measures
        intermediates = []
        for smiles in self.graph.nodes():
            if smiles in self.network.molecules:
                molecule = self.network.molecules[smiles]
                
                # Skip if it's a seed molecule (generation 0) or sink (no outgoing edges)
                generation = self.network.generation_info.get(smiles, 0)
                out_degree = self.graph.out_degree(smiles)
                in_degree = self.graph.in_degree(smiles)
                
                if generation > 0 and out_degree > 0 and in_degree > 0:
                    intermediates.append({
                        'molecule': molecule,
                        'smiles': smiles,
                        'generation': generation,
                        'betweenness_centrality': betweenness.get(smiles, 0),
                        'closeness_centrality': closeness.get(smiles, 0),
                        'degree': degree.get(smiles, 0),
                        'in_degree': in_degree,
                        'out_degree': out_degree,
                        'centrality_score': (
                            betweenness.get(smiles, 0) * 0.4 +
                            closeness.get(smiles, 0) * 0.3 +
                            degree.get(smiles, 0) / max(degree.values()) * 0.3
                        )
                    })
        
        # Sort by centrality score
        intermediates.sort(key=lambda x: x['centrality_score'], reverse=True)
        
        return intermediates[:top_n]
    
    def analyze_network_topology(self) -> Dict[str, Any]:
        """
        Analyze the overall network topology.
        
        Returns:
            Dictionary with topology metrics
        """
        graph = self.graph
        
        # Basic metrics
        num_nodes = graph.number_of_nodes()
        num_edges = graph.number_of_edges()
        
        # Degree distribution
        degrees = [d for n, d in graph.degree()]
        in_degrees = [d for n, d in graph.in_degree()]
        out_degrees = [d for n, d in graph.out_degree()]
        
        # Connectivity
        is_connected = nx.is_weakly_connected(graph)
        num_components = nx.number_weakly_connected_components(graph)
        
        # Path lengths (for connected components)
        avg_path_length = None
        diameter = None
        
        if is_connected:
            try:
                avg_path_length = nx.average_shortest_path_length(graph)
                diameter = nx.diameter(graph.to_undirected())
            except:
                pass
        
        # Clustering (convert to undirected for clustering)
        undirected = graph.to_undirected()
        clustering_coefficient = nx.average_clustering(undirected)
        
        # Generation analysis
        generation_counts = Counter(self.network.generation_info.values())
        
        return {
            'num_molecules': num_nodes,
            'num_reactions': num_edges,
            'is_connected': is_connected,
            'num_components': num_components,
            'average_path_length': avg_path_length,
            'diameter': diameter,
            'clustering_coefficient': clustering_coefficient,
            'degree_stats': {
                'mean': np.mean(degrees),
                'std': np.std(degrees),
                'max': max(degrees) if degrees else 0,
                'min': min(degrees) if degrees else 0
            },
            'in_degree_stats': {
                'mean': np.mean(in_degrees),
                'std': np.std(in_degrees),
                'max': max(in_degrees) if in_degrees else 0
            },
            'out_degree_stats': {
                'mean': np.mean(out_degrees),
                'std': np.std(out_degrees),
                'max': max(out_degrees) if out_degrees else 0
            },
            'generation_distribution': dict(generation_counts),
            'max_generation': max(generation_counts.keys()) if generation_counts else 0
        }
    
    def analyze_thermodynamics(self) -> Dict[str, Any]:
        """
        Analyze thermodynamic properties of the network.
        
        Returns:
            Dictionary with thermodynamic analysis
        """
        reaction_energies = []
        activation_energies = []
        exothermic_count = 0
        endothermic_count = 0
        
        for reaction in self.network.reactions.values():
            if reaction.reaction_energy is not None:
                reaction_energies.append(reaction.reaction_energy)
                if reaction.reaction_energy < 0:
                    exothermic_count += 1
                else:
                    endothermic_count += 1
            
            if reaction.activation_energy is not None:
                activation_energies.append(reaction.activation_energy)
        
        return {
            'num_reactions_with_energies': len(reaction_energies),
            'reaction_energy_stats': {
                'mean': np.mean(reaction_energies) if reaction_energies else None,
                'std': np.std(reaction_energies) if reaction_energies else None,
                'min': min(reaction_energies) if reaction_energies else None,
                'max': max(reaction_energies) if reaction_energies else None,
                'median': np.median(reaction_energies) if reaction_energies else None
            },
            'activation_energy_stats': {
                'mean': np.mean(activation_energies) if activation_energies else None,
                'std': np.std(activation_energies) if activation_energies else None,
                'min': min(activation_energies) if activation_energies else None,
                'max': max(activation_energies) if activation_energies else None,
                'median': np.median(activation_energies) if activation_energies else None
            },
            'exothermic_reactions': exothermic_count,
            'endothermic_reactions': endothermic_count,
            'thermodynamic_favorability': exothermic_count / len(reaction_energies) if reaction_energies else 0
        }
    
    def export_network_data(self, format: str = 'json') -> Dict[str, Any]:
        """
        Export network data in various formats.
        
        Args:
            format: Export format ('json', 'graphml', 'gexf')
            
        Returns:
            Network data dictionary
        """
        if format == 'json':
            return {
                'molecules': {
                    smiles: {
                        'name': mol.name,
                        'formula': mol.formula,
                        'molecular_weight': mol.molecular_weight,
                        'generation': self.network.generation_info.get(smiles, 0)
                    }
                    for smiles, mol in self.network.molecules.items()
                },
                'reactions': {
                    reaction.id: {
                        'name': reaction.name,
                        'reactants': [r.smiles for r in reaction.reactants],
                        'products': [p.smiles for p in reaction.products],
                        'reaction_energy': reaction.reaction_energy,
                        'activation_energy': reaction.activation_energy,
                        'operator': getattr(reaction, 'operator', None),
                        'template': getattr(reaction, 'template', None)
                    }
                    for reaction in self.network.reactions.values()
                },
                'network_stats': self.network.get_statistics(),
                'topology_analysis': self.analyze_network_topology(),
                'thermodynamic_analysis': self.analyze_thermodynamics()
            }
        else:
            raise ValueError(f"Unsupported export format: {format}")
    
    def find_reaction_cycles(self, max_cycle_length: int = 10) -> List[List[str]]:
        """
        Find reaction cycles in the network.
        
        Args:
            max_cycle_length: Maximum cycle length to search for
            
        Returns:
            List of cycles (each cycle is a list of molecule SMILES)
        """
        try:
            cycles = list(nx.simple_cycles(self.graph))
            # Filter by length
            cycles = [cycle for cycle in cycles if len(cycle) <= max_cycle_length]
            return cycles
        except Exception as e:
            logger.warning(f"Error finding cycles: {e}")
            return []
