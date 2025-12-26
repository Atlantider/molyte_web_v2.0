"""
Visualization utilities for RSNet.

This module provides functions for visualizing reaction networks,
molecular structures, and analysis results.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import networkx as nx
import numpy as np
from typing import Dict, List, Optional, Tuple, Any
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import Draw, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
import io
import base64

from ..core.molecule import Molecule
from ..core.reaction import Reaction

if False:  # TYPE_CHECKING hack to avoid import at runtime but keep linting happy if possible, or just avoid Circular
    from ..network.generator import ReactionNetwork

# For runtime usage if needed, import inside function or rely on duck typing
# But better: use string forward reference for type hints



def plot_network_graph(network: "ReactionNetwork", 
                      layout: str = 'spring',
                      figsize: Tuple[int, int] = (12, 8),
                      node_size_scale: float = 1.0,
                      edge_width_scale: float = 1.0,
                      show_labels: bool = True,
                      save_path: Optional[str] = None) -> plt.Figure:
    """
    Plot the reaction network as a graph.
    
    Args:
        network: ReactionNetwork object
        layout: Layout algorithm ('spring', 'circular', 'kamada_kawai', 'shell')
        figsize: Figure size (width, height)
        node_size_scale: Scale factor for node sizes
        edge_width_scale: Scale factor for edge widths
        show_labels: Whether to show molecule labels
        save_path: Path to save the figure (optional)
        
    Returns:
        matplotlib Figure object
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create NetworkX graph
    G = network.to_networkx()
    
    # Choose layout
    if layout == 'spring':
        pos = nx.spring_layout(G, k=1, iterations=50)
    elif layout == 'circular':
        pos = nx.circular_layout(G)
    elif layout == 'kamada_kawai':
        pos = nx.kamada_kawai_layout(G)
    elif layout == 'shell':
        pos = nx.shell_layout(G)
    else:
        pos = nx.spring_layout(G)
    
    # Get node properties
    node_colors = []
    node_sizes = []
    
    for node in G.nodes():
        mol = network.molecules[node]
        
        # Color by generation
        generation = getattr(mol, 'generation', 0)
        color = plt.cm.viridis(generation / max(1, network.get_statistics()['max_generation']))
        node_colors.append(color)
        
        # Size by molecular weight
        size = max(100, mol.molecular_weight * node_size_scale * 5)
        node_sizes.append(size)
    
    # Get edge properties
    edge_colors = []
    edge_widths = []
    
    for edge in G.edges():
        # Color by reaction energy (if available)
        reaction_id = G.edges[edge].get('reaction_id')
        if reaction_id and reaction_id in network.reactions:
            reaction = network.reactions[reaction_id]
            energy = getattr(reaction, 'reaction_energy', 0)
            
            # Color: blue for exothermic, red for endothermic
            if energy < 0:
                edge_colors.append('blue')
            elif energy > 0:
                edge_colors.append('red')
            else:
                edge_colors.append('gray')
            
            # Width by energy magnitude
            width = max(0.5, min(5.0, abs(energy) / 10.0 * edge_width_scale))
            edge_widths.append(width)
        else:
            edge_colors.append('gray')
            edge_widths.append(1.0)
    
    # Draw network
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, 
                          alpha=0.8, ax=ax)
    nx.draw_networkx_edges(G, pos, edge_color=edge_colors, width=edge_widths, 
                          alpha=0.6, ax=ax)
    
    if show_labels:
        # Create simplified labels
        labels = {}
        for node in G.nodes():
            mol = network.molecules[node]
            if hasattr(mol, 'name') and mol.name:
                labels[node] = mol.name
            else:
                # Use simplified SMILES
                smiles = mol.smiles
                if len(smiles) > 8:
                    labels[node] = smiles[:8] + '...'
                else:
                    labels[node] = smiles
        
        nx.draw_networkx_labels(G, pos, labels, font_size=8, ax=ax)
    
    ax.set_title(f'Reaction Network ({len(network.molecules)} molecules, {len(network.reactions)} reactions)')
    ax.axis('off')
    
    # Add colorbar for generations
    sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, 
                              norm=plt.Normalize(vmin=0, vmax=network.get_statistics()['max_generation']))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.8)
    cbar.set_label('Generation')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    return fig


def plot_reaction_pathway(network: "ReactionNetwork",
                         start_smiles: str,
                         end_smiles: str,
                         max_pathways: int = 3,
                         figsize: Tuple[int, int] = (14, 8),
                         save_path: Optional[str] = None) -> plt.Figure:
    """
    Plot reaction pathways between two molecules.
    
    Args:
        network: ReactionNetwork object
        start_smiles: SMILES of starting molecule
        end_smiles: SMILES of target molecule
        max_pathways: Maximum number of pathways to show
        figsize: Figure size
        save_path: Path to save the figure
        
    Returns:
        matplotlib Figure object
    """
    from ..network.analyzer import NetworkAnalyzer
    
    analyzer = NetworkAnalyzer(network)
    pathways = analyzer.find_pathways(start_smiles, end_smiles, max_pathways=max_pathways)
    
    if not pathways:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, f'No pathway found from {start_smiles} to {end_smiles}', 
                ha='center', va='center', transform=ax.transAxes, fontsize=14)
        ax.axis('off')
        return fig
    
    fig, axes = plt.subplots(len(pathways), 1, figsize=figsize, squeeze=False)
    
    for i, pathway in enumerate(pathways):
        ax = axes[i, 0]
        
        # Extract molecules and reactions from pathway
        molecules = pathway['molecules']
        reactions = pathway['reactions']
        
        # Create positions for molecules
        x_positions = np.linspace(0, 10, len(molecules))
        y_position = 0
        
        # Plot molecules
        for j, mol_smiles in enumerate(molecules):
            mol = network.molecules[mol_smiles]
            
            # Draw molecule circle
            circle = patches.Circle((x_positions[j], y_position), 0.3, 
                                  facecolor='lightblue', edgecolor='black')
            ax.add_patch(circle)
            
            # Add molecule label
            label = mol.name if hasattr(mol, 'name') and mol.name else mol_smiles[:6]
            ax.text(x_positions[j], y_position-0.6, label, ha='center', fontsize=8)
        
        # Plot reactions (arrows)
        for j in range(len(reactions)):
            reaction = network.reactions[reactions[j]]
            
            # Arrow from molecule j to molecule j+1
            arrow_start = x_positions[j] + 0.3
            arrow_end = x_positions[j+1] - 0.3
            
            ax.annotate('', xy=(arrow_end, y_position), xytext=(arrow_start, y_position),
                       arrowprops=dict(arrowstyle='->', lw=2, color='red'))
            
            # Add reaction energy if available
            if hasattr(reaction, 'reaction_energy') and reaction.reaction_energy is not None:
                energy_text = f'{reaction.reaction_energy:.1f}'
                ax.text((arrow_start + arrow_end) / 2, y_position + 0.2, energy_text, 
                       ha='center', fontsize=8, color='red')
        
        ax.set_xlim(-0.5, 10.5)
        ax.set_ylim(-1, 1)
        ax.set_title(f'Pathway {i+1} (Length: {len(molecules)-1} steps)')
        ax.axis('off')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    return fig


def plot_energy_distribution(network: "ReactionNetwork",
                           figsize: Tuple[int, int] = (10, 6),
                           save_path: Optional[str] = None) -> plt.Figure:
    """
    Plot distribution of reaction energies in the network.
    
    Args:
        network: ReactionNetwork object
        figsize: Figure size
        save_path: Path to save the figure
        
    Returns:
        matplotlib Figure object
    """
    # Collect reaction energies
    energies = []
    for reaction in network.reactions.values():
        if hasattr(reaction, 'reaction_energy') and reaction.reaction_energy is not None:
            energies.append(reaction.reaction_energy)
    
    if not energies:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, 'No reaction energy data available', 
                ha='center', va='center', transform=ax.transAxes, fontsize=14)
        ax.axis('off')
        return fig
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    
    # Histogram
    ax1.hist(energies, bins=20, alpha=0.7, color='skyblue', edgecolor='black')
    ax1.axvline(0, color='red', linestyle='--', label='Thermoneutral')
    ax1.set_xlabel('Reaction Energy (kcal/mol)')
    ax1.set_ylabel('Number of Reactions')
    ax1.set_title('Reaction Energy Distribution')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Box plot
    ax2.boxplot(energies, vert=True)
    ax2.axhline(0, color='red', linestyle='--', label='Thermoneutral')
    ax2.set_ylabel('Reaction Energy (kcal/mol)')
    ax2.set_title('Energy Distribution Summary')
    ax2.grid(True, alpha=0.3)
    
    # Add statistics
    stats_text = f'Mean: {np.mean(energies):.1f} kcal/mol\n'
    stats_text += f'Std: {np.std(energies):.1f} kcal/mol\n'
    stats_text += f'Min: {np.min(energies):.1f} kcal/mol\n'
    stats_text += f'Max: {np.max(energies):.1f} kcal/mol'
    
    ax2.text(0.02, 0.98, stats_text, transform=ax2.transAxes, 
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    return fig


def plot_molecule_structure(molecule: Molecule,
                          figsize: Tuple[int, int] = (6, 6),
                          save_path: Optional[str] = None) -> plt.Figure:
    """
    Plot 2D structure of a molecule.
    
    Args:
        molecule: Molecule object
        figsize: Figure size
        save_path: Path to save the figure
        
    Returns:
        matplotlib Figure object
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    try:
        # Generate 2D coordinates
        mol = molecule.rdkit_mol
        rdDepictor.Compute2DCoords(mol)
        
        # Create molecule image
        drawer = rdMolDraw2D.MolDraw2DCairo(400, 400)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        
        # Convert to matplotlib image
        img_data = drawer.GetDrawingText()
        img = plt.imread(io.BytesIO(img_data))
        
        ax.imshow(img)
        ax.axis('off')
        
        # Add title with molecule info
        title = f'{molecule.smiles}'
        if hasattr(molecule, 'name') and molecule.name:
            title = f'{molecule.name}\n{title}'
        title += f'\nMW: {molecule.molecular_weight:.1f} g/mol'
        
        ax.set_title(title, fontsize=12)
        
    except Exception as e:
        ax.text(0.5, 0.5, f'Error drawing molecule:\n{str(e)}', 
                ha='center', va='center', transform=ax.transAxes)
        ax.axis('off')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    return fig


def plot_generation_statistics(network: "ReactionNetwork",
                             figsize: Tuple[int, int] = (12, 5),
                             save_path: Optional[str] = None) -> plt.Figure:
    """
    Plot statistics by generation.
    
    Args:
        network: ReactionNetwork object
        figsize: Figure size
        save_path: Path to save the figure
        
    Returns:
        matplotlib Figure object
    """
    stats = network.get_statistics()
    molecules_by_gen = stats['molecules_by_generation']
    
    generations = sorted(molecules_by_gen.keys())
    molecule_counts = [molecules_by_gen[gen] for gen in generations]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    
    # Molecules by generation
    ax1.bar(generations, molecule_counts, alpha=0.7, color='skyblue', edgecolor='black')
    ax1.set_xlabel('Generation')
    ax1.set_ylabel('Number of Molecules')
    ax1.set_title('Molecules by Generation')
    ax1.grid(True, alpha=0.3)
    
    # Cumulative molecules
    cumulative = np.cumsum(molecule_counts)
    ax2.plot(generations, cumulative, marker='o', linewidth=2, markersize=6)
    ax2.set_xlabel('Generation')
    ax2.set_ylabel('Cumulative Molecules')
    ax2.set_title('Cumulative Network Growth')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    return fig
