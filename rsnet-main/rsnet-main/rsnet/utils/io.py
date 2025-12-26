"""
Input/Output utilities for RSNet.

This module provides functions for loading and saving molecules,
networks, and analysis results.
"""

import json
import pickle
import yaml
import csv
from pathlib import Path
from typing import Dict, List, Optional, Any, Union
import networkx as nx

from ..core.molecule import Molecule
from ..core.environment import Environment
from ..core.reaction import Reaction
# Import ReactionNetwork at runtime to avoid circular imports


def load_molecules_from_smiles_list(smiles_list: List[str],
                                   names: Optional[List[str]] = None) -> List[Molecule]:
    """
    从SMILES字符串列表加载分子。

    Args:
        smiles_list: SMILES字符串列表
        names: 分子名称列表（可选）

    Returns:
        分子对象列表
    """
    molecules = []

    for i, smiles in enumerate(smiles_list):
        try:
            # 确定分子名称
            if names and i < len(names):
                name = names[i]
            else:
                name = f'mol_{i+1}'

            # 创建分子
            mol = Molecule.from_smiles(smiles, name=name)
            molecules.append(mol)

        except Exception as e:
            print(f"Warning: Could not load SMILES '{smiles}': {e}")
            continue

    return molecules


def load_molecules_from_smiles(file_path: Union[str, Path],
                              name_column: Optional[str] = None) -> List[Molecule]:
    """
    Load molecules from a file containing SMILES strings.
    
    Args:
        file_path: Path to file (CSV, TXT, or JSON)
        name_column: Column name for molecule names (CSV only)
        
    Returns:
        List of Molecule objects
    """
    file_path = Path(file_path)
    molecules = []
    
    if file_path.suffix.lower() == '.csv':
        # CSV file
        import pandas as pd
        df = pd.read_csv(file_path)
        
        # Assume first column is SMILES if not specified
        smiles_column = df.columns[0]
        
        for idx, row in df.iterrows():
            smiles = row[smiles_column]
            name = row[name_column] if name_column and name_column in df.columns else f"mol_{idx}"
            
            try:
                mol = Molecule.from_smiles(smiles, name=name)
                molecules.append(mol)
            except Exception as e:
                print(f"Warning: Could not create molecule from SMILES '{smiles}': {e}")
    
    elif file_path.suffix.lower() == '.json':
        # JSON file
        with open(file_path, 'r') as f:
            data = json.load(f)
        
        if isinstance(data, list):
            for i, item in enumerate(data):
                if isinstance(item, str):
                    # Simple SMILES list
                    smiles = item
                    name = f"mol_{i}"
                elif isinstance(item, dict):
                    # Dictionary with smiles and name
                    smiles = item.get('smiles', item.get('SMILES'))
                    name = item.get('name', item.get('Name', f"mol_{i}"))
                else:
                    continue
                
                try:
                    mol = Molecule.from_smiles(smiles, name=name)
                    molecules.append(mol)
                except Exception as e:
                    print(f"Warning: Could not create molecule from SMILES '{smiles}': {e}")
    
    elif file_path.suffix.lower() == '.txt':
        # Text file (one SMILES per line)
        with open(file_path, 'r') as f:
            for i, line in enumerate(f):
                smiles = line.strip()
                if smiles and not smiles.startswith('#'):  # Skip empty lines and comments
                    try:
                        mol = Molecule.from_smiles(smiles, name=f"mol_{i}")
                        molecules.append(mol)
                    except Exception as e:
                        print(f"Warning: Could not create molecule from SMILES '{smiles}': {e}")
    
    else:
        raise ValueError(f"Unsupported file format: {file_path.suffix}")
    
    return molecules


def save_molecules_to_smiles(molecules: List[Molecule], 
                           file_path: Union[str, Path],
                           format: str = 'csv') -> None:
    """
    Save molecules to a file as SMILES strings.
    
    Args:
        molecules: List of Molecule objects
        file_path: Output file path
        format: Output format ('csv', 'json', 'txt')
    """
    file_path = Path(file_path)
    
    if format.lower() == 'csv':
        import pandas as pd
        data = []
        for mol in molecules:
            data.append({
                'smiles': mol.smiles,
                'name': getattr(mol, 'name', ''),
                'molecular_weight': mol.molecular_weight,
                'formula': mol.molecular_formula
            })
        
        df = pd.DataFrame(data)
        df.to_csv(file_path, index=False)
    
    elif format.lower() == 'json':
        data = []
        for mol in molecules:
            data.append({
                'smiles': mol.smiles,
                'name': getattr(mol, 'name', ''),
                'molecular_weight': mol.molecular_weight,
                'formula': mol.molecular_formula
            })
        
        with open(file_path, 'w') as f:
            json.dump(data, f, indent=2)
    
    elif format.lower() == 'txt':
        with open(file_path, 'w') as f:
            for mol in molecules:
                f.write(f"{mol.smiles}\n")
    
    else:
        raise ValueError(f"Unsupported format: {format}")


def save_network(network,
                file_path: Union[str, Path],
                format: str = 'json',
                include_coordinates: bool = False) -> None:
    """
    Save reaction network to file.
    
    Args:
        network: ReactionNetwork object
        file_path: Output file path
        format: Output format ('json', 'pickle', 'graphml')
        include_coordinates: Whether to include 3D coordinates
    """
    file_path = Path(file_path)
    
    if format.lower() == 'json':
        # Export to JSON
        data = network.to_dict(include_coordinates=include_coordinates)
        with open(file_path, 'w') as f:
            json.dump(data, f, indent=2)
    
    elif format.lower() == 'pickle':
        # Export to pickle (preserves all Python objects)
        with open(file_path, 'wb') as f:
            pickle.dump(network, f)
    
    elif format.lower() == 'graphml':
        # Export to GraphML (for network analysis tools)
        G = network.to_networkx()
        nx.write_graphml(G, file_path)
    
    else:
        raise ValueError(f"Unsupported format: {format}")


def load_network(file_path: Union[str, Path],
                format: Optional[str] = None):
    """
    Load reaction network from file.
    
    Args:
        file_path: Input file path
        format: Input format (auto-detected if None)
        
    Returns:
        ReactionNetwork object
    """
    file_path = Path(file_path)
    
    if format is None:
        format = file_path.suffix.lower().lstrip('.')
    
    if format == 'pickle':
        with open(file_path, 'rb') as f:
            return pickle.load(f)

    else:
        raise ValueError(f"Unsupported format for loading: {format}. Currently only 'pickle' is supported.")


def export_reaction_table(network,
                         file_path: Union[str, Path],
                         include_energies: bool = True) -> None:
    """
    Export reactions to a CSV table.
    
    Args:
        network: ReactionNetwork object
        file_path: Output CSV file path
        include_energies: Whether to include energy data
    """
    import pandas as pd
    
    data = []
    for reaction_id, reaction in network.reactions.items():
        row = {
            'reaction_id': reaction_id,
            'reactants': ' + '.join([mol.smiles for mol in reaction.reactants]),
            'products': ' + '.join([mol.smiles for mol in reaction.products]),
            'operator': reaction.operator_name,
            'reactant_names': ' + '.join([getattr(mol, 'name', mol.smiles) for mol in reaction.reactants]),
            'product_names': ' + '.join([getattr(mol, 'name', mol.smiles) for mol in reaction.products])
        }
        
        if include_energies:
            row['reaction_energy'] = getattr(reaction, 'reaction_energy', None)
            row['activation_energy'] = getattr(reaction, 'activation_energy', None)
            row['driving_force_score'] = getattr(reaction, 'driving_force_score', None)
        
        data.append(row)
    
    df = pd.DataFrame(data)
    df.to_csv(file_path, index=False)


def export_molecule_table(network,
                         file_path: Union[str, Path]) -> None:
    """
    Export molecules to a CSV table.
    
    Args:
        network: ReactionNetwork object
        file_path: Output CSV file path
    """
    import pandas as pd
    
    data = []
    for smiles, mol in network.molecules.items():
        row = {
            'smiles': smiles,
            'name': getattr(mol, 'name', ''),
            'molecular_weight': mol.molecular_weight,
            'formula': mol.molecular_formula,
            'generation': getattr(mol, 'generation', 0)
        }
        data.append(row)
    
    df = pd.DataFrame(data)
    df.to_csv(file_path, index=False)


def load_config(file_path: Union[str, Path]) -> Dict[str, Any]:
    """
    Load configuration from YAML file.
    
    Args:
        file_path: Path to YAML config file
        
    Returns:
        Configuration dictionary
    """
    with open(file_path, 'r') as f:
        return yaml.safe_load(f)


def save_config(config: Dict[str, Any], file_path: Union[str, Path]) -> None:
    """
    Save configuration to YAML file.
    
    Args:
        config: Configuration dictionary
        file_path: Output YAML file path
    """
    with open(file_path, 'w') as f:
        yaml.dump(config, f, default_flow_style=False, indent=2)


def create_analysis_report(network,
                          output_dir: Union[str, Path],
                          include_plots: bool = True) -> None:
    """
    Create a comprehensive analysis report.
    
    Args:
        network: ReactionNetwork object
        output_dir: Output directory for report files
        include_plots: Whether to generate plots
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Export data tables
    export_reaction_table(network, output_dir / 'reactions.csv')
    export_molecule_table(network, output_dir / 'molecules.csv')
    
    # Export network
    save_network(network, output_dir / 'network.json', format='json')
    save_network(network, output_dir / 'network.graphml', format='graphml')
    
    # Generate analysis
    from ..network.analyzer import NetworkAnalyzer
    analyzer = NetworkAnalyzer(network)
    
    # Network statistics
    stats = network.get_statistics()
    with open(output_dir / 'statistics.json', 'w') as f:
        json.dump(stats, f, indent=2)
    
    # Topology analysis
    topology = analyzer.analyze_network_topology()
    with open(output_dir / 'topology.json', 'w') as f:
        json.dump(topology, f, indent=2)
    
    # Generate plots if requested
    if include_plots:
        from .visualization import (plot_network_graph, plot_energy_distribution, 
                                   plot_generation_statistics)
        
        try:
            # Network graph
            fig = plot_network_graph(network)
            fig.savefig(output_dir / 'network_graph.png', dpi=300, bbox_inches='tight')
            plt.close(fig)
            
            # Energy distribution
            fig = plot_energy_distribution(network)
            fig.savefig(output_dir / 'energy_distribution.png', dpi=300, bbox_inches='tight')
            plt.close(fig)
            
            # Generation statistics
            fig = plot_generation_statistics(network)
            fig.savefig(output_dir / 'generation_stats.png', dpi=300, bbox_inches='tight')
            plt.close(fig)
            
        except Exception as e:
            print(f"Warning: Could not generate plots: {e}")
    
    # Create summary report
    create_summary_report(network, output_dir / 'summary.txt')


def create_summary_report(network, file_path: Union[str, Path]) -> None:
    """
    Create a text summary report.
    
    Args:
        network: ReactionNetwork object
        file_path: Output text file path
    """
    stats = network.get_statistics()
    
    with open(file_path, 'w') as f:
        f.write("RSNet Reaction Network Analysis Report\n")
        f.write("=" * 50 + "\n\n")
        
        f.write("Network Statistics:\n")
        f.write(f"  Total molecules: {stats['num_molecules']}\n")
        f.write(f"  Total reactions: {stats['num_reactions']}\n")
        f.write(f"  Generations: {stats['max_generation'] + 1}\n")
        f.write(f"  Average degree: {stats.get('avg_degree', 'N/A'):.2f}\n\n")
        
        f.write("Molecules by Generation:\n")
        for gen, count in sorted(stats['molecules_by_generation'].items()):
            f.write(f"  Generation {gen}: {count} molecules\n")
        f.write("\n")
        
        # Energy statistics
        energies = []
        for reaction in network.reactions.values():
            if hasattr(reaction, 'reaction_energy') and reaction.reaction_energy is not None:
                energies.append(reaction.reaction_energy)
        
        if energies:
            import numpy as np
            f.write("Energy Statistics:\n")
            f.write(f"  Mean reaction energy: {np.mean(energies):.2f} kcal/mol\n")
            f.write(f"  Energy range: {np.min(energies):.2f} to {np.max(energies):.2f} kcal/mol\n")
            f.write(f"  Exothermic reactions: {sum(1 for e in energies if e < 0)}\n")
            f.write(f"  Endothermic reactions: {sum(1 for e in energies if e > 0)}\n\n")
        
        f.write("Files Generated:\n")
        f.write("  - reactions.csv: Detailed reaction table\n")
        f.write("  - molecules.csv: Molecule properties\n")
        f.write("  - network.json: Complete network data\n")
        f.write("  - network.graphml: Network for external analysis\n")
        f.write("  - statistics.json: Detailed statistics\n")
        f.write("  - topology.json: Network topology analysis\n")
        if Path(file_path).parent.joinpath('network_graph.png').exists():
            f.write("  - network_graph.png: Network visualization\n")
            f.write("  - energy_distribution.png: Energy analysis\n")
            f.write("  - generation_stats.png: Growth statistics\n")
