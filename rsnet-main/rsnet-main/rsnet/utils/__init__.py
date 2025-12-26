"""
Utility functions for RSNet.
"""

from .structure_analysis import (
    detect_heteroatoms,
    detect_polar_bonds,
    detect_acidic_hydrogens,
    detect_hydrogen_bond_acceptors,
    detect_weak_bonds,
    detect_small_rings,
    detect_pi_systems,
    get_structure_tags,
    find_hydrogen_transfer_sites
)

from .io import (
    load_molecules_from_smiles_list,
    load_molecules_from_smiles,
    save_molecules_to_smiles,
    save_network,
    load_network,
    export_reaction_table,
    export_molecule_table,
    create_analysis_report
)

from .visualization import (
    plot_network_graph,
    plot_reaction_pathway,
    plot_energy_distribution,
    plot_molecule_structure,
    plot_generation_statistics
)

from .config import (
    RSNetConfig,
    get_default_config,
    load_config,
    save_config,
    get_preset_config,
    print_config_summary
)

__all__ = [
    # Structure analysis
    "detect_heteroatoms",
    "detect_polar_bonds",
    "detect_acidic_hydrogens",
    "detect_hydrogen_bond_acceptors",
    "detect_weak_bonds",
    "detect_small_rings",
    "detect_pi_systems",
    "get_structure_tags",
    "find_hydrogen_transfer_sites",

    # I/O
    "load_molecules_from_smiles_list",
    "load_molecules_from_smiles",
    "save_molecules_to_smiles",
    "save_network",
    "load_network",
    "export_reaction_table",
    "export_molecule_table",
    "create_analysis_report",

    # Visualization
    "plot_network_graph",
    "plot_reaction_pathway",
    "plot_energy_distribution",
    "plot_molecule_structure",
    "plot_generation_statistics",

    # Configuration
    "RSNetConfig",
    "get_default_config",
    "load_config",
    "save_config",
    "get_preset_config",
    "print_config_summary",
]
