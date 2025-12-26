"""
Network generation configuration for RSNet.

This module provides configuration classes and presets for different
types of reaction network generation scenarios.
"""

from typing import Dict, Any, Optional, List
from dataclasses import dataclass, field
from enum import Enum


class NetworkGenerationStrategy(Enum):
    """Network generation strategies."""
    EXHAUSTIVE = "exhaustive"  # Apply all applicable operators
    INTELLIGENT = "intelligent"  # Use intelligent operator selection
    FOCUSED = "focused"  # Focus on specific reaction types
    ADAPTIVE = "adaptive"  # Adapt strategy based on conditions


class OperatorSelectionStrategy(Enum):
    """Operator selection strategies."""
    ALL = "all"  # Use all active operators
    TOP_N = "top_n"  # Use top N operators by priority
    ADAPTIVE = "adaptive"  # Adaptive selection based on conditions
    THRESHOLD = "threshold"  # Use operators above priority threshold


@dataclass
class NetworkGenerationConfig:
    """
    Configuration for network generation.
    
    This class provides comprehensive configuration options for
    controlling reaction network generation behavior.
    """
    
    # Basic generation parameters
    max_generations: int = 5
    max_species: int = 100
    max_time: Optional[float] = None  # seconds
    
    # Energy and feasibility criteria
    energy_cutoff: float = 50.0  # kcal/mol
    activation_energy_cutoff: float = 80.0  # kcal/mol
    activation_energy_cutoff: float = 80.0  # kcal/mol
    min_new_reactions_per_generation: int = 1
    max_reactions_per_generation: int = 1000  # Safety limit per generation
    
    # Intelligent selection parameters
    generation_strategy: NetworkGenerationStrategy = NetworkGenerationStrategy.INTELLIGENT
    operator_selection_strategy: OperatorSelectionStrategy = OperatorSelectionStrategy.ADAPTIVE
    max_operators_per_generation: int = 5
    driving_force_threshold: float = 0.3
    operator_priority_threshold: float = 0.5
    
    # Structure-based filtering
    use_structure_based_filtering: bool = True
    max_molecular_weight: float = 500.0  # g/mol
    min_molecular_weight: float = 16.0   # g/mol
    max_ring_size: int = 8
    allow_unstable_intermediates: bool = False
    
    # Performance optimization
    parallel_screening: bool = True
    max_screening_workers: int = 4
    cache_molecular_properties: bool = True
    
    # Reaction filtering
    filter_duplicate_reactions: bool = True
    filter_trivial_reactions: bool = True
    min_structural_change: float = 0.1  # Minimum Tanimoto distance change
    
    # Advanced options
    prioritize_novel_structures: bool = True
    bias_towards_stable_products: bool = True
    consider_stereochemistry: bool = False
    
    # Logging and debugging
    verbose_logging: bool = False
    save_intermediate_results: bool = False
    debug_operator_selection: bool = False
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert configuration to dictionary."""
        result = {}
        for key, value in self.__dict__.items():
            if isinstance(value, Enum):
                result[key] = value.value
            else:
                result[key] = value
        return result
    
    @classmethod
    def from_dict(cls, config_dict: Dict[str, Any]) -> 'NetworkGenerationConfig':
        """Create configuration from dictionary."""
        # Handle enum conversions
        if 'generation_strategy' in config_dict:
            config_dict['generation_strategy'] = NetworkGenerationStrategy(
                config_dict['generation_strategy']
            )
        if 'operator_selection_strategy' in config_dict:
            config_dict['operator_selection_strategy'] = OperatorSelectionStrategy(
                config_dict['operator_selection_strategy']
            )
        
        return cls(**config_dict)


class NetworkGenerationPresets:
    """
    Predefined configuration presets for common scenarios.
    """
    
    @staticmethod
    def battery_chemistry() -> NetworkGenerationConfig:
        """Configuration optimized for battery chemistry applications."""
        return NetworkGenerationConfig(
            max_generations=7,
            max_species=5000, # Increased for larger SEI networks
            energy_cutoff=60.0,
            generation_strategy=NetworkGenerationStrategy.INTELLIGENT,
            operator_selection_strategy=OperatorSelectionStrategy.ADAPTIVE,
            driving_force_threshold=0.2,
            max_operators_per_generation=6,
            use_structure_based_filtering=True,
            prioritize_novel_structures=True,
            bias_towards_stable_products=False,
            allow_unstable_intermediates=True,
            verbose_logging=True,
            max_reactions_per_generation=2000 # Optimized safety limit
        )
    
    @staticmethod
    def organic_synthesis() -> NetworkGenerationConfig:
        """Configuration optimized for organic synthesis pathways."""
        return NetworkGenerationConfig(
            max_generations=6,
            max_species=150,
            energy_cutoff=45.0,
            activation_energy_cutoff=70.0,
            generation_strategy=NetworkGenerationStrategy.INTELLIGENT,
            operator_selection_strategy=OperatorSelectionStrategy.TOP_N,
            driving_force_threshold=0.4,
            max_operators_per_generation=4,
            use_structure_based_filtering=True,
            max_molecular_weight=400.0,
            filter_trivial_reactions=True,
            min_structural_change=0.2,
            prioritize_novel_structures=True,
            consider_stereochemistry=True
        )
    
    @staticmethod
    def combustion_chemistry() -> NetworkGenerationConfig:
        """Configuration optimized for combustion chemistry."""
        return NetworkGenerationConfig(
            max_generations=8,
            max_species=300,
            energy_cutoff=70.0,
            generation_strategy=NetworkGenerationStrategy.EXHAUSTIVE,
            operator_selection_strategy=OperatorSelectionStrategy.ALL,
            driving_force_threshold=0.5,  # High temperature favors many reactions
            use_structure_based_filtering=False,  # Allow all structures
            max_molecular_weight=200.0,  # Smaller fragments
            allow_unstable_intermediates=True,
            bias_towards_stable_products=False,
            parallel_screening=True,
            max_screening_workers=8
        )
    
    @staticmethod
    def catalysis_screening() -> NetworkGenerationConfig:
        """Configuration optimized for catalysis screening."""
        return NetworkGenerationConfig(
            max_generations=4,
            max_species=80,
            energy_cutoff=40.0,
            activation_energy_cutoff=60.0,
            generation_strategy=NetworkGenerationStrategy.FOCUSED,
            operator_selection_strategy=OperatorSelectionStrategy.THRESHOLD,
            operator_priority_threshold=0.6,
            driving_force_threshold=0.3,
            max_operators_per_generation=3,
            use_structure_based_filtering=True,
            filter_trivial_reactions=True,
            min_structural_change=0.15,
            prioritize_novel_structures=False,  # Focus on known pathways
            bias_towards_stable_products=True
        )
    
    @staticmethod
    def exploratory() -> NetworkGenerationConfig:
        """Configuration for exploratory network generation."""
        return NetworkGenerationConfig(
            max_generations=10,
            max_species=500,
            energy_cutoff=80.0,
            generation_strategy=NetworkGenerationStrategy.EXHAUSTIVE,
            operator_selection_strategy=OperatorSelectionStrategy.ALL,
            driving_force_threshold=0.1,  # Very permissive
            use_structure_based_filtering=False,
            allow_unstable_intermediates=True,
            filter_duplicate_reactions=True,
            filter_trivial_reactions=False,
            prioritize_novel_structures=True,
            parallel_screening=True,
            verbose_logging=True,
            save_intermediate_results=True
        )
    
    @staticmethod
    def fast_screening() -> NetworkGenerationConfig:
        """Configuration for fast screening with limited depth."""
        return NetworkGenerationConfig(
            max_generations=3,
            max_species=50,
            max_time=300.0,  # 5 minutes
            energy_cutoff=35.0,
            generation_strategy=NetworkGenerationStrategy.INTELLIGENT,
            operator_selection_strategy=OperatorSelectionStrategy.TOP_N,
            max_operators_per_generation=2,
            driving_force_threshold=0.5,
            use_structure_based_filtering=True,
            filter_trivial_reactions=True,
            min_structural_change=0.3,
            parallel_screening=True,
            cache_molecular_properties=True
        )


def get_preset_config(preset_name: str) -> NetworkGenerationConfig:
    """
    Get a preset configuration by name.
    
    Args:
        preset_name: Name of the preset
        
    Returns:
        Configuration object
        
    Raises:
        ValueError: If preset name is not recognized
    """
    presets = {
        'battery_chemistry': NetworkGenerationPresets.battery_chemistry,
        'organic_synthesis': NetworkGenerationPresets.organic_synthesis,
        'combustion_chemistry': NetworkGenerationPresets.combustion_chemistry,
        'catalysis_screening': NetworkGenerationPresets.catalysis_screening,
        'exploratory': NetworkGenerationPresets.exploratory,
        'fast_screening': NetworkGenerationPresets.fast_screening
    }
    
    if preset_name not in presets:
        available = ', '.join(presets.keys())
        raise ValueError(f"Unknown preset '{preset_name}'. Available presets: {available}")
    
    return presets[preset_name]()
