"""
Configuration management for RSNet.

This module provides utilities for loading, validating, and managing
configuration files and parameters.
"""

import yaml
import json
from pathlib import Path
from typing import Dict, Any, Optional, Union, List
from dataclasses import dataclass, asdict
import os

from ..network.config import NetworkGenerationConfig, NetworkGenerationPresets


@dataclass
class RSNetConfig:
    """Main RSNet configuration class."""
    
    # Network generation
    network: NetworkGenerationConfig = None
    
    # Environment settings
    environment: Dict[str, Any] = None
    
    # Operator settings
    operators: Dict[str, Any] = None
    
    # Computational settings
    computation: Dict[str, Any] = None
    
    # I/O settings
    io: Dict[str, Any] = None
    
    # Visualization settings
    visualization: Dict[str, Any] = None
    
    def __post_init__(self):
        """Initialize default values."""
        if self.network is None:
            self.network = NetworkGenerationConfig()
        
        if self.environment is None:
            self.environment = {
                'temperature': 298.15,
                'pressure': 1.0,
                'solvent': None,
                'electrode_type': None,
                'voltage': None
            }
        
        if self.operators is None:
            self.operators = {
                'use_registry': True,
                'selection_strategy': 'adaptive',
                'max_operators_per_generation': 5,
                'operator_timeout': 30.0
            }
        
        if self.computation is None:
            self.computation = {
                'xtb_path': None,  # Auto-detect
                'parallel': False,
                'max_workers': 4,
                'memory_limit': '2GB',
                'temp_dir': None
            }
        
        if self.io is None:
            self.io = {
                'output_dir': './rsnet_output',
                'save_intermediates': True,
                'export_formats': ['json', 'csv'],
                'compression': False
            }
        
        if self.visualization is None:
            self.visualization = {
                'plot_network': True,
                'plot_energies': True,
                'plot_pathways': True,
                'figure_format': 'png',
                'figure_dpi': 300,
                'figure_size': [12, 8]
            }
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        result = {}
        
        # Handle network config specially
        if self.network:
            result['network'] = self.network.to_dict()
        
        # Add other fields
        for field in ['environment', 'operators', 'computation', 'io', 'visualization']:
            value = getattr(self, field)
            if value is not None:
                result[field] = value
        
        return result
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'RSNetConfig':
        """Create from dictionary."""
        config = cls()
        
        # Handle network config specially
        if 'network' in data:
            config.network = NetworkGenerationConfig.from_dict(data['network'])
        
        # Set other fields
        for field in ['environment', 'operators', 'computation', 'io', 'visualization']:
            if field in data:
                setattr(config, field, data[field])
        
        return config


def get_default_config() -> RSNetConfig:
    """Get default RSNet configuration."""
    return RSNetConfig()


def load_config(file_path: Union[str, Path]) -> RSNetConfig:
    """
    Load configuration from file.
    
    Args:
        file_path: Path to configuration file (YAML or JSON)
        
    Returns:
        RSNetConfig object
    """
    file_path = Path(file_path)
    
    if not file_path.exists():
        raise FileNotFoundError(f"Configuration file not found: {file_path}")
    
    with open(file_path, 'r') as f:
        if file_path.suffix.lower() in ['.yaml', '.yml']:
            data = yaml.safe_load(f)
        elif file_path.suffix.lower() == '.json':
            data = json.load(f)
        else:
            raise ValueError(f"Unsupported configuration file format: {file_path.suffix}")
    
    return RSNetConfig.from_dict(data)


def save_config(config: RSNetConfig, file_path: Union[str, Path]) -> None:
    """
    Save configuration to file.
    
    Args:
        config: RSNetConfig object
        file_path: Output file path
    """
    file_path = Path(file_path)
    data = config.to_dict()
    
    with open(file_path, 'w') as f:
        if file_path.suffix.lower() in ['.yaml', '.yml']:
            yaml.dump(data, f, default_flow_style=False, indent=2)
        elif file_path.suffix.lower() == '.json':
            json.dump(data, f, indent=2)
        else:
            raise ValueError(f"Unsupported configuration file format: {file_path.suffix}")


def validate_config(config: RSNetConfig) -> List[str]:
    """
    Validate configuration and return list of issues.
    
    Args:
        config: RSNetConfig object
        
    Returns:
        List of validation error messages
    """
    issues = []
    
    # Validate network config
    if config.network:
        if config.network.max_generations <= 0:
            issues.append("max_generations must be positive")
        
        if config.network.max_species <= 0:
            issues.append("max_species must be positive")
        
        if config.network.energy_cutoff <= 0:
            issues.append("energy_cutoff must be positive")
    
    # Validate environment
    if config.environment:
        temp = config.environment.get('temperature', 298.15)
        if temp <= 0:
            issues.append("temperature must be positive")
        
        pressure = config.environment.get('pressure', 1.0)
        if pressure <= 0:
            issues.append("pressure must be positive")
    
    # Validate computation settings
    if config.computation:
        max_workers = config.computation.get('max_workers', 4)
        if max_workers <= 0:
            issues.append("max_workers must be positive")
        
        # Check xTB path if specified
        xtb_path = config.computation.get('xtb_path')
        if xtb_path and not Path(xtb_path).exists():
            issues.append(f"xTB executable not found: {xtb_path}")
    
    # Validate I/O settings
    if config.io:
        output_dir = config.io.get('output_dir')
        if output_dir:
            try:
                Path(output_dir).mkdir(parents=True, exist_ok=True)
            except Exception as e:
                issues.append(f"Cannot create output directory {output_dir}: {e}")
    
    return issues


def create_config_template(file_path: Union[str, Path], 
                          preset: Optional[str] = None) -> None:
    """
    Create a configuration template file.
    
    Args:
        file_path: Output file path
        preset: Preset name (e.g., 'battery_chemistry', 'organic_synthesis')
    """
    if preset:
        # Use preset network configuration
        if hasattr(NetworkGenerationPresets, preset):
            network_config = getattr(NetworkGenerationPresets, preset)()
        else:
            raise ValueError(f"Unknown preset: {preset}")
        
        config = RSNetConfig(network=network_config)
    else:
        config = get_default_config()
    
    save_config(config, file_path)


def get_config_from_env() -> Dict[str, Any]:
    """
    Get configuration overrides from environment variables.
    
    Returns:
        Dictionary of configuration overrides
    """
    overrides = {}
    
    # Network settings
    if 'RSNET_MAX_GENERATIONS' in os.environ:
        overrides.setdefault('network', {})['max_generations'] = int(os.environ['RSNET_MAX_GENERATIONS'])
    
    if 'RSNET_MAX_SPECIES' in os.environ:
        overrides.setdefault('network', {})['max_species'] = int(os.environ['RSNET_MAX_SPECIES'])
    
    if 'RSNET_ENERGY_CUTOFF' in os.environ:
        overrides.setdefault('network', {})['energy_cutoff'] = float(os.environ['RSNET_ENERGY_CUTOFF'])
    
    # Environment settings
    if 'RSNET_TEMPERATURE' in os.environ:
        overrides.setdefault('environment', {})['temperature'] = float(os.environ['RSNET_TEMPERATURE'])
    
    if 'RSNET_PRESSURE' in os.environ:
        overrides.setdefault('environment', {})['pressure'] = float(os.environ['RSNET_PRESSURE'])
    
    # Computation settings
    if 'RSNET_XTB_PATH' in os.environ:
        overrides.setdefault('computation', {})['xtb_path'] = os.environ['RSNET_XTB_PATH']
    
    if 'RSNET_MAX_WORKERS' in os.environ:
        overrides.setdefault('computation', {})['max_workers'] = int(os.environ['RSNET_MAX_WORKERS'])
    
    # I/O settings
    if 'RSNET_OUTPUT_DIR' in os.environ:
        overrides.setdefault('io', {})['output_dir'] = os.environ['RSNET_OUTPUT_DIR']
    
    return overrides


def merge_configs(base_config: RSNetConfig, overrides: Dict[str, Any]) -> RSNetConfig:
    """
    Merge configuration overrides into base configuration.
    
    Args:
        base_config: Base RSNetConfig object
        overrides: Dictionary of overrides
        
    Returns:
        Merged RSNetConfig object
    """
    # Convert base config to dict
    base_dict = base_config.to_dict()
    
    # Deep merge overrides
    def deep_merge(base: Dict, override: Dict) -> Dict:
        result = base.copy()
        for key, value in override.items():
            if key in result and isinstance(result[key], dict) and isinstance(value, dict):
                result[key] = deep_merge(result[key], value)
            else:
                result[key] = value
        return result
    
    merged_dict = deep_merge(base_dict, overrides)
    
    return RSNetConfig.from_dict(merged_dict)


def load_config_with_overrides(file_path: Optional[Union[str, Path]] = None,
                              env_overrides: bool = True) -> RSNetConfig:
    """
    Load configuration with environment variable overrides.
    
    Args:
        file_path: Path to configuration file (optional)
        env_overrides: Whether to apply environment variable overrides
        
    Returns:
        RSNetConfig object
    """
    # Start with default config
    if file_path and Path(file_path).exists():
        config = load_config(file_path)
    else:
        config = get_default_config()
    
    # Apply environment overrides
    if env_overrides:
        env_config = get_config_from_env()
        if env_config:
            config = merge_configs(config, env_config)
    
    return config


def print_config_summary(config: RSNetConfig) -> None:
    """
    Print a summary of the configuration.
    
    Args:
        config: RSNetConfig object
    """
    print("RSNet Configuration Summary")
    print("=" * 30)
    
    # Network settings
    print(f"Network Generation:")
    print(f"  Max generations: {config.network.max_generations}")
    print(f"  Max species: {config.network.max_species}")
    print(f"  Energy cutoff: {config.network.energy_cutoff} kcal/mol")
    print(f"  Strategy: {config.network.generation_strategy.value}")
    
    # Environment
    print(f"\nEnvironment:")
    print(f"  Temperature: {config.environment['temperature']} K")
    print(f"  Pressure: {config.environment['pressure']} atm")
    
    # Operators
    print(f"\nOperators:")
    print(f"  Use registry: {config.operators['use_registry']}")
    print(f"  Selection strategy: {config.operators['selection_strategy']}")
    print(f"  Max operators/gen: {config.operators['max_operators_per_generation']}")
    
    # Computation
    print(f"\nComputation:")
    print(f"  Parallel: {config.computation['parallel']}")
    print(f"  Max workers: {config.computation['max_workers']}")
    
    # I/O
    print(f"\nI/O:")
    print(f"  Output directory: {config.io['output_dir']}")
    print(f"  Export formats: {config.io['export_formats']}")
    
    print()


# Preset configurations
PRESET_CONFIGS = {
    'battery_chemistry': lambda: RSNetConfig(
        network=NetworkGenerationPresets.battery_chemistry(),
        environment={
            'temperature': 500.0,
            'electrode_type': 'cathode',
            'voltage': 4.2,
            'li_activity': 0.5
        }
    ),
    
    'organic_synthesis': lambda: RSNetConfig(
        network=NetworkGenerationPresets.organic_synthesis(),
        environment={
            'temperature': 350.0,
            'solvent': 'water'
        }
    ),
    
    'fast_screening': lambda: RSNetConfig(
        network=NetworkGenerationPresets.fast_screening(),
        computation={
            'parallel': True,
            'max_workers': 8
        }
    )
}


def get_preset_config(preset_name: str) -> RSNetConfig:
    """
    Get a preset configuration.
    
    Args:
        preset_name: Name of the preset
        
    Returns:
        RSNetConfig object
    """
    if preset_name not in PRESET_CONFIGS:
        available = list(PRESET_CONFIGS.keys())
        raise ValueError(f"Unknown preset '{preset_name}'. Available: {available}")
    
    return PRESET_CONFIGS[preset_name]()
