"""
Physical driving forces evaluation for RSNet.

This module provides functions to evaluate different types of driving forces
that can activate chemical reactions in various environments.
"""

from typing import Dict, List, Optional, Tuple
import math
from ..core.environment import Environment
from ..core.molecule import Molecule


class DrivingForceAnalyzer:
    """驱动力分析器 - 分析环境条件并确定激活的驱动力"""

    def __init__(self):
        """初始化驱动力分析器"""
        pass

    def analyze_environment(self, environment: Environment) -> Dict[str, bool]:
        """
        分析环境条件，返回激活的驱动力

        Args:
            environment: 反应环境

        Returns:
            激活的驱动力字典 {驱动力名称: 是否激活}
        """
        return environment.get_active_drives()


def get_driving_forces(molecules: List[Molecule], environment: Environment) -> Dict[str, float]:
    """
    Evaluate all driving forces for a set of molecules in an environment.
    
    Args:
        molecules: List of molecules to evaluate
        environment: Reaction environment
        
    Returns:
        Dictionary of driving force names and their strengths (0.0 to 1.0)
    """
    forces = {}
    
    # Get environment-based driving forces
    env_forces = environment.get_active_drives()
    
    # Convert boolean drives to strength values
    for drive_name, is_active in env_forces.items():
        if is_active:
            forces[drive_name] = evaluate_drive_strength(drive_name, environment)
        else:
            forces[drive_name] = 0.0
    
    # Evaluate molecule-specific driving forces
    for mol in molecules:
        mol_forces = evaluate_molecular_driving_forces(mol, environment)
        for force_name, strength in mol_forces.items():
            # Take maximum strength if multiple molecules contribute
            forces[force_name] = max(forces.get(force_name, 0.0), strength)
    
    return forces


def evaluate_drive_strength(drive_name: str, environment: Environment) -> float:
    """
    Evaluate the strength of a specific driving force.
    
    Args:
        drive_name: Name of the driving force
        environment: Reaction environment
        
    Returns:
        Driving force strength (0.0 to 1.0)
    """
    if drive_name == 'thermal':
        return evaluate_thermal_driving_force(environment)
    elif drive_name == 'high_temperature':
        return evaluate_high_temperature_driving_force(environment)
    elif drive_name == 'oxidation':
        return evaluate_oxidation_driving_force(environment)
    elif drive_name == 'reduction':
        return evaluate_reduction_driving_force(environment)
    elif drive_name == 'electrochemical':
        return evaluate_electrochemical_driving_force(environment)
    elif drive_name == 'li_coordination':
        return evaluate_li_coordination_driving_force(environment)
    elif drive_name == 'surface_reaction':
        return evaluate_surface_reaction_driving_force(environment)
    elif drive_name == 'radical_environment':
        return evaluate_radical_environment_driving_force(environment)
    else:
        return 0.5  # Default moderate strength for unknown drives


def evaluate_thermal_driving_force(environment: Environment) -> float:
    """Evaluate thermal driving force strength."""
    T = environment.temperature
    
    # Sigmoid function to map temperature to driving force strength
    # Room temperature (298K) -> 0.0, High temperature (1000K) -> 1.0
    normalized_temp = (T - 298.15) / (1000.0 - 298.15)
    return max(0.0, min(1.0, normalized_temp))


def evaluate_high_temperature_driving_force(environment: Environment) -> float:
    """Evaluate high temperature driving force strength."""
    T = environment.temperature
    
    if T < 500.0:
        return 0.0
    elif T > 1200.0:
        return 1.0
    else:
        # Linear scaling between 500K and 1200K
        return (T - 500.0) / (1200.0 - 500.0)


def evaluate_oxidation_driving_force(environment: Environment) -> float:
    """Evaluate oxidation driving force strength."""
    if environment.voltage is None:
        return 0.0
    
    V = environment.voltage
    
    # Higher voltage = stronger oxidation driving force
    if environment.electrode_type == "cathode":
        # Cathode: high voltage vs Li/Li+ indicates strong oxidation
        if V > 4.5:
            return 1.0
        elif V > 3.5:
            return (V - 3.5) / (4.5 - 3.5)
        else:
            return 0.0
    elif environment.electrode_type == "anode":
        # Anode: voltage above 1V indicates oxidation
        if V > 2.0:
            return 1.0
        elif V > 1.0:
            return (V - 1.0) / (2.0 - 1.0)
        else:
            return 0.0
    
    return 0.0


def evaluate_reduction_driving_force(environment: Environment) -> float:
    """Evaluate reduction driving force strength."""
    if environment.voltage is None:
        return 0.0
    
    V = environment.voltage
    
    # Lower voltage = stronger reduction driving force
    if environment.electrode_type == "anode":
        # Anode: low voltage vs Li/Li+ indicates strong reduction
        if V < 0.5:
            return 1.0
        elif V < 1.0:
            return (1.0 - V) / (1.0 - 0.5)
        else:
            return 0.0
    elif environment.electrode_type == "cathode":
        # Cathode: voltage below 3V indicates reduction
        if V < 2.0:
            return 1.0
        elif V < 3.0:
            return (3.0 - V) / (3.0 - 2.0)
        else:
            return 0.0
    
    return 0.0


def evaluate_electrochemical_driving_force(environment: Environment) -> float:
    """Evaluate general electrochemical driving force strength."""
    if environment.electrode_type is None or environment.voltage is None:
        return 0.0
    
    # Electrochemical driving force is stronger at extreme voltages
    V = environment.voltage
    
    # Distance from neutral voltage (around 2.5V for typical battery systems)
    neutral_voltage = 2.5
    voltage_deviation = abs(V - neutral_voltage)
    
    # Maximum deviation of 2V gives full strength
    return min(1.0, voltage_deviation / 2.0)


def evaluate_li_coordination_driving_force(environment: Environment) -> float:
    """Evaluate lithium coordination driving force strength."""
    activity = environment.li_activity
    
    # Higher Li+ activity = stronger coordination driving force
    # Logarithmic scaling to handle wide range of activities
    if activity <= 0.001:
        return 0.0
    elif activity >= 1.0:
        return 1.0
    else:
        # Logarithmic scaling
        return (math.log10(activity) + 3) / 3  # -3 to 0 -> 0 to 1


def evaluate_surface_reaction_driving_force(environment: Environment) -> float:
    """Evaluate surface reaction driving force strength."""
    if environment.interface_type == "bulk":
        return 0.0
    elif environment.interface_type in ["SEI", "CEI"]:
        return 0.8  # Strong driving force at interfaces
    else:
        return 0.3  # Moderate for other interface types


def evaluate_radical_environment_driving_force(environment: Environment) -> float:
    """Evaluate radical environment driving force strength."""
    strength = 0.0
    
    # Temperature contribution
    if environment.temperature > 400.0:
        temp_contribution = min(1.0, (environment.temperature - 400.0) / 600.0)
        strength = max(strength, temp_contribution)
    
    # Electrochemical contribution
    if environment.voltage is not None:
        if environment.voltage > 4.0 or environment.voltage < 0.5:
            electrochemical_contribution = 0.7
            strength = max(strength, electrochemical_contribution)
    
    # Interface contribution
    if environment.interface_type in ["SEI", "CEI"]:
        interface_contribution = 0.6
        strength = max(strength, interface_contribution)
    
    return strength


def evaluate_molecular_driving_forces(molecule: Molecule, environment: Environment) -> Dict[str, float]:
    """
    Evaluate molecule-specific driving forces.
    
    Args:
        molecule: Molecule to evaluate
        environment: Reaction environment
        
    Returns:
        Dictionary of driving force names and strengths
    """
    forces = {}
    
    # Get molecular structure tags
    try:
        from .structure_tags import get_structure_tags
        tags = get_structure_tags(molecule)
    except ImportError:
        # Fallback if structure_tags not available
        tags = {}
    
    # Ring strain driving force
    if 'small_rings' in tags:
        ring_strain = 0.0
        for ring in tags['small_rings']:
            if ring['size'] <= 4:
                ring_strain = max(ring_strain, ring['strain_energy'] / 30.0)  # Normalize
        forces['ring_strain'] = min(1.0, ring_strain)
    
    # Polar bond driving force
    if 'polar_bonds' in tags:
        polar_strength = 0.0
        for bond in tags['polar_bonds']:
            if bond['electronegativity_diff'] > 1.0:
                polar_strength = max(polar_strength, bond['electronegativity_diff'] / 3.0)
        forces['polar_bond_cleavage'] = min(1.0, polar_strength)
    
    # Weak bond driving force
    if 'weak_bonds' in tags:
        weak_bond_strength = len(tags['weak_bonds']) / 5.0  # Normalize by typical number
        forces['weak_bond_cleavage'] = min(1.0, weak_bond_strength)
    
    # π-system driving force
    if 'pi_systems' in tags:
        pi_strength = len(tags['pi_systems']) / 3.0  # Normalize
        forces['pi_system_reaction'] = min(1.0, pi_strength)
    
    # Functional group driving force
    if 'functional_groups' in tags:
        fg_strength = len(tags['functional_groups']) / 5.0  # Normalize
        forces['functional_group_reaction'] = min(1.0, fg_strength)
    
    return forces


def evaluate_thermodynamic_driving_force(
    reactants: List[Molecule], 
    products: List[Molecule], 
    environment: Environment
) -> float:
    """
    Evaluate thermodynamic driving force for a reaction.
    
    Args:
        reactants: List of reactant molecules
        products: List of product molecules
        environment: Reaction environment
        
    Returns:
        Thermodynamic driving force strength (0.0 to 1.0)
    """
    # This would require energy calculations
    # For now, return a placeholder based on environment
    
    # Higher temperature generally increases thermodynamic driving force
    temp_factor = min(1.0, environment.temperature / 1000.0)
    
    # Electrochemical environments can provide additional driving force
    if environment.electrode_type is not None:
        electrochemical_factor = 0.3
    else:
        electrochemical_factor = 0.0
    
    return min(1.0, temp_factor + electrochemical_factor)


def evaluate_kinetic_driving_force(
    reactants: List[Molecule], 
    products: List[Molecule], 
    environment: Environment
) -> float:
    """
    Evaluate kinetic driving force for a reaction.
    
    Args:
        reactants: List of reactant molecules
        products: List of product molecules
        environment: Reaction environment
        
    Returns:
        Kinetic driving force strength (0.0 to 1.0)
    """
    # Temperature strongly affects kinetics (Arrhenius equation)
    T = environment.temperature
    
    # Exponential dependence on temperature
    # Room temperature (298K) -> 0.1, High temperature (1000K) -> 1.0
    kinetic_strength = 1.0 - math.exp(-(T - 298.15) / 200.0)
    
    return max(0.0, min(1.0, kinetic_strength))


def combine_driving_forces(forces: Dict[str, float], weights: Optional[Dict[str, float]] = None) -> float:
    """
    Combine multiple driving forces into a single strength value.
    
    Args:
        forces: Dictionary of driving force strengths
        weights: Optional weights for different forces
        
    Returns:
        Combined driving force strength (0.0 to 1.0)
    """
    if not forces:
        return 0.0
    
    if weights is None:
        # Default equal weights
        weights = {name: 1.0 for name in forces.keys()}
    
    # Weighted average
    total_weight = 0.0
    weighted_sum = 0.0
    
    for force_name, strength in forces.items():
        weight = weights.get(force_name, 1.0)
        weighted_sum += strength * weight
        total_weight += weight
    
    if total_weight == 0.0:
        return 0.0
    
    return weighted_sum / total_weight
