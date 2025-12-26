"""
Environment class for RSNet.

This module provides the Environment class that represents reaction conditions
such as temperature, pressure, and solvent, along with physical driving forces.
"""

from typing import Optional, Dict, Any, List, Set
import yaml
import math


class Environment:
    """
    Represents the reaction environment conditions.
    
    This class encapsulates environmental parameters that affect reaction
    thermodynamics and kinetics, such as temperature, pressure, and solvent.
    """
    
    def __init__(
        self,
        temperature: float = 298.15,
        pressure: float = 1.0,
        solvent: str = "gas",
        # Electrochemical parameters
        electrode_type: Optional[str] = None,  # "cathode", "anode", or None
        voltage: Optional[float] = None,  # V vs Li/Li+
        li_activity: float = 1.0,  # Li+ activity
        interface_type: str = "bulk",  # "bulk", "SEI", "CEI"
        # Additional parameters
        **kwargs
    ):
        """
        Initialize reaction environment.

        Args:
            temperature: Temperature in Kelvin (default: 298.15 K)
            pressure: Pressure in atm (default: 1.0 atm)
            solvent: Solvent name (default: "gas" for gas phase)
            electrode_type: Type of electrode ("cathode", "anode", or None)
            voltage: Electrode potential in V vs Li/Li+
            li_activity: Lithium ion activity (default: 1.0)
            interface_type: Interface type ("bulk", "SEI", "CEI")
            **kwargs: Additional environment parameters
        """
        self.temperature = temperature
        self.pressure = pressure
        self.solvent = solvent
        self.electrode_type = electrode_type
        self.voltage = voltage
        self.li_activity = li_activity
        self.interface_type = interface_type
        self.extra_params = kwargs
    
    @classmethod
    def from_config(cls, config: Dict[str, Any]) -> 'Environment':
        """
        Create Environment from configuration dictionary.
        
        Args:
            config: Configuration dictionary
            
        Returns:
            Environment object
        """
        env_config = config.get('environment', {})
        return cls(
            temperature=env_config.get('temperature', 298.15),
            pressure=env_config.get('pressure', 1.0),
            solvent=env_config.get('solvent', 'gas'),
            **{k: v for k, v in env_config.items() 
               if k not in ['temperature', 'pressure', 'solvent']}
        )
    
    @classmethod
    def from_yaml(cls, filepath: str) -> 'Environment':
        """
        Create Environment from YAML configuration file.
        
        Args:
            filepath: Path to YAML configuration file
            
        Returns:
            Environment object
        """
        with open(filepath, 'r') as f:
            config = yaml.safe_load(f)
        return cls.from_config(config)
    
    @property
    def is_gas_phase(self) -> bool:
        """Check if this is a gas phase environment."""
        return self.solvent.lower() == "gas"
    
    @property
    def is_solution(self) -> bool:
        """Check if this is a solution phase environment."""
        return not self.is_gas_phase
    
    def get_thermal_energy(self) -> float:
        """
        Get thermal energy kT in kcal/mol.

        Returns:
            Thermal energy in kcal/mol
        """
        # kT = k_B * T, where k_B = 1.987e-3 kcal/(mol·K)
        return 1.987e-3 * self.temperature

    def get_active_drives(self) -> Dict[str, bool]:
        """
        Determine active physical driving forces based on environment conditions.

        Returns:
            Dictionary of driving force names and their activation status
        """
        drives = {}

        # Thermal driving forces
        drives['thermal'] = self.temperature > 298.15  # Above room temperature
        drives['high_temperature'] = self.temperature > 500.0  # High temperature reactions
        drives['cryogenic'] = self.temperature < 200.0  # Low temperature effects

        # Electrochemical driving forces
        if self.electrode_type and self.voltage is not None:
            drives['oxidation'] = self._is_oxidizing_condition()
            drives['reduction'] = self._is_reducing_condition()
            drives['li_coordination'] = self.li_activity > 0.1  # Li+ coordination possible
            drives['electrochemical'] = True
        else:
            drives['oxidation'] = False
            drives['reduction'] = False
            drives['li_coordination'] = False
            drives['electrochemical'] = False

        # Interface-specific driving forces
        drives['surface_reaction'] = self.interface_type in ['SEI', 'CEI']
        drives['bulk_reaction'] = self.interface_type == 'bulk'
        drives['sei_formation'] = self.interface_type == 'SEI'
        drives['cei_formation'] = self.interface_type == 'CEI'

        # Pressure-related driving forces
        drives['high_pressure'] = self.pressure > 10.0  # High pressure effects
        drives['gas_phase'] = self.is_gas_phase
        drives['solution_phase'] = self.is_solution

        # Radical environment (estimated from conditions)
        drives['radical_environment'] = self._estimate_radical_environment()

        return drives

    def _is_oxidizing_condition(self) -> bool:
        """Check if environment provides oxidizing conditions."""
        if self.voltage is None:
            return False

        # High voltage vs Li/Li+ indicates oxidizing conditions
        if self.electrode_type == "cathode":
            return self.voltage > 3.5  # Above typical cathode potentials
        elif self.electrode_type == "anode":
            return self.voltage > 1.0  # Above typical anode potentials
        return False

    def _is_reducing_condition(self) -> bool:
        """Check if environment provides reducing conditions."""
        if self.voltage is None:
            return False

        # Low voltage vs Li/Li+ indicates reducing conditions
        if self.electrode_type == "anode":
            return self.voltage < 1.0  # Below typical anode potentials
        elif self.electrode_type == "cathode":
            return self.voltage < 3.0  # Below typical cathode potentials
        return False

    def _estimate_radical_environment(self) -> bool:
        """Estimate if radical reactions are likely based on conditions."""
        # High temperature promotes radical formation
        thermal_radicals = self.temperature > 400.0

        # Electrochemical conditions can generate radicals
        electrochemical_radicals = (
            self.electrode_type is not None and
            self.voltage is not None and
            (self.voltage > 4.0 or self.voltage < 0.5)
        )

        # Interface reactions often involve radicals
        interface_radicals = self.interface_type in ['SEI', 'CEI']

        return thermal_radicals or electrochemical_radicals or interface_radicals
    
    def to_dict(self) -> Dict[str, Any]:
        """
        Convert environment to dictionary.

        Returns:
            Dictionary representation
        """
        result = {
            'temperature': self.temperature,
            'pressure': self.pressure,
            'solvent': self.solvent,
            'electrode_type': self.electrode_type,
            'voltage': self.voltage,
            'li_activity': self.li_activity,
            'interface_type': self.interface_type,
        }
        result.update(self.extra_params)
        return result
    
    def copy(self) -> 'Environment':
        """
        Create a copy of this environment.

        Returns:
            New Environment object
        """
        return Environment(
            temperature=self.temperature,
            pressure=self.pressure,
            solvent=self.solvent,
            electrode_type=self.electrode_type,
            voltage=self.voltage,
            li_activity=self.li_activity,
            interface_type=self.interface_type,
            **self.extra_params
        )
    
    def __eq__(self, other) -> bool:
        """Check if two environments are equal."""
        if not isinstance(other, Environment):
            return False
        return (
            self.temperature == other.temperature and
            self.pressure == other.pressure and
            self.solvent == other.solvent and
            self.electrode_type == other.electrode_type and
            self.voltage == other.voltage and
            self.li_activity == other.li_activity and
            self.interface_type == other.interface_type and
            self.extra_params == other.extra_params
        )
    
    def __str__(self) -> str:
        """String representation of the environment."""
        return (f"Environment(T={self.temperature:.1f}K, "
                f"P={self.pressure:.1f}atm, solvent='{self.solvent}')")
    
    def __repr__(self) -> str:
        """Detailed string representation."""
        return self.__str__()
