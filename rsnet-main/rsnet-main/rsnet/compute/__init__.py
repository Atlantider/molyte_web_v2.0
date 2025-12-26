"""
Computational modules for RSNet.

Contains:
- XTB energy calculation
- Reaction screening
- CEI/SEI chemistry (cathode/anode interface)
- Electrode species injection
- Oxidation operators
"""

from .xtb_calculator import XTBCalculator
from .reaction_screener import ReactionScreener

# CEI Chemistry modules
from .cathode_materials import (
    CathodeMaterial,
    CATHODE_MATERIALS,
    get_cathode_material,
    list_available_materials
)
from .electrode_species_injector import (
    ElectrodeSpeciesInjector,
    auto_inject_electrode_species
)
from .interface_components import (
    InterfaceComponentRecognizer,
    SEI_COMPONENTS,
    CEI_COMPONENTS,
    recognize_sei_component,
    recognize_cei_component
)
from .oxidation_operators import (
    OxidationOperatorManager,
    OXIDATION_OPERATORS,
    get_oxidation_operators_for_cathode
)
from .cei_chemistry import (
    CEIChemistryIntegrator,
    integrate_cei_chemistry
)

__all__ = [
    # Core
    'XTBCalculator',
    'ReactionScreener',
    
    # CEI Chemistry
    'CEIChemistryIntegrator',
    'integrate_cei_chemistry',
    
    # Cathode Materials
    'CathodeMaterial',
    'CATHODE_MATERIALS',
    'get_cathode_material',
    'list_available_materials',
    
    # Electrode Species
    'ElectrodeSpeciesInjector',
    'auto_inject_electrode_species',
    
    # Interface Components
    'InterfaceComponentRecognizer',
    'SEI_COMPONENTS',
    'CEI_COMPONENTS',
    'recognize_sei_component',
    'recognize_cei_component',
    
    # Oxidation Operators
    'OxidationOperatorManager',
    'OXIDATION_OPERATORS',
    'get_oxidation_operators_for_cathode',
]

