"""
RSNet - Reaction Network Generation Engine

A Python package for automated generation and analysis of chemical reaction networks.
"""

__version__ = "0.1.0"
__author__ = "RSNet Development Team"

# Import main classes for easy access
from .core.molecule import Molecule
from .core.environment import Environment
from .core.reaction import Reaction, ReactionTemplate
from .core.operator import Operator

__all__ = [
    "Molecule",
    "Environment", 
    "Reaction",
    "ReactionTemplate",
    "Operator",
]
