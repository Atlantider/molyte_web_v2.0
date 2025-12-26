"""
Base Operator class for RSNet.

This module provides the base class for all reaction operators.
"""

from abc import ABC, abstractmethod
from typing import List, Dict, Any, Optional
from .molecule import Molecule
from .reaction import Reaction, ReactionTemplate
from .environment import Environment


class Operator(ABC):
    """
    Abstract base class for reaction operators.
    
    Operators implement specific reaction types (e.g., hydrogen transfer,
    bond breaking) and can generate reactions from input molecules.
    """
    
    def __init__(self, name: str, description: str, config: Optional[Dict[str, Any]] = None):
        """
        Initialize the operator.
        
        Args:
            name: Operator name
            description: Operator description
            config: Configuration parameters
        """
        self.name = name
        self.description = description
        self.config = config or {}
        self.enabled = self.config.get('enabled', True)
        self.energy_threshold = self.config.get('energy_threshold', 50.0)
    
    @abstractmethod
    def get_templates(self) -> List[ReactionTemplate]:
        """
        Get reaction templates for this operator.
        
        Returns:
            List of reaction templates
        """
        pass
    
    @abstractmethod
    def apply(
        self, 
        molecules: List[Molecule], 
        environment: Environment
    ) -> List[Reaction]:
        """
        Apply this operator to generate reactions.
        
        Args:
            molecules: Input molecules
            environment: Reaction environment
            
        Returns:
            List of generated reactions
        """
        pass
    
    def is_applicable(self, molecules: List[Molecule]) -> bool:
        """
        Check if this operator can be applied to the given molecules.
        
        Args:
            molecules: Input molecules
            
        Returns:
            True if operator is applicable
        """
        return self.enabled and len(molecules) > 0
    
    def filter_reactions(self, reactions: List[Reaction]) -> List[Reaction]:
        """
        Filter reactions based on energy thresholds and other criteria.
        
        Args:
            reactions: List of reactions to filter
            
        Returns:
            Filtered list of reactions
        """
        filtered = []
        for reaction in reactions:
            if (reaction.activation_energy is None or 
                reaction.activation_energy <= self.energy_threshold):
                filtered.append(reaction)
        return filtered
    
    def __str__(self) -> str:
        """String representation of the operator."""
        return f"Operator('{self.name}')"
    
    def __repr__(self) -> str:
        """Detailed string representation."""
        return f"Operator(name='{self.name}', enabled={self.enabled})"
