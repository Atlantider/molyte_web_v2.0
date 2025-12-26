"""
Reaction and ReactionTemplate classes for RSNet.

This module provides classes for representing chemical reactions and reaction templates.
"""

from typing import List, Optional, Dict, Any, Tuple
import uuid
from .molecule import Molecule


class Reaction:
    """
    Represents a chemical reaction with reactants, products, and energetics.
    """
    
    def __init__(
        self,
        reactants: List[Molecule],
        products: List[Molecule],
        reaction_energy: Optional[float] = None,
        activation_energy: Optional[float] = None,
        name: Optional[str] = None,
        operator_name: Optional[str] = None
    ):
        """
        Initialize a reaction.
        
        Args:
            reactants: List of reactant molecules
            products: List of product molecules
            reaction_energy: Reaction energy in kcal/mol (ΔE)
            activation_energy: Activation energy in kcal/mol (Ea)
            name: Optional name for the reaction
        """
        self.reactants = reactants
        self.products = products
        self.reaction_energy = reaction_energy
        self.activation_energy = activation_energy
        self.name = name or self._generate_name()
        self.operator_name = operator_name
        self.id = str(uuid.uuid4())
    
    def _generate_name(self) -> str:
        """Generate a default name for the reaction."""
        reactant_names = " + ".join([mol.name for mol in self.reactants])
        product_names = " + ".join([mol.name for mol in self.products])
        return f"{reactant_names} -> {product_names}"
    
    @property
    def is_thermodynamically_favorable(self) -> bool:
        """Check if reaction is thermodynamically favorable (ΔE < 0)."""
        return self.reaction_energy is not None and self.reaction_energy < 0
    
    @property
    def is_kinetically_accessible(self, threshold: float = 50.0) -> bool:
        """Check if reaction is kinetically accessible (Ea < threshold)."""
        return (self.activation_energy is not None and 
                self.activation_energy < threshold)
    
    def __str__(self) -> str:
        """String representation of the reaction."""
        return self.name
    
    def __repr__(self) -> str:
        """Detailed string representation."""
        return f"Reaction('{self.name}', ΔE={self.reaction_energy}, Ea={self.activation_energy})"


class ReactionTemplate:
    """
    Represents a reaction template that can be applied to generate reactions.
    """
    
    def __init__(
        self,
        name: str,
        description: str,
        pattern: Optional[str] = None
    ):
        """
        Initialize a reaction template.
        
        Args:
            name: Template name
            description: Template description
            pattern: SMARTS pattern (if applicable)
        """
        self.name = name
        self.description = description
        self.pattern = pattern
        self.id = str(uuid.uuid4())
    
    def apply(self, molecules: List[Molecule]) -> List[Reaction]:
        """
        Apply this template to generate reactions.
        
        Args:
            molecules: List of molecules to apply template to
            
        Returns:
            List of generated reactions
        """
        # This is a placeholder - will be implemented by specific operators
        raise NotImplementedError("Template application must be implemented by subclasses")
    
    def __str__(self) -> str:
        """String representation of the template."""
        return f"ReactionTemplate('{self.name}')"
    
    def __repr__(self) -> str:
        """Detailed string representation."""
        return f"ReactionTemplate(name='{self.name}', description='{self.description}')"
