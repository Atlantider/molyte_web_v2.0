"""
Base operator implementation for RSNet.

This module provides the BaseOperator class that extends the core Operator
with common functionality for reaction operators.
"""

from typing import List, Dict, Any, Optional
from ..core.operator import Operator
from ..core.molecule import Molecule
from ..core.reaction import Reaction, ReactionTemplate
from ..core.environment import Environment


class BaseOperator(Operator):
    """
    Base implementation of reaction operator with common functionality.
    """
    
    def __init__(
        self, 
        name: str, 
        description: str, 
        config: Optional[Dict[str, Any]] = None
    ):
        """
        Initialize the base operator.
        
        Args:
            name: Operator name
            description: Operator description
            config: Configuration parameters
        """
        super().__init__(name, description, config)
        self.arity = 1  # Default to unimolecular
        self.reaction_type = "unspecified"
        self._templates = []
    
    def get_templates(self) -> List[ReactionTemplate]:
        """
        Get reaction templates for this operator.
        
        Returns:
            List of reaction templates
        """
        return self._templates
    
    def add_template(self, template: ReactionTemplate):
        """
        Add a reaction template to this operator.
        
        Args:
            template: Reaction template to add
        """
        self._templates.append(template)
    
    def apply(
        self, 
        molecules: List[Molecule], 
        environment: Environment
    ) -> List[Reaction]:
        """
        Apply this operator to generate reactions.
        
        This base implementation provides a framework that can be extended
        by specific operators.
        
        Args:
            molecules: Input molecules
            environment: Reaction environment
            
        Returns:
            List of generated reactions
        """
        if not self.is_applicable(molecules):
            return []
        
        reactions = []
        
        # Apply each template
        for template in self._templates:
            try:
                template_reactions = template.apply(molecules)
                reactions.extend(template_reactions)
            except NotImplementedError:
                # Template doesn't implement apply method yet
                continue
            except Exception as e:
                # Log error and continue with other templates
                print(f"Warning: Template {template.name} failed: {e}")
                continue
        
        # Filter reactions based on criteria
        filtered_reactions = self.filter_reactions(reactions)
        
        return filtered_reactions
    
    def is_applicable(self, molecules: List[Molecule]) -> bool:
        """
        Check if this operator can be applied to the given molecules.

        Args:
            molecules: Input molecules

        Returns:
            True if operator is applicable
        """
        if not super().is_applicable(molecules):
            return False

        # Additional checks can be implemented by subclasses
        return True

    def can_apply(self, molecules: List[Molecule], environment: Environment) -> bool:
        """
        Check if this operator can be applied to molecules in the given environment.

        Args:
            molecules: Input molecules
            environment: Reaction environment

        Returns:
            True if operator can be applied
        """
        return self.is_applicable(molecules)
