"""
Reaction screening using quantum chemical calculations.

This module provides tools to screen reactions based on thermodynamic
and kinetic criteria using xTB calculations.
"""

import logging
from typing import Dict, List, Optional, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed
import time

from ..core.reaction import Reaction
from ..core.molecule import Molecule
from ..core.environment import Environment
from .xtb_calculator import XTBCalculator


logger = logging.getLogger(__name__)


class ReactionScreener:
    """
    Screen reactions using quantum chemical calculations.
    
    This class evaluates reactions based on:
    - Thermodynamic feasibility (ΔG)
    - Kinetic accessibility (estimated barriers)
    - Stability of products
    """
    
    def __init__(self, 
                 calculator: Optional[XTBCalculator] = None,
                 max_workers: int = 4,
                 energy_threshold: float = 50.0,  # kcal/mol
                 optimize_geometries: bool = True):
        """
        Initialize reaction screener.
        
        Args:
            calculator: xTB calculator instance
            max_workers: Maximum number of parallel calculations
            energy_threshold: Maximum allowed reaction energy (kcal/mol)
            optimize_geometries: Whether to optimize geometries before energy calculation
        """
        self.calculator = calculator or XTBCalculator()
        self.max_workers = max_workers
        self.energy_threshold = energy_threshold
        self.optimize_geometries = optimize_geometries
        
        # Cache for calculated energies
        self._energy_cache = {}
    
    def _get_molecule_key(self, molecule: Molecule) -> str:
        """Generate cache key for molecule."""
        return f"{molecule.smiles}_{molecule.formula}_{hash(molecule.smiles)}"
    
    def _calculate_molecule_energy(self, molecule: Molecule, environment: Environment) -> Optional[float]:
        """
        Calculate energy for a single molecule with caching.
        
        Args:
            molecule: Molecule to calculate
            environment: Environmental conditions
            
        Returns:
            Energy in kcal/mol, or None if calculation failed
        """
        # Check cache first
        cache_key = f"{self._get_molecule_key(molecule)}_{environment.temperature}_{environment.solvent}"
        
        if cache_key in self._energy_cache:
            return self._energy_cache[cache_key]
        
        try:
            if self.optimize_geometries:
                result = self.calculator.optimize(molecule, environment)
                energy = result.get('energy_kcal_mol')
                
                # Store optimized molecule if available
                if result.get('optimized_molecule'):
                    # Update the molecule with optimized geometry
                    # (In a full implementation, we might want to return the optimized molecule)
                    pass
            else:
                result = self.calculator.single_point(molecule, environment)
                energy = result.get('energy_kcal_mol')
            
            # Cache the result
            if energy is not None:
                self._energy_cache[cache_key] = energy
            
            return energy
            
        except Exception as e:
            logger.warning(f"Energy calculation failed for {molecule.name}: {e}")
            return None
    
    def calculate_reaction_energy(self, reaction: Reaction, environment: Environment) -> Dict:
        """
        Calculate reaction energy (ΔE) for a reaction.
        
        Args:
            reaction: Reaction to evaluate
            environment: Environmental conditions
            
        Returns:
            Dictionary with reaction energetics
        """
        # Calculate energies for reactants
        reactant_energies = []
        for reactant in reaction.reactants:
            energy = self._calculate_molecule_energy(reactant, environment)
            if energy is None:
                return {
                    'success': False,
                    'error': f'Failed to calculate energy for reactant {reactant.name}'
                }
            reactant_energies.append(energy)
        
        # Calculate energies for products
        product_energies = []
        for product in reaction.products:
            energy = self._calculate_molecule_energy(product, environment)
            if energy is None:
                return {
                    'success': False,
                    'error': f'Failed to calculate energy for product {product.name}'
                }
            product_energies.append(energy)
        
        # Calculate reaction energy
        total_reactant_energy = sum(reactant_energies)
        total_product_energy = sum(product_energies)
        reaction_energy = total_product_energy - total_reactant_energy
        
        return {
            'success': True,
            'reaction_energy': reaction_energy,  # kcal/mol
            'reactant_energies': reactant_energies,
            'product_energies': product_energies,
            'total_reactant_energy': total_reactant_energy,
            'total_product_energy': total_product_energy,
            'exothermic': reaction_energy < 0,
            'endothermic': reaction_energy > 0
        }
    
    def estimate_activation_energy(self, reaction: Reaction, environment: Environment) -> Optional[float]:
        """
        Estimate activation energy for a reaction.
        
        This is a simplified estimation based on:
        - Bond breaking/forming energies
        - Reaction energy
        - Empirical correlations
        
        Args:
            reaction: Reaction to evaluate
            environment: Environmental conditions
            
        Returns:
            Estimated activation energy in kcal/mol, or None if estimation failed
        """
        # Get reaction energy
        energetics = self.calculate_reaction_energy(reaction, environment)
        
        if not energetics['success']:
            return None
        
        reaction_energy = energetics['reaction_energy']
        
        # 改进的活化能估算 - 基于反应类型的Bell-Evans-Polanyi参数
        # 参考: Bell-Evans-Polanyi原理, Marcus理论, Hammond假说
        #
        # 反应类型 -> (α, β) 参数
        # α: Brønsted系数, 描述过渡态与产物的相似程度
        # β: 本征势垒 (kcal/mol)
        
        reaction_name = getattr(reaction, 'name', '').lower()
        
        # 根据反应类型选择BEP参数
        if 'electron' in reaction_name or 'redox' in reaction_name:
            # 电子转移 - Marcus理论
            alpha, beta = 0.35, 5.0
        elif 'h_transfer' in reaction_name or 'hydrogen' in reaction_name:
            # 氢转移
            alpha, beta = 0.45, 12.0
        elif 'ring_open' in reaction_name:
            # 环开裂 - 考虑环应变
            alpha, beta = 0.35, 18.0
        elif 'li' in reaction_name:
            # 锂相关反应 - 低势垒
            alpha, beta = 0.25, 5.0
        elif 'radical' in reaction_name:
            # 自由基反应
            alpha, beta = 0.30, 20.0
        elif 'polymerization' in reaction_name:
            # 聚合反应
            alpha, beta = 0.35, 8.0
        else:
            # 默认值 - 保守估计
            alpha, beta = 0.45, 15.0
        
        # 计算活化能 (考虑Hammond假说)
        if reaction_energy > 0:  # 吸能反应
            # 过渡态更接近产物
            activation_energy = alpha * reaction_energy + beta
        else:  # 放能反应
            # 过渡态更接近反应物 (Hammond假说)
            activation_energy = beta - (1 - alpha) * abs(reaction_energy)
            activation_energy = max(activation_energy, 2.0)  # 最小势垒2 kcal/mol
        
        # 电化学校正 (如果有电极信息)
        voltage = getattr(reaction, 'voltage', 0.0)
        if voltage != 0:
            # Butler-Volmer: Ea = Ea0 - 0.5 * n * F * η
            F = 23.06  # kcal/mol/V
            correction = 0.5 * 1 * F * abs(voltage)
            activation_energy = max(2.0, activation_energy - correction)
        
        return max(0, activation_energy)  # Ensure non-negative
    
    def screen_reaction(self, reaction: Reaction, environment: Environment) -> Dict:
        """
        Comprehensive screening of a single reaction.
        
        Args:
            reaction: Reaction to screen
            environment: Environmental conditions
            
        Returns:
            Dictionary with screening results
        """
        start_time = time.time()
        
        # Calculate reaction energy
        energetics = self.calculate_reaction_energy(reaction, environment)
        
        if not energetics['success']:
            return {
                'reaction': reaction,
                'success': False,
                'error': energetics.get('error', 'Unknown error'),
                'calculation_time': time.time() - start_time
            }
        
        reaction_energy = energetics['reaction_energy']
        
        # Estimate activation energy
        activation_energy = self.estimate_activation_energy(reaction, environment)
        
        # Apply screening criteria
        thermodynamically_feasible = abs(reaction_energy) <= self.energy_threshold
        kinetically_accessible = activation_energy is None or activation_energy <= self.energy_threshold
        
        # Overall feasibility
        feasible = thermodynamically_feasible and kinetically_accessible
        
        # Calculate rate constant estimate (Arrhenius equation)
        rate_constant = None
        if activation_energy is not None:
            # k = A * exp(-Ea/RT)
            # Using typical pre-exponential factor A = 10^13 s^-1
            import math
            R = 1.987e-3  # kcal/(mol·K)
            T = environment.temperature
            rate_constant = 1e13 * math.exp(-activation_energy / (R * T))
        
        # Update reaction object with calculated values
        reaction.reaction_energy = reaction_energy
        reaction.activation_energy = activation_energy
        
        return {
            'reaction': reaction,
            'success': True,
            'feasible': feasible,
            'thermodynamically_feasible': thermodynamically_feasible,
            'kinetically_accessible': kinetically_accessible,
            'reaction_energy': reaction_energy,
            'activation_energy': activation_energy,
            'rate_constant': rate_constant,
            'exothermic': energetics['exothermic'],
            'endothermic': energetics['endothermic'],
            'reactant_energies': energetics['reactant_energies'],
            'product_energies': energetics['product_energies'],
            'calculation_time': time.time() - start_time
        }
    
    def screen_reactions(self, reactions: List[Reaction], environment: Environment) -> List[Dict]:
        """
        Screen multiple reactions in parallel.
        
        Args:
            reactions: List of reactions to screen
            environment: Environmental conditions
            
        Returns:
            List of screening results
        """
        if not reactions:
            return []
        
        logger.info(f"Screening {len(reactions)} reactions with {self.max_workers} workers")
        
        results = []
        
        if self.max_workers == 1:
            # Sequential processing
            for reaction in reactions:
                result = self.screen_reaction(reaction, environment)
                results.append(result)
        else:
            # Parallel processing
            with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
                # Submit all tasks
                future_to_reaction = {
                    executor.submit(self.screen_reaction, reaction, environment): reaction
                    for reaction in reactions
                }
                
                # Collect results as they complete
                for future in as_completed(future_to_reaction):
                    reaction = future_to_reaction[future]
                    try:
                        result = future.result()
                        results.append(result)
                    except Exception as e:
                        logger.error(f"Error screening reaction {reaction.name}: {e}")
                        results.append({
                            'reaction': reaction,
                            'success': False,
                            'error': str(e),
                            'calculation_time': 0
                        })
        
        # Sort results by feasibility and reaction energy
        results.sort(key=lambda x: (
            not x.get('feasible', False),  # Feasible reactions first
            abs(x.get('reaction_energy', float('inf')))  # Then by energy magnitude
        ))
        
        return results
    
    def filter_feasible_reactions(self, reactions: List[Reaction], environment: Environment) -> List[Reaction]:
        """
        Filter reactions to keep only feasible ones.
        
        Args:
            reactions: List of reactions to filter
            environment: Environmental conditions
            
        Returns:
            List of feasible reactions
        """
        screening_results = self.screen_reactions(reactions, environment)
        
        feasible_reactions = []
        for result in screening_results:
            if result.get('success', False) and result.get('feasible', False):
                feasible_reactions.append(result['reaction'])
        
        logger.info(f"Filtered {len(feasible_reactions)} feasible reactions from {len(reactions)} total")
        
        return feasible_reactions
    
    def get_energy_statistics(self) -> Dict:
        """Get statistics about cached energy calculations."""
        if not self._energy_cache:
            return {'cache_size': 0}
        
        energies = list(self._energy_cache.values())
        
        return {
            'cache_size': len(self._energy_cache),
            'min_energy': min(energies),
            'max_energy': max(energies),
            'mean_energy': sum(energies) / len(energies)
        }
