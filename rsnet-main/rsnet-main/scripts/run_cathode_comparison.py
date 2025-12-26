"""
Comparative Cathode Simulation Framework
=========================================
Compares network formation under different salt conditions at cathode.
"""

import logging
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from rdkit import Chem
from rsnet.core.molecule import Molecule
from rsnet.core.environment import Environment
from rsnet.network.evolution import NetworkEvolver
from rsnet.network.config import NetworkGenerationPresets
from rsnet.compute.reaction_screener import ReactionScreener
import json

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

class MockCalculator:
    """Mock calculator for quick testing."""
    def optimize(self, molecule, environment):
        import random
        n_atoms = molecule.rdkit_mol.GetNumAtoms()
        base_energy = n_atoms * -100.0
        charge = Chem.GetFormalCharge(molecule.rdkit_mol)
        charge_penalty = (charge ** 2) * 50.0
        return {
            'energy_kcal_mol': base_energy + charge_penalty + random.uniform(-10, 10),
            'success': True,
            'optimized_molecule': molecule
        }
    
    def single_point(self, molecule, environment):
        return self.optimize(molecule, environment)

def run_cathode_simulation(salt_type='LiPF6', output_prefix='cathode'):
    """
    Run cathode simulation with specified salt.
    
    Args:
        salt_type: 'LiPF6', 'LiTFSI', or 'mixed'
        output_prefix: Prefix for output files
    """
    
    print("="*60)
    print(f"CATHODE SIMULATION: {salt_type}")
    print("="*60)
    print()
    
    # Setup environment (cathode at 4.2V)
    env = Environment(
        temperature=298.0,
        pressure=1.0,
        solvent="mixture",
        voltage=4.2,  # Cathode voltage
        electrode_type="cathode",
        electrode_material="NMC"
    )
    
    # Setup initial molecules
    mols = []
    
    # Solvents
    mols.append(Molecule(Chem.MolFromSmiles("O=C1OCCO1"), name="EC"))
    mols.append(Molecule(Chem.MolFromSmiles("COC(=O)OC"), name="DMC"))
    mols.append(Molecule(Chem.MolFromSmiles("[Li+]"), name="Li+"))
    
    # Salts
    if salt_type == 'LiPF6' or salt_type == 'mixed':
        mols.append(Molecule(Chem.MolFromSmiles("F[P-](F)(F)(F)(F)F"), name="PF6-"))
    
    if salt_type == 'LiTFSI' or salt_type == 'mixed':
        # LiTFSI anion: (CF3SO2)2N-
        try:
            tfsi = Molecule(Chem.MolFromSmiles("FC(F)(F)S(=O)(=O)[N-]S(=O)(=O)C(F)(F)F"), name="TFSI-")
            mols.append(tfsi)
        except:
            logger.warning("Could not create TFSI-, using simplified version")
            mols.append(Molecule(Chem.MolFromSmiles("[N-]"), name="TFSI-"))
    
    print(f"Initial species: {[m.name for m in mols]}")
    print()
    
    # Configure simulation
    config = NetworkGenerationPresets.battery_chemistry()
    config.max_generations = 3  # Shorter for comparison
    config.max_species = 1000
    config.max_reactions_per_generation = 1000
    
    # Run simulation
    mock_calc = MockCalculator()
    screener = ReactionScreener(calculator=mock_calc)
    evolver = NetworkEvolver(environment=env, screener=screener, config=config)
    
    print(f"Running simulation (max {config.max_generations} generations)...")
    final_species = evolver.evolve(mols, max_generations=config.max_generations)
    
    # Save results
    summary = {
        'salt_type': salt_type,
        'voltage': 4.2,
        'electrode': 'cathode',
        'generations': config.max_generations,
        'species_count': len(final_species),
        'reaction_count': len(evolver.reactions),
        'species_list': [s.molecule.smiles for s in final_species],
        'reaction_list': [r.name for r in evolver.reactions]
    }
    
    summary_file = f'{output_prefix}_{salt_type}_summary.json'
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    print()
    print("="*60)
    print("RESULTS")
    print("="*60)
    print(f"Total species: {len(final_species)}")
    print(f"Total reactions: {len(evolver.reactions)}")
    print(f"Summary saved to: {summary_file}")
    print()
    
    return summary

def compare_simulations():
    """Run all three simulations and compare."""
    
    print("\n" + "="*60)
    print("COMPARATIVE CATHODE ANALYSIS")
    print("="*60)
    print()
    
    results = {}
    
    # Run simulations
    for salt_type in ['LiPF6', 'LiTFSI', 'mixed']:
        results[salt_type] = run_cathode_simulation(salt_type, 'cathode')
    
    # Compare
    print("\n" + "="*60)
    print("COMPARISON SUMMARY")
    print("="*60)
    print()
    
    print(f"{'Condition':<15} {'Species':<10} {'Reactions':<12} {'Difference'}")
    print("-" * 60)
    
    for salt_type, data in results.items():
        diff = ""
        if salt_type == 'LiPF6':
            base_species = data['species_count']
        else:
            diff = f"({data['species_count'] - base_species:+d} species)"
        
        print(f"{salt_type:<15} {data['species_count']:<10} {data['reaction_count']:<12} {diff}")
    
    print()
    print("Next step: Run visualize_comparison.py to create visual comparison")
    print()

if __name__ == "__main__":
    compare_simulations()
