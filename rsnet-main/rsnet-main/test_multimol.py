"""
Multi-Molecular Clustering Test
===============================
Verify that the system can generate solvation shells (Li+ + 4EC -> Li(EC)4)
using the new Anchor-Based Clustering strategy.
"""

import logging
from rsnet.core.molecule import Molecule
from rsnet.core.environment import Environment
from rsnet.network.evolution import NetworkEvolver

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Mock xTB Calculator to bypass external dependency
from unittest.mock import MagicMock
import rsnet.compute.xtb_calculator
import rsnet.compute.reaction_screener

class MockXTBCalculator:
    def __init__(self, *args, **kwargs):
        pass
    def calculate_energy(self, mol, *args, **kwargs):
        return -100.0  # Dummy energy
    def optimize_geometry(self, mol, *args, **kwargs):
        return mol  # Return same mol
    def _check_xtb_availability(self):
        pass  # Bypass check

# Patch the real class
rsnet.compute.xtb_calculator.XTBCalculator = MockXTBCalculator
rsnet.compute.reaction_screener.XTBCalculator = MockXTBCalculator


def test_solvation_clustering():
    print("\n" + "="*70)
    print("锚点聚合策略验证 (Solvation Clustering)")
    print("="*70)
    
    # 1. Setup Species
    ec = Molecule.from_smiles('C1COC(=O)O1', name='EC_Solvent')
    li = Molecule.from_smiles('[Li+]', name='Lithium_Ion')
    
    # 2. Environment (Bulk Electrolyte)
    env = Environment(
        temperature=298.15,
        solvent='EC',
        interface_type='bulk' # Bulk solution favours solvation
    )
    
    print(f"Injecting: {li.name} + {ec.name} (into SpeciesPool)")
    
    # 3. Evolve (Generation 0 -> 1 should form cluster)
    evolver = NetworkEvolver(environment=env)
    
    # We run for just 1 generation to see if it catches the cluster immediately
    products = evolver.evolve([li, ec], max_generations=1)
    
    print("\nEvolution Result:")
    found_complex = False
    for p in products:
        name = p.molecule.name
        smiles = p.molecule.smiles
        print(f"- {name}: {smiles}")
        
        if 'Li' in name and 'EC' in name and ('4' in name or 'complex' in name.lower()):
            found_complex = True
            print("  >>> FOUND SOLVATION COMPLEX! <<<")
            
    if found_complex:
        print("\n✅ Verification SUCCESS: Anchor-Based Clustering worked!")
    else:
        print("\n❌ Verification FAILED: No complex formed.")

if __name__ == "__main__":
    test_solvation_clustering()
