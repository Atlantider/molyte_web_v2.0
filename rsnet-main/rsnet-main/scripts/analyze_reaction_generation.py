
import logging
import sys
import os
import time
from collections import Counter
from rdkit import Chem

# Ensure we can import rsnet
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from rsnet.core.molecule import Molecule
from rsnet.core.environment import Environment
from rsnet.network.generator import NetworkGenerator
from rsnet.network.config import NetworkGenerationConfig, NetworkGenerationPresets, NetworkGenerationStrategy, OperatorSelectionStrategy
from rsnet.operators.registry import OPERATOR_REGISTRY
from rsnet.compute.reaction_screener import ReactionScreener

class MockCalculator:
    """Mock Calculator when xTB is missing."""
    def __init__(self):
        pass
        
    def optimize(self, molecule, environment):
        import random
        # Improved Mock Physics:
        n_atoms = molecule.rdkit_mol.GetNumAtoms()
        base_energy = n_atoms * -100.0
        
        # Add a "stability bonus"
        smi = molecule.smiles
        stability_bonus = 0.0
        if "Li" in smi and "F" in smi: stability_bonus = -50.0
        if "C(=O)O" in smi: stability_bonus = -20.0
        
        # Add electrostatic penalty for high charge
        try:
            charge = Chem.GetFormalCharge(molecule.rdkit_mol)
        except:
            charge = 0
        charge_penalty = (charge ** 2) * 50.0
        
        final_energy = base_energy + stability_bonus + charge_penalty + random.uniform(-10, 10)
        
        return {
            'energy_kcal_mol': final_energy,
            'success': True,
            'optimized_molecule': molecule
        }
        
    def single_point(self, molecule, environment):
        return self.optimize(molecule, environment)

# Configure Logging to capture details
logging.basicConfig(
    level=logging.DEBUG,
    format='%(message)s',
    handlers=[logging.StreamHandler()]
)
logger = logging.getLogger("analysis")

def analyze_reaction_generation():
    print("="*60)
    print("REACTION GENERATION ANALYSIS")
    print("="*60)

    # 1. Setup Environment matching run_full_simulation.py
    env = Environment(
        temperature=298.0,
        pressure=1.0,
        solvent="mixture", 
        voltage=0.1,  # Anode
        electrode_type="anode",
        electrode_material="graphite"
    )

    # 2. Setup Initial Molecules
    mols = []
    # EC
    mols.append(Molecule(Chem.MolFromSmiles("O=C1OCCO1"), name="EC"))
    # DMC
    mols.append(Molecule(Chem.MolFromSmiles("COC(=O)OC"), name="DMC"))
    # Li+
    mols.append(Molecule(Chem.MolFromSmiles("[Li+]"), name="Li+"))
    # PF6- (simplified)
    try:
        mols.append(Molecule(Chem.MolFromSmiles("F[P-](F)(F)(F)(F)F"), name="PF6-"))
    except:
        mols.append(Molecule(Chem.MolFromSmiles("[F-]"), name="F-"))

    print(f"Initial Species: {[m.name for m in mols]}")

    # 3. Configure Generator
    # Use battery chemistry preset but limit to 3 generations for deep analysis
    config = NetworkGenerationPresets.battery_chemistry()
    config.max_generations = 3
    config.verbose_logging = True
    # Ensure we get many candidates to analyze the "explosion"
    config.max_species = 10000 
    
    config.generation_strategy = NetworkGenerationStrategy.EXHAUSTIVE
    config.operator_selection_strategy = OperatorSelectionStrategy.ALL
    config.filter_duplicate_reactions = False # See raw count
    config.filter_trivial_reactions = False
    config.max_reactions_per_generation = 5000 # Allow enough for analysis but prevent crash
    
    print("\nConfiguration (Modified):")
    print(f"  Strategy: {config.generation_strategy}")
    print(f"  Operator Selection: {config.operator_selection_strategy}")
    print(f"  Driving Force Threshold: {config.driving_force_threshold}")

    # 4. Initialize Generator with Mock Calculator
    mock_calc = MockCalculator()
    screener = ReactionScreener(calculator=mock_calc)
    generator = NetworkGenerator(config=config, screener=screener)

    # 5. Run Generation
    print("\nRunning Generation 0...")
    start_time = time.time()
    network = generator.generate_network(mols, env)
    duration = time.time() - start_time

    # 6. Analyze Results
    print(f"\nGeneration completed in {duration:.2f}s")
    print(f"Total Molecules: {len(network.molecules)}")
    print(f"Total Reactions: {len(network.reactions)}")

    if len(network.reactions) == 0:
        print("WARNING: No reactions generated. Check connectivity or operator logic.")
        return

    # Analyze by Operator
    operator_counts = Counter()
    operator_energy_stats = {} # operator -> [energies]

    for reaction in network.reactions.values():
        op_name = reaction.operator_name
        operator_counts[op_name] += 1
        
        if op_name not in operator_energy_stats:
            operator_energy_stats[op_name] = []
        
        if reaction.reaction_energy is not None:
            operator_energy_stats[op_name].append(reaction.reaction_energy)

    print("\nReactions by Operator:")
    print("-" * 40)
    print(f"{'Operator':<25} | {'Count':<10} | {'Avg Energy':<10} | {'% of Total':<10}")
    print("-" * 40)
    
    total_rxns = len(network.reactions)
    
    for op, count in operator_counts.most_common():
        energies = operator_energy_stats.get(op, [])
        avg_energy = sum(energies)/len(energies) if energies else 0.0
        percent = (count / total_rxns) * 100
        print(f"{op:<25} | {count:<10} | {avg_energy:<10.2f} | {percent:<9.1f}%")

    # Analyze Top Products
    print("\nTop 5 Most Common Product Molecules (by connectivity):")
    # Just listing random products for now as "frequency" in network is connectivity
    # A better metric might be "which products are created by the most reactions"
    
    product_counts = Counter()
    for rxn in network.reactions.values():
        for p in rxn.products:
            product_counts[p.name] += 1
            
    for name, count in product_counts.most_common(5):
        print(f"  {name}: generated by {count} reactions")

if __name__ == "__main__":
    analyze_reaction_generation()
