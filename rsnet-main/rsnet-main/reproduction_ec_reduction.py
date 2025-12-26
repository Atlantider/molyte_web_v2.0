
from rdkit import Chem
from rsnet.core.molecule import Molecule
from rsnet.core.environment import Environment
from rsnet.operators.electron_transfer import ElectronTransferOperator

def test_ec_reduction():
    # 1. Setup EC
    ec = Molecule(Chem.MolFromSmiles("O=C1OCCO1"), name="EC")
    print(f"Testing EC Reduction. Initial SMILES: {ec.smiles}")
    
    # 2. Setup Environment (Reductive)
    env = Environment(voltage=0.1, solvent="mixture", electrode_type="anode")
    print(f"Environment Voltage: {env.voltage} V")
    
    # 3. Setup Operator
    op = ElectronTransferOperator()
    print(f"Operator Threshold: {op.reduction_threshold} V")
    
    # 4. Apply
    if op.can_apply([ec], env):
        print("Operator can apply.")
        reactions = op.apply([ec], env)
        print(f"Generated {len(reactions)} reactions.")
        for r in reactions:
            p = r.products[0]
            print(f"Product: {p.name} - {p.smiles}")
            print(f"Charge: {Chem.GetFormalCharge(p.rdkit_mol)}")
            print(f"Radicals: {sum(a.GetNumRadicalElectrons() for a in p.rdkit_mol.GetAtoms())}")
    else:
        print("Operator cannot apply (conditions not met).")
        
if __name__ == "__main__":
    test_ec_reduction()
