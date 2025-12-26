from rdkit import Chem

smiles_list = [
    "O=C1OCCO1.[Li+]", 
    "[Li+]", 
    "COC(=O)OC.[Li+]",
    "F[P-](F)(F)(F)(F)F.[Li+]"
]

print("Checking Formal Charges:")
for smi in smiles_list:
    mol = Chem.MolFromSmiles(smi)
    charge = Chem.GetFormalCharge(mol)
    print(f"SMILES: {smi} -> Charge: {charge}")

print("\nChecking Explicit Valence for problematic atom logs:")
# The log had "Explicit valence for atom # 1 C, 4, is greater than permitted"
# This usually happens when adding Hydrogens or sanitizing bad structures.
# Let's try to reproduce with a reduced carbonyl
smi_red = "[O-]C1OCCO1" # EC radical anion opened?
mol_red = Chem.MolFromSmiles(smi_red)
if mol_red:
    print(f"Valid Mol: {smi_red}")
else:
    print(f"Invalid Mol: {smi_red}")
