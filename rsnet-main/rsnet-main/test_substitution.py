
import unittest
from rdkit import Chem
from rsnet.operators.substitution import SubstitutionOperator
from rsnet.core.molecule import Molecule
from rsnet.core.environment import Environment

class TestSubstitutionOperator(unittest.TestCase):
    def setUp(self):
        self.operator = SubstitutionOperator()
        self.env = Environment(temperature=298.0, solvent="DMC")
        
    def test_transesterification(self):
        """Test Methoxide + DMC (Dimethyl Carbonate) transesterification."""
        # Methoxide [O-]C
        nuc = Molecule(Chem.MolFromSmiles("[O-]C"), name="Methoxide")
        # DMC COC(=O)OC
        sub = Molecule(Chem.MolFromSmiles("COC(=O)OC"), name="DMC")
        
        # Expectation: CH3O- + CH3OCOOCH3 -> CH3OCOOCH3 + CH3O- (Identity reaction actually, let's use EMC)
        # EMC: COC(=O)OCC
        sub = Molecule(Chem.MolFromSmiles("COC(=O)OCC"), name="EMC")
        
        molecules = [nuc, sub]
        reactions = self.operator.apply(molecules, self.env)
        
        found = False
        for rxn in reactions:
            print(f"Reaction: {rxn.name}")
            for p in rxn.products:
                # Remove explicit hydrogens for SMILES comparison to be robust
                mol_no_h = Chem.RemoveHs(p.rdkit_mol)
                smi = Chem.MolToSmiles(mol_no_h)
                print(f"  Product: {smi}")
                
                # Check for DMC (dimethyl carbonate)
                # SMILES: COC(=O)OC
                if smi == "COC(=O)OC": # DMC formed
                    found = True
                    
        self.assertTrue(found, "Should find DMC product from EMC + Methoxide transesterification")

    def test_sn2_alkylation(self):
        """Test Methoxide + DMC Sn2 alkylation (attack on methyl)."""
        # Methoxide [O-]C
        nuc = Molecule(Chem.MolFromSmiles("[O-]C"), name="Methoxide")
        # DMC COC(=O)OC
        sub = Molecule(Chem.MolFromSmiles("COC(=O)OC"), name="DMC")
        
        molecules = [nuc, sub]
        reactions = self.operator.apply(molecules, self.env)
        
        found_ether = False
        for rxn in reactions:
            for p in rxn.products:
                mol_no_h = Chem.RemoveHs(p.rdkit_mol)
                smi = Chem.MolToSmiles(mol_no_h)
                # Sn2: MeO- + Me-O-C(=O)OMe -> Me-O-Me + (-)O-C(=O)OMe
                if smi == "COC": # Dimethyl ether
                    found_ether = True
                    
        self.assertTrue(found_ether, "Should find Dimethyl Ether from Sn2 attack")

if __name__ == "__main__":
    unittest.main()
