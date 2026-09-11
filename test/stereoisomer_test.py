import pathlib
from rdkit import Chem
from molscrub.protonate import enumerate_stereoisomers


test_smiles = 'BrC(Cl)(F)C=CC(O)C'
test_mol = Chem.MolFromSmiles(test_smiles)
expected_smiles = ["C[C@H](O)/C=C/[C@](F)(Cl)Br",
"C[C@@H](O)/C=C/[C@](F)(Cl)Br",
"C[C@H](O)/C=C/[C@@](F)(Cl)Br",
"C[C@@H](O)/C=C/[C@@](F)(Cl)Br",
"C[C@H](O)/C=C\\[C@](F)(Cl)Br",
"C[C@@H](O)/C=C\\[C@](F)(Cl)Br",
"C[C@H](O)/C=C\\[C@@](F)(Cl)Br",
"C[C@@H](O)/C=C\\[C@@](F)(Cl)Br"]

def test_stereoisomers():
    isomers = enumerate_stereoisomers(test_mol)
    isomers_smiles = [Chem.MolToSmiles(x) for x in isomers]

    assert set(isomers_smiles) == set(expected_smiles)


