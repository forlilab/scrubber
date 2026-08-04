import pathlib
from rdkit import Chem
from molscrub import AcidBaseConjugator

from molscrub.core import load_model

model = load_model()

conjugator = AcidBaseConjugator.from_default_data_files(model) 

test_smiles = 'N1CC2C3C(C1)CNC4CCNC(C43)CN2'
test_mol = Chem.MolFromSmiles(test_smiles)

def test_protonation_generator():
    mols = conjugator.generate_all_protonation_states(test_mol)
    assert len(mols) == 16


def test_pairs():
    mols = conjugator.generate_all_protonation_states(test_mol)
    pairs = conjugator.build_pairs()
    assert len(pairs) == 32

