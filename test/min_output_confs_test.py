from rdkit import Chem
from molscrub import Scrub


def multiple_output_confs(scrub):
    mol = Chem.MolFromSmiles("OCCC")
    mols = scrub(mol)
    x_coord_1 = mols[0].GetConformer(0).GetPositions()[0][0]
    x_coord_2 = mols[0].GetConformer(1).GetPositions()[0][0]
    assert x_coord_1 != x_coord_2

def test_multiple_output_confs_no_rng_seed():
    scrub = Scrub(min_output_confs=2)
    multiple_output_confs(scrub)

def test_multiple_output_confs_minus_one_rng_seed():
    scrub = Scrub(min_output_confs=2, etkdg_rng_seed=-1)
    multiple_output_confs(scrub)

def test_multiple_output_confs_specified_rng_seed():
    scrub = Scrub(min_output_confs=2, etkdg_rng_seed=12345)
    multiple_output_confs(scrub)
