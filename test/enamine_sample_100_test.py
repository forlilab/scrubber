import pathlib
from rdkit import Chem
from molscrub import Scrub
import multiprocessing


workdir = pathlib.Path(__file__)
datadir = workdir.parents[0] / "data"
mols_filename = str(datadir / "enamine-sample-100-scrubbed.sdf")
mol_supplier = Chem.SDMolSupplier(mols_filename)
nr_mols = len(mol_supplier)

def run(scrub):
    nr_cumulative_states = 0
    with multiprocessing.Pool() as pool:
        for i, output in enumerate(pool.imap(scrub, mol_supplier)):
            msg = f"no output for mol with {i=} from {mols_filename}"
            assert len(output), msg
    return

def test_enamine_sample_defaults():
    run(Scrub())

def test_enamine_sample_pka_etr1():
    run(Scrub(pka_model="etr1"))
