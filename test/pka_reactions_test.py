from rdkit import Chem
from molscrub import AcidBaseConjugator

acids = Chem.SmilesMolSupplier("test/data/acids.smi", titleLine=False)
bases = Chem.SmilesMolSupplier("test/data/bases.smi", titleLine=False)

def test_acidbase_reactions():
    conjugator = AcidBaseConjugator.from_default_data_files()

    acids_mols = [conjugator.generate_all_protonation_states(m) for m in acids]
    bases_mols = [conjugator.generate_all_protonation_states(m) for m in bases]

    assert len(acids_mols) == len(acids)
    assert len(bases_mols) == len(bases)







