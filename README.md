# Scrubber
Process large numbers of small molecules for docking with AutoDock.
May be useful for structure-based modeling in general.

What happens:
 - generate 3D coordinates using RDKit's ETKDGv3 and UFF minimization
 - enumerate tautomers (aiming at low energy states only)
 - enumerate pH corrections
 - convert boats to chairs (6-member rings) and enumerate both chair states
 - enumerate chiral centers (not implemented right now)


# Installation
```sh
conda activate <desired-environment>    # if you are using conda environments

git clone git@github.com:forlilab/scrubber.git
cd scrubber
pip install -e .
```

Depends on the RDKit, which can be installed from conda-forge in the desired environment:
```sh
conda activate <desired-environment>
conda install rdkit -c conda-forge
```

## Python scripting
```python
from rdkit import Chem
from scrubber import Scrub

scrub = Scrub(
    ph_low=7.4,
    ph_high=7.4,
)

mol = Chem.MolFromSmiles("Clc1c(OCCC3)c3ccc1C(=O)Nc2nc[nH]c2")

# each state (e.g. tautomer) an rdkit mol and may have multiple conformers
for mol_state in scrub(mol):
    print(Chem.MolToSmiles(mol_state), "nr conformers: %d" % mol_state.GetNumConformers())
```

## Command line tool examples
```sh
scrub.py "c1cc[nH]c(=O)c1" -o scrubbed.sdf --pH 5 --skip_gen3d
scrub.py input_mols.sdf -o scrubbed.sdf
scrub.py input_mols.smi -o scrubbed.sdf
```

Other options described in the help message:
```sh
scrub.py -h
```

Where "input\_mols.smi" can look like this:
```
CC(=O)O aceticacid
CN(C)C trimethylamine 
Clc1cc(O)ccc1C(=O)Nc2nc[nH]c2 hello_mol
c1cccc1 rdkit_will_cry
CCC good4bbq
CCO alsogood4bbq
c1cccnc1CC(=O)C a_ketone
```

## Using a .toml file

Another way to use `scrub.py` is to place all command line options inside 
a .toml input file:

```toml
input = "c1cc[nH]c(=O)c1"
out_fname = "output.sdf"
ph = 5
skip_gen3d = true
```

This can be executed with `scrub.py input.toml` as usual. 

This simplifies the reusability of the same command line options. Additionally, 
multiple SMILES strings can be placed as an array inside the .toml file: 

```toml
input = [
    "c1cc[nH]c(=O)c1",
    "CC(=O)O",
    "c1cccnc1CC(=O)C"
]
out_fname = "output.sdf"
ph = 5
```

Of course, you may still pass additional command line arguments after the toml input file.
The command line options always take precedence over any options included in the toml file. 