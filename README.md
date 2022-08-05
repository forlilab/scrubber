# Scrubber
Process large numbers of small molecules for docking with AutoDock.
May be useful for structure-based modeling in general.

What happens:
 - generate 3D coordinates using RDKit's ETKDGv3 and UFF minimization
 - enumerate tautomers (aiming at low energy states only)
 - enumerate pH corrections
 - convert boats to chairs (6-member rings) and enumerate both chair states
 - enumerate chiral centers


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
from scrubber.core import ScrubberCore

scrub = ScrubberCore()

mol = Chem.MolFromSmiles("Clc1ccccc1C(=O)Nc2nc[nH]c2")
microstates = []
microstate_generator = scrub.process(mol)

# each state is an rdkit mol
for state in microstate_generator:
    microstates.append(state)
```

## Command line tool
```sh
scrub.py --in_fname test.smi --out_fname test.sdf
```

Where "test.smi" can look like this:
```
CC(=O)O aceticacid
CN(C)C trimethylamine 
Clc1cc(O)ccc1C(=O)Nc2nc[nH]c2 hello_mol
c1cccc1 rdkit_will_cry
CCC good4bbq
CCO alsogood4bbq
c1cccnc1CC(=O)C a_ketone
```
