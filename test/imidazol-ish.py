#!/usr/bin/env python

import pathlib

import scrubber
from scrubber.transform import MoleculeIsomers
from scrubber.transform import exhaustive_reaction
from scrubber.transform import parse_reaction_file
from scrubber.transform import enumerate_tautomers

from rdkit import Chem

p = pathlib.Path(scrubber.__file__).parents[0]/"data"/"tautomers.txt"

mol = Chem.MolFromSmiles("Cc1nc[nH]c1")
miso = MoleculeIsomers()

reactions = parse_reaction_file(p)
reactions = [(r[3], r[0]) for r in reactions]

exhaust_pool = exhaustive_reaction(mol, reactions)

taut_pool = enumerate_tautomers(mol, reactions, verbose=True)
