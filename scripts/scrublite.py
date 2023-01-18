#!/usr/bin/env python

import sys
import pathlib

from scrubber import enumerate_pka
from scrubber import parse_reaction_file
from scrubber import build_pka_reactions
from scrubber import parse_tautomers_config_file
from scrubber import enumerate_tautomers
import scrubber

from rdkit import Chem

scrubber_dir = pathlib.Path(scrubber.__file__).parents[0]
pka_file = scrubber_dir/"data"/"pka_reactions.txt"
tauto_file = scrubber_dir/"data"/"tautomers.txt"

if __name__ == "__main__":
    smiles = sys.argv[1]    
    ph_low = float(sys.argv[2])
    if len(sys.argv) > 3:
        ph_high = float(sys.argv[3])
    else:
        ph_high = ph_low
    mol = Chem.MolFromSmiles(smiles)
    reactions = parse_reaction_file(pka_file)
    pka_reactions = build_pka_reactions(reactions)
    output = enumerate_pka(mol, pka_reactions, ph_low, ph_high)
    print("Num acid-base states:", len(output))
    for o in output:
        print(Chem.MolToSmiles(o))
    print()

    # tautomers
    tauto_reactions, keepmax_smarts = parse_tautomers_config_file(tauto_file)
    tautomers = enumerate_tautomers(mol, tauto_reactions, keepmax_smarts)
    print("Num tautomers:", len(tautomers))
    for o in tautomers:
        print(Chem.MolToSmiles(o))
