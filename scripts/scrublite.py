#!/usr/bin/env python

import argparse
import pathlib
import sys

from scrubber import enumerate_pka
from scrubber import parse_reaction_file
from scrubber import build_pka_reactions
from scrubber import parse_tautomers_config_file
from scrubber import enumerate_tautomers
from scrubber import Scrub
import scrubber

from rdkit import Chem

scrubber_dir = pathlib.Path(scrubber.__file__).parents[0]
pka_file = scrubber_dir/"data"/"pka_reactions.txt"
tauto_file = scrubber_dir/"data"/"tautomers.txt"

parser = argparse.ArgumentParser(description="Protonate molecules and add 3D coordinates")
parser.add_argument("input", help="input filename (.sdf/.mol/.smi) or SMILES string")
parser.add_argument("-o", "--out_fname", help="output filename (.sdf) default: 'scrubbed.sdf'")
parser.add_argument("--ph", help="pH value for acid/base transformations", default=7.4, type=float)
parser.add_argument("--ph_low", help="low end of pH range (superseeds --ph)", type=float)
parser.add_argument("--ph_high", help="high end of pH range (superseeds --ph)", type=float)
args = parser.parse_args()

if args.ph_low is None and args.ph_high is None:
    ph_low = args.ph
    ph_high = args.ph
elif args.ph_low is not None and args.ph_high is not None:
    ph_low = args.ph_low
    ph_high = args.ph_high
else:
    print("--ph_low and --ph_high work together, either use both or none.")
    sys.exit()

extension = pathlib.Path(args.input).suffix
if extension == ".sdf":
    supplier = Chem.SDMolSupplier(args.input)
elif extension == ".mol":
    supplier = [Chem.MolFromMolFile(args.input)]
else:
    mol = Chem.MolFromSmiles(args.input)
    if mol is None:
        print("Input parsed as SMILES string, but conversion to RDKit mol failed.")
        print("Extension must be .sdf/.mol/.smi to parse input as filename.")
        sys.exit()
    supplier = [mol]

scrub = Scrub(ph_low, ph_high)

counter_input_mols = 0
counter_input_none = 0
for input_mol in supplier:
    print(input_mol.GetProp("_Name"))
    if input_mol is None:
        counter_input_none += 1
        continue
    counter_input_mols += 1
    for mol in scrub(input_mol):
        print("%2d" % counter_input_mols, mol.GetProp("_Name"), mol.GetNumConformers())
