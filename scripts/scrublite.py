#!/usr/bin/env python

import argparse
import json
import pathlib
import sys

from scrubber import Scrub
from scrubber import SMIMolSupplierWrapper

from rdkit import Chem
from rdkit.Chem import rdMolInterchange

try:
    import h5py
    _got_h5py = True
except ImportError as e:
    _h5py_import_error = e
    _got_h5py = False

class SDWriter:
    """support Python's `with` statement and always write all conformers"""

    def __init__(self, filename):
        self.filename = filename

    def __enter__(self):
        self.rdkit_sdwriter = Chem.SDWriter(self.filename)
        self.counter_mol_group = 0
        return self

    def __exit__(self, *args):
        self.rdkit_sdwriter.close()

    def write_mols(self, mol_group):
        for mol in mol_group:
            for i, conf in enumerate(mol.GetConformers()):
                mol.SetProp("ScrubInfo", json.dumps({
                    "group_id": self.counter_mol_group,
                    "conformer": i,
                    "nr_conformers:":  mol.GetNumConformers(),
                }))
                self.rdkit_sdwriter.write(mol, confId=i)
        self.counter_mol_group += 1

class HDF5Writer:

    def __init__(self, filename):
        self.filename = filename

    def __enter__(self):
        self.h5file = h5py.File(self.filename, "w")
        dt = h5py.string_dtype()
        self.mols = self.h5file.create_dataset("mols", (0,), maxshape=(None,), dtype=dt, chunks=True)
        self.group_id = self.h5file.create_dataset("group_id", (0,), maxshape=(None,), dtype="i8", chunks=True)
        self.counter_mol_group = 0
        return self

    def __exit__(self, *args):
        self.h5file.close()

    def write_mols(self, mol_group):
        for mol in mol_group:
            index = self.mols.shape[0]
            self.mols.resize((index + 1, ))
            self.mols[index] = rdMolInterchange.MolToJSON(mol)
            self.group_id.resize((index + 1, ))
            self.group_id[index] = self.counter_mol_group
        self.counter_mol_group += 1


def get_info_str():
    s = ""
    s += "Input molecules processed: %d, skipped by rdkit: %d\n" % (
            counter_input_mols, counter_supplied - counter_input_mols)
    s += "nr isomers (tautomers and acid/base conjugates): %d (avg. %.3f per mol)\n" % (
            counter_isomers, counter_isomers/counter_input_mols)
    s += "nr conformers:  %d (avg. %.3f per isomer, %.3f per mol)\n" % (
            counter_conformers,
            counter_conformers/counter_isomers,
            counter_conformers/counter_input_mols)
    return s


parser = argparse.ArgumentParser(description="Protonate molecules and add 3D coordinates")
parser.add_argument("input", help="input filename (.sdf/.mol/.smi) or SMILES string")
parser.add_argument("-o", "--out_fname", help="output filename (.sdf/.hdf5)", required=True)
parser.add_argument("--ph", help="pH value for acid/base transformations", default=7.4, type=float)
parser.add_argument("--ph_low", help="low end of pH range (superseeds --ph)", type=float)
parser.add_argument("--ph_high", help="high end of pH range (superseeds --ph)", type=float)
parser.add_argument("--no_acidbase", help="skip enumeration of acid/base conjugates", action="store_true")
parser.add_argument("--no_tautomers", help="skip enumeration of tautomers", action="store_true")
parser.add_argument("--no_ringfix", help="skip fixes of six-member rings", action="store_true")
parser.add_argument("--no_gen3d", help="skip generation of 3D coordinates (also skips ring fixes)", action="store_true")
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

# input
extension = pathlib.Path(args.input).suffix
if extension == ".sdf":
    supplier = Chem.SDMolSupplier(args.input)
elif extension == ".mol":
    supplier = [Chem.MolFromMolFile(args.input)]
elif extension == ".smi":
    supplier = SMIMolSupplierWrapper(args.input)
else:
    mol = Chem.MolFromSmiles(args.input)
    if mol is None:
        print("Input parsed as SMILES string, but conversion to RDKit mol failed.")
        print("The SMILES might be incorrect.")
        print("If you want to pass a filename, its extension must be .sdf/.mol/.smi.")
        sys.exit()
    supplier = [mol]

# output
extension = pathlib.Path(args.out_fname).suffix
if extension == ".sdf":
    Writer = SDWriter
elif extension == ".hdf5":
    if _got_h5py:
        Writer = HDF5Writer
    else:
        print(_h5py_import_error, file=sys.stderr)
        print("Could not import h5py. Install h5py to write .hdf5")
        sys.exit()
else:
    print("output file extension must be .sdf/.hdf5")
    sys.exit()

scrub = Scrub(
    ph_low,
    ph_high,
    pka_fname=None,
    tauto_fname=None,
    no_acidbase=args.no_acidbase,
    no_tautomers=args.no_tautomers,
    no_ringfix=args.no_ringfix,
    no_gen3d=args.no_gen3d,
)

counter_supplied = 0
counter_input_mols = 0
counter_isomers = 0
counter_conformers = 0
with Writer(args.out_fname) as w:
    for input_mol in supplier:
        counter_supplied += 1
        if input_mol is None:
            continue
        counter_input_mols += 1
        isomer_list = scrub(input_mol)
        w.write_mols(isomer_list)
        counter_isomers += len(isomer_list)
        counter_conformers += sum([mol.GetNumConformers() for mol in isomer_list])
        if counter_supplied % 100 == 0:
            print("Scrub in progress. Here's how things are going:")
            print(get_info_str())

print("Scrub completed.\nSummary of what happened:")
print(get_info_str(), end="")
