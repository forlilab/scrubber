#!/usr/bin/env python

import argparse
import json
import multiprocessing
import pathlib
import sys

from scrubber import Scrub
from scrubber import SMIMolSupplierWrapper

from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import rdMolInterchange

Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.MolProps |
                                Chem.PropertyPickleOptions.PrivateProps)
RDLogger.DisableLog("rdApp.*")

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

    def write_mols(self, mol_group, add_suffix=False):

        if add_suffix and len(mol_group) > 0:
            if mol_group[0].HasProp("_Name"):
                name = mol_group[0].GetProp("_Name") # assumes all mols have same name, which they should
            else:
                name = ""
        nr_isomers = len(isomer_list)
        for i, mol in enumerate(mol_group):
            nr_confs = mol.GetNumConformers()
            for j, conf in enumerate(mol.GetConformers()):
                mol.SetProp("ScrubInfo", json.dumps({
                    "isomerGroup": self.counter_mol_group,
                    "isomerId": i,
                    "confId": j,
                    "nr_conformers:":  mol.GetNumConformers(),
                    "nr_isomers:":  len(mol_group),
                }))
                if add_suffix:
                    if nr_isomers > 1 and nr_confs > 1:
                        suffix = "_i%d-c%d" % (i, j)
                    elif nr_isomers > 1 and nr_confs <= 1:
                        suffix = "_i%d" % i
                    elif nr_isomers <= 1 and nr_confs > 1:
                        suffix = "_c%d" % j
                    else:
                        suffix = ""
                    mol.SetProp("_Name", name + suffix)
                self.rdkit_sdwriter.write(mol, confId=j)
        self.counter_mol_group += 1


def suffix_iter(isomer_list):
    if self.do_name_suffixes:
        for i, mol in enumerate(isomer_list):
            mol.SetProp("_Name", name + suffix)


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

    def write_mols(self, mol_group, add_suffix=False):
        for mol in mol_group:
            if add_suffix:
                if mol.HasProp("_Name"):
                    name = mol.GetProp("_Name")
                else:
                    name = ""
                if len(mol_group) > 1:
                    suffix = "_i%d" % len(mol_group)
                else:
                    suffix = ""
                mol.SetProp("_Name", name + suffix)
            index = self.mols.shape[0]
            self.mols.resize((index + 1, ))
            self.mols[index] = rdMolInterchange.MolToJSON(mol)
            self.group_id.resize((index + 1, ))
            self.group_id[index] = self.counter_mol_group
        self.counter_mol_group += 1


def get_info_str(counter=None):
    if counter is not None:
        counter_supplied   = counter["supplied"]
        counter_rdkit_nope = counter["rdkit_nope"]
        counter_ok_mols    = counter["ok_mols"]
        counter_isomers    = counter["isomers"]
        counter_conformers = counter["conformers"]
        counter_failed     = counter["failed"]

    s = ""
    s += "Input molecules supplied: %d\n" % counter_supplied
    s += "mols processed: %d, skipped by rdkit: %d, failed: %d\n" % (
            counter_ok_mols, counter_rdkit_nope, counter_failed)
    if counter_ok_mols == 0:
        return s
    s += "nr isomers (tautomers and acid/base conjugates): %d (avg. %.3f per mol)\n" % (
            counter_isomers, counter_isomers/counter_ok_mols)
    s += "nr conformers:  %d (avg. %.3f per isomer, %.3f per mol)\n" % (
            counter_conformers,
            counter_conformers/counter_isomers,
            counter_conformers/counter_ok_mols)
    return s


parser = argparse.ArgumentParser(description="Protonate molecules and add 3D coordinates")
parser.add_argument("input", help="input filename (.sdf/.mol/.smi) or SMILES string")
parser.add_argument("-o", "--out_fname", help="output filename (.sdf/.hdf5)", required=True)
parser.add_argument("--name_from_prop", help="set molecule name from RDKit/SDF property")
parser.add_argument("--ph", help="pH value for acid/base transformations", default=7.4, type=float)
parser.add_argument("--ph_low", help="low end of pH range (superseeds --ph)", type=float)
parser.add_argument("--ph_high", help="high end of pH range (superseeds --ph)", type=float)
parser.add_argument("--skip_acidbase", help="skip enumeration of acid/base conjugates", action="store_true")
parser.add_argument("--skip_tautomers", help="skip enumeration of tautomers", action="store_true")
parser.add_argument("--skip_ringfix", help="skip fixes of six-member rings", action="store_true")
parser.add_argument("--skip_gen3d", help="skip generation of 3D coordinates (also skips ring fixes)", action="store_true")
parser.add_argument("--cpu", help="number of processes to run in parallel", default=0, type=int)
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
do_gen2d = False # if output SDF and skip_gen3d, we will need 2D conformers
extension = pathlib.Path(args.out_fname).suffix
if extension == ".sdf":
    Writer = SDWriter
    if args.skip_gen3d:
        do_gen2d = True
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
    skip_acidbase=args.skip_acidbase,
    skip_tautomers=args.skip_tautomers,
    skip_ringfix=args.skip_ringfix,
    skip_gen3d=args.skip_gen3d,
    do_gen2d=do_gen2d,
    name_from_prop=args.name_from_prop,
)

counter_supplied = 0
counter_ok_mols = 0
counter_isomers = 0
counter_conformers = 0
counter_failed = 0

counter = {
    "supplied": 0,
    "rdkit_nope": 0,
    "ok_mols": 0,
    "isomers": 0,
    "conformers": 0,
    "failed": 0,
}

def write_and_log(isomer_list, log, counter):
    counter["supplied"] += 1
    if log["input_mol_none"]:
        counter["rdkit_nope"] += 1
    elif len(isomer_list):
        try:
            w.write_mols(isomer_list, add_suffix=True)
            counter["ok_mols"] += 1
        except Exception as e:
            print(e, file=sys.stderr)
            counter["failed"] += 1
            return
        counter["isomers"] += len(isomer_list)
        counter["conformers"] += sum([mol.GetNumConformers() for mol in isomer_list])
        if counter["supplied"] % 100 == 0:
            print("Scrub in progress. Here's how things are going:")
            print(get_info_str(counter))
    else:
        counter["failed"] += 1
        if "exception" in log:
            print(log["exception"], file=sys.stderr)

with Writer(args.out_fname) as w:
    if args.cpu == 1:
        for input_mol in supplier:
            isomer_list, log = scrub(input_mol)
            write_and_log(isomer_list, log, counter)
    else:
        if args.cpu < 1:
            nr_proc = multiprocessing.cpu_count()
        else:
            nr_proc = args.cpu
        p = multiprocessing.Pool(nr_proc - 1) # leave 1 for main process
        for (isomer_list, log) in p.imap_unordered(scrub, supplier):
            write_and_log(isomer_list, log, counter)

print("Scrub completed.\nSummary of what happened:")
print(get_info_str(counter), end="")
