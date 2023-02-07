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

    def write_mols(self, mol_group, add_suffix=False, add_serial_suffix=False):

        add_suffix |= add_serial_suffix
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
                if add_serial_suffix:
                    if nr_isomers > 1 or nr_confs > 1:
                        suffix = "_%d" % (i + 1) 
                        mol.SetProp("_Name", name + suffix)
                elif add_suffix:
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
        nr_isomers = len(isomer_list)
        for i, mol in enumerate(mol_group):
            nr_confs = mol.GetNumConformers()
            if add_serial_suffix:
                if nr_isomers > 1 or nr_confs > 1:
                    suffix = "_%d" % (i + 1) 
            elif add_suffix:
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


class MolSupplier:
    """wraps other suppliers (e.g. Chem.SDMolSupplier) to change non-integer
        molecule names to integers, and to set rdkit mol names from properties
    """

    def __init__(self, supplier, name_from_prop=None, rename_to_int=False, nr_digits=10):
        self.supplier = supplier
        self.name_from_prop = name_from_prop
        self.rename_to_int = rename_to_int
        self.nr_digits = nr_digits
        self.names = {}
        self.counter = 0
        
    def __iter__(self):
        self.supplier.reset()
        return self

    def __next__(self):
        mol = self.supplier.__next__()
        if mol is None:
            return mol
        if self.name_from_prop:
            name = mol.GetProp(self.name_from_prop)
            mol.SetProp("_Name", name)
        if self.rename_to_int:
            name = mol.GetProp("_Name")
            newname = self._rename(name)
            mol.SetProp("_Name", newname)
        return mol
        
    def _rename(self, name):
        """rename if name is not an integer, or a sequence of alphabet chars
            followed by an integer."""

        # special case for Enamine's molecules
        if name.startswith("PV-") and name[3:].isdigit():
            return "PV" + name[3:] # remove dash from Enamine's PV-000000000000
        is_good = False
        if name.isalnum():
            # make sure all letters preceed the decimals, no mix
            is_good = True
            num_started = False
            for c in name:
                num_started |= c.isdecimal()
                if num_started and not c.isdecimal():
                    is_good = False
                    break
        if is_good:
            return name
        
        self.counter += 1
        #if name in self.names:
        #    raise RuntimeError("repeated molecule name: %s" % name)
        #self.names[name] = self.counter
        self.names[self.counter] = name
        tmp = "RN%0" + "%d" % self.nr_digits + "d"
        return tmp % self.counter


def get_info_str(counter):
    c = counter
    s = ""
    s += "Input molecules supplied: %d\n" % c["supplied"]
    s += "mols processed: %d, skipped by rdkit: %d, failed: %d\n" % (
            c["ok_mols"], c["rdkit_nope"], c["failed"])
    if c["ok_mols"] == 0:
        return s
    s += "nr isomers (tautomers and acid/base conjugates): %d (avg. %.3f per mol)\n" % (
            c["isomers"], c["isomers"]/c["ok_mols"])
    s += "nr conformers:  %d (avg. %.3f per isomer, %.3f per mol)\n" % (
            c["conformers"],
            c["conformers"]/c["isomers"],
            c["conformers"]/c["ok_mols"])
    return s


parser = argparse.ArgumentParser(description="Protonate molecules and add 3D coordinates")
parser.add_argument("input", help="input filename (.sdf/.mol/.smi/.cxsmiles) or SMILES string")
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
parser.add_argument("--debug", help="errors are raised", action="store_true")
parser.add_argument("--wcg", help="make sure mol names and suffixes are integers", action="store_true")
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
elif extension == ".cxsmiles":
    supplier = SMIMolSupplierWrapper(args.input, is_enamine_cxsmiles=True, titleLine=True)
else:
    mol = Chem.MolFromSmiles(args.input)
    if mol is None:
        print("Input parsed as SMILES string, but conversion to RDKit mol failed.")
        print("The SMILES might be incorrect.")
        print("If you want to pass a filename, its extension must be .sdf/.mol/.smi.")
        sys.exit()
    supplier = [mol]

if args.wcg or args.name_from_prop:
    supplier = MolSupplier(
        supplier,
        name_from_prop=args.name_from_prop,
        rename_to_int=args.wcg
    )

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
)

counter = {
    "supplied": 0,
    "rdkit_nope": 0,
    "ok_mols": 0,
    "isomers": 0,
    "conformers": 0,
    "failed": 0,
}

def scrub_and_catch_errors(input_mol):
    log = {}
    if input_mol is None:
        log["input_mol_none"] = True
        isomer_list = []
    else:
        log["input_mol_none"] = False
        try:
            isomer_list = scrub(input_mol)
        except Exception as e:
            log["exception"] = e
            isomer_list = []
    return (isomer_list, log)

def scrub_and_debug(input_mol):
    log = {"input_mol_none": input_mol is None}
    isomer_list = scrub(input_mol)
    return (isomer_list, log)

def write_and_log(isomer_list, log, counter):
    counter["supplied"] += 1
    if log["input_mol_none"]:
        counter["rdkit_nope"] += 1
    elif len(isomer_list):
        try:
            w.write_mols(isomer_list, add_suffix=True, add_serial_suffix=args.wcg)
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

if args.debug:
    scrub_fn = scrub_and_debug
else:
    scrub_fn = scrub_and_catch_errors

with Writer(args.out_fname) as w:
    if args.cpu == 1:
        for input_mol in supplier:
            isomer_list, log = scrub_fn(input_mol)
            write_and_log(isomer_list, log, counter)
    else:
        if args.cpu < 1:
            nr_proc = multiprocessing.cpu_count()
        else:
            nr_proc = args.cpu
        p = multiprocessing.Pool(nr_proc - 1) # leave 1 for main process
        for (isomer_list, log) in p.imap_unordered(scrub_fn, supplier):
            write_and_log(isomer_list, log, counter)

print("Scrub completed.\nSummary of what happened:")
print(get_info_str(counter), end="")

if args.wcg:
    fname = pathlib.Path(args.out_fname).with_suffix(".renaming.json")
    print("Writing %s" % (fname))
    with open(fname, "w") as f:
        json.dump(supplier.names, f)
    print("Done.")
