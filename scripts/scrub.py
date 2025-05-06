#!/usr/bin/env python

import argparse
import io
import json
import multiprocessing
from os import linesep
import pathlib
import sys

from molscrub import Scrub
from molscrub import SMIMolSupplierWrapper

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
        nr_isomers = len(mol_group)
        serial_suffix = 0
        for i, mol in enumerate(mol_group):
            nr_confs = mol.GetNumConformers()
            for j in range(mol.GetNumConformers()):
                mol.SetProp("ScrubInfo", json.dumps({
                    "isomerGroup": self.counter_mol_group,
                    "isomerId": i,
                    "confId": j,
                    "nr_conformers:":  mol.GetNumConformers(),
                    "nr_isomers:":  nr_isomers,
                }))
                if add_serial_suffix:
                    serial_suffix += 1
                    if nr_isomers > 1 or nr_confs > 1:
                        mol.SetProp("_Name", name + "_%d" % serial_suffix)
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

    def write_mols(self, mol_group, add_suffix=False, add_serial_suffix=False):
        nr_isomers = len(isomer_list)
        serial_suffix = 0
        for i, mol in enumerate(mol_group):
            nr_confs = mol.GetNumConformers()
            if add_serial_suffix:
                if nr_isomers > 1 or nr_confs > 1:
                    serial_suffix += 1
                    mol.SetProp("_Name", name + "_%d" % serial_suffix)
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

parser_essential = argparse.ArgumentParser(description="Protonate molecules and add 3D coordinates", add_help=False)

parser_essential.add_argument("input", help="input filename (.sdf/.mol/.smi/.cxsmiles) or SMILES string")

basic = parser_essential.add_argument_group("options")
basic.add_argument("-o", "--out_fname", help="output filename (.sdf/.hdf5)", required=True)
basic.add_argument("--write_failed_mols", help="filename for failed molecules (.sdf)")
basic.add_argument("--name_from_prop", help="set molecule name from RDKit/SDF property")
basic.add_argument("--ph", help="pH value for acid/base transformations", default=7.4, type=float)
basic.add_argument("--skip_acidbase", help="skip enumeration of acid/base conjugates", action="store_true")
basic.add_argument("--skip_tautomers", help="skip enumeration of tautomers", action="store_true")
basic.add_argument("--skip_ringfix", help="skip fixes of six-member rings", action="store_true")
basic.add_argument("--skip_gen3d", help="skip generation of 3D coordinates (also skips ring fixes)", action="store_true")

misc = parser_essential.add_argument_group("miscellaneous")
misc.add_argument("--cpu", help="number of processes to run in parallel", default=0, type=int)
misc.add_argument("--debug", help="errors are raised", action="store_true")
misc.add_argument("-h", "--help", help="show this help message and exit", action="help")
misc.add_argument("--help_advanced", help="show advanced options and exit", action="store_true")

parser_advanced = argparse.ArgumentParser() # for --help_advanced

acidbase = parser_advanced.add_argument_group("acid base enumeration")
acidbase.add_argument("--ph_low", help="low end of pH range (superseeds --ph)", type=float)
acidbase.add_argument("--ph_high", help="high end of pH range (superseeds --ph)", type=float)

geom = parser_advanced.add_argument_group("3D coordinates")
geom.add_argument("--max_ff_iter", help="maximum number of force field optimization steps", type=int, default=400)
geom.add_argument("--numconfs", help="Number of conformers to generate", type=int)
geom.add_argument("--etkdg_rng_seed", help="seed for random number generator used in ETKDG", type=int)
geom.add_argument("--ff", help="uff, mmff94, mmff94s, espaloma", choices=["uff", "mmff94", "mmff94s","espaloma"], default="mmff94s")
geom.add_argument("--template", help="Template molecule for 3D embedding with constraints")
geom.add_argument("--template_smarts", help="SMARTs patter matching atoms of template and query molecules for 3D embedding")
geom.add_argument("--ring_minimize", help="use FF energy minimization to determine optimal ring conformer", action="store_true")
geom.add_argument("--energy_threshold", help="energy threshold for conformer distinction", default=0.5, type=float)

misc2 = parser_advanced.add_argument_group("more miscellaneous options")
misc2.add_argument("--wcg", help="make sure mol names and suffixes are integers", action="store_true")

if "--help_advanced" in sys.argv:
    parser_essential.print_help()
    f = io.StringIO()
    parser_advanced.print_help(f)
    f.seek(0)
    advanced_help = f.read()
    advanced_help = linesep + linesep.join(advanced_help.split(linesep)[5:-1])
    print(advanced_help)
    sys.exit()

args_essential, remaining_args = parser_essential.parse_known_args()
args_advanced = parser_advanced.parse_args(remaining_args)
args = argparse.Namespace(**vars(args_essential), **vars(args_advanced))

if args.ph_low is None and args.ph_high is None:
    ph_low = args.ph
    ph_high = args.ph
elif args.ph_low is not None and args.ph_high is not None:
    ph_low = args.ph_low
    ph_high = args.ph_high
else:
    print("--ph_low and --ph_high work together, either use both or none.")
    sys.exit()

force_single_process = False
if args.ff == "espaloma":
    if args.cpu > 1:  # default is zero
        print("--ff espaloma can't be used with multiprocessing")
        sys.exit(2)
    if args.cpu == 0:
        print("will use only one process because of espaloma")
        force_single_process = True

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

if args.template is not None:
    extension_template = pathlib.Path(args.template).suffix
    if extension_template == ".sdf":
        template_mol = Chem.SDMolSupplier(args.template, removeHs=True)[0]
    elif extension_template == ".mol":
        template_mol = Chem.MolFromMolFile(args.template, removeHs=True)
    else:
        print("You must provide a template with 3D coordinates in .sdf or .mol format")
        sys.exit()
else:
    template_mol = None

if args.template_smarts is not None:
    if args.template is None:
        print('If passing a template SMARTs you must provide a template molecule first.')
        sys.exit()
    if not isinstance(args.template_smarts, str):
        print("Template SMARTs must be a string")
    else:
        template_smarts = Chem.MolFromSmarts(args.template_smarts)
    if template_smarts is None:
        print('The provided SMARTs could not be converted into a molecule. Check your inputs.')
        sys.exit()
else:
    template_smarts = None

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


# ring_minimize and espaloma incompatible for now. 
if args.ring_minimize and args.ff == "espaloma":
    # use colors
    error_message = "\x1b[31mring_minimize and ff=espaloma are incompatible... for now.\033[0m"
    raise ValueError(error_message)
 
# if ring_minimize is chosen, then numconfs is automatically 3
if args.ring_minimize and args.numconfs is None:
    nconfs = 3
else:
    if args.numconfs is None:
        nconfs = 1
    else:
        nconfs = args.numconfs



scrub = Scrub(
    ph_low,
    ph_high,
    pka_fname=None,
    tauto_fname=None,
    skip_acidbase=args.skip_acidbase,
    skip_tautomers=args.skip_tautomers,
    skip_ringfix=args.skip_ringfix,
    skip_gen3d=args.skip_gen3d,
    template=template_mol,
    template_smarts=template_smarts,
    do_gen2d=do_gen2d,
    max_ff_iter=args.max_ff_iter,
    numconfs = nconfs,
    etkdg_rng_seed=args.etkdg_rng_seed,
    ff=args.ff,
    ring_minimize=args.ring_minimize,
    energy_threshold=args.energy_threshold,
    debug=args.debug
)

counter = {
    "supplied": 0,
    "rdkit_nope": 0,
    "ok_mols": 0,
    "isomers": 0,
    "conformers": 0,
    "failed": 0,
}

def scrub_and_catch_errors(input_mol, sdwriter_failed_mols=None):
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
            if sdwriter_failed_mols is not None:
                input_mol.SetProp("exception", str(e))
                sdwriter_failed_mols.write(input_mol)
                sdwriter_failed_mols.flush() # slow?
    return (isomer_list, log)

def scrub_and_debug(input_mol, _=None):
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

if args.debug and args.write_failed_mols:
    print("--write_failed_mols does not work with --debug, exiting", file=sys.stderr)
    sys.exit(2)
    scrub_fn = scrub_and_debug
    sdwriter_failures = None
elif args.debug:
    scrub_fn = scrub_and_debug
    sdwriter_failures = None
elif args.write_failed_mols is not None:
    if args.cpu != 1:
        print("--write_failed_mols does not work with multiprocessing, needs --cpu 1, exiting", file=sys.stderr)
        sys.exit(2)
    scrub_fn = scrub_and_catch_errors
    sdwriter_failures = Chem.SDWriter(args.write_failed_mols)
else:
    scrub_fn = scrub_and_catch_errors
    sdwriter_failures = None

if __name__ == '__main__':
    with Writer(args.out_fname) as w:
        if args.cpu == 1 or force_single_process:
            for input_mol in supplier:
                isomer_list, log = scrub_fn(input_mol, sdwriter_failures)
                write_and_log(isomer_list, log, counter)
        else:
            if args.cpu < 1:
                nr_proc = multiprocessing.cpu_count()
            else:
                nr_proc = args.cpu
            p = multiprocessing.Pool(nr_proc - 1) # leave 1 for main process
            for (isomer_list, log) in p.imap_unordered(scrub_fn, supplier):
                write_and_log(isomer_list, log, counter)


    if sdwriter_failures is not None:
        sdwriter_failures.close()

    print("Scrub completed.\nSummary of what happened:")
    print(get_info_str(counter), end="")

    if args.wcg:
        fname = pathlib.Path(args.out_fname).with_suffix(".renaming.json")
        print("Writing %s" % (fname))
        with open(fname, "w") as f:
            json.dump(supplier.names, f)
        print("Done.")
