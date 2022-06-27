import sys
from io import BytesIO
import pickle
import re

# TODO use cPickle
import threading
import multiprocessing
import queue
import os
from rdkit import Chem, RDLogger
from rdkit.Chem.PropertyMol import PropertyMol

from .common import ScrubberClass

""" this file contains all the  molecule providers
    - files
    - pipe streams
    - web services
    - ...

should be implemented here

"""
VALID_FORMATS = ["smi", "sdf"]


class MoleculeProvider(ScrubberClass):
    """This is the class of Molecule Providers iterator, which behave all the same:
    - instanciated optionally with a file
    - can be iterated returning a molecule at each step

    - an optional property can be specified to set the molecule name (often the vendor catalog id)
    - if the molecule has no name ("_Name" property), a default ("MOL") will be assigned

    """

    # TODO add multithreaded mol supplier
    #  - https://iwatobipen.wordpress.com/2021/05/04/read-sdf-with-multi-thread-rdkit-memo-chemoinformatics/
    # TODO add gz support (custom func for extension)
    default_mol_name = "MOL"

    def __init__(
        self,
        fname=None,
        ftype=None,
        pipe: bool = False,
        sanitize: bool = True,
        removeHs: bool = False,
        strictParsing: bool = True,
        preserve_properties: bool = True,
        name_property:str= None,
        wrapped: bool = True,
        discarded_data_file: str = None,
        # use_PropertyMol: bool = True,
        start_count: int = 0,
        end_count: int = -1,
        _stop_at_defaults: bool = False,
    ):
        self.fname = fname
        self.ftype = ftype
        self.pipe = pipe
        self.sanitize = sanitize
        self.removeHs = removeHs
        self.strictParsing = strictParsing
        self.preserve_properties = preserve_properties
        self.name_property = name_property
        self.wrapped = wrapped
        self.discarded_data_file = discarded_data_file
        # self.use_PropertyMol = use_PropertyMol
        self.start_count = start_count
        self.end_count = end_count
        if _stop_at_defaults:
            return
        self._counter = 0
        self._counter_problematic = 0
        # self._build_opts_dict()
        if not int(self.fname is not None) + (self.pipe) == 1:
            msg = (
                "Error: either filename or pipe mode must be "
                "specified (current: fname:%s, pipe:%s)" % (self.fname, str(self.pipe))
            )
            raise ValueError(msg)
        # check extension and activate proper source
        if pipe:
            self._source = PipeMolSupplier2()
        elif fname is not None:
            print("[ FILE MODE ]")
            name, ext = os.path.splitext(fname)
            ext = ext[1:].lower()
            if self.ftype == None:
                self.ftype = ext
            if not self.ftype in VALID_FORMATS:
                msg = (
                    "Error: the specified format [ %s ] is not valid. "
                    "Allowed values: [ %s ]" % (self.ftype, ",".join(VALID_FORMATS))
                )
                raise ValueError(msg)
            if self.ftype == "sdf":
                if self.wrapped:
                    self._source = SDFMolSupplierWrapper(
                        fname,
                        sanitize=self.sanitize,
                        removeHs=self.removeHs,
                        strictParsing=self.strictParsing,
                        discarded_data_file=self.discarded_data_file,
                    )
                else:
                    self._source = Chem.SDMolSupplier(
                        fname,
                        sanitize=self.sanitize,
                        removeHs=self.removeHs,
                        strictParsing=self.strictParsing,
                    )
            elif self.ftype == "smi":
                if self.wrapped:
                    self._source = SMIMolSupplierWrapper(
                        fname, sanitize=self.sanitize, titleLine=False
                    )
                else:
                    self._source = Chem.SmilesMolSupplier(
                        fname,
                        sanitize=self.sanitize,
                        titleLine=False,
                        discarded_data_file=self.discarded_data_file,
                    )

    def __iter__(self):
        return self

    def __next__(self):
        """return the next molecule"""
        try:
            while True:
                self._counter += 1
                next_mol = self._source.__next__()
                if self._counter < self.start_count:
                    # print("SKIPPING START")
                    continue
                if (not self.end_count == -1) and (self._counter > self.end_count):
                    # print("SKIPPING END")
                    if self.wrapped:
                        self._source._close_fp()
                    raise StopIteration
                if not next_mol is None:
                    next_mol = PropertyMol(next_mol)
                    # name priority will be given to the requested property
                    if not self.name_property is None:
                        if next_mol.HasProp(self.name_property):
                            # preserve original name
                            if next_mol.HasProp("_Name"):
                                next_mol.SetProp("Original_name_from_input", mol.GetProp("_Name"))
                            # rename the molecule using the requested field
                            next_mol.SetProp("_Name", mol.GetProp(self.name_property))
                        # else:
                        #     # set default molecule name
                        #     next_mol.SetProp("_Name", self.default_mol_name)
                    # if the molecule has no name, a default one will be assigned
                    if not next_mol.HasProp("_Name"):
                        next_mol.SetProp("_Name", self.default_mol_name)
                    next_mol.SetProp("_scrubber_source", "input")
                    return (self._counter, next_mol)
                else:
                    self._counter_problematic += 1
        except StopIteration:
            if self.wrapped:
                self._source._close_fp()
            raise StopIteration


class MoleculeStorage(ScrubberClass, multiprocessing.Process):
    """Class to write molecules processed;

    destinations can be : file, stdout, dbaseo

    # write molecules to a SMILES file
    >>> ms = MoleculeStorage(out_fname = 'output.txt', out_fname_format = 'smi')

    # write molecules in a SDF file
    >>> ms = MoleculeStorage(out_fname = 'output.sdf', out_fname_format = 'auto')

    # write molecules in to a database (must have a compatible write() method)
    >>> ms = MoleculeStorage(out_container = db_object)

    # write to STOUD
    >>> ms = MoleculeStorage(pipe=True)

    """

    # TODO this should be split into an abstract class and MoleculeFileStorage, MoleculePipeStorage, etc...

    out_format_opts_default = {
        "single": {"smi": {}, "sdf": {}},
        "split": {"smi": {}, "sdf": {}},
    }
    out_format_file_writers = {
        "single": {"smi": Chem.SmilesWriter, "sdf": Chem.SDWriter},
        "split": {"smi": Chem.MolToSmiles, "sdf": Chem.MolToMolBlock},
    }

    def __init__(
        self,
        mode: str = "single",  # "single", "split", "pipe"
        out_fname: str = None,
        out_format: str = None,  # smi, sdf, None (auto)
        out_format_opts: dict = out_format_opts_default,
        naming: str = "auto",  # auto, (progressive), name (molname)
        naming_field: str = "",
        max_lig_per_dir=0,  # it not 0, create automatically subdirectories each [max_lig_per_dir] ligands
        sanitize_name: bool = True,  # remove all
        queue: multiprocessing.Queue = None,
        preserve_properties: bool = True,  # preserve any extra properties found in the molecule
        comm_pipe: multiprocessing.Pipe=None,
        workers_count: int = 1,
        # overwrite:bool=False,
        disable_rdkit_warnings: bool = True,
        _stop_at_defaults=False,
    ):
        """
        alternative output modes:

            - PIPE
                done
            - SINGLE (FILE)
                - OUTFILE = filename.ext
                - OUTFORMAT = file format | .ext
                - WRITER = persistent (OUTFILE)
            - SPLIT (FILE)
                 (case 1)
                    OUTFILE = None
                    OUTFORMAT = file format
                    OUTFILE = mol_# | mol_[PROPERTY]
                    WRITER = disposable (OUTFILE)
                 (case 2)
                    OUTFILE = filename.ext
                    OUTFORMAT = file format | .ext
                    OUTFILE = filename_#.ext | filename_[PROPERTY].ext
                    WRITER = disposable (OUTFILE)

        """
        self.mode = mode
        self.out_fname = out_fname
        self.out_format = out_format
        self.out_format_opts = out_format_opts
        self.naming = naming
        # self.naming_field = naming_field
        self.max_lig_per_dir = max_lig_per_dir
        self.sanitize_name = sanitize_name
        self.queue = queue
        self.preserve_properties = preserve_properties
        self.comm_pipe = comm_pipe
        self.workers_count = workers_count
        # self.overwrite = overwrite
        self.disable_rdkit_warnings = disable_rdkit_warnings
        if _stop_at_defaults:
            return
        if self.disable_rdkit_warnings:
            RDLogger.DisableLog("rdApp.*")
        multiprocessing.Process.__init__(self)
        self._counter = 0
        self._dir_counter = 0
        self.writer = None
        # in single mode, filename is mandatory
        if self.mode == "single":
            if not self.naming == "auto":
                print(
                    "Warning: naming scheme requested [%s] will be ignored in 'single' mode"
                    % self.naming
                )
            if self.out_fname is None:
                raise ValueError("Filename must be specified in 'single' mode")
            self._basename, self._ext = os.path.splitext(self.out_fname)
            self._ext = self._ext[1:].lower()
            if self.out_format is None:
                self.out_format = self._ext
            if not self.out_format in VALID_FORMATS:
                raise ValueError(
                    "Invalid outpuf file format: current(%s), "
                    "accepted (%s)" % (self._ext, ",".join(VALID_FORMATS))
                )
            self.writer = self.out_format_file_writers[self.mode][self.out_format](
                self.out_fname, **self.out_format_opts[self.mode][self.out_format]
            )
        # in single mode many things can happen...
        elif self.mode == "split":
            if self.out_fname is None:
                self._basedir = os.path.curdir
                self._basename = []
                self._ext = None
            else:
                self._basedir = os.path.dirname(self.out_fname)
                name, ext = os.path.splitext(self.out_fname)
                self._basename = [os.path.basename(name)]
                self._ext = ext[1:].lower()
            if self.out_format is None:
                self.out_format = self._ext
            if not self.out_format in VALID_FORMATS:
                raise ValueError(
                    "Invalid outpuf file format: current(%s), "
                    "accepted (%s)" % (self._ext, ",".join(VALID_FORMATS))
                )

        elif self.mode == "pipe":
            # write to pipe
            pass

    def run(self):
        """multithreading default function with listening loop that waits for
        molecules to be written"""
        # TODO handle these:
        # - subdir creations
        # - duplicate fnames
        # print("THE WRITER IS IN BUSINESS!!!")
        while True:
            # TODO addd try/except for thread safety
            try:
                package = self.queue.get()
            except KeyboardInterrupt:
                if not self.writer is None:
                    self.writer.close()
                print("writer done (Ctrl-C)")
                return
            if package is None:
                # poison pill
                self.workers_count -= 1
                if self.workers_count == 0:
                    # TODO this might be superfluouus
                    if not self.writer is None:
                        self.writer.close()
                    if not self.comm_pipe is None:
                        self.comm_pipe.send(self._counter)
                    return
                else:
                    continue
            try:
                mol = package["mol"]
            except:
                print("PROBLEMATIC PACKAGE!", package)
                sys.exit(1)
            if self.mode == "single":
                # save all non-private properties in the current molecule
                if self.preserve_properties and self.out_format == "sdf":
                    # this triggers a warning if done after the initialization
                    # of the writer therefore it requires the suppression of
                    # all the RDKit warnings
                    if not self.disable_rdkit_warnings:
                        RDLogger.DisableLog("rdApp.*")
                    self.writer.SetProps(mol.GetPropNames())
                    if not self.disable_rdkit_warnings:
                        RDLogger.EnableLog("rdApp.*")
                self._counter+=1
                self.writer.write(mol)
            elif self.mode == "split":
                outfname = self._get_outfname(mol)
                writer = self.out_format_file_writers["single"][self.out_format]
                # this triggers a warning if done after the initialization
                # of the writer therefore it requires the suppression of
                # all the RDKit warnings
                # if not self.disable_rdkit_warnings:
                #     RDLogger.DisableLog("rdApp.*")
                # writer.SetProps(mol.GetPropNames())
                # if not self.disable_rdkit_warnings:
                #     RDLogger.EnableLog("rdApp.*")
                self._counter+=1
                with writer(outfname) as fp:
                    fp.write(mol)
                # with open(outfname, "w") as fp:
                #     fp.write(writer(mol, **self.out_format_opts[self.mode][self._ext]))
            elif self.mode == "pipe":
                pickle.dump(mol, sys.stdout.buffer)
            # self._counter += 1

    def _get_outfname(self, mol,):
        """function to automate the output molecule name, and generate
        progressive subdirectories, if required"""

        # TODO check for Scrubber properties here, if they can be used for naming, e.g.: p1_t2_s4

        default_name = "MOL"
        dirname = [self._basedir]
        ext = self.out_format
        # get the name
        basename = []
        # generate the file name
        if self.naming == "name":
            # extract the field from RDKit molecule
            # print("MOLSZZZ", mol.GetPropsAsDict())
            basename = mol.GetProp("_Name")
        # elif self.naming == "field":
        #     # extract the name from RDKit molecule
        #     try:
        #         basename = mol.GetProp(self.naming_field).strip()
        #     except KeyError:
        #         basename = ""
        if (self.naming == "auto") or (len(basename) == 0):
            if len(self._basename):
                basename = "%s_%s" % (self._basename[0], self._counter)
            else:
                basename = "%s_%s" % (default_name, self._counter)
        if self.sanitize_name:
            basename = self._sanitize_string(basename)
        # generate the directory name
        if self.max_lig_per_dir > 0:
            if self._counter % self.max_lig_per_dir == 0:
                self._dir_counter += 1
            new_dir = "{:0>8}".format(self._dir_counter - 1)
            dirname.append(new_dir)
            dirname = os.path.sep.join(dirname)
        else:
            dirname = dirname[0]
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        # else:
        #     if not self.overwrite:
        #         raise FileExistsError
        fullpath = "%s%s%s.%s" % (dirname, os.path.sep, basename, ext)
        if os.path.exists(fullpath):
            attempt = 0
            while os.path.exists(fullpath):
                attempt += 1
                fullpath = "%s%s%s_v%d.%s" % (
                    dirname,
                    os.path.sep,
                    basename,
                    attempt,
                    ext,
                )
                # fullpath = "%s_v%d.%s" % (os.path.sep.join(dirname + basename), attempt, ext)
        return fullpath

    def _sanitize_string(self, string):
        """function to apply rules to generate a valid filename from a
        molecular name"""
        # replace spaces and tabs with underscore
        # print("----------------")
        # print("RECEIVED", string)
        string = re.sub("\s+","_", string)
        string = re.sub("\t+","_", string)
        string = string.replace("{", "(")
        string = string.replace("[", "(")
        string = string.replace("]", ")")
        string = string.replace("}", ")")
        # print("PRODUCED", string)
        # print("----------------")
        return string


class SDFMolSupplierWrapper(object):
    """RDKit SDF molecule wapper to provide molecules in a robust way.
    If an error is encountered when parsing the molecule, the raw text is dumped in a log file

    """

    def __init__(
        self,
        filename: str,
        sanitize: bool = True,
        removeHs: bool = False,
        strictParsing: bool = True,
        discarded_data_file: str = None,
        _stop_at_defaults: bool = False,
    ):
        self.filename = filename
        self.sanitize = sanitize
        self.removeHs = removeHs
        self.strictParsing = strictParsing
        self.discarded_data_file = discarded_data_file
        if _stop_at_defaults:
            return
        self.fp_input = open(filename, "r")
        self.fp_errors = None
        self._buff = []

    def __iter__(self):
        return self

    def _close_fp(self):
        """close all potential file pointers"""
        self.fp_input.close()
        if not self.fp_errors is None:
            self.fp_errors.close()

    def __next__(self):
        """iterator step"""
        while True:
            line = self.fp_input.readline()
            # empty line
            if not line:
                # last molecule in the file
                if len(self._buff):
                    try:
                        # buffer full, returning last molecule
                        mol = Chem.MolFromMolBlock(
                            "".join(self._buff),
                            sanitize=self.sanitize,
                            removeHs=self.removeHs,
                            strictParsing=self.strictParsing,
                        )
                        self._buff = []
                        break  # go to last try/except
                    except:
                        self._counter_problematic += 1
                        if self.discarded_data_file is None:
                            continue
                        if self.fp_errors is None:
                            self.fp_errors = open(self.discarded_data_file, "w")
                        buff = "\n".join(self._buff)
                        self.fp_errors.write(buff + "\n")
                        return None
                else:
                    # buffer empty, stopping the iteration
                    self._close_fp()
                    raise StopIteration
            else:
                self._buff.append(line)
                if "$$$$" in line:
                    # molecule completed
                    complete_buff = self._buff[:]
                    self._buff = []
                    break  # go to last try/except...
        try:
            # found a complete molecule, returning it
            mol = Chem.MolFromMolBlock(
                "".join(complete_buff),
                sanitize=self.sanitize,
                removeHs=self.removeHs,
                strictParsing=self.strictParsing,
            )
            return mol
        except:
            if not self.discarded_data_file is None:
                if self.fp_errors is None:
                    self.fp_errors = open(self.discarded_data_file)
                buff = "\n".join(self._buff)
                self.fp_errors.write(buff + "\n")
            return None


class SMIMolSupplierWrapper(object):
    """RDKit SMI molecule supplier wrapper."""

    def __init__(
        self,
        filename: str,
        sanitize: bool = True,
        titleLine: bool = False,
        discarded_input_fname: str = None,
        _stop_at_defaults: bool = False,
    ):
        self.filename = filename
        self.sanitize = sanitize
        self.titleLine = titleLine
        self.discarded_input_fname = discarded_input_fname
        if _stop_at_defaults:
            return
        self.fp_input = open(filename, "r")
        self.fp_errors = None
        self._buff = []
        self._first_line = False
        print("INITIALIZED", self.titleLine)

    def _close_fp(self):
        """close all potential file pointers"""
        self.fp_input.close()
        if not self.fp_errors is None:
            self.fp_errors.close()

    def __iter__(self):
        return self

    def __next__(self):
        """iterator step"""
        while True:
            line = self.fp_input.readline()
            if self.titleLine and not self._first_line:
                self._first_line = True
                continue
            # empty line
            if not line:
                self.fp_input.close()
                if not self.fp_errors is None:
                    self.fp_errors.close()
                raise StopIteration
            try:
                mol = Chem.MolFromSmiles(line, sanitize=self.sanitize)
                return mol
            except:
                if self.discarded_input_fname is None:
                    continue
                if self.fp_errors is None:
                    self.fp_errors = open(self.discarded_input_fname)
                self.fp_errors.write(line + "\n")


class PipeMolSupplier(threading.Thread):
    """Threded pipe reader
    https://www.askpython.com/python/oops/threading-with-classes

    """

    PICKLE_HEADER = pickle.dumps(None)[0]

    def __init__(self, max_proc=None):
        threading.Thread.__init__(self)
        if max_proc is None:
            max_proc = multiprocessing.cpu_count()
        self.buffer = queue.Queue(maxsize=max_proc)
        self._stop = False
        # start the thread
        # self.start()

    def run(self):
        unread_bytes = sys.stdin.buffer.read()
        while len(unread_bytes):
            unread_buffer = BytesIO(unread_bytes)
            header_byte = unread_buffer.getbuffer()[0]
            if header_byte == PipeMolSupplier.PICKLE_HEADER:
                # try:
                io = pickle.load(unread_buffer)
                # if not isinstance(io, str ):
                #     break
                self.buffer.put(io)
                # print(">>>SENDING THIS:", io)
            # except pickle.UnpicklingError:
            #     print("FATAL ERROR UNPICKLING THIS MOLECULE!")
            # sys.exit(1)
            else:
                # just passing through...
                print("UNREADEX", unread_bytes.decode())
                unread_bytes = unread_buffer.read()
                # break
            unread_bytes = unread_buffer.read()
        print("DONE reading...")
        # raise StopIteration

    def __next__(self):
        """overloading this function"""
        # print("> called get buffer")
        print("LEN QUEUE", self.buffer.qsize())
        out = self.buffer.get()
        if out == "terminate":
            raise StopIteration
        else:
            return self.buffer.get()


# TODO: implement a SQLite tmp-file to store all molecules seen so far?
# - remove duplicates
# - prevent file names collisions and overwriting


class PipeMolSupplier2:
    """ """

    PICKLE_HEADER = pickle.dumps(None)[0]

    def __init__(self, sanitize=True, removeHs=False, cleanupSubstructure=True):
        self._opts = {
            "sanitize": sanitize,
            "removeHs": removeHs,
            "cleanupSubstructures": cleanupSubstructure,
        }
        self._counter = 0

    def __iter__(self):
        """standard iterator definition"""
        return self

    def __next__(self):
        """iterator reader"""
        print("CALLED NEXT", self._counter)
        while True:
            unread_bytes = sys.stdin.buffer.read()
            self._counter += 1
            # if self._counter == 6:
            #     print("*************  TOO MANY")
            #     raise StopIteration
            # while unread_bytes:
            print("ITERATING...")
            # self._counter += 1
            if not len(unread_bytes):
                print("DONE READONG", self._counter)
                raise StopIteration
            unread_buffer = BytesIO(unread_bytes)
            header_byte = unread_buffer.getbuffer()[0]
            if header_byte == PipeMolSupplier.PICKLE_HEADER:
                try:
                    mol = pickle.load(unread_buffer)
                    # print("THIS IS PICKLED!", io)
                    # unread_bytes = sys.stdin.buffer.read()
                    print("returning mol to the queue", mol)
                    return mol

                    # return mol
                except pickle.UnpicklingError:
                    print("ERROR UNPICKLING THIS MOLECULE!", self._counter)
            else:
                # just passing through...
                # unread_bytes = sys.stdin.buffer.read()
                print(q)
                # print(">>PRINT", unread_bytes.decode())
                return unread_bytes.decode()


# class SimpleSmilesReader(object):
#     """simple SMILES file iterator"""
#     def __init__(self, fname, sanitize=True):
#         self._fp = open(fname, 'w')
#         return

#     def __iter__(self):
#         return self

#     def __next__(self):
#         try:
#         # line = self._fp.readline()
#         # if not line == "":
#             smiles, *name , *_ = line.split(None, 3)


if __name__ == "__main__":
    pr = PipeReader(None)
    pr.start()
