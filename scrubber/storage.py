import sys
from io import BytesIO
import pickle
import re

import threading
import multiprocessing
import queue
import os
from rdkit import Chem, RDLogger
from rdkit.Chem.PropertyMol import PropertyMol

from .common import ScrubberBase

""" this file contains all the  molecule providers
    - files
    - pipe streams
    - web services
    - ...

should be implemented here

"""
VALID_FORMATS = ["smi", "sdf"]

# TODO use cPickle for pipes
# TODO implement a non-strict parsing to use molvis-like fixing operations (i.e.: nitro-fixes?)


class MoleculeProvider(ScrubberBase):
    """This is the class of Molecule Providers iterator, which behave all the same:
    - instanciated optionally with a file
    - can be iterated returning a molecule at each step

    - an optional property can be specified to set the molecule name (often the vendor catalog id)
    - if the molecule has no name ("_Name" property), a default ("MOL") will be assigned

    """

    # TODO add multithreaded mol supplier
    #  - https://iwatobipen.wordpress.com/2021/05/04/read-sdf-with-multi-thread-rdkit-memo-chemoinformatics/
    # TODO add gz support (custom func for extension)
    # TODO add desalting method/options
    default_mol_name = "MOL"

    def __init__(
        self,
        fname=None,
        ftype=None,
        use_pipe: bool = False,  # not yet supported...
        sanitize: bool = True,
        removeHs: bool = False,
        strictParsing: bool = True,
        # preserve_properties: bool = True,
        name_property: str = None,
        safeparsing: bool = True,
        discarded_datafile: str = None,
        queue_err: multiprocessing.Queue = None,
        pipe_comm: multiprocessing.Pipe = None,
        # use_PropertyMol: bool = True,
        start_count: int = 0,
        end_count: int = -1,
        _stop_at_defaults: bool = False,
    ):
        self.fname = fname
        self.ftype = ftype
        self.use_pipe = use_pipe
        self.sanitize = sanitize
        self.removeHs = removeHs
        self.strictParsing = strictParsing
        # self.preserve_properties = preserve_properties
        self.name_property = name_property
        self.safeparsing = safeparsing
        self.discarded_datafile = discarded_datafile
        self.queue_err = queue_err
        self.pipe_comm = pipe_comm
        # self.use_PropertyMol = use_PropertyMol
        self.start_count = start_count
        self.end_count = end_count
        if _stop_at_defaults:
            return
        self._counter = 0
        self._counter_problematic = 0
        # self._build_opts_dict()
        if not int(self.fname is not None) + (self.use_pipe) == 1:
            msg = (
                "Error: either filename or pipe mode must be "
                "specified (current: fname:%s, pipe:%s)" % (self.fname, str(self.use_pipe))
            )
            raise ValueError(msg)
        # check extension and activate proper source
        if self.use_pipe:
            self._source = PipeMolSupplier2()
        elif fname is not None:
            #print("[ FILE MODE ]")
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
                if self.safeparsing:
                    self._source = SDFMolSupplierWrapper(
                        fname,
                        sanitize=self.sanitize,
                        removeHs=self.removeHs,
                        strictParsing=self.strictParsing,
                        discarded_datafile=self.discarded_datafile,
                        queue_err=self.queue_err,
                    )
                else:
                    self._source = Chem.SDMolSupplier(
                        fname,
                        sanitize=self.sanitize,
                        removeHs=self.removeHs,
                        strictParsing=self.strictParsing,
                    )
            elif self.ftype == "smi":
                if self.safeparsing:
                    self._source = SMIMolSupplierWrapper(
                        fname,
                        sanitize=self.sanitize,
                        titleLine=False,
                        queue_err=self.queue_err,
                    )
                else:
                    self._source = Chem.SmilesMolSupplier(
                        fname,
                        sanitize=self.sanitize,
                        titleLine=False,
                        discarded_datafile=self.discarded_datafile,
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
                if (not self.end_count == -1) and (self._counter == self.end_count):
                    # print("SKIPPING END")
                    if self.safeparsing:
                        self._source._close_fp()
                        if not self.pipe_comm is None:
                            self.pipe_comm.send(
                                "input:%d" % (self._counter - self.start_count)
                            )
                    raise StopIteration
                if not next_mol is None:
                    next_mol = PropertyMol(next_mol)
                    # name priority will be given to the requested property
                    if not self.name_property is None:
                        if next_mol.HasProp(self.name_property):
                            # preserve original name
                            if next_mol.HasProp("_Name"):
                                next_mol.SetProp(
                                    "Original_name_from_input", mol.GetProp("_Name")
                                )
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
            if self.safeparsing:
                self._source._close_fp()
            if not self.pipe_comm is None:
                self.pipe_comm.send("input:%d" % (self._counter - self.start_count))
            raise StopIteration


class MoleculeStorage(ScrubberBase, multiprocessing.Process):
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
        fname: str = None,
        ftype: str = None,  # smi, sdf, None (auto)
        mode: str = "single",  # "single", "split", "pipe"
        format_opts: dict = out_format_opts_default,
        naming: str = "auto",  # auto, (progressive), name (molname)
        naming_field: str = "",
        max_lig_per_dir=0,  # it not 0, create automatically subdirectories each [max_lig_per_dir] ligands
        disable_name_sanitize: bool = False,  # disable sanitizing output filename (based on mol name)
        disable_preserve_properties: bool = False,  # disable preserving any extra properties found in the molecule
        num_poison_pills_needed: int = 1,
        # disable_rdkit_warnings: bool = True,
        queue: multiprocessing.Queue = None,
        comm_pipe: multiprocessing.Pipe = None,
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
        self.fname = fname
        self.ftype = ftype
        self.mode = mode
        self.format_opts = format_opts
        self.naming = naming
        self.max_lig_per_dir = max_lig_per_dir
        self.disable_name_sanitize = disable_name_sanitize
        self.disable_preserve_properties = disable_preserve_properties
        self.num_poison_pills_needed = num_poison_pills_needed
        self.queue = queue
        self.comm_pipe = comm_pipe
        if _stop_at_defaults:
            return
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
            if self.fname is None:
                raise ValueError("Filename must be specified in 'single' mode")
            self._basename, self._ext = os.path.splitext(self.fname)
            self._ext = self._ext[1:].lower()
            if self.ftype is None:
                self.ftype = self._ext
            if not self.ftype in VALID_FORMATS:
                raise ValueError(
                    "Invalid outpuf file format: current(%s), "
                    "accepted (%s)" % (self._ext, ",".join(VALID_FORMATS))
                )
            self.writer = self.out_format_file_writers[self.mode][self.ftype](
                self.fname, **self.format_opts[self.mode][self.ftype]
            )
        # in single mode many things can happen...
        elif self.mode == "split":
            if self.fname is None:
                self._basedir = os.path.curdir
                self._basename = []
                self._ext = None
            else:
                self._basedir = os.path.dirname(self.fname)
                name, ext = os.path.splitext(self.fname)
                self._basename = [os.path.basename(name)]
                self._ext = ext[1:].lower()
            if self.ftype is None:
                self.ftype = self._ext
            if not self.ftype in VALID_FORMATS:
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
                print("writer done (Ctrl-C)")
                if not self.writer is None:
                    self.writer.close()
                if not self.comm_pipe is None:
                    self.comm_pipe.send("writer:%d" % self._counter)
                return
            if package is None: # poison pill
                self.num_poison_pills_needed -= 1
                print("number of lives left: %d" % self.num_poison_pills_needed)
                if self.num_poison_pills_needed == 0:
                    # TODO this might be superfluous
                    if not self.writer is None:
                        self.writer.close()
                    if not self.comm_pipe is None:
                        # self.comm_pipe.send(self._counter)
                        self.comm_pipe.send("writer:%d" % self._counter)
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
                if not self.disable_preserve_properties and self.ftype == "sdf":
                    self.writer.SetProps(mol.GetPropNames())
                self._counter += 1
                self.writer.write(mol)
            elif self.mode == "split":
                outfname = self._get_outfname(mol)
                writer = self.out_format_file_writers["single"][self.ftype]
                self._counter += 1
                with writer(outfname) as fp:
                    fp.write(mol)
            elif self.mode == "pipe":
                pickle.dump(mol, sys.stdout.buffer)
            # self._counter += 1

    def _get_outfname(
        self,
        mol,
    ):
        """function to automate the output molecule name, and generate
        progressive subdirectories, if required"""

        # TODO check for Scrubber properties here, if they can be used for naming, e.g.: p1_t2_s4
        default_name = "MOL"
        dirname = [self._basedir]
        ext = self.ftype
        # get the name
        basename = []
        # generate the file name
        if self.naming == "name":
            # extract the field from RDKit molecule
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
        if not self.disable_name_sanitize or True:
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
        # TODO this is inaccurate... does not seem to work...
        # replace spaces and tabs with underscore
        # print("----------------")
        # print(">RECEIVED|%s|" % string)
        string = re.sub("\s+", "_", string)
        string = re.sub("\t+", "_", string)
        string = string.replace("{", "(")
        string = string.replace("[", "(")
        string = string.replace("]", ")")
        string = string.replace("}", ")")
        string = string.replace("/", "=")
        # print("PRODUCED", string)
        # print("----------------")
        # trunkate filenames that are too long
        if len(string) > 100:
            string = string[:100]
        # print(">RETURNNIG|%s|" % string)
        return string


class SDFMolSupplierWrapper(object):
    """RDKit SDF molecule wapper to provide molecules in a robust way.
    If an error is encountered when parsing the molecule, the raw text is
    dumped in a log file and/or passed in the error queue,
    """

    def __init__(
        self,
        filename: str,
        sanitize: bool = True,
        removeHs: bool = False,
        strictParsing: bool = True,
        discarded_datafile: str = None,
        queue_err: multiprocessing.Queue = None,
        _stop_at_defaults: bool = False,
    ):
        self.filename = filename
        self.sanitize = sanitize
        self.removeHs = removeHs
        self.strictParsing = strictParsing
        self.discarded_datafile = discarded_datafile
        self.queue_err = queue_err
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
                        if not self.discarded_datafile is None:
                            if self.fp_errors is None:
                                self.fp_errors = open(self.discarded_datafile, "w")
                            buff = "\n".join(self._buff)
                            self.fp_errors.write(buff + "\n")
                        if not self.queue_err is None:
                            self.queue_err.put(("input", buff))
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
            if not self.discarded_datafile is None:
                if self.fp_errors is None:
                    self.fp_errors = open(self.discarded_datafile)
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
        queue_err: multiprocessing.Queue = None,
        discarded_input_fname: str = None,
        _stop_at_defaults: bool = False,
    ):
        self.filename = filename
        self.sanitize = sanitize
        self.titleLine = titleLine
        self.queue_err = queue_err
        self.discarded_input_fname = discarded_input_fname
        if _stop_at_defaults:
            return
        self.fp_input = open(filename, "r")
        self.fp_errors = None
        self._buff = []
        self._first_line = False
        # print("INITIALIZED", self.titleLine)
        # print("INITIALIZED", self.queue_err)

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
                if mol is None:
                    if not self.queue_err is None:
                        self.queue_err.put(("input", line))
                    if not self.discarded_input_fname is None:
                        if self.fp_errors is None:
                            self.fp_errors = open(self.discarded_input_fname)
                        self.fp_errors.write(line + "\n")
                return mol
            except:
                if not self.queue_err is None:
                    self.queue_err.put(("input", line))
                if not self.discarded_input_fname is None:
                    if self.fp_errors is None:
                        self.fp_errors = open(self.discarded_input_fname)
                    self.fp_errors.write(line + "\n")
                return None


class MoleculeIssueStorage(ScrubberBase, multiprocessing.Process):
    """Multiprocessing class to save problematic data encountered during the processing.
    Problematic data encountered when parsing the input will be written with a
    "raw", writer that will dump the text the wrapped parsers (i.e., non-RDKit
    parsers) tried to process.

    Problematic encountered during the molecule processing will be attempted to
    be saved with RDKit writers

    To prevent the creation of empty files in case of successful processing of
    all data, output files will be written only if at least one molecule will
    be passed to each of the writers
    """

    source_to_writer = {"geom": "format", "input": "raw"}

    def __init__(
        self,
        log_from_input: str = None,
        log_from_process: str = None,
        queue: multiprocessing.Queue = None,
        comm_pipe: multiprocessing.Pipe = None,
        _stop_at_defaults: bool = False,
    ):
        self.log_from_input = log_from_input
        self.log_from_process = log_from_process
        self.queue = queue
        self.comm_pipe = comm_pipe
        if _stop_at_defaults:
            return
        multiprocessing.Process.__init__(self)
        self._counter = 0
        self._input_err_fp = None
        self._process_err_fp_smi = None
        self._process_err_fp_sdf = None
        # early checks: all errors need to be raised as early as possible to
        # prevent spending time on long calculations and failing when
        # encountering errors
        #
        if (self.input_err_fname is None) and (self.process_err_fname is None):
            msg = "Either input_errof_fname or process_err_fname must be specified"
            raise ValueError(msg)
        # check that none of the output files exists
        fname_base, _ = os.path.splitext(self.process_err_fname)
        self._process_err_fname_smi = "%s.smi" % fname_base
        self._process_err_fname_sdf = "%s.sdf" % fname_base
        for fname in self._process_err_fname_smi, self._process_err_fname_smi:
            if os.path.exists(fname):
                msg = "Filename [%s] already exists." % fname
                raise FileExistsError(msg)

    def _close_fp(self):
        """close all open file pointers"""
        if self._process_err_fp_smi is not None:
            self._process_err_fp_smi.close()
        if self._process_err_fp_sdf is not None:
            self._process_err_fp_sdf.close()
        if self._input_err_fp is not None:
            self._input_err_fp.close()

    def run(self):
        """start waiting for input in the queue"""
        try:
            while True:
                package = self.queue.get()
                if package == None:
                    self._close_fp()
                    break
                self._write_package(package)
        except KeyboardInterrupt:
            self._close_fp()

    def _write_package(self, package):
        """Write the problematic molecule according to the data type. Raw input
        data tha couldn't be parsed properly will be written using a raw writer
        that will dump the problematic buffer in the output file. Molecules
        that failed at any other processing step will be written depending on
        their dimensionality: 1D molecules in SMILES files, >1D in SDF files
        """
        source, mol = package
        self._counter += 1
        if isinstance(mol, Chem.rdchem.Mol):
            if self._process_err_fp_sdf is None:
                self._process_err_fp_sdf = Chem.SDWriter(self._process_err_fname_sdf)
            self._process_err_fp_sdf.write(mol)
            if not self.comm_pipe is None:
                self.comm_pipe.send("err_process:%d" % self._counter)
        elif isinstance(mol, str):
            if self._input_err_fp is None:
                self._input_err_fp = open(self.input_err_fname, "w")
            self._input_err_fp.write(mol)
            self._input_err_fp.flush()
            if not self.comm_pipe is None:
                self.comm_pipe.send("err_input:%d" % self._counter)
        else:
            print("THIS SHOULDN'T HAPPEN")
            print(q)


# TODO: implement a SQLite tmp-file to store all molecules seen so far?
# - remove duplicates
# - prevent file names collisions and overwriting


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


if __name__ == "__main__":
    pr = PipeReader(None)
    pr.start()
