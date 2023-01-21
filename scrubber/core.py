import multiprocessing


# IMPORTANT RESOURCES, INVESTIGATE
# https://sefiks.com/2021/07/05/handling-hang-in-python-multiprocessing/
# https://pythonspeed.com/articles/python-multiprocessing/ (alterantive to fork() )

# from time import sleep
# import sys
import random
import pathlib
import time
from rdkit.Chem.PropertyMol import PropertyMol
from rdkit import Chem
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdForceFieldHelpers
from rdkit.Geometry import Point3D

from .storage import MoleculeProvider
from .storage import MoleculeStorage
from .storage import MoleculeIssueStorage
from .geom.geometry import ParallelGeometryGenerator
from .geom.geometry import GeometryGenerator
from .transform.isomer import MoleculeIsomers
from .protonate import AcidBaseConjugator
from .protonate import Tautomerizer
from .common import UniqueMoleculeContainer
from .ringfix import fix_rings


"""
This file contains the core scrubber object
INSPIRATION: https://xkcd.com/1343/
"""

class ScrubberCore(object):
    """read input files, manipulate them to perform the following operations:

        - read input files
            - convert from OB?
        - elaborate structures ( * EXTERNAL * )
            - tautomers, protomers, chiral centers
        - pass molecules to 3D processors ( * EXTERNAL * )
        - write the results

    for mol in input_lib:
        for mol_n in transform(mol):
            gen_3d(mol_n)
            write(mol_n)
    """

    __default_options = {
        "input": {
            "values": MoleculeProvider.get_defaults(),
            "ignore": ["use_pipe", "queue_err", "pipe_comm", "handbrake"],
        },
        "output": {
            "values": MoleculeStorage.get_defaults(),
            "ignore": ["queue", "comm_pipe", "workers_count", "handbrake"],
        },
        "isomers": {
            "active": True,
            "values": MoleculeIsomers.get_defaults(),
            "ignore": [],
        },
        "geometry": {
            "active": True,
            "values": ParallelGeometryGenerator.get_defaults(),
            "ignore": [
                "queue_in",
                "queue_out",
                "queue_err",
                "nice_level",
                "max_proc",
                "handbrake",
            ],
        },
        "general": {
            "values": {
                "max_proc": multiprocessing.cpu_count(),
                "nice_level": None,
            },
            "ignore": [],
        },
        "errors": {
            "values": MoleculeIssueStorage.get_defaults(),
            "ignore": [
                "queue",
                "comm_pipe",
                "handbrake",
            ],
        },
        # -- REMOVE SALTS
    }

    def __init__(self, options: dict = None):
        """this is where the classes doing all operations will do"""
        if options is None:
            options = self.get_defaults()
        self.options = self._conditional_activation_isomers_geometry(options)
        self.geometry_optimize = None
        self.isomer = None
        self._success_cutoff = 0
        self._registered_workers = []
        self._pipe_listener, self._pipe_remote = None, None
        self.queue_out = None
        self.queue_err = None
        self.handbrake = None
        self.mol_provider = None
        self.mol_writer = None
        self.mol_issues = None
        self.max_proc = self.options["general"]["values"]["max_proc"]
        self.counter_data = {
            "input": 0,
            "err_input": 0,
            "err_process": 0,
            "writer": 0,
        }

    def _initialize_mp_pipeline(self):
        """Initialize all objects required to run a parallel pipeline.
        If the optional flag 'queue_err' is provided, an error queue is
        generated even if a log file is not specified"""
        print("[ initializing multiprocessing pipeline ]")
        self.handbrake = multiprocessing.Event()
        self._registered_workers = []
        # initialize pipes for multiprocessing communication
        self._pipe_listener, self._pipe_remote = multiprocessing.Pipe(False)
        # queue out is where processed molecules are sent
        nice = self.options["general"]["values"]["nice_level"]
        # queue that receives the final products of the pipeline
        # if the geometry optimization is not active, it serves also as
        # self._target_queue
        self.queue_out = multiprocessing.Queue(maxsize=self.max_proc * 3)
        ###########################
        # initialize error logging, if necessary
        if self.options["errors"]["values"]["log_basename"] is not None:
            print("[ setting up errors storage ]")
            self.queue_err = multiprocessing.Queue(maxsize=-1)
            self.options["errors"]["values"]["queue"] = self.queue_err
            self.options["errors"]["values"]["comm_pipe"] = self._pipe_remote
            self.options["errors"]["values"]["handbrake"] = self.handbrake
            self.options["errors"]["values"]["handbrake"] = self.handbrake
            self.mol_issues = MoleculeIssueStorage(**self.options["errors"]["values"])
            self.mol_issues.start()
            # self._registered_workers.append(self.mol_issues)
        else:
            self.queue_err = None
            self.mol_issues = None
        ###########################
        # initialize molecular provider and molecular storage, if files are defined
        #
        if not self.options["input"]["values"]["fname"] is None:
            self.options["input"]["values"]["queue_err"] = self.queue_err
            self.options["input"]["values"]["pipe_comm"] = self._pipe_remote
            self.options["input"]["values"]["handbrake"] = self.handbrake
            self.mol_provider = MoleculeProvider(**self.options["input"]["values"])
        if not self.options["output"]["values"]["fname"] is None:
            self.options["output"]["values"]["queue"] = self.queue_out
            self.options["output"]["values"]["comm_pipe"] = self._pipe_remote
            self.options["output"]["values"]["handbrake"] = self.handbrake
            self.options["output"]["values"]["workers_count"] = self.max_proc
            self.mol_writer = MoleculeStorage(**self.options["output"]["values"])
            self.mol_writer.start()
            # self._registered_workers.append(self.mol_writer)
        ###########################
        # geometry
        #
        if self.options["geometry"]["active"]:
            # define acceptance threshold for geometry optimization
            if self.options["geometry"]["values"]["strict"]:
                self._success_cutoff = 0
            else:
                self._success_cutoff = -1
            # queue in is the source of molecules to process
            self.queue_in = multiprocessing.JoinableQueue(self.max_proc)
            self._target_queue = self.queue_in
            self.options["geometry"]["values"]["queue_in"] = self.queue_in
            self.options["geometry"]["values"]["queue_out"] = self.queue_out
            self.options["geometry"]["values"]["queue_err"] = self.queue_err
            self.options["geometry"]["values"]["handbrake"] = self.handbrake
            self.options["geometry"]["values"]["nice_level"] = nice
            self.options["geometry"]["values"]["max_proc"] = self.max_proc
            self.geometry_optimize = ParallelGeometryGenerator(
                **self.options["geometry"]["values"]
            )
        else:
            self.geometry_optimize = None
            self._target_queue = self.queue_out

        ###########################
        # isomeric transformations
        #
        if self.options["isomers"]["active"]:
            self.isomer = MoleculeIsomers(**self.options["isomers"]["values"])
        else:
            self.isomer = None

        # populate the list of workers to wait for completion
        # the order of workers is sorted by the order in which
        # they're expected to complete their job.
        if not self.geometry_optimize is None:
            self._registered_workers.append(
                ("geometry optimization", self.geometry_optimize)
            )
        if self.mol_writer is not None:
            self._registered_workers.append(("results writer", self.mol_writer))
        if self.mol_issues is not None:
            self._registered_workers.append(("problematic writer", self.mol_issues))

    def _check_still_alive(self):
        """check that all pending workers have terminated their job"""
        print(
            "REGISTER WORKERS", len(self._registered_workers), self._registered_workers
        )
        alive = False
        for worker in self._registered_workers:
            print(worker.is_alive(), "WORKER>", worker)
            if worker.is_alive():
                alive = True
        if alive:
            return True
        return False

    @classmethod
    def get_defaults(cls, terse=False):
        """return the defaults of the Scrubber Core, except for the 'ignore' group"""
        defaults = cls.__default_options.copy()
        if terse:
            for group, opt in defaults.items():
                defaults[group] = {k: v for k, v in opt.items() if not k == "ignore"}
        return defaults

    def _conditional_activation_isomers_geometry(self, options: dict) -> None:
        """check for the isomers/geometry and flag
        the class as active if at least one of the enumerations is requested"""
        # isomers
        options["isomers"]["active"] = any(
            options["isomers"]["values"][x]
            for x in ["stereo_enum", "proto_enum", "tauto_enum"]
        )
        # geometry operations
        options["geometry"]["active"] = any(
           options["geometry"]["values"][x] for x in ["gen3d", "fix_ring_corners"]
        )
        return options

    def process(self, mol):
        """process a molecule and return one or more valid molecules"""
        if self.max_proc > 1:
            self._initialize_mp_pipeline()
            # parallel processing pripeline
            if not self.isomer is None:
                isomer_report = self.isomer.process(mol)
                print("ISOMER REPORT", isomer_report)
                mol_pool = self.isomer.mol_pool
            else:
                mol_pool = [mol]
            for mol_raw in mol_pool:
                self._target_queue.put(PropertyMol(mol_raw), block=True)
            workers_count = self.max_proc
            self._send_poison_pills()
            while True:
                report = self.queue_out.get()
                if report is None:
                    workers_count -= 1
                    if workers_count == 0:
                        return
                else:
                    print("[ DEBUG> parallel queue packet received : ", report, "]")
                    if report["accepted"] >= self._success_cutoff:
                        yield report["mol"]
                    else:
                        self.queue_err.put(("geom", report["mol"]), block=True)
        else:
            # serial processing pripeline
            self._initialize_pipeline()
            if not self.isomer is None:
                isomer_report = self.isomer.process(mol)
                mol_pool = self.isomer.mol_pool
            else:
                mol_pool = [mol]
            for mol_raw in mol_pool:
                report = self.geometry_optimize.process(mol_raw)
                print("[ DEBUG> serial queue packet received : ", report, "]")
                if report["accepted"] >= self._success_cutoff:
                    yield report["mol"]
                else:
                    print("IDISCCARDED!!!!", report)
                    self.queue_err.put(("geom", report["mol"]), block=True)
            self.queue_err.put(None, block=True)

    def _initialize_pipeline(self):
        """initialize the serial processing pipeline"""
        print("[ initializing serial pipeline ]")
        # isomer enumerator
        self.isomer = MoleculeIsomers(**self.options["isomers"]["values"])
        # geometry optimization
        g_opts = self.options["geometry"]["values"]
        self.geometry_optimize = GeometryGenerator(
            add_h=g_opts["add_h"],
            force_trans_amide=g_opts["force_trans_amide"],
            force_field=g_opts["force_field"],
            max_iterations=g_opts["max_iterations"],
            auto_iter_cycles=g_opts["auto_iter_cycles"],
            gen3d=g_opts["gen3d"],
            gen3d_max_attempts=g_opts["gen3d_max_attempts"],
            fix_ring_corners=g_opts["fix_ring_corners"],
            preserve_mol_properties=g_opts["preserve_mol_properties"],
        )
        # define acceptance threshold for geometry optimization
        if g_opts["strict"]:
            self._success_cutoff = 0
        else:
            self._success_cutoff = -1
        # error queue
        self.queue_err = multiprocessing.Queue(maxsize=self.max_proc * 3)

    def get_problematic(self):
        """function to retrieve problematic molecules generated using the
        single-molecule process()"""
        if self.queue_err is None:
            return []
        errors = []
        while True:
            packet = self.queue_err.get()
            if packet is None:
                break
            else:
                errors.append(packet)
        return errors

    def process_file(self, quiet=False):
        """function where all file processing happens
        ideas:
              isomers? ( 1 core )
            - out queue: writer?
        """
        if self.max_proc > 1:
            self._initialize_mp_pipeline()
        else:
            self._initialize_pipeline()
        counter = 0
        skipped = 0
        # time =a TIME GOES HERE
        progress = "⣾⣽⣻⢿⡿⣟⣯⣷"
        progress = ["▁", "▃", "▄", "▅", "▆", "▇", "█", "▇", "▆", "▅", "▄", "▃"]
        # progress = "▉▊▋▌▍▎▏▎▍▌▋▊▉"
        # progress = [ "◜", "◝", "◞", "◟" ]
        t_start = time.time()
        mol_sec = -1
        mol_step = 10
        try:
            for counter, mol in self.mol_provider:
                if counter % mol_step == 0:
                    mol_sec = mol_step / (time.time() - t_start)
                    t_start = time.time()
                if mol_sec != -1:
                    timing = "( %2.3f mol/sec.)" % mol_sec
                else:
                    timing = ""
                if mol is None:
                    skipped += 1
                    continue
                if not quiet:
                    bars = "|".join(
                        random.choice(progress) for x in range(self.max_proc)
                    )
                    print(
                        "\r[ %s processing input mol. %d ]  %s"
                        % (bars, counter, timing),
                        end="",
                    )
                if not self.isomer is None:
                    self.isomer.process(mol)
                    mol_pool = self.isomer.mol_pool
                else:
                    mol_pool = [mol]
                for mol_raw in mol_pool:
                    mol_raw.SetProp("Scrubber_was_here", "Yes!")
                    self._target_queue.put(PropertyMol(mol_raw), block=True)

            print(
                "\r ----- COMPLETED -----                                                         "
            )
            self.counter_data["input"] = counter
            self._send_poison_pills()
        except KeyboardInterrupt:
            print(
                "\n\n *** Keyboard interruption captured. Attempting to exit gracefully... ***\n\n"
            )
            self.handbrake.set()
            self.counter_data["input"] = counter
        self.wait_pending(skipped, quiet)
        if not quiet:
            self.print_summary()

    def wait_pending(self, skipped, quiet=False):
        """ wait for all pending operations: join registered workers;
        retrieve information from the comm pipes; populate summary inforation """
        for name, worker in self._registered_workers:
            if not quiet:
                t_start = time.time()
                print("[ waiting for pending task: %s ..." % name, end="")
            worker.join()
            if not quiet:
                print(" DONE (%2.3f s) ]" % (time.time() - t_start))
        if not quiet:
            print("[ cleaning pipes", end="")
            t_start = time.time()
        while self._pipe_listener.poll():
            if not quiet:
                print(".", end="")
            packet = self._pipe_listener.recv()
            label, value = packet.split(":")
            self.counter_data[label] = int(value)
        if not quiet:
            print(" DONE (%2.3f s) ]" % (time.time() - t_start))
        self._pipe_listener.close()
        for k, v in self.counter_data.items():
            if v is None:
                self.counter_data[k] = "n/a"
        try:
            pc = self.counter_data["writer"] / (self.counter_data["input"] - skipped)
            if pc > 1.0:
                net_result = "+%2.2f%%" % ((pc - 1.0) * 100)
            else:
                net_result = "-%2.2f%%" % ((1.0 - pc) * 100)
        except ZeroDivisionError:
            net_result = "n/a"
        self.counter_data["net_result"] = net_result

    def print_summary(self):
        """print the summary of the calculation"""
        print("\n==============================================")
        print("Summary")
        print("----------------------------------------------")
        print(" Input parsed   : %s " % (self.counter_data["input"]))
        print(" Input errors   : %s" % (self.counter_data["err_input"]))
        print(" Process errors : %s" % (self.counter_data["err_process"]))
        print(
            " Written        : %s  | %s"
            % (self.counter_data["writer"], self.counter_data["net_result"])
        )
        print("==============================================")

    def _send_poison_pills(self):
        """send poison pills to fill the output queue"""
        # print("FILLING TARGET QUEUE", self.max_proc)
        for _ in range(self.max_proc):
            # send poison pills
            # print("sebd poison pills...")
            self._target_queue.put(None, block=True)
        if not self.queue_err is None:
            # print("QUEUE ERROR STATUS IS", self.queue_err.qsize())
            # print("send pill to queue")
            # print("FILLING ERROR QUEUE")
            self.queue_err.put(None, timeout=60, block=False)

    def __recursive_dict_match(self, source: dict, target: dict) -> None:
        """function to use an arbitrarily nested dict to change values in an
        arbitrarily nested dict

        Warning: it should be used only to assign values as standard Python
        types, no guarantee it is going to work with other types (i.e.,
        classes)
        """
        curr_source = source
        curr_target = target
        for k, v in curr_source.items():
            if not k in curr_target:
                print("Warning: unrecognized keyword:", k, curr_target)
                continue
            if isinstance(v, dict):
                self.__recursive_dict_match(v, curr_target[k])
            else:
                curr_target[k] = v

class Scrub:

    def __init__(
        self,
        ph_low=7.4,
        ph_high=None,
        pka_fname=None,
        tauto_fname=None,
        no_acidbase=False,
        no_tautomers=False,
        no_ringfix=False,
        no_gen3d=False,
    ):
        self.acid_base_conjugator = AcidBaseConjugator()
        self.tautomerizer = Tautomerizer()
        self.ph_low = ph_low
        if ph_high is None:
            ph_high = ph_low
        self.ph_high = ph_high
        self.do_acidbase = not no_acidbase
        self.do_tautomers = not no_tautomers
        self.no_ringfix = no_ringfix # not avoid double negative to pass directly to gen3d
        self.do_gen3d = not no_gen3d

    def __call__(self, input_mol):
        mol = Chem.RemoveHs(input_mol) 
        pool = [input_mol]

        if self.do_acidbase:
            molset = UniqueMoleculeContainer()
            for mol in pool:
                for mol_out in self.acid_base_conjugator(mol, self.ph_low, self.ph_high):
                    molset.add(mol_out)
            pool = list(molset)

        if self.do_tautomers:
            molset = UniqueMoleculeContainer()
            for mol in pool:
                for mol_out in self.tautomerizer(mol): 
                    molset.add(mol_out)
            pool = list(molset)
        print("Nr mols to gen3d", len(pool))
        output_mol_list = []
        if self.do_gen3d:
            for mol in pool:
                mol_out = gen3d(mol, no_ringfix=self.no_ringfix)
                output_mol_list.append(mol_out)

        return output_mol_list

etkdg_config = rdDistGeom.ETKDGv3()

def gen3d(mol, no_ringfix=False):
    mol.RemoveAllConformers()
    mol = Chem.AddHs(mol)
    rdDistGeom.EmbedMolecule(mol, etkdg_config)
    etkdg_coords = mol.GetConformer().GetPositions()
    mol.RemoveAllConformers()
    if no_ringfix:
        coords_list = [etkdg_coords]
    else:
        coords_list = fix_rings(mol, etkdg_coords)
        print("nr coords", len(coords_list))
    for coords in coords_list:
        c = Chem.Conformer(mol.GetNumAtoms())
        for i, (x, y, z) in enumerate(coords):
            c.SetAtomPosition(i, Point3D(x, y, z))
        mol.AddConformer(c, assignId=True)
    rdForceFieldHelpers.UFFOptimizeMoleculeConfs(mol)
    return mol



if __name__ == "__main__":
    # import sys
    import json
    config = ScrubberCore.get_defaults()
    print("CONFIG!")
