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
from rdkit.Chem import AllChem, rdMolAlign
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
from .espaloma_minim import EspalomaMinimizer


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
        """wait for all pending operations: join registered workers;
        retrieve information from the comm pipes; populate summary inforation"""
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

################################
# code from here up is not used. 
##############################3
class Scrub:

    def __init__(
        self,
        ph_low=7.4,
        ph_high=None,
        pka_fname=None,
        tauto_fname=None,
        skip_acidbase=False,
        skip_tautomers=False,
        skip_ringfix=False,
        skip_gen3d=False,
        template=None,
        template_smarts=None,
        do_gen2d=False,
        max_ff_iter=200,
        numconfs=1,
        etkdg_rng_seed=None,
        ff="mmff94s",
    ):
        self.acid_base_conjugator = AcidBaseConjugator.from_default_data_files()
        self.tautomerizer = Tautomerizer.from_default_data_files()
        self.ph_low = ph_low
        if ph_high is None:
            ph_high = ph_low
        self.ph_high = ph_high
        self.do_acidbase = not skip_acidbase
        self.do_tautomers = not skip_tautomers
        self.skip_ringfix = (
            skip_ringfix  # not avoiding negative to pass directly to gen3d
        )
        self.do_gen3d = not skip_gen3d
        self.template = template
        self.template_smarts = template_smarts
        self.do_gen2d = do_gen2d
        self.max_ff_iter = max_ff_iter
        self.numconfs = numconfs
        self.etkdg_rng_seed = (
            etkdg_rng_seed if etkdg_rng_seed else random.randint(0, 1000000)
        )
        self.ff = ff

        if ff == "espaloma":
            self.espaloma = EspalomaMinimizer()
        else:
            self.espaloma = None

    def __call__(self, input_mol):

        mol = Chem.RemoveHs(input_mol)
        pool = [input_mol]

        if self.do_acidbase:
            molset = UniqueMoleculeContainer()
            for mol in pool:
                for mol_out in self.acid_base_conjugator(
                    mol, self.ph_low, self.ph_high
                ):
                    molset.add(mol_out)
            pool = list(molset)

        if self.do_tautomers:
            molset = UniqueMoleculeContainer()
            for mol in pool:
                for mol_out in self.tautomerizer(mol):
                    molset.add(mol_out)
            pool = list(molset)

        if self.do_gen3d:
            output_mol_list = []
            for mol in pool:
                mol_out = gen3d(
                    mol,
                    skip_ringfix=self.skip_ringfix,
                    max_ff_iter=self.max_ff_iter,
                    etkdg_rng_seed=self.etkdg_rng_seed,
                    numconfs=self.numconfs,
                    ff=self.ff,
                    espaloma=self.espaloma,
                    template=self.template,
                    template_smarts=self.template_smarts,
                )
                output_mol_list.append(mol_out)
        elif self.do_gen2d:  # useful to write SD files
            output_mol_list = []
            for mol in pool:
                AllChem.Compute2DCoords(mol)
                output_mol_list.append(mol)
        else:
            output_mol_list = pool

        return output_mol_list


def constrained_embeding(
    query_mol,
    core_mol,
    template_smarts: str = None,
    numconfs: int = 1,
    ff: str = "mmff94s",
    ps=None,
):
    """Generate an embedding of a query molecule where part of the molecule
    is constrained to have particular coordinates derived from a core.
    Alternatively, a SMARTs pattern can be provided to match specific atoms from both molecules
    """

    force_constant = 1000
    confId = -1

    if ff == "uff":
        getForceField = AllChem.UFFGetMoleculeForceField
    elif ff == "mmff94":
        getForceField = lambda x: AllChem.MMFFGetMoleculeForceField(
            x, AllChem.MMFFGetMoleculeProperties(x), confId=confId
        )
    # This is just is a minimization with restraints to force the querry mol to match the template. 
    # If you chose espaloma as ff it will still minimize it at the end with that forcefield.
    elif ff == "mmff94s" or ff == "espaloma":
        getForceField = lambda x: AllChem.MMFFGetMoleculeForceField(
            x,
            AllChem.MMFFGetMoleculeProperties(x, mmffVariant="MMFF94s"),
            confId=confId,
        )

    if template_smarts == None:
        query_match = query_mol.GetSubstructMatches(core_mol)
        if not query_match:
            raise ValueError("molecule doesn't match the core")
        elif len(query_match) > 1:
            raise RuntimeError("Expected one match but multiple matches were found.")
        else:
            query_match = query_match[0]

        algMap = [(j, i) for i, j in enumerate(query_match)]

    else:
        query_match = query_mol.GetSubstructMatches(template_smarts)

        if not query_match:
            raise ValueError("SMARTs doesn't match the molecule")
        elif len(query_match) > 1:
            raise RuntimeError("Expected one match but multiple matches were found.")
        else:
            query_match = query_match[0]

        core_match = core_mol.GetSubstructMatches(template_smarts)

        if not core_match:
            raise ValueError("SMARTs doesn't match the template")
        elif len(core_match) > 1:
            raise RuntimeError("Expected one match but multiple matches were found.")
        else:
            core_match = core_match[0]

        algMap = [(j, i) for j, i in zip(query_match, core_match)]

    cids = rdDistGeom.EmbedMultipleConfs(query_mol, numconfs, ps)

    # rotate the embedded conformation onto the core:
    rms = rdMolAlign.AlignMol(query_mol, core_mol, atomMap=algMap)
    ff = getForceField(query_mol, confId=confId)
    conf = core_mol.GetConformer()
    for atom in algMap:
        p = conf.GetAtomPosition(atom[1])
        pIdx = ff.AddExtraPoint(p.x, p.y, p.z, fixed=True) - 1
        ff.AddDistanceConstraint(pIdx, atom[0], 0, 0, force_constant)
    ff.Initialize()
    n = 4
    more = ff.Minimize(energyTol=1e-4, forceTol=1e-3)
    while more and n:
        more = ff.Minimize(energyTol=1e-4, forceTol=1e-3)
        n -= 1

    # realign
    rms = rdMolAlign.AlignMol(query_mol, core_mol, atomMap=algMap)

    query_mol.SetProp("EmbedRMS", str(rms))

    return query_mol, cids


def _ConfToMol(mol, conf_id):
    conf = mol.GetConformer(conf_id)
    new_mol = Chem.Mol(mol)
    new_mol.RemoveAllConformers()
    new_mol.AddConformer(Chem.Conformer(conf), assignId=True)
    return new_mol


def translate_failures(failure_counts):
    """A function to catch embeding failures and translate the codes to meaninful error messages.
    An error will all failures is raised, catched by scrub_and_catch_errors function, and printed to the SDF property of the failed SDF file.

    Here is what the individual failures mean:

    INITIAL_COORDS: generation of the initial coordinates from the random distance matrix (default) or from a set of random coordinates (when using random coordinate embedding) failed.
    FIRST_MINIMIZATION: the initial optimization of the atom positions using the distance-geometry force field failed to produce a low-enough energy conformer. The check here has thresholds for both average energy per atom and the individual atom energies. I’m not providing the threshold values here since the energies from the distance-geometry force field are not physically meaningful - the threshold values are not interpretable.
    CHECK_TETRAHEDRAL_CENTERS: at least one tetrahedral C and N centers either has a volume around it which is too small or is outside the volume defined by its neighbors
    MINIMIZE_FOURTH_DIMENSION: the minmization to force the values of the fourth-dimensional component of each atom position failed
    ETK_MINIMIZATION: after the minimization with the ET and/or K terms, at least one atom which should have been planar was not
    FINAL_CHIRAL_BOUNDS: the neighborhood of an atom with specified chirality was too distorted (it violated distance constraints)
    FINAL_CENTER_IN_VOLUME: an atom with specified chirality was outside of the volume defined by its neighbors
    LINEAR_DOUBLE_BOND: one of the end atoms of a double bond had a linear geometry
    BAD_DOUBLE_BOND_STEREO: the stereochemistry of a double bond with specified stereochemistry was wrong in the generated conformer

    """

    if sum(failure_counts) != 0:
        failure_msgs = {}
        for i, k in enumerate(rdDistGeom.EmbedFailureCauses.names):
            if failure_counts[i] != 0:
                failure_msgs[k] = failure_counts[i]

        raise RuntimeError(failure_msgs)
    else:
        return None


def gen3d(
    mol,
    skip_ringfix: bool = False,
    max_ff_iter: int = 200,
    etkdg_rng_seed: int = 42,
    numconfs: int = 1,
    ff: str = "mmff94s",
    espaloma=None,
    template=None,
    template_smarts=None,
):
    mol.RemoveAllConformers()
    mol = Chem.AddHs(mol)

    # Set up the ETKDG parameters
    ps = rdDistGeom.ETKDGv3()
    ps.randomSeed = etkdg_rng_seed
    ps.trackFailures = True
    ps.enforceChirality = True
    ps.useSmallRingTorsions = True
    ps.useMacrocycleTorsions = True
    ps.clearConfs = True

    if template is not None:
        mol, cids = constrained_embeding(
            query_mol=mol,
            core_mol=template,
            template_smarts=template_smarts,
            numconfs=numconfs,
            ff=ff,
            ps=ps,
        )

    else:
        cids = rdDistGeom.EmbedMultipleConfs(mol, numconfs, ps)

    if len(cids) == 0:
        translate_failures(ps.GetFailureCounts())

    etkdg_coords = [c.GetPositions() for c in mol.GetConformers()]

    mol.RemoveAllConformers()  # to be added back after ringfix

    if skip_ringfix:
        coords_list = etkdg_coords
    else:
        coords_list = []
        [coords_list.extend(fix_rings(mol, c)) for c in etkdg_coords]

    for coords in coords_list:
        c = Chem.Conformer(mol.GetNumAtoms())
        for i, (x, y, z) in enumerate(coords):
            c.SetAtomPosition(i, Point3D(x, y, z))
        mol.AddConformer(c, assignId=True)

    if ff not in ["uff", "mmff94", "mmff94s", "espaloma"]:
        raise RuntimeError(
            f"ff is {ff} but must be 'uff', 'mmff94', 'mmff94s', or 'espaloma'"
        )

    if ff == "espaloma":
        if espaloma is None:
            raise ValueError("espaloma minimizer needs to be passed")
        mol, energies = espaloma.minim_espaloma(mol)
    else:
        optimize_func = {
            "uff": rdForceFieldHelpers.UFFOptimizeMoleculeConfs,
            "mmff94": rdForceFieldHelpers.MMFFOptimizeMoleculeConfs,
            "mmff94s": lambda mol, maxIters: rdForceFieldHelpers.MMFFOptimizeMoleculeConfs(
                mol, maxIters=maxIters, mmffVariant="mmff94s"
            ),
        }[ff]
        _energies = optimize_func(mol, maxIters=max_ff_iter)
        energies = [e[1] for e in _energies]

    best_energy_index = min(zip(cids, energies), key=lambda x: x[1])[0]
    final_mol = _ConfToMol(mol, best_energy_index)

    return final_mol


if __name__ == "__main__":
    config = ScrubberCore.get_defaults()
