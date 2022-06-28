import multiprocessing
from time import sleep

from rdkit.Chem.PropertyMol import PropertyMol

# storage import . import storage
from .storage import (
    MoleculeProvider,
    MoleculeStorage,
    SDFMolSupplierWrapper,
    SMIMolSupplierWrapper,
)
from .geom.geometry import ParallelGeometryGenerator  # , GeometryGenerator
from .filters import MoleculeFilter
from .transform.isomer import MoleculeIsomers
from .transform.reaction import Reactor


"""
This file contains the core scrubber object


# TODO add a dump for json format
    - dump all defaults
    - dump active only (stuff that's been modified)

"""


class ScrubberCore(object):
    """read input files, manipulate them to perform the following operations:

        - read input files
            - convert from OB?
        - apply filters
            - generic SMARTS wanted/unwanted
            - PAINS
            - mw, logP, others...? (which can be used to apply a series of controversial filters, like Lipinski)
        - elaborate structures ( * EXTERNAL * )
            - tautomers, protomers, chiral centers ( chem modofications/warheads? )
        - pass molecules to 3D processors ( * EXTERNAL * )
        - write the results

    for mol in input_lib:
        for mol_n in transform(mol):
            if not filter(mol_n):
                continue
            gen_3d(mol_n)
            write(mol_n)
    """

    __default_options = {
        "input": {
            "values": MoleculeProvider.get_defaults(),
            "ignore": ["pipe"],
        },
        "output": {
            "values": MoleculeStorage.get_defaults(),
            "ignore": ["queue", "comm_pipe", "workers_count"],
        },
        "filter_pre": {
            "active": False,
            "values": MoleculeFilter.get_defaults(),
            "ignore": [],
        },
        "filter_post": {
            "active": False,
            "values": MoleculeFilter.get_defaults(),
            "ignore": [],
        },
        "reaction": {
            "active": False,
            "values": Reactor.get_defaults(),
            "desc": "Options for chemical modifications to be performed on the input",
            "ignore": [],
        },
        "isomers": {
            "active": False,
            "values": MoleculeIsomers.get_defaults(),
            "ignore": [],
        },
        "geometry": {
            "active": True,
            "values": ParallelGeometryGenerator.get_defaults(),
            "ignore": ["queue_in", "queue_out", "nice_level"],
            # "values": GeometryGenerator.get_default_options(),
        },
        "general": {
            "values": {
                "max_proc": multiprocessing.cpu_count(),
                "nice_level": None,
            },
            "ignore": [],
        },
        # TODO STUFF BELOW HERERE
        "log": {
            "values": {
                "problematic_logfile": None,
                "other": None,
            },
            "ignore": [],
        },
        # -- CALC PROPERTIES
        # calc_properties : {}
        # -- REMOVE SALTS
        #  "clean" : {}
    }

    def __init__(self, options: dict = None):
        """this is where the classes doing all operations will do"""
        self.options = options
        # self._parse_init(self.options)
        ###########################
        # initialize molecular provider and molecular storage
        #
        # print("MOL PROV INIT", self.options["input"]["values"])
        self.mol_provider = MoleculeProvider(**self.options["input"]["values"])

        # queue out is where processed molecules are sent
        # if more than one processor is used, then use 4 as many for writing (low load)
        self.max_proc = self.options["general"]["values"]["max_proc"]
        nice = self.options["general"]["values"]["nice_level"]
        # if self.max_proc > 1:
        #     self.max_proc *= 4
        self._pipe_listener, self._pipe_remote = multiprocessing.Pipe(False)
        self.queue_out = multiprocessing.Queue(maxsize=-1)  # self.max_proc)
        self.options["output"]["values"]["queue"] = self.queue_out
        self.options["output"]["values"]["workers_count"] = self.max_proc
        self.options["output"]["values"]["comm_pipe"] = self._pipe_remote
        self.mol_writer = MoleculeStorage(**self.options["output"]["values"])
        # self.mol_writer._counter+=1000
        self.mol_writer.start()
        # self.mol_writer.join()

        ###########################
        # geometry
        #
        if self.options["geometry"]["active"]:
            # queue in is the source of molecules to process
            # nice = self.options["general"]["nice_level"]
            self.queue_in = multiprocessing.JoinableQueue(self.max_proc)
            self._target_queue = self.queue_in
            self.options["geometry"]["values"]["queue_in"] = self.queue_in
            self.options["geometry"]["values"]["queue_out"] = self.queue_out
            self.options["geometry"]["values"]["nice_level"] = nice
            self.geometry_optimize = ParallelGeometryGenerator(
                **self.options["geometry"]["values"]
            )
            # TODO move the multiproc queues here?
        else:
            self.geometry_optimize = None
            self._target_queue = self.queue_out

        ###########################
        # filters
        #
        if self.options["filter_pre"]["active"]:
            self.filter_pre = MoleculeFilter(**self.options["filter_pre"]["values"])
        else:
            self.filter_pre = None
        if self.options["filter_post"]["active"]:
            self.filter_post = MoleculeFilter(**self.options["filter_post"]["values"])
        else:
            self.filter_post = None

        ###########################
        # reactions
        #
        if self.options["reaction"]["active"]:
            self.reactor = Reactor(**self.options["reaction"]["values"])
        else:
            self.reactor = None

        ###########################
        # isomeric transformations
        #
        if self.options["isomers"]["active"]:
            print("XXX", self.options["isomers"]["values"])
            self.isomer = MoleculeIsomers(**self.options["isomers"]["values"])
        else:
            self.isomer = None

    @classmethod
    def get_defaults(cls):
        return cls.__default_options.copy()

    def _parse_init(self, options: dict) -> None:
        """extend the parent parse_init to add a check for the isomers and flag
        the class as active if at least one of the enumerations is requested"""
        # isomers
        self.options["isomers"]["active"] = any(
            self.options["isomers"]["values"][x]
            for x in ["stereoisomer_enum", "protomer_enum", "tautomer_enum"]
        )
        # geometry operations
        self.options["geometry"]["active"] = any(
            self.options["geometry"]["general"][x]
            for x in ["gen3d", "fix_ring_corners"]
        )

    def process(self):
        """function where all processing happens
        ideas:

            - pre-processing queue to parallelize: filtering, reactions,
              isomers? ( 1 core )
            - post-processing queue: 3D geometry, rings? ( max_cpu-1 cores )
            - out queue: writer?

        TODO: explore the possibility of writing a file with all the generated
        proto/stereo/tautomers/reactions??
        """
        # this is required here to prevent exceptions in case of early failures
        counter = 0
        try:
            for counter, mol_org in self.mol_provider:
                if mol_org is None:
                    continue
                print("Processing mol:", counter)
                # print("MOL COUNTER (pre-PILL)", self.mol_writer._RANDOM)
                # prefilter
                if not self.filter_pre is None:
                    if not self.filter_pre.filter(mol_org):
                        continue
                # chemical modifications/transformations
                if not self.reactor is None:
                    mol_pool = self.reactor.react(mol_org)
                else:
                    mol_pool = [mol_org]
                # isomers
                for mol_out in mol_pool:
                    if not self.isomer is None:
                        isomer_report = self.isomer.process(mol_out)
                        mol_isomers = self.isomer.mol_pool
                    else:
                        mol_isomers = [mol_out]
                    for mol_raw in mol_isomers:
                        if not self.filter_post is None:
                            if not self.filter_post.filter(mol_raw):
                                continue
                        mol_raw.SetProp("Scrubber_was_here", "Yes!")
                        self._target_queue.put(PropertyMol(mol_raw))
            self._send_poison_pills()
        except KeyboardInterrupt:
            print("[ Keyboard interruption requested ]")
            if not self.geometry_optimize is None:
                self.geometry_optimize.terminate()
            else:
                self._send_poison_pills()
        # except Exception as err:
        #     print("GENERIC ERROR",err )
        #     if not self.geometry_optimize is None:
        #         # self.geometry_optimize.halt_workers()
        #         self.geometry_optimize.terminate()
        #         print("\n\n**** WE'RE DONE ****\n\n")
        #     else:
        #         self._send_poison_pills()

        # except ValueError:
        #     if not self.geometry_optimize is None:
        #         # self.geometry_optimize.halt_workers()
        #         self.geometry_optimize.terminate()
        #         print("\n\n**** WE'RE DONE ****\n\n")
        #     else:
        #         self._send_poison_pills()
        print("[ saving all molecules... ]")
        if not self.mol_writer.is_alive():
            print("WRITER DEAD... trying to get the number...")
            if self._pipe_listener.poll():
                written_count = self._pipe_listener.recv()
        else:
            while self.mol_writer.is_alive():
                sleep(0.1)
                if self._pipe_listener.poll():
                    written_count = self._pipe_listener.recv()
        print("FINAL: %d input molecules processed" % counter)
        print(" total molecules written:: %d" % written_count)
        # x = self._target_queue.get()
        # print("XXX", x)
        # self._target_queue.

    def _send_poison_pills(self):
        """send poison pills to fill the output queue"""
        for i in range(self.max_proc):
            # send poison pills
            self._target_queue.put(None)

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


if __name__ == "__main__":
    import sys
    from pprint import pprint as pp

    config = ScrubberCore.get_defaults()
    print("CONFIG!")
    pp(config)
