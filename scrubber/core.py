import multiprocessing
from time import sleep
import sys


from rdkit import Chem
from rdkit.Chem.PropertyMol import PropertyMol

from .storage import (
    MoleculeProvider,
    MoleculeStorage,
    MoleculeIssueStorage,
)

from .geom.geometry import ParallelGeometryGenerator
from .geom.geometry import GeometryGenerator
from .filters import MoleculeFilter
from .transform.isomer import MoleculeIsomers
from .transform.reaction import Reactor


"""
This file contains the core scrubber object


# TODO add a dump for json format
    - dump all defaults
    - dump active only (stuff that's been modified)

INSPIRATION: https://xkcd.com/1343/

"""


class MolPrep:
    
    def __init__(self, isomerizer, geom_opt):

        assert(isomerizer is not None) # TODO
        assert(geom_opt is not None) # TODO
        self.isomerizer = isomerizer
        self.geom_opt = geom_opt

    def __call__(self, input_mol):
        mol_list = []
        isomer_report = self.isomerizer.process(input_mol)
        for mol in self.isomerizer.mol_pool: # TODO mol_pool gonna fail with concurrency
            report = self.geom_opt.process(mol)
            #optmol = report["mol"]
            mol_list.append(report)
        return mol_list


class MolPrepWorker(multiprocessing.Process):

    def __init__(self, molprep, queue_in, queue_out):
        multiprocessing.Process.__init__(self)
        self.molprep = molprep
        self.queue_in = queue_in
        self.queue_out = queue_out

    def run(self):
        while True:
            try:
                input_mol = self.queue_in.get()
                if input_mol is None:
                    self.queue_out.put(None)
                    break
                mol_list = self.molprep(input_mol)
                for mol in mol_list:
                    self.queue_out.put(mol)
            except KeyboardInterrupt:
                print("Keyboard Interruption at MolPrepWorker")
                break




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
            "ignore": ["use_pipe", "queue_err", "pipe_comm"],
        },
        "output": {
            "values": MoleculeStorage.get_defaults(),
            "ignore": ["queue", "comm_pipe", "num_poison_pills_needed"],
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
            "ignore": ["queue_in", "queue_out", "queue_err", "nice_level", "max_proc"],
        },
        "general": {
            "values": {
                "max_proc": multiprocessing.cpu_count(),
                "nice_level": None,
            },
            "ignore": [],
        },
        # TODO STUFF BELOW HERERE
        "errors": {
            "values": MoleculeIssueStorage.get_defaults(),
            "ignore": ["queue", "disable_name_sanitize", "comm_pipe"],
        },
        # -- CALC PROPERTIES
        # calc_properties : {}
        # -- REMOVE SALTS
        #  "clean" : {}
    }

    def __init__(self, options: dict = None):
        """this is where the classes doing all operations will do"""
        self.options = options
        self._parse_init(self.options)
        self._registered_workers = []
        # initialize pipes for multiprocessing communication
        self._pipe_listener, self._pipe_remote = multiprocessing.Pipe(False)
        ###########################
        # initialize error logging, if requested
        #
        if (self.options["errors"]["values"]["log_from_input"] is not None) or (
            self.options["errors"]["values"]["log_from_process"] is not None
        ):
            self.queue_err = multiprocessing.Queue(maxsize=-1)
            self.options["errors"]["values"]["queue"] = self.queue_err
            self.options["errors"]["values"]["comm_pipe"] = self._pipe_remote
            self.mol_issues = MoleculeIssueStorage(**self.options["errors"]["values"])
            self.mol_issues.start()
            self._registered_workers.append(self.mol_issues)
        else:
            self.queue_err = None
            self.mol_issues = None
        ###########################
        # initialize molecular provider and molecular storage
        #
        self.options["input"]["values"]["queue_err"] = self.queue_err
        self.options["input"]["values"]["pipe_comm"] = self._pipe_remote
        self.mol_provider = MoleculeProvider(**self.options["input"]["values"])
        # queue out is where processed molecules are sent
        # if more than one processor is used, then use 4 as many for writing (low load)
        self.max_proc = self.options["general"]["values"]["max_proc"]
        nice = self.options["general"]["values"]["nice_level"]
        self.queue_out = multiprocessing.Queue()
        self.options["output"]["values"]["queue"] = self.queue_out
        # self.options["output"]["values"]["queue_err"] = self.queue_err
        self.options["output"]["values"]["num_poison_pills_needed"] = self.max_proc-1 # one for writer
        self.options["output"]["values"]["comm_pipe"] = self._pipe_remote
        self.mol_writer = MoleculeStorage(**self.options["output"]["values"])
        self.mol_writer.start()
        self._registered_workers.append(self.mol_writer)
        ###########################
        # geometry
        #
        if self.options["geometry"]["active"]:
            g_opts = self.options["geometry"]["values"]
            geometry_optimize = GeometryGenerator(
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
            # TODO move the multiproc queues here?
        else:
            geometry_optimize = None
            #self._target_queue = self.queue_out

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
            isomer = MoleculeIsomers(**self.options["isomers"]["values"])
        else:
            isomer = None

        print(self.max_proc)
        self.queue_in = multiprocessing.JoinableQueue(self.max_proc)
        self.molprep = MolPrep(isomer, geometry_optimize)

    def _check_still_alive(self):
        """check that all pending workers have terminated their job"""
        for worker in self._registered_workers:
            if worker.is_alive():
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

    def _parse_init(self, options: dict) -> None:
        """extend the parent parse_init to add a check for the isomers and flag
        the class as active if at least one of the enumerations is requested"""
        # isomers
        self.options["isomers"]["active"] = any(
            self.options["isomers"]["values"][x]
            for x in ["stereo_enum", "proto_enum", "tauto_enum"]
        )
        # geometry operations
        self.options["geometry"]["active"] = any(
            self.options["geometry"]["values"][x] for x in ["gen3d", "fix_ring_corners"]
        )

    def process(self):
        """function where all processing happens
        ideas:

            - pre-processing queue to parallelize: filtering, reactions,
              isomers? ( 1 core )
            - out queue: writer?

        TODO: explore the possibility of writing a file with all the generated
        proto/stereo/tautomers/reactions??
        """
        # this is required here to prevent exceptions in case of early failures
        counter = 0
        skipped = 0
        
        for i in range(self.max_proc - 1): # one of the processes is the writer
            worker = MolPrepWorker(self.molprep, self.queue_in, self.queue_out)
            worker.start()
            self._registered_workers.append(worker)

        try:
            for counter, mol_org in self.mol_provider:
                if mol_org is None:
                    skipped += 1
                    continue
                print("\rProcessing mol: % 10d" % counter, end="")
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
                    self.queue_in.put(mol_out)
            self._send_poison_pills()
        except KeyboardInterrupt:
            print("[ Keyboard interruption requested ]")
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

        print("\r[ saving all molecules... ]                           ")
        counter_data = {
            "input": -1,
            "err_input": -1,
            "err_process": -1,
            "writer": -1,
        }
        print("\r[ waiting for the writers", end="")
        # while (self.mol_writer.is_alive()) or (self.mol_issues.is_alive()):
        while self._check_still_alive():
            print(".", end="")
            sys.stdout.flush()
            sleep(0.5)
        print(" ]")
        while self._pipe_listener.poll():
            packet = self._pipe_listener.recv()
            label, value = packet.split(":")
            counter_data[label] = int(value)
        self._pipe_listener.close()
        for k, v in counter_data.items():
            if v is None:
                counter_data[k] = "n/a"
            # else:
            #     counter_data[k] = '%d' % v
        try:
            pc = counter_data["writer"] / (counter_data["input"] - skipped)
            if pc > 1.0:
                net_result = "+%2.2f%%" % ((pc - 1.0) * 100)
            else:
                net_result = "-%2.2f%%" % ((1.0 - pc) * 100)
        except ZeroDivisionError:
            net_result = "n/a"

        print("\n==============================================")
        print("Summary")
        print("----------------------------------------------")
        print(" Input parsed   : %s " % (counter_data["input"]))
        print(" Input errors   : %s" % (counter_data["err_input"]))
        print(" Process errors : %s" % (counter_data["err_process"]))
        print(" Written        : %s  | %s" % (counter_data["writer"], net_result))
        print("==============================================")
        # print("%d input molecules processed" % counter)
        # # calculate percentages
        # print("%d molecules written (%s)" % (written_count, net_result))
        # print("%d molecules skipped (errors)" % skipped)

    def _send_poison_pills(self):
        """send poison pills to fill the output queue"""
        for i in range(self.max_proc):
            # send poison pills
            self.queue_in.put(None)
        if not self.queue_err is None:
            self.queue_err.put(None)

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
    # TODO add an option to create run an instance of this object and run it with a JSON file input?
    # import sys
    import json
    # from pprint import pprint as pp
    config = ScrubberCore.get_defaults()
    print("CONFIG!")
    # pp(config)
