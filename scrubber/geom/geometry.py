import multiprocessing
import os

from rdkit import Chem
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdForceFieldHelpers
from rdkit.Chem import rdMolTransforms
from rdkit.Chem.PropertyMol import PropertyMol

from scrubber import ScrubberClass, copy_mol_properties


class GeometryGenerator(ScrubberClass):
    """Generate and optimize

    manage the generation of 3D coordinates of molecules and optimize them using force fields

    The class is meant to be used by initializing it with a given configuration
    then processing molecules (i.e.: no parameter changes after initialization)

    NOTE: This class should be used for managing molecules formed by a single fragment only

    """

    # TODO add ring corner flipping here?
    # TODO implement dedicated macrocycle minimization?
    FORCE_FIELD_LIST = ["uff", "mmff94", "mmff94s"]

    def __init__(
        self,
        add_h: bool = True,
        force_trans_amide: bool = True,
        force_field: str = "uff",
        max_iterations: int = 200,
        auto_iterations: int = 10,
        gen3d: bool = True,
        gen3d_max_attempts: int = 3,
        fix_ring_corners: bool = False,
        preserve_mol_properties:bool=True,
        _stop_at_defaults: bool = False,
    ):
        """forcefield       :   (uff, MMFF94, MMFF94s)
        max_iterations  :   ( int, 200 )
        auto            :   automatic minimization: how many times the
                            minimization with max_iterations will be repeated
        gen3d           :   generate 3D coordinates (using the ETKDGv3 method)
        gen3d_max_attempts  : how many times distance geometry attempts in case
                            of failure in generating 3D coords
        """
        self.add_h = add_h
        self.force_trans_amide = force_trans_amide
        self.force_field = force_field
        self.max_iterations = max_iterations
        self.auto_iterations = auto_iterations
        self.gen3d = gen3d
        self.gen3d_max_attempts = gen3d_max_attempts
        self.fix_ring_corners = fix_ring_corners
        self.preserve_mol_properties = preserve_mol_properties
        if _stop_at_defaults:
            return
        # this flag is used to control if at any time the minimization has been asked to stop
        self._handbrake = False
        if not self.force_field in self.FORCE_FIELD_LIST:
            msg = "*** ERROR *** invalid force field type |%s|, allowed: %s" % (
                force_field,
                ",".join(self.FORCE_FIELD_LIST),
            )
            raise ValueError(msg)
        self.ff_parms = {"maxIters": self.max_iterations}
        if self.force_field in ["mmff94", "mmff94s"]:
            # print("MMFF", force_field)
            self.ff_parms["mmffVariant"] = self.force_field
            self.ff_optimize = rdForceFieldHelpers.MMFFOptimizeMolecule
        elif self.force_field == "uff":
            self.ff_optimize = rdForceFieldHelpers.UFFOptimizeMolecule
        self.gen3d_max_attempts = self.gen3d_max_attempts
        if self.gen3d:
            # TODO use v2 or v3 if there are rings?
            self.gen3d_engine = rdDistGeom.ETKDGv3()
        else:
            self.gen3d_engine = None
        self.auto_iterations = max(1, self.auto_iterations)
        self.opt_add_h = self.add_h
        self.opt_force_trans_amide = self.force_trans_amide
        # TODO add support for flexible ring corners

    def process(
        self,
        mol_input,
    ):
        """process the molecule"""
        # print("MOLSS", mol_input)
        # print("MOLSS", mol_input.GetProp("Scrubber_was_here"))
        # print("MOLINPUT", mol_input.GetProp("_Name"))
        # add hydrogens if necessary
        if self.opt_add_h:
            mol = Chem.AddHs(mol_input)
            copy_mol_properties(mol_input, mol)
        else:
            mol = Chem.Mol(mol_input)
        report = {
            "state": None,
            "cycles": 0,
            "mol": mol,
            "iterations": self.ff_parms["maxIters"],
            "accepted" : 1  # 1: full, 0: partial, -1: rejected
        }
        # run the cycle once (better than try/except?)
        while True:
            if self._handbrake:
                report["state"] = "interrupted (handbrake)"
                report['accepted'] = -1
                break
            # generate 3D coordinates if requested
            if not self.gen3d_engine is None:
                report["gen3d_max_attempts"] = self.gen3d_max_attempts
                try:
                    rdDistGeom.EmbedMolecule(
                        mol, self.gen3d_engine
                    )  # , maxAttempts=self.gen3d_max_attempts)
                    # TODO check why disabled?
                except Exception as err:
                    report['state'] = err.__str__()
                    report['accepted'] = -1
                    # return report
                    break
                if mol.GetNumConformers() == 0:
                    report["state"] = "fail_3d"
                    report['accepted'] = -1
                    break
            if self.opt_force_trans_amide:
                self._fix_amide(mol)
            for steps in range(self.auto_iterations):
                report["cycles"] += 1
                # try:
                if True:
                    out = self.ff_optimize(mol, **self.ff_parms)
                # TODO this part of the code should be uncommented only at the verey end...
                # except Exception as err:
                #     report['state'] = err.__str__()
                #     report['accepted'] = -1
                #     # print("CAPTURED EXOTIC ERROR!", report['state'])
                #     # return report
                #     break
                if out == -1:
                    report["state"] = "ff_fail"
                    report['accepted'] = -1
                    break
                elif out == 0:
                    report["state"] = "mini_converged"
                    report['accepted'] = 1
                    break
                elif out == 1:
                    report["state"] = "mini_not_converged"
                    report['accepted'] = 0
            # return report
            break
        # print("RETURN", report, report['mol'].GetPropsAsDict())
        report['mol'] = PropertyMol(report['mol'])
        return report

    def _fix_amide(self, mol):
        """check that secondary amides are in trans (~180 deg) or trans-like ()
        configuration"""
        pattern = "[H][NX3;R0]([!#1])[CX3;R0](=[OX1])[#6]"
        found_amides = mol.GetSubstructMatches(Chem.MolFromSmarts(pattern))
        if not found_amides:
            return
        conf = mol.GetConformer(0)
        for a in found_amides:
            dihe = rdMolTransforms.GetDihedralDeg(conf, a[0], a[1], a[3], a[4])
            # print("DIHE", a[0], a[1], a[3], a[4], " IS ", dihe)
            if -90 > dihe < 90:
                rdMolTransforms.SetDihedralDeg(conf, a[0], a[1], a[3], a[4], 180)
                # print("\n\n\n FIXED AMIDE FIXED! \n\tTEST IT AND SEE IF IT WORKED!")


class GeometryGeneratorMPWorker(multiprocessing.Process, GeometryGenerator):
    """MP worker version of the GeometryGenerator based on multiprocessing queues

    self.queue_in   :  source of molecules to process
    self.queue_out  :  destination of processed molecules

    """

    def __init__(self, queue_in, queue_out, nice_level: int = None, strict:bool=False, **args):
        """initialize"""
        multiprocessing.Process.__init__(self)
        # set nice level if requested
        # if "nice_level" in args:
        if not nice_level is None:
            # print("NICENESS!!!")
            os.nice(nice_level)
        # initialize the parent class
        GeometryGenerator.__init__(self, **args)
        self.queue_in = queue_in
        self.queue_out = queue_out
        self.strict = strict
        if self.strict:
            self._success_cutoff = 0
        else:
            self._success_cutoff = -1

    def run(self):
        """overload of multiprocessing run method"""
        while True:
            if self._handbrake:
                print("trying to exit gracefully...")
            mol = self.queue_in.get()
            if mol is None or self._handbrake:
                # print("GEOM RECEIVED POISON PILL")
                self.queue_out.put(None)
                break
            # print("MOL", mol.GetPropsAsDict())
            report = self.process(mol)
            if report['accepted'] >= self._success_cutoff:
            # report["name"] = mol_name
                self.queue_out.put(report)
            # except:
            #     print("PROBLEMATIC MOLECULE captured...")
            #     continue


        return


class ParallelGeometryGenerator(object):
    """Parallelized (multiprocessing) 3D geometry builder
    instanciate multiple workers and connect them with the in/out queues

         queue_in   : queue from which molecules to be processed are pulled
         queue_out  : queue in which processed molecules are pushed
         geom_opts  : dictionary containig parmeters for the optimization to be passed to GeometryGenerator
         max_proc   : max number of parallel processes
         nice_level : nce level
         strict     : flag to define which molecules are accepted; if True,
                      only converged molecules are accepted, otherwise any minimized
                      molecule is accepted
    """
    def __init__(
        self,
        queue_in: multiprocessing.Queue = None,
        queue_out: multiprocessing.Queue = None,
        geom_opts: dict = GeometryGenerator.get_defaults(),
        max_proc: int = None,
        nice_level: int = None,
        strict:bool = False,
        _stop_at_defaults=False,
    ):
        self.queue_in = queue_in
        self.queue_out = queue_out
        self.geom_opts = geom_opts
        self.max_proc = max_proc
        self.nice_level = nice_level
        self.strict = strict
        if _stop_at_defaults:
            return
        print("=============================")
        print("self.queue_in", self.queue_in)
        print("self.queue_out", self.queue_out)
        print("self.geom_opts", self.geom_opts)
        print("self.max_proc ", self.max_proc)
        print("self.nice_level", self.nice_level)
        print("self.strict", self.strict)
        print("=============================")
        self.__workers = []
        if self.max_proc is None:
            self.max_proc = multiprocessing.cpu_count()
        if self.queue_in is None:
            self.queue_in = multiprocessing.JoinableQueue(maxsize=self.max_proc)
        if self.queue_out is None:
            self._queue_size = self.max_proc * 4
            self.queue_out = multiprocessing.Queue(maxsize=self.max_proc * 4)
        else:
            self._queue_size = self.max_proc
        for i in range(self.max_proc):
            # print("DOING THIS", self.queue_in, s)
            w = GeometryGeneratorMPWorker(
                self.queue_in,
                self.queue_out,
                self.nice_level,
                self.strict,
                **geom_opts,
            )
            self.__workers.append(w)
            w.start()
        print("%d workers initialized" % len(self.__workers))

    @classmethod
    def get_defaults(cls):
        """method to return the default values of init options"""
        return cls(_stop_at_defaults=True).__dict__

    def halt_workers(self):
        """function to stop all pending calculations"""
        for w in self.__workers:
            w._handbrake = True

    def terminate(self):
        """multiprocessing method override"""
        # ask all processes to terminate gracefully
        for idx, w in enumerate(self.__workers):
            w.terminate()
        for i in range(self._queue_size):
            self.queue_out.put(None)

    # def stop_calculations(self):
    #     """  """
    #     # self.queue_in.put(q, block=True)
    #     self.queue_in.put(None)
