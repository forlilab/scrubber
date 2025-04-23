import multiprocessing
from multiprocessing.synchronize import Event as EventClass
import os

from rdkit import Chem
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdForceFieldHelpers
from rdkit.Chem import rdMolTransforms
from rdkit.Chem.PropertyMol import PropertyMol
from rdkit.Chem import rdMolAlign
from rdkit.Geometry import Point3D



from rdkit.Chem.rdchem import Mol 
from rdkit.Chem import AllChem

from .common import ScrubberBase, copy_mol_properties
from .ringfix import fix_rings
from .utils import find_best_conformer

# from .ringcorners import RingManager


# TODO
# print add proper reporting for each molecule that gets put into the queue_err, so a Scrubber_error field can be added in the output molecule


class GeometryGenerator(ScrubberBase):
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
        auto_iter_cycles: int = 10,
        gen3d: bool = True,
        gen3d_max_attempts: int = 3,
        fix_ring_corners: bool = False,
        preserve_mol_properties: bool = True,
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
        self.auto_iter_cycles = auto_iter_cycles
        self.gen3d = gen3d
        self.gen3d_max_attempts = gen3d_max_attempts
        self.fix_ring_corners = fix_ring_corners
        self.preserve_mol_properties = preserve_mol_properties
        if _stop_at_defaults:
            return
        # this flag is used to control if at any time the minimization has been asked to stop
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
        self.auto_iter_cycles = max(1, self.auto_iter_cycles)
        self.opt_add_h = self.add_h
        self.opt_force_trans_amide = self.force_trans_amide
        # TODO add support for flexible ring corners

    def process(
        self,
        mol_input,
    ):
        """process the molecule"""
        # add hydrogens if necessary
        if self.opt_add_h:
            mol = Chem.AddHs(mol_input)
            copy_mol_properties(mol_input, mol)
        else:
            mol = Chem.Mol(mol_input)
        report = {
            "state": None,
            "iter_cycles": 0,
            "mol": mol,
            "iterations": self.ff_parms["maxIters"],
            "accepted": 1,  # 1:full, 0:partial, -1:rejected
        }
        # run the cycle once (better than try/except?)
        while True:
            try:
                # generate 3D coordinates if requested
                if not self.gen3d_engine is None:
                    report["gen3d_max_attempts"] = self.gen3d_max_attempts
                    try:
                        rdDistGeom.EmbedMolecule(
                            mol, self.gen3d_engine)
                    except Exception as err:
                        report["state"] = err.__str__()
                        report["accepted"] = -1
                        # return report
                        break
                    if mol.GetNumConformers() == 0:
                        report["state"] = "fail_3d"
                        report["accepted"] = -1
                        break
                if self.opt_force_trans_amide:
                    self._fix_amide(mol)
                for steps in range(self.auto_iter_cycles):
                    report["iter_cycles"] += 1
                    try:
                        out = self.ff_optimize(mol, **self.ff_parms)
                    except Exception as err:
                        report['state'] = err.__str__()
                        report['accepted'] = -1
                        print("CAPTURED EXOTIC ERROR!", report['state'])
                        return report
                        break
                    if out == -1:
                        report["state"] = "ff_fail"
                        report["accepted"] = -1
                        break
                    elif out == 0:
                        report["state"] = "mini_converged"
                        report["accepted"] = 1
                        break
                    elif out == 1:
                        report["state"] = "mini_not_converged"
                        report["accepted"] = 0
            except RuntimeError as err:
                report["state"] =  err.__str__()  #"ff_fail"
                report["accepted"] = -1
                break
            break
        report["mol"] = PropertyMol(report["mol"])
        if self.fix_ring_corners:
            print("RING CORNERS FLIPPED HERE")
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

geom_default = GeometryGenerator.get_defaults()

class GeometryGeneratorMPWorker(multiprocessing.Process, GeometryGenerator):
    """MP worker version of the GeometryGenerator based on multiprocessing queues

    self.queue_in   :  source of molecules to process
    self.queue_out  :  destination of processed molecules

    """

    def __init__(
        self,
        queue_in: multiprocessing.Queue,
        queue_out: multiprocessing.Queue,
        queue_err: multiprocessing.Queue=None,
        nice_level: int = None,
        strict: bool = False,
        handbrake: EventClass = None,
        # geom_opts: dict = GeometryGenerator.get_defaults(),
        add_h: bool = geom_default["add_h"],
        force_trans_amide: bool = geom_default["force_trans_amide"],
        force_field: str = geom_default["force_field"],
        max_iterations: int = geom_default["max_iterations"],
        auto_iter_cycles: int = geom_default["auto_iter_cycles"],
        gen3d: bool = geom_default["gen3d"],
        gen3d_max_attempts: int = geom_default["gen3d_max_attempts"],
        fix_ring_corners: bool = geom_default["fix_ring_corners"],
        preserve_mol_properties: bool = geom_default["preserve_mol_properties"],
    ):
        """initialize"""
        multiprocessing.Process.__init__(self)
        # set nice level if requested
        # if "nice_level" in args:
        if not nice_level is None:
            # print("NICENESS!!!")
            os.nice(nice_level)
        # initialize the parent class
        GeometryGenerator.__init__(
            self,
            add_h,
            force_trans_amide,
            force_field,
            max_iterations,
            auto_iter_cycles,
            gen3d,
            gen3d_max_attempts,
            fix_ring_corners,
            preserve_mol_properties,
        )
        self.queue_in = queue_in
        self.queue_out = queue_out
        self.queue_err = queue_err
        self.strict = strict
        self.handbrake = handbrake
        if self.strict:
            self._success_cutoff = 0
        else:
            self._success_cutoff = -1

    def run(self):
        """overload of multiprocessing run method"""
        while True:
            try:
                if self.handbrake.is_set():
                    # print("WORKER__NAME: trying to exit gracefully...")
                    self.queue_out.put(None, block=True)
                    break
                mol = self.queue_in.get()
                if mol is None: # or self.handbrake.is_set():
                    # print("FOUND POISON PILL INGEOMETRY")
                    self.queue_out.put(None, block=True)
                    break
                # print("MOL", mol.GetPropsAsDict())
                report = self.process(mol)
                if report["accepted"] > self._success_cutoff:
                    # report["name"] = mol_name
                    self.queue_out.put(report, block=True)
                else:
                    if self.queue_err is None:
                        continue
                    self.queue_err.put( ("geom_" + report['state'], report['mol']), block=True)
                            # , report['mol'].GetProp("_Name")) )
                # except:
                #     print("PROBLEMATIC MOLECULE captured...")
                #     continue
            except KeyboardInterrupt:
                # print("[geom] Caught Ctrl-C...")
                self.queue_out.put(None, block=True)
                return
        return

class ParallelGeometryGenerator():
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
        queue_err: multiprocessing.Queue = None,
        # geom_opts: dict = GeometryGenerator.get_defaults(),
        add_h: bool = geom_default["add_h"],
        force_trans_amide: bool = geom_default["force_trans_amide"],
        force_field: str = geom_default["force_field"],
        max_iterations: int = geom_default["max_iterations"],
        auto_iter_cycles: int = geom_default["auto_iter_cycles"],
        gen3d: bool = geom_default["gen3d"],
        gen3d_max_attempts: int = geom_default["gen3d_max_attempts"],
        fix_ring_corners: bool = geom_default["fix_ring_corners"],
        preserve_mol_properties: bool = geom_default["preserve_mol_properties"],
        max_proc: int = None,
        nice_level: int = None,
        handbrake: EventClass = None,
        strict: bool = False,
        _stop_at_defaults=False,
    ):
        self.queue_in = queue_in
        self.queue_out = queue_out
        self.queue_err = queue_err
        self.add_h = add_h
        self.force_trans_amide = force_trans_amide
        self.force_field = force_field
        self.max_iterations = max_iterations
        self.auto_iter_cycles = auto_iter_cycles
        self.gen3d = gen3d
        self.gen3d_max_attempts = gen3d_max_attempts
        self.fix_ring_corners = fix_ring_corners
        self.preserve_mol_properties = preserve_mol_properties
        self.max_proc = max_proc
        self.nice_level = nice_level
        self.strict = strict
        self.handbrake = handbrake
        if _stop_at_defaults:
            return
        # print("============================= PARELL GEOM")
        # print("self.queue_in", self.queue_in)
        # print("self.queue_out", self.queue_out)
        # # print("self.geom_opts", self.geom_opts)
        # print("self.max_proc ", self.max_proc)
        # print("self.nice_level", self.nice_level)
        # print("self.strict", self.strict)
        # print("=============================")
        self.__workers = []
        if self.max_proc is None:
            self.max_proc = multiprocessing.cpu_count()
            print(
                "> parallel geometry initialized with independent multiprocessors: %d"
                % self.max_proc
            )
        if self.queue_in is None:
            self.queue_in = multiprocessing.JoinableQueue(maxsize=self.max_proc)
            print("> parallel geometry initialized with independent queue in")
        if self.queue_out is None:
            self._queue_size = self.max_proc * 4
            self.queue_out = multiprocessing.Queue(maxsize=self.max_proc * 4)
            print(
                "> parallel geometry initialized with independent queue out with size:%d"
                % self._queue_size
            )
        else:
            self._queue_size = self.max_proc
        for i in range(self.max_proc):
            w = GeometryGeneratorMPWorker(
                queue_in = self.queue_in,
                queue_out = self.queue_out,
                queue_err = self.queue_err,
                nice_level = self.nice_level,
                strict = self.strict,
                add_h = self.add_h,
                force_trans_amide = self.force_trans_amide,
                force_field = self.force_field,
                max_iterations = self.max_iterations,
                auto_iter_cycles = self.auto_iter_cycles,
                gen3d = self.gen3d,
                gen3d_max_attempts = self.gen3d_max_attempts,
                fix_ring_corners = self.fix_ring_corners,
                preserve_mol_properties = self.preserve_mol_properties,
                handbrake = self.handbrake,
            )
            self.__workers.append(w)
            # w.daemon = True
            w.start()
        print("[ %d geometry workers initialized ]" % len(self.__workers))

    def join(self):
        """function to wrap the join functions of the workers """
        for w in self.__workers:
            # print("ParallelGeometryGenerator> sub-join:", w)
            w.join()

    @classmethod
    def get_defaults(cls):
        """method to return the default values of init options"""
        return cls(_stop_at_defaults=True).__dict__

    def is_alive(self):
        """function to implement a similar behavior of multiprocessing
        is_alive() method; return true if at least one of the workers
        initialized is still alive"""
        for w in self.__workers:
            if w.is_alive():
                return True
        return False

    def halt_workers(self):
        """function to stop all pending calculations"""
        for w in self.__workers:
            w._handbrake = True

    def terminate(self):
        """multiprocessing method override"""
        # ask all processes to terminate gracefully
        for i in range(self.max_proc):
            print("GEOW: PILLL IN PARALLEL")
            self.queue_in.put(None, block=True)
        self.queue_in.close()

        for idx, w in enumerate(self.__workers):
            print("GEOW: Terminating idx", w)
            print("GEOW: Is alive?", w.is_alive())
            # print("DIR WORKERS", dir(w))
            w.terminate()
            print("GEOW: Is alive?", w.is_alive())
            w.join()
            print("GEOW: Is alive?", w.is_alive())
            # w.terminate()
            # w.kill()
            # w.close()
            # print("NOW->", w)
        # for i in range(self._queue_size):
        #     self.queue_out.put(None)


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
    use_energy=False,
    energy_threshold=0.5,
    debug=False
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
        # if ring is minimized, take best of numconfs = 3
        if use_energy:
            mol, cids = find_best_conformer(mol, ps, numconfs)
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
        [coords_list.extend(fix_rings(mol, c, use_energy, energy_threshold, debug=debug)) for c in etkdg_coords]
   
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
