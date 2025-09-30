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

from .ringfix import fix_rings
from .utils import find_best_conformer


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

        getForceField = lambda mol, confId: AllChem.MMFFGetMoleculeForceField(
            mol, AllChem.MMFFGetMoleculeProperties(mol), confId=confId

        )
    # This is just is a minimization with restraints to force the querry mol to match the template. 
    # If you chose espaloma as ff it will still minimize it at the end with that forcefield.
    elif ff == "mmff94s" or ff == "espaloma":

        getForceField = lambda mol, confId: AllChem.MMFFGetMoleculeForceField(
            mol,
            AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s"),

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
    max_ff_iter: int = 400,
    etkdg_rng_seed: int = 42,
    numconfs: int = 1,
    ff: str = "mmff94s",
    espaloma=None,
    template=None,
    template_smarts=None,
    ring_minimize=False,
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
        if ring_minimize:
            mol, cids = find_best_conformer(mol, ps, numconfs, max_ff_iter, ff)
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
        [coords_list.extend(fix_rings(mol, c, ring_minimize, energy_threshold, debug=debug, max_ff_iter=max_ff_iter, ff=ff)) for c in etkdg_coords]
   
    for coords in coords_list:
        c = Chem.Conformer(mol.GetNumAtoms())
        for i, (x, y, z) in enumerate(coords):
            c.SetAtomPosition(i, Point3D(x, y, z))
        mol.AddConformer(c, assignId=True)
    if ff not in ["uff", "mmff94", "mmff94s", "espaloma"]:
        raise RuntimeError(
            f"ff is {ff} but must be 'uff', 'mmff94', 'mmff94s', or 'espaloma'"
        )

    # if ring_minimize is used, then don't do the optimization again
    # instead just return what was given. 


    if (debug):
        from . import amine_flip as am
        print("XYZ coordinates produced by ringfix")
        for coords in coords_list:
            print(am.molToXYZ(mol, coords))

    if ring_minimize:
        final_mol = mol
    else:
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
