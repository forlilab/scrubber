#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
import numpy as np
import warnings
warnings.simplefilter("ignore", category=UserWarning)
from rdkit import Chem

class EspalomaMinimizer:
    
    def __init__(self,
                version='latest'
        ):
        try:
            print("importing espaloma...")
            from espaloma import get_model, Graph 
            from espaloma.graphs.deploy import openmm_system_from_graph
            print("imported espaloma.")
        except ImportError as err:
            raise ImportError("Espaloma is required") from err
        try:
            print("importing openff toolkit...")
            from openff.toolkit.topology import Molecule
            from openff.units.openmm import to_openmm
            print("imported openff toolkit.")
        except ImportError as err:
            raise ImportError("OpenFF is required") from err
        try:
            print("importing openmm...")
            import openmm.unit as Unit
            from openmm import VerletIntegrator
            from openmm.app import Simulation
            print("imported openmm")
        except ImportError as err:
            raise ImportError("OpenMM is required") from err
   
        # load pretrained espaloma model
        self.espaloma_model = get_model(version)

        # store methods in instance, otherwise they are out of scope
        self.CreateMolecule = Molecule
        self.EspalomaGraph = Graph
        self.System_from_Graph = openmm_system_from_graph
        self.Integrator = VerletIntegrator
        self.Simulator = Simulation
        self.Unit = Unit
        self.offquantity_to_openmm = to_openmm

    def minimize(self, mol):

        molecule = self.CreateMolecule.from_rdkit(mol, allow_undefined_stereo=True)

        # Make an OpenFF Topology so we can parameterize the system
        off_top = molecule.to_topology()

        # Convert the OpenFF Topology to OpenMM Topology
        omm_top = off_top.to_openmm()

        # create an Espaloma Graph object to represent the molecule of interest
        molecule_graph = self.EspalomaGraph(molecule)

        # apply a trained espaloma model to assign parameters
        self.espaloma_model(molecule_graph.heterograph)

        # create an OpenMM System for the specified molecule
        system = self.System_from_Graph(molecule_graph)

        # Set up integrator and simulation.
        integrator = self.Integrator(1*self.Unit.femtoseconds)
        simulation = self.Simulator(omm_top, system, integrator)
        mol_copy = self.CreateMolecule(molecule)
        mol_copy._conformers = None
        energies = []
        for conformer in molecule.conformers:

            conformer = self.offquantity_to_openmm(conformer)
            simulation.context.setPositions(conformer)

            simulation.minimizeEnergy()

            min_state = simulation.context.getState(getEnergy=True, getPositions=True)
            min_coords = min_state.getPositions()
            min_coords = np.array([ [atom.x, atom.y, atom.z] for atom in min_coords]) * self.Unit.nanometer

            mol_copy.add_conformer(min_coords)
            energies.append(min_state.getPotentialEnergy()/self.Unit.kilojoules_per_mole)
        rdmol = mol_copy.to_rdkit()

        return rdmol, energies


def _snap_to_int(value, tolerance=0.12):
    for inc in [-1, 0, 1]:
        if abs(value - int(value) - inc) <= tolerance:
            return int(value) + inc
    return None

def divide_int_gracefully(integer, weights):
    for weight in weights:
        if type(weight) not in [int, float, np.float32, np.float64] or weight < 0:
            raise ValueError("weights must be numeric and non-negative")
    if type(integer) is not int:
        raise ValueError("integer must be integer")
    inv_total_weight = 1.0 / sum(weights)
    shares = [w * inv_total_weight for w in weights]  # normalize
    result = [_snap_to_int(integer * s, tolerance=0.5) for s in shares]
    surplus = integer - sum(result)
    if surplus == 0:
        return result
    data = [(i, w) for (i, w) in enumerate(weights)]
    data = sorted(data, key=lambda x: x[1], reverse=True)
    idxs = [i for (i, _) in data]
    groups = []
    last_weight = None
    for i in idxs:
        if weights[i] == last_weight:
            groups[-1] += 1
        else:
            groups.append(1)
        last_weight = weights[i]

    # iterate over all possible combinations of groups
    # this is potentially very slow
    nr_groups = len(groups)
    for j in range(1, 2**nr_groups):
        n_changes = 0
        combo = []
        for grpidx in range(nr_groups):
            is_changed = bool(j & 2**grpidx)
            combo.append(is_changed)
            n_changes += is_changed * groups[grpidx]
        if n_changes == abs(surplus):
            break

    # add or subtract 1 to distribute surplus
    increment = surplus / abs(surplus)
    index = 0
    for i, is_changed in enumerate(combo):
        if is_changed:
            for j in range(groups[i]):
                result[idxs[index]] += increment
                index += 1

    return result


def rectify_charges(q_list, net_charge=None, decimals=3) -> list[float]:

    if net_charge is None:
        net_charge = _snap_to_int(sum(q_list), tolerance=0.15)
        if net_charge is None:
            msg = "net charge could not be predicted from input q_list. (residual is beyond tolerance) "
            msg = "Please set the net_charge argument directly"
            raise RuntimeError(msg)
    elif type(net_charge) != int:
        raise TypeError("net charge must be an integer")

    fstr = "%%.%df" % decimals
    charges_dec = [float(fstr % q) for q in q_list]
    surplus = net_charge - sum(charges_dec)
    surplus_int = _snap_to_int(10**decimals * surplus)

    if surplus_int == 0:
        return charges_dec

    weights = [abs(q) for q in q_list]
    surplus_int_splits = divide_int_gracefully(surplus_int, weights)
    for i, increment in enumerate(surplus_int_splits):
        charges_dec[i] += 10**-decimals * increment

    return charges_dec


class EspalomaCharger:
    def __init__(self, version="latest"):
        print("importing espaloma and openff toolkit...")
        import espaloma
        try:
            from openff.toolkit import Molecule
        except ImportError:
            print("A recent version of OpenFF is required for Espaloma charges")

        print("imported espaloma and openff toolkit.")
        self.espaloma_model = espaloma.get_model(version)
        self.Molecule = Molecule
        self.espaloma = espaloma

    def get_espaloma_charges(self, rdkit_mol: Chem.Mol):
        '''
            compute espaloma charges from rdkit mol,
            return espaloma charges in form of array
        '''
        openff_mol = self.Molecule.from_rdkit(
            rdkit_mol,
            hydrogens_are_explicit=True,
            allow_undefined_stereo=True,
        )
        molgraph = self.espaloma.Graph(openff_mol)
        self.espaloma_model(molgraph.heterograph)
        charges = [float(q) for q in molgraph.nodes["n1"].data["q"]]
        return charges

    def mol_with_charges(self, 
                    rdkit_mol: Chem.Mol , 
                    decimals: int = 3, 
                    prop: str ="atom.dprop.PartialCharge"):
        '''
            compute the Espaloma charges of the rdkit mol, return a new mol with charges
            as properties. 
        '''

        mol = Chem.Mol(rdkit_mol) # new mol object
        charges = self.get_espaloma_charges(mol)
        charges = rectify_charges(charges, decimals=decimals)
        fstr = "%%.%df" % decimals
        charges = " ".join([fstr % q for q in charges]) 
        mol.SetProp(prop, charges)
        return mol

class NaglCharger:
    def __init__(self, version="latest"):
        print("importing the openff toolkit...")
        try:
            from openff.toolkit import Molecule
        except ImportError:
            print("A recent version of OpenFF is required for NAGL charges")
        self.Molecule = Molecule



    def get_nagl_charges(self, rdkit_mol: Chem.Mol):
        '''
            compute nagl charges from rdkit mol,
            return nagl charges in form of array
        '''
        openff_mol = self.Molecule.from_rdkit(
            rdkit_mol,
            hydrogens_are_explicit=True,
            allow_undefined_stereo=True,
        )

        try:
            openff_mol.assign_partial_charges(
                partial_charge_method="openff-gnn-am1bcc-1.0.0.pt"
            )
            charges = openff_mol.partial_charges.magnitude
        except Exception as e:
            print("NAGL charge computation failed with with exception:")
            print(e)
            print("Make sure you've installed the latest version of openff")
        return charges

    def mol_with_charges(self, 
                         rdkit_mol: Chem.Mol, 
                         decimals: int = 3, 
                         prop: str = "atom.dprop.PartialCharge"):
        '''
            compute the NAGL charges of the rdkit mol, return a new mol with charges
            as properties. 
        '''
        
        #
        mol = Chem.Mol(rdkit_mol)
        charges = self.get_nagl_charges(rdkit_mol)
        charges = rectify_charges(charges, decimals=decimals)
        fstr = "%%.%df" % decimals
        charges = " ".join([fstr % q for q in charges]) 
        mol.SetProp(prop, charges)
        return mol