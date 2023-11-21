#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
import numpy as np
import warnings
warnings.simplefilter("ignore", category=UserWarning)

class EspalomaMinimizer:
    
    def __init__(self,
                version='latest'
        ):
        try:
            from espaloma import get_model, Graph 
            from espaloma.graphs.deploy import openmm_system_from_graph
        except ImportError:
            raise ImportError("Espaloma is required")
        try:
            from openff.toolkit.topology import Molecule
            from openff.units.openmm import to_openmm
        except ImportError:
            raise ImportError("OpenFF is required")
        try:
            import openmm.unit as Unit
            from openmm import VerletIntegrator
            from openmm.app import Simulation
        except ImportError:
            raise ImportError("OpenMM is required")
   
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

    def minim_espaloma(self, mol):

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


