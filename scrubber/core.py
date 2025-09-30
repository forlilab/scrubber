import multiprocessing


# IMPORTANT RESOURCES, INVESTIGATE
# https://sefiks.com/2021/07/05/handling-hang-in-python-multiprocessing/
# https://pythonspeed.com/articles/python-multiprocessing/ (alterantive to fork() )

import random

from rdkit import Chem
from rdkit.Chem import AllChem

from .geometry import find_best_conformer

from .protonate import AcidBaseConjugator
from .protonate import Tautomerizer
from .common import UniqueMoleculeContainer
from .espaloma_minim import EspalomaMinimizer
from .geometry import gen3d

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
        max_ff_iter=400,
        numconfs=1,
        etkdg_rng_seed=None,
        ff="mmff94s",
        ring_minimize=False,
        energy_threshold=0.5,
        debug=False
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
        self.ring_minimize = ring_minimize
        self.energy_threshold = energy_threshold
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
        self.debug = debug

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
                    ring_minimize = self.ring_minimize,
                    energy_threshold = self.energy_threshold,
                    debug=self.debug
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

