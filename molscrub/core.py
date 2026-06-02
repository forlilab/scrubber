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
from .openff_toolkit import EspalomaMinimizer
from .openff_toolkit import EspalomaCharger, NaglCharger
from .geometry import gen3d

import joblib
import tempfile
import urllib.request
import shutil
from pathlib import Path


#handle pka model downloading. 
MODEL_URL = "https://github.com/forlilab/pkaPrediction/raw/refs/heads/main/models/ETR_latest_compressed.joblib"
MODEL_FILENAME = "ETR_latest_compressed.joblib"

CACHE_DIR = Path.home() / ".cache" / "molscrub"
MODEL_PATH = CACHE_DIR / MODEL_FILENAME


def download_model():
    """
    Download pKa model and store in chace directory
    `$HOME/.cache/molscrub`
    """
    CACHE_DIR.mkdir(parents=True, exist_ok=True)

    print()
    print("Downloading pickled ML model (one time only) from:")
    print(MODEL_URL)
    print("This will be saved in $HOME/.cache/molscrub/")
    print()

    with tempfile.NamedTemporaryFile(delete=False) as tmp_file:
        urllib.request.urlretrieve(MODEL_URL, tmp_file.name)
        tmp_path = Path(tmp_file.name)

    shutil.move(str(tmp_path), MODEL_PATH)


def verify_model():
    """
    Check if model is in cach directory. If not, download. 
    """
    if not MODEL_PATH.exists():
        download_model()
    else:
        print()
        print("pKa ML model found: ", MODEL_PATH)
        print("download skipped")
        print()


def load_model():
    """
    Load and return etr1 model. Checks if it's available and 
    downloads into .cache directory of it's not. 
    """
    verify_model()
    model = joblib.load(MODEL_PATH)
    return model


class Scrub:

    def __init__(
        self,
        ph_low=7.4,
        ph_high=None,
        pka_fname=None,
        pka_model="rules",
        model_file=None,
        tauto_fname=None,
        skip_acidbase=False,
        skip_tautomers=False,
        skip_ringfix=False,
        skip_gen3d=False,
        template=None,
        template_smarts=None,
        do_gen2d=False,
        max_ff_iter=400,
        skip_etkdg=False,
        num_internal_confs=3,
        min_output_confs=1,
        etkdg_rng_seed=None,
        use_random_coords=False,
        ff="mmff94s",
        ring_minimize=False,
        energy_threshold=0.5,
        keep_all_frags=False,
        charge_model=None,
        debug=False,
        num_etkdg_attempts=1,
    ):
        
        # this is needed if using API instead of CLI
        if pka_model != "rules" and model_file == None:
            model_file = load_model()

        if pka_fname is None:
            self.acid_base_conjugator = AcidBaseConjugator.from_default_data_files(model=model_file)
        else:
            reactions = AcidBaseConjugator.parse_reaction_file(pka_fname)
            self.acid_base_conjugator = AcidBaseConjugator(reactions)
        if tauto_fname is None:
            self.tautomerizer = Tautomerizer.from_default_data_files()
        else:
            reactions, keepmax_smarts = Tautomerizer.parse_tautomers_config_file(tauto_fname)
            self.tautomerizer = Tautomerizer(reactions, keepmax_smarts)
        self.ph_low = ph_low
        if ph_high is None:
            ph_high = ph_low
        self.ph_high = ph_high
        self.do_acidbase = not skip_acidbase
        self.pka_model = pka_model
        self.model_file = model_file
        self.do_tautomers = not skip_tautomers
        self.skip_ringfix = (
            skip_ringfix  # not avoiding negative to pass directly to gen3d
        )
        self.ring_minimize = ring_minimize
        self.skip_etkdg = skip_etkdg
        self.energy_threshold = energy_threshold
        self.do_gen3d = not skip_gen3d
        self.template = template
        self.template_smarts = template_smarts
        self.do_gen2d = do_gen2d
        self.max_ff_iter = max_ff_iter
        self.num_internal_confs = num_internal_confs
        if min_output_confs < 1:
            raise ValueError(f"min_output_confs must be greater than zero. Consider skip_gen3d=True instead.")
        self.min_output_confs = min_output_confs
        self.etkdg_rng_seed = (
            etkdg_rng_seed if etkdg_rng_seed else random.randint(0, 1000000)
        )
        self.use_random_coords = use_random_coords
        self.ff = ff
        self.keep_all_frags = keep_all_frags
        self.charge_model = charge_model
        self.debug = debug
        self.num_etkdg_attempts = num_etkdg_attempts

        if ff == "espaloma":
            self.espaloma = EspalomaMinimizer()
        else:
            self.espaloma = None

        if charge_model == "espaloma":
            self.charger = EspalomaCharger()
        elif charge_model == "nagl":
            self.charger = NaglCharger()
        elif charge_model is None:
            self.charger = None
        else:
            raise ValueError(f"{charge_model=} not supported")

    def __call__(self, input_mol: Chem.Mol) -> Chem.Mol:
        """
        the main scrubbing function, produces a list
        of scrubbed molecules
        """

        #check for fragments and keep the largest. 
        frags = Chem.GetMolFrags(input_mol, asMols=True)
        if len(frags) > 1 and not self.keep_all_frags:
            if input_mol.HasProp("_Name") and input_mol.GetProp("_Name"):
                name = input_mol.GetProp("_Name") 
                print(f"Molecule {name} contains {len(frags)} fragments")
            else:
                print(f"Unnamed molecule contains {len(frags)} fragments")
            print("Only the largest fragment will be processed")
            input_mol = max(frags, key=lambda x: x.GetNumAtoms())

        ref_mol = Chem.Mol(input_mol) # keep a copy of the original mol
        input_mol = Chem.RemoveHs(input_mol)
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
            
            if self.skip_etkdg:
                from .geometry import copy_mcs_coordinates
                # constrained embedding from input mol
                print("skip_etkdg choosen, using constrained embedding with reference coordinates.")
                pool = copy_mcs_coordinates(ref_mol, pool)

            for mol in pool:
                confs = []
                loop_count = 0
                while len(confs) < self.min_output_confs:
                    loop_count += 1
                    mol_out = None
                    last_exc = None
                    try:
                        mol_out = gen3d(
                            mol,
                            skip_ringfix=self.skip_ringfix,
                            max_ff_iter=self.max_ff_iter,
                            skip_etkdg=self.skip_etkdg,
                            etkdg_rng_seed=self.etkdg_rng_seed + loop_count if self.etkdg_rng_seed != -1 else -1,
                            use_random_coords=self.use_random_coords,
                            num_internal_confs=self.num_internal_confs,
                            ff=self.ff,
                            espaloma=self.espaloma,
                            template=self.template,
                            template_smarts=self.template_smarts,
                            ring_minimize=self.ring_minimize,
                            energy_threshold=self.energy_threshold,
                            debug=self.debug,
                            num_etkdg_attempts = self.num_etkdg_attempts
                        )
                    except Exception as e:
                        last_exc = e

                    if mol_out is None:
                        raise last_exc
                    confs += [conf for conf in mol_out.GetConformers()]
                mol_out = Chem.Mol(mol_out)  # got MemoryError without this
                mol_out.RemoveAllConformers()
                for conf in confs:
                    mol_out.AddConformer(conf, assignId=True)
                output_mol_list.append(mol_out)
        elif self.do_gen2d:  # useful to write SD files
            output_mol_list = []
            for mol in pool:
                AllChem.Compute2DCoords(mol)
                output_mol_list.append(mol)
        else:
            output_mol_list = pool


        if self.charge_model is not None:
            output_mol_list = list(map(self.charger.mol_with_charges, output_mol_list))

        return output_mol_list

    def scrub_and_catch_errors(self, input_mol):
        log = {}
        if input_mol is None:
            log["input_mol_none"] = True
            isomer_list_if_ok_else_input = []
            return isomer_list_if_ok_else_input, log
        log["input_mol_none"] = False
        try:
            isomer_list_if_ok_else_input = self(input_mol)
        except Exception as e:
            log["exception"] = e
            isomer_list_if_ok_else_input = input_mol
            if self.debug:
                raise e
        return isomer_list_if_ok_else_input, log