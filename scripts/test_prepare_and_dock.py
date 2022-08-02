import sys, os
from rdkit import Chem
from time import sleep
sys.path.append("../")
try:
    from scrubber.core import ScrubberCore
except Exception as e:
    print("The script must be executed from the script/ directory")
    raise e

def dock(mol):
    print("Docking molecule: |-",  end="")
    for i in range(10):
        # 10328**99*3
        print("-", end="")
        sys.stdout.flush()
        # sleep(0.05)
    print("| DONE!")

def test_scrubber(config, smiles_list):
    scrub = ScrubberCore(config)
    errors = []
    for smi in smiles_list:
        print("\n>>> Processing SMILES", smi.strip())
        if 'errorium' in smi:
            print("SCIPPING")
            continue
        mol = Chem.MolFromSmiles(smi)
        if not mol is None:
            for m in scrub.process(mol):
                dock(mol)
        else:
            print("failed at start...")
        # errors += scrub.errors
    print("RESCUING PROBLEMATICS")
    problematic = scrub.get_problematic()
    print("Problematic molecules encountered [%d]:" % len(problematic), problematic)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        smiles_list =  ["O=C([C@H](CC1=CNC=N1)N)O histidine", 'CcCCcC problematic' ]
    else:
        with open(sys.argv[1]) as fp:
            smiles_list = fp.readlines()
    # get and customize configuration (these should be the defaults)
    print("PROCESSING %d MOLS..." % len(smiles_list))
    config = ScrubberCore.get_defaults()
    config["isomers"]["values"]["ph_datafile"] = "..//scrubber/data/test_model.txt"
    config["isomers"]["values"]["stereo_enum"] = "all"
    config["geometry"]["values"]["strict"] = True
    config["geometry"]["values"]["auto_iter_cycles"] = 1

    # initialize parallel scrubber
    config['general']['values']["max_proc"] = 8
    print("\n\n=============================")
    print("test PARALLEL scrubber")
    print("=============================")
    # scrub = ScrubberCore(config)

    test_scrubber(config, smiles_list)

    # initialize serial scrubber
    config['general']['values']["max_proc"] = 1
    print("\n\n=============================")
    print("test SERIAL scrubber")
    print("=============================")
    test_scrubber(config, smiles_list)

# HERE IT HANDS TODO
