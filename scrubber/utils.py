from rdkit import Chem
from rdkit.Chem.rdchem import Mol
from rdkit.Chem import rdDistGeom
from rdkit.Chem import AllChem

from rdkit.ForceField.rdForceField import ForceField


def optimize_conformers(mol: Mol, use_mmff: bool=True):
    optimized_energies = []
    
    for conf_id in range(mol.GetNumConformers()):
        if use_mmff and AllChem.MMFFHasAllMoleculeParams(mol):
            ff : ForceField = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol), confId=conf_id)
        else:
            ff : ForceField = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
        
        success = ff.Minimize(maxIts=400)
        # print(f"Success: {success}")
        # Check if minimization is successful. 
        if success == 0:
            energy = ff.CalcEnergy()
            optimized_energies.append((conf_id, energy))
        else:
            ff.Minimize(maxIts=800)
            energy == ff.CalcEnergy()
            optimized_energies.append((conf_id, energy))
    
    return optimized_energies

def add_conformers_to_mol(mol: Mol, conf_coords_list):
    mol = Chem.Mol(mol)  # Make a copy to avoid modifying the original mol
    mol.RemoveAllConformers()  # Clear any existing conformers

    for conf_id, coords in enumerate(conf_coords_list):
        conf = Chem.Conformer(mol.GetNumAtoms())
        for i, (x, y, z) in enumerate(coords):
            conf.SetAtomPosition(i, (float(x), float(y), float(z)))
        conf.SetId(conf_id)
        mol.AddConformer(conf, assignId=True)

    return mol

def find_best_conformer(mol: Mol, ps, num_confs=3):
    cids = rdDistGeom.EmbedMultipleConfs(mol, num_confs, ps)
    energies = optimize_conformers(mol, False)
    if not energies:
        raise ValueError("No conformers could be optimized during initial generation.")
    
    best_conf_id, best_energy = min(energies, key=lambda x: x[1])


    # Create a new molecule with only the best conformer
    best_mol = Chem.Mol(mol)
    conf = mol.GetConformer(best_conf_id)

    best_mol.RemoveAllConformers()

    best_mol.AddConformer(conf, assignId=True)

    return best_mol, [0]