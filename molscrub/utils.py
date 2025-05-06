from rdkit import Chem
from rdkit.Chem.rdchem import Mol
from rdkit.Chem import rdDistGeom
from rdkit.Chem import AllChem

from rdkit.ForceField.rdForceField import ForceField


def optimize_conformers(mol: Mol, use_mmff: bool=True, max_ff_iter: int = 400):
    optimized_energies = []
    
    for conf_id in range(mol.GetNumConformers()):
        if use_mmff and AllChem.MMFFHasAllMoleculeParams(mol):
            ff : ForceField = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol), confId=conf_id)
        else:
            ff : ForceField = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
        
        success = ff.Minimize(maxIts=max_ff_iter)
        # print(f"Success: {success}")
        # Check if minimization is successful. 
        if success == 0:
            energy = ff.CalcEnergy()
            optimized_energies.append((conf_id, energy))
        else:
            ff.Minimize(maxIts=2*max_ff_iter)
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

def find_best_conformer(mol: Mol, ps, num_confs=3, max_ff_iter=400):
    cids = rdDistGeom.EmbedMultipleConfs(mol, num_confs, ps)
    energies = optimize_conformers(mol, False, max_ff_iter)
    if not energies:
        raise ValueError("No conformers could be optimized during initial generation.")
    
    best_conf_id, best_energy = min(energies, key=lambda x: x[1])


    # Create a new molecule with only the best conformer
    best_mol = Chem.Mol(mol)
    conf = mol.GetConformer(best_conf_id)

    best_mol.RemoveAllConformers()

    best_mol.AddConformer(conf, assignId=True)

    return best_mol, [0]

#debug
def write_conformers_to_sdf(mol, filename="test.sdf"):
    writer = Chem.SDWriter(filename)
    
    for conf_id in range(mol.GetNumConformers()):
        mol.SetProp("_Name", f"Conformer {conf_id}")  # Optional: Label conformers
        writer.write(mol, confId=conf_id)
    
    writer.close()
    print(f"All conformers written to {filename}")
