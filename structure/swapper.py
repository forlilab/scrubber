import itertools
try:
    from openbabel import openbabel as ob
except:
    import openbabel as ob  
    from molcopy import MoleculeDuplicator


# TODO add protonated nitrogen only option
# TODO convert to yeld

class ChiralEnumerator3D:
    """ class to enumerate chiral centers in 3D structures 
        
        max_iter    :   maximum number of enantiomers to be enumerated (default: 32, i.e., 5 chiral centers)
        strict      :   do not generate any enantiomers if a single intra-ring chiral center is found
                        (default: False, generate chiral centers that can be enumerated)
        nitro_protinated_only:  enumerate only protomers of quaternary amines
    
    """

    def __init__(self, mol, max_iter=32, strict=False, refine=False, nitro_protonated_only=False):
        self.mol = mol
        if not self.mol.Has3D():
            print("Error: the molecule does not have 3D coordinates")
            raise Exception
        self._mask = None
        self._MAX_ITER=max_iter
        self._strict = strict
        self.valid = True
        self.chiral = False
        self.perceive_chirals()
        self._refine = refine

    def perceive_chirals(self):
        """ enumerate chiral centers in the molecule """
        self._chiral_atoms = {}
        for a in ob.OBMolAtomIter(mol):
            if a.IsChiral():
                self.chiral = True
                idx = a.GetIdx()
                neigh = self.get_neighbors(a)
                if neigh==None:
                    #print("No swappable atoms for this chiral center! %d" % idx)
                    if self._strict:
                        self.valid = False
                    continue
                self._chiral_atoms[idx] = neigh
        if len(self._chiral_atoms)==0:
            self.valid = False
            return
        # enumerate all the possible permutations of neighbors swapping, 
        # except the input conformation; 1=swap, 0=keep
        # 0...00 - 1...11
        self._mask = list(map(list, itertools.product([0, 1], repeat=len(self._chiral_atoms))))[1:]
        if len(self._mask)>self._MAX_ITER:
            print("WARNING: permutations [%d] exceed maximum number of iterations [%d]" % (len(self._mask), self._MAX_ITER))
            self._mask = self._mask[self._MAX_ITER:]

    def get_neighbors(self, atom):
        """ retrieve neighbors to be swapped 
            only neighbors that are not connected by ring bonds are accepted
        """
        neigh = []
        idx = atom.GetIdx()
        for b in ob.OBAtomBondIter(atom):
            if b.IsInRing():
                continue
            bond_atoms = (b.GetBeginAtomIdx(), b.GetEndAtomIdx())
            if bond_atoms[0] == idx:
                neigh.append(bond_atoms[1])
            else:
                neigh.append(bond_atoms[0])
            if len(neigh)==2:
                return neigh
        return None

    def get_enantiomers(self,unique=True):
        """ generate enantiomers for the molecule """
        if not self.valid:
            print("No enantiomers are possible for this molecule")
            return None
        if not self.chiral:
            print("NO CHIRAL")
            return
        if self._mask == None:
            print("NO MASK")
            return
        out = []
        builder = ob.OBBuilder()
        for enantiomer in self._mask:
            #print("enumerating configuration", enantiomer)
            mol = ob.OBMol(self.mol)
            chiral_list = list(self._chiral_atoms.keys())
            for idx, chiral_center in enumerate(chiral_list):
                enantiomer_id = enantiomer[idx]
                if enantiomer_id == 0:
                    continue
                neigh = self._chiral_atoms[chiral_center]
                # print(">swapping atom %d:   %d <-> %d" % (chiral_center, neigh[0], neigh[1]))
                result = builder.Swap(mol, chiral_center, neigh[0], chiral_center, neigh[1] )
            out.append(mol)
        if unique:
            out = self._unique_enantiomers(out)
        if self._refine:
            self._refine_rotamers(out)
        return out
            

    def _unique_enantiomers(self, mol_list):
        """ save unique enantiomers """
        unique = {}
        canonical = ob.OBConversion()
        canonical.SetOutFormat('can')
        for mol in mol_list:
            can = canonical.WriteString(mol).split()[0]
            unique[can] = mol
        return list(unique.values())

    def _refine_rotamers(self, mol_list, num_conf = 15, geom_steps=0):
        """ perform rotamer minimization """
        mmff = 'mmff94s'
        uff = 'uff'
        for mol in mol_list:
            ff= ob.OBForceField.FindType(mmff)
            ff.SetLogLevel=ob.OBFF_LOGLVL_NONE
            setup_status = ff.Setup(mol)
            if not setup_status:
                print("Warning no MMFF95s parameters for", mol)
                ff= ob.OBForceField.FindType(uff)
                ff.SetLogLevel=ob.OBFF_LOGLVL_NONE
            ff.WeightedRotorSearch(num_conf, geom_steps)
            ff.UpdateCoordinates(mol)





if __name__=='__main__':
    import os
    import sys

    def getNameExt(fname):
        """ extract name and extension from the input file, removing the dot
            filename.ext -> [filename, ext]
        """
        name, ext = os.path.splitext(fname)
        return name, ext[1:] #.lower()

    def loadMolecule(fname=None, ftype=None, string=None, setFormalCharge=True):
        """ load molecule with openbabel"""
        if ftype == None:
            n, ftype = getNameExt(fname)
            ftype = ftype.lower()
        if fname == None and (string == None):
            print("NONE NONE")
            raise Exception

        # print("CALLED", fname)
        mol = ob.OBMol()
        conv = ob.OBConversion()
        conv.SetInFormat(ftype)
        if fname:
            conv.ReadFile(mol, fname)
        else:
            conv.ReadString(mol,string)
        return mol

    def writeMolecule(mol, fname=None, ftype=None):
        """ save a molecule with openbabel"""
        if ftype == None:
            n, ftype = getNameExt(fname)
            ftype = ftype.lower()
        conv = ob.OBConversion()
        conv.SetOutFormat(ftype)
        if not fname == None:
            conv.WriteFile(mol, fname)
        else:
            return conv.WriteString(mol)

    infile = sys.argv[1]
    name,ext = getNameExt(infile)
    mol = loadMolecule(infile)
    chiral = ChiralEnumerator3D(mol, max_iter=32, strict=False, refine=True)
    if not chiral.chiral:
        print("- molecule is not chiral")
        sys.exit(0)

    if not chiral.valid:
        print("- molecule is chiral but not VALID! (intra-ring chiral center)")
        sys.exit(1)


    unique=True
    enantiomers = chiral.get_enantiomers(unique=True)
    print("%d enantiomers found." % len(enantiomers))
    c = 1
    for e in enantiomers:
        writeMolecule(e, "%s_%d.%s" % (name, c, ext))
        c+=1
    sys.exit()



