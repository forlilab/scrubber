import itertools


from openbabel import openbabel as ob


# TODO XXX IMPORTANT BUG
# - scan also non-protonated amines, to flip the chirality!?
# https://en.wikipedia.org/wiki/Amine#Chirality

class ChiralEnumerator(object):
    """ chirality generator
        initialized with a molecule provides a generator
        that can be used to enumerate isomers of a molecule


        mol :   obmol
        chiralityType   : 'all' = generate all enantiomers for any chiral center
                          'undefined' = generate enantiomers for undefined chiral centers
                          'protomer' = generate enantiomers for protonated atoms, which could switch chirality because of proton flip
                                       NOTE: this option overrides 'undefined'
                          'off' = disable enantiomer generation
        maxIter         : maximum number of enantiomers to generate
                 : print status of generation
    """
    def __init__(self, mol, chiralityType='all', maxIter=50, debug=False):
        """ """
        self.debug = debug
        #self. = verbose
        self.mol = mol
        self.maxIter = maxIter
        self.chiralityType = chiralityType
        self.active = False
        self._stop = False
        self._iterations = 0
        self.chiralSet = []
        self.uniqueSMILES = []
        self._canonical = ob.OBConversion()
        self._canonical.SetInAndOutFormats('smi', 'can')
        self._protomers = [] # chiral atoms that are protonated (e.g., protonated tertiary amine)
        if self._is_chiral(mol) and (not self.chiralityType=='off'):
            self.chiral = []
            self.chiralDefined = []
            self.blacklist = self.findBlacklistCenters()
            self.facade = ob.OBStereoFacade(self.mol)
            # find which atoms are to be considered chiral
            for a in ob.OBMolAtomIter(self.mol):
                idx = a.GetIdx()
                # TODO this function here should be replaced with a custom function/class
                #      that checks not only the classical tetrahedral molecules but also 
                #      non-planar nitrogen, sulfur etc,
                if self.facade.HasTetrahedralStereo(idx) and not idx in self.blacklist:
                    self.chiral.append( idx )
                    ts = self.facade.GetTetrahedralStereo(idx)
                    config = ts.GetConfig()
                    if config.specified:
                        self.chiralDefined.append(idx)
                    else:
                        config.specified = True
                        ts.SetConfig(config)
                if self._isProtonatedChiral(a):
                    self._protomers.append(idx)
            if len(self.chiral) > 0:
                self.active = True
                msg = ('[ChiralEnumerator] chirality Implicit[%d] Explicit[%d]'
                        '' % (len(self.chiral), len(self.chiralDefined)))
                self.dprint(msg)
                self.chiralAccepted = set([])
                if self.chiralityType == 'all':
                    self.chiralAccepted = self.chiral
                    msg = ('[ChiralEnumerator] chiral set ALL [%d] atoms' % len(self.chiralAccepted))
                elif self.chiralityType == 'undefined':
                    self.chiralAccepted = set(self.chiral) - set(self.chiralDefined)
                    msg = ('[ChiralEnumerator] chiral set UNDEFINED [%d] atoms' % len(self.chiralAccepted))
                if self.chiralityType == 'protomer':
                    self.chiralAccepted = set(self.chiralAccepted) | set(self._protomers)
                    msg += ('\n[ChiralEnumerator] chiral set PROTOMER added [%d] atoms' % len(self.chiralAccepted))
                self.dprint(msg)
                self.chiralSet = []
                tmp = {}
                for idx in self.chiralAccepted:
                    ref1 = self.facade.GetTetrahedralStereo(idx).GetConfig().refs
                    # TODO use OBBuilder.Swap(OBMol, ref1[0], ref[1], ref[0], ref[2]
                    ref2 = ( ref1[1], ref1[0], ref1[2] )
                    self.chiralSet.append((ref1, ref2))
                self.chiralSet = list( itertools.product( *self.chiralSet) )
        self.generateUniqueEnantiomers()
        self.dprint("\n\n\n[ChiralEnumerator]: %d chiral centers accepted" % len(self.chiralSet) )

    def generateUniqueEnantiomers(self):
        """ generate the unique set of SMILES enantiomers"""
        for refs in self.chiralSet:
            msg = ('[enumchiral] |%s| enantiomer to be generated' % str(refs)) 
            self.dprint(msg)
            for idx, aIdx in enumerate(self.chiralAccepted):
                ts = self.facade.GetTetrahedralStereo(aIdx)
                config = ts.GetConfig()
                config.refs = ( refs[idx][0], refs[idx][1], refs[idx][2]) 
                x = ts.SetConfig(config)
            can = self._canonical.WriteString(self.mol)
            self.uniqueSMILES.append(can.strip())
            self.dprint("CANONICAL SMILE ENANTIOMER: %s"% can.strip())
        redundant = len(self.uniqueSMILES)
        self.uniqueSMILES = list(set(self.uniqueSMILES))
        unique = len(self.uniqueSMILES)
        self.dprint("\n\n\n[ChiralEnumerator]: %d unique enantiomers generated ( %d total)" % (unique, redundant))
        

    def dprint(self, *args, **kwargs):
        if self.debug:
            print(*args, **kwargs)

    def _is_chiral(self, mol):
        """ checks if molecule is chiral """
        for a in ob.OBMolAtomIter(mol):
            if a.IsChiral():
                return True
        return False

    def _isProtonatedChiral(self, atom):
        """ scan atom neighbors to see if there are protons that could flip"""
        #print "ATOM", atom.GetIdx(), atom.GetAtomicNum(), atom.IsChiral()
        if not atom.GetFormalCharge() > 0:
            return False
        for n in ob.OBAtomAtomIter(atom):
            if n.IsHydrogen():
                return True
        return False


    def __iter__(self):
        """" """
        return self

    def __next__(self):
        """ """
        if self._iterations >= self.maxIter:
            self._stop = True
        if self._stop:
            raise StopIteration
        if not self.active:
            # no (undefined) chiral centers are found in the molecule
            # so return the input
            self._stop = True
            return self.mol
        else:
            try:
                smi = self.uniqueSMILES.pop()
                msg = ('[enumchiral] |%s| enantiomer to be generated' % str(smi)) 
                #self.dprint(msg)
                self._iterations += 1
                enantiomer = ob.OBMol()
                self.dprint("CANONICAL SMILE SERVED: %s"% smi.strip())
                self._canonical.ReadString(enantiomer, smi)
                return enantiomer
            except IndexError:
                raise StopIteration

    def findBlacklistCenters(self):
        """ identify unwanted patterns (adamantane, for now) to avoid combinatorial 
            enumeration of its chiral centers
        """
        pattern = ['C1C3CC2CC(CC1C2)C3'] #, '[C@@H]', '[C@H]' ] 
        m = ob.OBSmartsPattern()
        for p in pattern:
            m.Init(p)
            found = m.Match(self.mol)
            if found:
                self.dprint("[findBlacklistCenters]: found unwanted! [%s]" % (pattern))
                return [ x for x in m.GetUMapList() ][0]
        return ()


if __name__ == '__main__':
    import sys, os
    name, ext = os.path.splitext(sys.argv[1])
    ext = ext[1:].lower()
    for pmol in pybel.readfile(ext, sys.argv[1]):
        mol = pmol.OBMol
        print("INPUT> ",pmol.__str__().strip())
        for chiral in ChiralEnumerator(mol, chiralityType='all',
                        maxIter = 50, debug = True):
            s = pybel.Molecule(chiral).__str__()
            print(s.strip())
            
