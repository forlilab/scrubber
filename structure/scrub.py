#!/usr/bin/env python


from openbabel import openbabel as ob
import numpy as np
import math, string
import sys
import multiprocessing as mp

OPTIMALAMIDEDIHE = math.pi
ZEROISH = 0.005

# suppress warning messages
ob.obErrorLog.SetOutputLevel(0)


###
### IMPORTANT: get example to suppress log from Dalke
#              self.lvl = ob.obErrorLog.GetOutputLevel()
#              ob.obErrorLog.SetOutputLevel(-1)

###

# XXX refactor this class to have multiple set methods
# XXX so objects are initialized once, performance will be better
class Scrub(object):
    """ This class is devoted to generate a minimum conformation
            for a ligand.

        [ http://www.merriam-webster.com/dictionary/scrub ]
        2scrub verb
                : to rub (something) hard with a rough object or
                  substance and often with soap in order to clean it
    """
    # TODO
    #   - add tautomeric search
    #   - add quick conformational search
    #   - include (as opposed to exclude)
    #   - add PAINS support
    #   - read from compressed formats
    def __init__(self, ff='MMFF94s', ff_extra='UFF',
        sdsteps=250, sdconv= '1e-6',
                cgsteps=100, cgconv='1e-8',
                sdsteps_extra=100, sdconv_extra= '1e-6',
                cgsteps_extra=100, cgconv_extra='1e-8',
                rotamer_conf = 10, rotamer_geom_steps=5,
                extra_3d_mini = False,
                flip_cis_amide = True,
                checkhydro = True,
                pH = 7.4, stripsalts=True,
                chargemodel='gasteiger',
                name = None,
                verbose=False,
                strict=False,
                inFormat=None,
                outFormat='mol2',
                debug=False,
                auto=False):
        self.debug = debug
        self.verbose = verbose
        self.stripsalts = stripsalts
        self.ready = True
        self.log = []
        self.mol = None
        self.mol_name = name
        self.ff = ff.lower()
        self.ff_extra = ff_extra.lower()
        if not pH == None:
            self.pH = self._double(pH)
        else:
            self.pH = None
        self.sdsteps = sdsteps
        self.sdconv = sdconv
        self.cgsteps = cgsteps
        self.cgconv= cgconv
        self.sdsteps_extra = sdsteps_extra
        self.sdconv_extra = sdconv_extra
        self.cgsteps_extra = cgsteps_extra
        self.cgconv_extra = cgconv_extra
        self.rotamer_conf = rotamer_conf
        self.rotamer_geom_steps = rotamer_geom_steps
        self.chargemodel = chargemodel
        self.strict = strict
        self.inFormat = inFormat
        self.outFormat = outFormat
        self.flip_cis_amide = flip_cis_amide
        self.check_hydro = checkhydro
        self._init_essentials() # XXX <- this should help when refactoring into reusable class

    def _init_essentials(self):
        """ initialize OB object used during the processing"""
        # 3D structure builder
        self.builder = ob.OBBuilder()
        # molecule writer to create string
        self.parser = ob.OBConversion()
        self.parser.SetOutFormat(self.outFormat)
        # check this? why none is even considered?
        if not self.inFormat == None:
            self.parser.SetInFormat(self.inFormat)

    def vprint(self, string, addnewline=True, flush=True, rewind=False):
        """ debugging printing"""
        if not self.verbose:
            return
        msg = ''
        if rewind:
            msg += '\r'
        msg += "VERBOSE: [SCRUB] %s" % string
        if not addnewline:
            print(msg, end=' ')
        else:
            print(msg)
        if flush:
            sys.stdout.flush()

    def _double(self, value):
        """ helper function to generate double floats"""
        return ob.double_array([value])[0]

    def process(self, mol_string, in_type, out_type, unique_atom_names=True):
        """ perform the operations requested"""
        # try:
        if 1:
            if not in_type == None:
                self.parser.SetInAndOutFormats(in_type,out_type)
            self.mol = ob.OBMol()
            self.parser.ReadString(self.mol, mol_string)
            # get a molecule name
            if self.mol_name == None:
                self.mol_name = self.mol.GetTitle()
            self._preprocess_mol()
            # perform minimization
            self.optimize()
            self.ready = True
            # assign unique atom names
            if unique_atom_names:
                self._set_unique_atom_names()
            out = self.parser.WriteString(self.mol)
            return ("%s\n" % out)
        else:
        # except:
            # problem somewhere!
            print("Error on mol: %s" % self.mol_name)
            return None


    def _preprocess_mol(self):
        """ initialize presettings"""
        # strip fragments, salts...
        if self.stripsalts:
            self.mol.StripSalts()
        # check if structure needs 3D coordinates
        self.checkStructure()
        if not self.ready:
            return
        # check hydrogens
        self.fix_ph()
        # calculate charges
        self.calc_charges()


    def checkStructure(self):
        """ check if the molecule has a 3D structure
            otherwise generate it
            if the structure is already 3D,
            check that hydrogens are properly added
        """
        # 3D + hydrogens
        dimension = self.mol.GetDimension()
        if dimension < 3:
            canonicize = False
            self.vprint("[checkStructure] Generating 3D coords...")
            if dimension < 2:
                canonicize=True
            self.make_3d(canonicize)
            self.vprint("[checkStructure] [ DONE ]")
        # hydrogens only
        else:
            if self.check_hydro:
                self._check_missing_h()

    def _check_missing_h(self):
        """ check that implicit and explicit valences
            match, otherwise add appropriate hydrogens
        """
        #self.mol.PerceiveBondOrders()
        #print "SKIPP CHECKISS"
        return
        for a in ob.OBMolAtomIter(self.mol):
            valence = a.GetValence()
            implicit = a.GetImplicitValence()
            if not (valence == implicit):
                self.mol.AddHydrogens(a)

    def fix_ph(self, pH=None):
        """ performs pH magics"""
        if pH == None:
            pH = self.pH
        if self.pH == None:
            return
        self.verbose=True
        # self.vprintCan('fix_ph: START')
        self.mol.DeleteHydrogens()
        # self.vprintCan('fix_ph: DEL_H   ')
        self.mol.UnsetFlag(ob.OB_PH_CORRECTED_MOL)
        # self.vprintCan('fix_ph: UNSETFLAG   ')
        for a in ob.OBMolAtomIter(self.mol):
            a.SetFormalCharge(0)
        # self.vprintCan('fix_ph: FORMALCHARGE')
        self.mol.SetAutomaticFormalCharge(True)
        # insanity fix for nitrogen for OB total disregard of chemistry...
        for a in ob.OBMolAtomIter(self.mol):
            if not a.GetAtomicNum() == 7:
                continue
            if (a.GetTotalDegree()==4) and (a.GetFormalCharge()==0):
                a.SetFormalCharge(1)
        # self.vprintCan('fix_ph: AUTOFORMAL  ')
        self.mol.AddHydrogens(False, True, pH)
        # self.vprintCan('fix_ph: ADDYDRO  ')
        self.verbose=False

    def calc_charges(self):
        """ calculate partial charges using selected charge model"""
        # XXX partial charges are queried to trigger something in the
        #     openbabel lib that fixes charge bug
        # TODO this function should be disabled if an option to use original
        #      charges is used
        #print "WARNING: destroying any input charges!"
        self.vprint("Charge assignment; model [ %s ]" % self.chargemodel )
        charger = ob.OBChargeModel.FindType(self.chargemodel)
        # print("WARNING: unset partial charge perceived disabled (OB 2->3) FIXME")
        for a in ob.OBMolAtomIter(self.mol):
            a.SetPartialCharge(0.0)
        self.mol.SetPartialChargesPerceived(False)
        self.mol.SetAutomaticPartialCharge(False)
        report = charger.ComputeCharges(self.mol)
        c = sum([a.GetPartialCharge() for a in ob.OBMolAtomIter(self.mol)])
        self.vprint("Total partial charges sum: %2.3f"  % c)
        if not report:
            msg = ('[calc_charges] WARNING: charge model is '
                     'missing parameters for some atoms')
            if self.verbose:
                print(msg)
            if self.strict:
                msg = ('MISSING CHARGES (--strict mode on) charge model '
                       'not available for molecule -> rejecting')
                print(msg)
                self.ready = False
                return
        partialCharges = charger.GetPartialCharges()
        if self.verbose:
            #self.vprint("[calc_charges] Charges assigned:")
            totalCharge = sum(partialCharges)
            buff = '[calc_charges] Charges assigned, total charge: %2.3f\n' % totalCharge
            indices = list(range(len(partialCharges)))
            elements = [ self.mol.GetAtom(i+1).GetAtomicNum() for i in indices]
            charges = ['%3.2f' % x for x in partialCharges]

            idxString = ",".join([str(x) for x in indices])
            elementString = ",".join([str(x) for x in elements])
            chargeString = ",".join(charges)

            buff += "  AtomIdx\t:" + idxString + "\n"
            buff += "  Element\t:" + elementString + "\n"
            buff += "  Charge\t:" + chargeString + "\n"
            self.vprint(buff)

    def _check_amide(self):
        """ check if amide is in cis conformation, and flips it if necessary"""
        # TODO CONVERT TO SMARTS AND MERGE WITH FORCE_TRANS_AMIDE?
        for b in ob.OBMolBondIter(self.mol):
            if b.IsPrimaryAmide() or b.IsSecondaryAmide():
                if not b.IsInRing():
                    self._force_trans_amide(b)

    def _force_trans_amide(self, bond):
        """ measure amide torsion angle and force it to be 180"""
        # TODO convert to SMARTS!!!!
        C, N, O, H = 6,7,8,1

        optimal = OPTIMALAMIDEDIHE
        tol = math.radians(36)
        r180 = math.radians(180)
        begin = bond.GetBeginAtom()
        end = bond.GetEndAtom()
        if begin.GetAtomicNum() == C:
            carbon = begin
            nitro = end
        else:
            carbon = end
            nitro = begin
        oxy = self._findAttached(carbon, O)
        if not oxy:
            print("WARNING! MISSING O FOR AMIDE!")
            return
        hydro = self._findAttached(nitro, H)
        if not hydro:
            #print "WARNING! MISSING H FOR AMIDE!"
            return
        dihe = self.calc_dihedral(oxy, carbon, nitro, hydro)
        if abs(dihe) <= r180:
            deviation = abs(dihe)
        else:
            deviation = abs(dihe) - r180
        if deviation <= tol:
            if self.verbose:
                print(('[VERBOSE] Amide TRANS conformation detected'
                        ' ( %3.2f deg, deviation: %3.3f): OK' % (dihe, deviation  )))
            return
        if self.verbose:
                print(('[VERBOSE] Amide CIS conformation '
                       'detected ( %2.2f deg): FIXING!' % math.degrees(dihe) ))
        angleFix = ob.double_array([-dihe])[0]
        self.mol.SetTorsion(oxy,carbon,nitro,hydro, angleFix)
        dihe2 = self.calc_dihedral(oxy, carbon, nitro, hydro)


    def _flipDihedral(self, atoms=[], bond=None, angle=None):
        """ rotate a bond by the required angle"""
        idx = [ x.GetIdx() for x in atoms]
        rotor = ob.OBRotor()
        rotor.SetBond(bond)
        rotor.SetDihedralAtoms(idx)
        return


    def calc_dihedral(self, a1, a2, a3, a4):
        """ given 4 OBAtom return the dihedral
            angle between them
        """
        v1 = self._obatom_to_vec(a1, a2)
        v2 = self._obatom_to_vec(a3, a2)
        v3 = self._obatom_to_vec(a3, a4)

        v4 = np.cross(v1, v2)
        v5 = np.cross(v2, v4)
        try:
            dihe = math.atan2(np.dot(v3,v4), np.dot(v3,v5) * math.sqrt(np.dot(v2,v2)))
        except ZeroDivisionError:
            dihe = 0.
        return dihe

    def _obatom_to_vec(self, a1, a2):
        """ return the vector between two atoms"""
        c1 = self._get_obatom_coord(a1)
        c2 = self._get_obatom_coord(a2)
        return c2-c1

    def _get_obatom_coord(self, a):
        """ given an atom return its coords as float"""
        coord = [ a.GetX(), a.GetY(), a.GetZ() ]
        return np.array(coord, 'f')

    def _findAttached(self, atom, anum):
        """ return the first neighbor of atom that is of type anum"""
        for n in ob.OBAtomAtomIter(atom):
            if n.GetAtomicNum() == anum:
                return n
        return None

    def set_ff_outlev(self):
        """ set reduced output for the forcefield"""
        pass


    # def _init_force_field(self):
    #     """ retrieve known forcefields"""
    #     # XXX change to INIT FORCEFIELD
    #     # forcefields =  ob.vectorString()
    #     # ob.OBPlugin.ListAsVector('forcefields', None, forcefields)
    #     # self.knownForcefields = [x.split()[0].lower() for x in forcefields]
    #     # if not self.ff in self.knownForcefields:
    #     #     print("ERROR: Unknown forcefield! [%s]" % self.ff)
    #     #     self.ready = False
    #     #     return False
    #     self.forcefield = ob.OBForceField.FindType(self.ff)
    #     # setup the outlevel
    #     self.set_ff_outlev()

    #     outlev = ob.OBFF_LOGLVL_NONE
    #     if self.verbose:
    #         outlev = ob.OBFF_LOGLVL_LOW
    #     self.forcefield.SetLogLevel(outlev)

    def set_force_field(self, ff_name):
        """ initialize the requested forcefield and set the lowest log levels"""
        ff= ob.OBForceField.FindType(ff_name)
        outlev = ob.OBFF_LOGLVL_NONE
        if self.verbose:
            outlev = ob.OBFF_LOGLVL_LOW
        ff.SetLogLevel(outlev)
        # XXX ff.SetLogFile() here!
        # setup molecule
        setup = ff.Setup(self.mol)
        if self.verbose:
            buff = '[setForceField] Atom types assigned:\n'
            idx, atype = [], []
            for a in ob.OBMolAtomIter(self.mol):
                idx.append(a.GetIdx())
                atype.append(a.GetType())
            idx = [str(x) for x in idx]
            buff += '  AtomIdx \t' + ','.join(idx) + '\n'
            buff += '  AtomType\t' +  ','.join(atype) + '\n'
            self.vprint(buff)

        if not setup:
            if self.strict:
                msg =  ('MISSING FF PARAMETERS (--strict mode on) forcefield parameters'
                       'not available for molecule -> rejecting')
                self.ready = False
                return False
            print("WARNING: [%s] Missing force field parameters for molecule [%s] !" % (ff_name, self.mol_name))
        return ff

    def vprintCan(self, msg='VPRINTCAN'):
        """ print canonical SMI string"""
        if not self.verbose:
            return
        canonicForm = ob.OBMol()
        obc = ob.OBConversion()
        obc.SetInAndOutFormats('smi', 'smi')
        can = obc.WriteString(self.mol)
        self.vprint('[VPRINTCAN] Generated canonical form: %s' %can.strip())
        print('[%s] %s' % (msg, can.strip()))


    def make_3d(self, canonicize=False, extra_3d_mini=False):
        """ make 3D on request
            by default, a canonical SMILES form is generated
            to avoid twisted/interconnected groups of
            atoms ? USEFUL?
        """
        self.vprint('[make_3d] Generating 3D structure')
        outcome = self.builder.CorrectStereoAtoms(self.mol)
        # print("TEST correct stereo atoms:", outcome)
        outcome = self.builder.CorrectStereoBonds(self.mol)
        # print("TEST correct stereo bonds:", outcome)
        outcome = self.builder.Build(self.mol)
        # print("TEST configuration 3D:", outcome)
        self.mol.SetDimension(3)
        self.mol.AddHydrogens(False, False)
        ff = self.set_force_field('mmff94s')
        setup_out = ff.Setup(self.mol)
        if not setup_out:
            if self.strict:
                print("ERROR SETTING UP THE MOLECULE! (strict, bailing out)")
                self.ready=False
                return
            # try the alternative force field
            ff = self.set_force_field('uff')
            setup_out = ff.Setup(self.mol)
            if not setup_out:
                print("ERROR SETTING UP THE MOLECULE! (non-strict, bailing out after alt_FF attempt)")

        ff.EnableCutOff(True);
        ff.SetVDWCutOff(10.0);
        ff.SetElectrostaticCutOff(20.0);
        # ff.SetUpdateFrequency(10) # slow?
        steps = 500
        if extra_3d_mini:
            steps = 1000
        ff.SetUpdateFrequency(steps) # fast?
        self.minimize(ff=ff, gradient='sd', steps=steps)
        #self.verbose = False
        if self.debug:
            obc = ob.OBConversion()
            obc.SetInAndOutFormats('mol2', 'mol2')
            can = obc.WriteFile(self.mol, 'builder.mol2')


    def optimize(self):
        """ generate optmized structure calling the basic and extra minimization steps"""
        # initialize the primary force field
        # self.set_force_field(self.ff)
        # if not self.ready:
        #     return
        # run the standard basic minimization (with the first FF, usually MMFF94), which can includes rotameric search
        if not self.sdsteps == self.cgsteps == 0:
            self._run_ff_optimization(self.ff,
                    self.sdsteps, self.sdconv,
                    self.cgsteps, self.cgconv,
                    (self.rotamer_conf>0))

        # check if amides are in cis-form
        if self.flip_cis_amide:
            self._check_amide()
        # run the extra minimization (UFF)
        # ff = self.set_force_field(self.ff_extra)
        # if not self.ready:
        #     return
        if not self.sdsteps_extra == self.cgsteps_extra == 0:
            self._run_ff_optimization(self.ff_extra,
                    self.sdsteps_extra, self.sdconv_extra,
                    self.cgsteps_extra, self.cgconv_extra,
                    False)

    def _run_ff_optimization(self, ff, sdsteps, sdconv, cgsteps, cgconv, rotamer_search=False):
        """ run a full minimization process, which includes SD, [rotameric search], and CG"""
        ff = self.set_force_field(self.ff)
        if not self.ready:
            return
        # XXX Add logging here
        #print "%s:" % self.name,
        # steepest descent
        if sdsteps:
            self.minimize(ff, 'sd', steps=sdsteps, convergence=self.sdconv)
        # rotameric search
        if rotamer_search: # and self.rotamer_geom_steps:
            self.rotamer_search(ff, num_conf=self.rotamer_conf, geom_steps=self.rotamer_geom_steps)
        # conjugated gradients
        if self.cgsteps:
            self.minimize(ff, 'cg', steps=self.cgsteps, convergence=self.cgconv)

    def rotamer_search(self, ff, num_conf=10, geom_steps=1):
        """ perform rotameric search"""
        self.vprint("[rotamer_search] numConf=%d, geomSteps=%d" % (num_conf, geom_steps))
        coord = self.mol.GetAtom(1).GetX()
        ff.WeightedRotorSearch(num_conf, geom_steps)
        ff.UpdateCoordinates(self.mol)

    def minimize(self, ff=None, gradient='sd', steps=None, convergence=1e-5):
        """ minimize energy by combining multiple minimizers """
        convergence = ob.double_array([float(convergence)])[0]
        if gradient == 'sd':
            gradientObj = ff.SteepestDescent
        else:
            gradientObj = ff.ConjugateGradients
        gradientObj(steps, convergence)
        ff.UpdateCoordinates(self.mol)

    def _set_unique_atom_names(self):
        """ create unique names for each atom in the molecule"""
        #etable = ob.OBElementTable()
        self.mol._atomNames = {}
        res_list = [x for x in ob.OBResidueIter(self.mol)]
        if len(res_list) == 0:
            res = self.mol.NewResidue()
            for a in ob.OBMolAtomIter(self.mol):
                res.AddAtom(a)
            res_list = [ res ]

        for res in res_list:
            for atom in ob.OBResidueAtomIter(res):
                name = res.GetAtomID(atom)
                if name == "":
                    num = atom.GetAtomicNum()
                    name = ob.GetSymbol(num)
                    #name = etable.GetSymbol( atom.GetAtomicNum() )
                # print "NAME", name,
                if not name in self.mol._atomNames:
                    self.mol._atomNames[name] = 1
                    #newName = name
                else:
                    self.mol._atomNames[name] += 1
                newName = "%s%d" % (name,  self.mol._atomNames[name])
                #print newName
                res.SetAtomID(atom, newName)

