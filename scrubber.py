from math import pi
from string import ascii_letters, digits
import sys
import argparse 
import os
import itertools
import psutil
#import multiprocessing as mp

from openbabel import openbabel as ob

from structure.scrub import Scrub
from structure.chiralenumerator import ChiralEnumerator
OPTIMALAMIDEDIHE = pi
ZEROISH = 0.005
# format used to generate string of molecule in the processing queue 
"""
ScrubMultiprocess and MolWriterMultiProcess work in separate threads:
 - ScrubMultiprocess receives data from queue_in, generates optimal coordinates and puts 
    (molString, fname) in queue_out. Multiple instances of ScrubMultiprocess are created
    in multiprocessing. Each time ScrubMultiprocess receives a poison pill, it will pass 
    it to the MolWriterMultiProcess, then will exit.

 - MolWriterMultiProcess reads 'molString' and writes it in fname; if MolWriterMultiProcess is initialized 
    with a fname in the constructor, then fhane specified in the queue is ignored 
    and all outputs are written in the constructor single fname (self.fp).
    When initialized, the number of ScrubMultiprocess instances is defined, so that when
    the same number of poison pills is received, MolWriterMultiProcess will flush the data
    in the output file, then will exit.
"""

### class ScrubMultiprocess(mp.Process, Scrub):
###     def __init__(self, queue_in, queue_out, nice=None, in_type='sdf', out_type='sdf', scrub_opts={}):
###         mp.Process.__init__(self)
###         self.in_type = in_type
###         self.out_type = out_type
###         if 'auto' in scrub_opts:
###             scrub_opts['auto'] = True
###         Scrub.__init__(self, **scrub_opts)
###         self.queue_in = queue_in
###         self.queue_out = queue_out
###         self.nice = nice
###         # total processed
###         self._c = 0
###         # failed
###         self._f = 0
### 
###     def run(self):
###         # nice from here: https://andrewbolster.info/2014/05/multiprocessing-niceness-in-python
###         # TODO add IOPRIO classes in both linux and windows
###         if not self.nice == None:
###             if not psutil.WINDOWS:
###                 p = psutil.Process(os.getpid())
###                 niceness = os.nice(0)
###                 os.nice(self.nice - niceness)
###             else:
###                 print("Warning: ignoring nice value on Windows")
###         while True:
###             # TODO
###             # this should become: name, string, fname
###             string, fname = self.queue_in.get()
###             if string == None:
###                 # poison pill!
###                 self.queue_in.task_done()
###                 self.queue_out.put((None, None))
###                 break
###             # this should become (name, string, self.in_type, self.out_type)
###             mol_string = self.process(string, self.in_type, self.out_type, unique_atom_names=True)
###             if not mol_string is None:
###                 self.queue_out.put((mol_string, fname))
###             else:
###                 self._f += 1
###             self.queue_in.task_done()
###             self._c +=1
###         return
### 
### class MolWriterMultiProcess(mp.Process):
###     def __init__(self, queue, threads, fname=None, ext=None):
###         """ molecular writer process, used to write 
###             molecules from the queue.
### 
###             by default, the queue provides tuples as
###             (mol_string, fname), which will be written in
###             the provided filename.
### 
###             If fname is provided at the constructor, 
###             a file pointer is created when the class is 
###             instanciated (self.fp). The 'fname' from the queue will
###             be ignored and all results will be appended in
###             self.fp. 
### 
###             if 'ext' is provided, it will be used as extension
###             for the output files from the queue
###         """
###         mp.Process.__init__(self)
###         self.queue = queue
###         self.threads = threads
###         self.fp = None
###         if not fname == None:
###             self.fp = open(fname, 'wb', 0)
###         self.ext = ext
###         self._c=0
### 
###     def run(self):
###         """ """
###         while True:
###             # TODO XXX handle error messages here 
###             string, fname = self.queue.get()
###             if string == None:
###                 # poison pills are coming!
###                 self.threads -= 1
###                 if self.threads == 0:
###                     break
###                 else:
###                     #print "Pill received...", self.threads
###                     pass
###             else:
###                 self._c+=1
###                 if not (self.fp == None):
###                     self.fp.write(str.encode(string))
###                 else:
###                     with open("%s.%s" % (fname, self.ext), 'wb', 0) as fp:
###                         fp.write(str.encode(string))
###                     #fp.close()
###         if not self.fp == None:
###             self.fp.close()
###         # DEBUG
###         #print "Structures written [%d]" % self._c
###         return 

class MolecularHound:
    """ a(n allegedly) efficient multi-SMARTS matcher"""
    def __init__(self, patterns):
        """ """
        self.patterns = patterns

    def sniff_mol(self, mol):
        """ """
        matchers = self._init_matchers()
        c = 1
        for m in matchers:
            if m.Match(mol):
                return c
            c += 1
        return False

    def _init_matchers(self):
        """ generator of OB matchers"""
        for p in self.patterns:
            m = ob.OBSmartsPattern()
            m.Init(p)
            yield m

class MiniMee(object):

    # XXX update documentation: usemolname + outfile creates prefix_name
    # XXX add a property filter?

    def __init__(self, debug=False):
        """ """
        self.debug = debug
        self.progname = 'Raccoon Scrub'
        self._set_opt_defaults()
        # init stuff
        self._init_ob_params()
        self._init_docs()
        self._init_options()
        self._init_opt_parser()
        self._parse_opt()
        self.start()

    def dprint(self, *args, **kwargs):
        if self.debug:
            print(*args, **kwargs)

    def _set_opt_defaults(self):
        """ set all predefined options"""
        # in/out
        self.default_out = 'mol2'
        self.default_informat = None
        self.default_outformat = None
        self.default_usemolname = None
        self.default_usemolnamesafe = None
        # slicing
        self.default_split = False
        self.default_only = False
        self.default_begin = 1
        self.default_end = None
        # mol processing
        self.default_pH = 7.4
        self.default_nopH = False
        self.default_charges = 'gasteiger'
        self.default_noflipamide = False
        self.default_nostripsalts = False
        self.default_nocheckhydro = False
        self.default_noprocess = False
        self.default_enumchiral = 'off' # all, undefined, protomer, off
        self.default_maxenumchiral = 10
        self.default_nomini = False
        self.default_noextra = False
        self.default_norotamer = True
        self.default_exclude = None
        self.default_fix_non_std_autodock = False
        # minimizer
        self.default_sdconv = 1e-5
        self.default_sdsteps = 300
        self.default_cgconv = 1e-6
        self.default_cgsteps = 300
        self.default_forcefield = 'mmff94s'
        self.default_heuristics = 'quick'
        self.default_heuristics_list = [ 'quick', 'accurate', 'long', 'extreme' ] 

        # extra steps
        #self.default_sdconv = 1e-9
        #self.default_sdsteps = 5000000
        #self.default_cgconv = 1e-10
        #self.default_cgsteps = 10000000

        # rotamer search
        self.default_rotamer_conf = 10
        self.default_rotamer_steps = 5

        # minimizer (post-processing)
        self.default_sdconv_extra = 1e-5
        self.default_sdsteps_extra = 1000
        self.default_cgconv_extra = 1e-6
        self.default_cgsteps_extra = 1000
        self.default_forcefield_extra = 'uff'

        self.default_strict = False # abort minimization if no parameters?
        # logging and info
        self.default_log = None
        self.default_verbose = False

        # multiprocessor
        self.default_multiproc_mode = 'all'
        #print("XXX URGENT: add psutil to the packages")
        self.default_multiproc_max = psutil.cpu_count()
        #self.default_multiproc_max = 1
        self.multiproc_max = self.default_multiproc_max
        # NICE XXX TODO 
        # psutil
        # source: https://stackoverflow.com/questions/1023038/change-process-priority-in-python-cross-platform
        #self.default_multiproc_nice = psutil.

        self._allowedOutFormats = ['mol2', 'pdb', 'pdbqt', 'sdf', 'ent', 
            'gpr', 'mdl', 'mol', 'mop', 'sd']
        self._validFnameChars = "-_.()%s%s" % (ascii_letters, digits)


    def _init_opt_parser(self):
        """ initialize the parser"""
        # http://stackoverflow.com/questions/3853722/python-argparse-
        # how-to-insert-newline-the-help-text
        class SmartFormatter(argparse.HelpFormatter):
            def _split_lines(self, text, width):
                # this is the RawTextHelpFormatter._split_lines
                if text.startswith('R|'):
                    return text[2:].splitlines()  
                return argparse.HelpFormatter._split_lines(self, text, width)
        self.parser = argparse.ArgumentParser(description=self.desc,
            usage=self.usage, epilog=self.epilog,
            formatter_class = argparse.RawDescriptionHelpFormatter)
        for group_name, group_desc in self._optparser_groups_order:
            group = self.parser.add_argument_group(group_name, group_desc)
            for name in self._optparser_args_order[group_name]:
                opts = self._args_dict[name]
                group.add_argument(name, **opts)
        # help advanced checking
        if '--help_advanced' in sys.argv:
            self.show_advanced_help()
            sys.exit(0)
        self.args = self.parser.parse_args()
        # activate verbose as soon as possible, if requested
        self.verbose = self.args.verbose
        self._data = []
        print("CKECKING", self.args)
        # sys.exit(0)

    def print_msg(self, msg, mtype = 'e', q=None):
        """ utility to report messages """
        if mtype == 'e':
            buff = '*** ERROR *** '
        elif mtype == 'w':
            buff = '*** WARNING *** '
        qs = self._get_string(q)
        print(buff + msg + qs) 
        if mtype == 'e':
            sys.exit(1)
            
    def _get_string(self,q):
        """ Base/64"""
        return ""
        return self._data[q]

    def _parse_opt(self):
        """ perform the parsing of options"""
        self._parse_opt_logging()
        self._parse_opt_mol_names()
        self._parse_opt_slicing()
        self._parse_opt_input_output()
        #self._opt_update_output_fname()
        self._parse_opt_minimizer()
        self._parse_opt_mol_processing()
        self._parse_opt_exclude()
        self._parse_opt_multiproc()
        
    def _parse_opt_logging(self):
        """ set verbose mode and logging"""
        #self.verbose = self.args.verbose
        if self.verbose:
            print("==========VERBOSE MODE================\nOPTS")
            print("Called with following settings:")
            for group_name, group_desc in self._optparser_groups_order:
                for i in self._optparser_args_order[group_name]:
                    print(i, "=>", getattr(self.args,i[2:]))
            print("-" * 30)
            print("Supported forcefields:", end=' ')
            print(self.ff_list_names)
            print("Supported charge models:", end=' ')
            print(self.charge_models_list_names)
            print("======================")
        # XXX INIT LOG FILE POINTER HERE
        # logging
        log = self.args.log
        #if log == None:
        #    self.log = self.outname + '.log'

    def _parse_opt_mol_names(self):
        """ check what's the policy for molecule names
        
            if using molecule names is requested, split is 
            switched on, and the dictionary to keep track
            of molecules is created
        """
        self.usemolname =  self.args.usemolname
        self.safenames =  self.args.usemolnamesafe
        if self.safenames:
            self.usemolname = True
        self.mol_names_dict = {}
        if self.usemolname:
            self.args.split = True
            self.dprint("[usemolname] output file names will use molecule name if found")
            if self.safenames:
                self.dprint("[safenames] molecule names will be sanitized")
            self.dprint("[split] split is forced to be ON by usemolname")

    def _parse_opt_input_output(self):
        """ parse input and output file management"""
        self._parse_opt_input()
        self._parse_opt_multi()
        self._parse_opt_output()

    def _parse_opt_input(self):
        """ manage the input file name options"""
        # infile
        self.infile = self.args.infile
        self.inname, self.inext = self.get_name_ext(self.infile)
        self.informat = self.args.informat
        self.intype = None

        if (self.inext == '') and (self.informat == None):
            msg = ('The input file has no extension and no input '
                   'format has been specified.\nUse \'--informat\' to '
                   'specify a supported format type\.\n'
                   '"One goes to the right, the other to the left; ' 
                   'both are wrong, but in different directions." '
                   '-- Quintus Horatius Flaccus, Satires (II,3,5)'
                   )
            print(msg)
            sys.exit(1)

        if self.informat:
            self.intype = self.informat
            self.dprint("[_parse_opt_input] input file type defined as: %s" % self.informat)
        else:
            self.intype = self.inext
            msg = ("[_parse_opt_input] input file name \"%s\" (guessed type: %s)") 
            self.dprint( msg %  (self.infile, self.inext))
            self.informat = self.intype

    def _parse_opt_multi(self):
        """ perform the multi-structure check"""
        # check for multi-structure 
        try:
            self.multi = self.is_multi(self.infile)
            if self.multi == True:
                msg = ("[_parse_opt_multi] Input file is a multi-structure file")
                self.dprint(msg)
            elif self.multi == False:
                self.args.split = False
                msg = ('(_parse_opt_multi) Input file is not a multi-structure '
                       'file (forcing --split off)')
                self.dprint(msg)
            elif self.multi == None:
                print("*** WARNING ***\nand now?")
                    
        except IOError as e:
            msg =  'An error occurred while reading file [%s]:' % self.infile
            msg += 'I/O error (%d) : %s' % (e.errno, e.strerror)
            self.print_msg(msg)
    
    def _parse_opt_output(self):
        """ parse output information to determine output filename
            and format.

        Priorities in the assignments are:

        Filename:
            outname and usemolname specified ->  %outname_%usemolname
            outname only specified           -> %outname_(auto_progressive naming)
            no outname specified             -> %inname 

        File type:
            outformat specified -> %outformat
            outfile specified   -> %outfile_extension
            no outfile specified-> %informat (if valid); default otherwise

        File extension:
            (non-empty) -> use specified file extension
            (empty) + outformat -> keep empty
            (empty)             -> use informat if valid; default otherwise
        """
        # check if the output file is actually a directory
        self.outdir = None
        self.outfile = self.args.outfile
        if os.path.isdir(self.outfile):
            self.outdir = self.outfile
            self.outfile = ''
        self.outformat = self.args.outformat.lower()

        self.outtype = ''
        self.outname, self.outext = self.get_name_ext(self.outfile)

        #### file name
        if not len(self.outname) and not self.usemolname:
                # "%inname_processed"
                self.outname = '%s_processed' % self.inname 
        ### file type
        if self.outformat:
            # requested format
            if self.outformat in self._allowedOutFormats:
                self.outtype = self.outformat
                #self.outext = self.outformat
            else:
                msg = ('Specified format [%s] is not an acceptable output format' % self.outformat)
                self.print_msg(msg)
        else:
            if len(self.outext):
                # format from output file extension
                if self.outext.lower() in self._allowedOutFormats:
                    self.outtype = self.outext.lower()
                # format from input type (if acceptable)
            elif self.informat in self._allowedOutFormats:
                    self.outtype = self.informat
            else:
                # default format type
                self.dprint('[_parse_opt_output] output set to default (%s)' % self.default_out)
                self.outtype = self.default_out

        # file extension
        if len(self.outext) == 0:
            if self.outformat:
                self.outext = self.outformat
                pass
            elif self.informat in self._allowedOutFormats:
                self.outext = self.informat
            else:
                self.outext = self.default_out
        self._opt_update_output_fname()


    def _parse_opt_slicing(self):
        """ check that all splicing options are well formed"""
        self.single = single = self.args.single
        begin = self.args.begin
        end = self.args.end
        self.split = self.args.split
        # single
        if self.single:
            if begin or end:
                msg =  "Use either --single or --begin/--end"
                self.print_msg(msg, q=1)
            self.begin = single
            self.end = single
            msg = "[_parse_opt_slicing] processing only molecule #%d (%s,%s)"
            self.dprint(msg % (single,self.begin,self.end))
        # begin/end
        else:
            if begin or end:
                self.dprint("[_parse_opt_slicing] begin/end detected (%s,%s)" % (begin, end))
            self.begin = ( begin or self.default_begin)
            self.end = end
            if self.end and (self.begin > self.end):
                msg = ('The value specified in --begin should be '
                         'smaller than the one used in --end')
                self.print_msg(msg, q=2)
            self.dprint("[_parse_opt_slicing] processing molecule range:(%s,%s)" % (begin,end))

    
    def _opt_update_output_fname(self):
        """ check various option combinations to provide
            the appropriate file name
        """
        self.outfile = self.outname
        # updating naming scheme (redundant 'usemolname', but good reminder)
        if self.split or self.single or self.usemolname: 
            if len(self.outfile):
                self.outfile += "_%s"
            else:
                self.outfile = "%s"
        # create output filename
        if len(self.outext):
            self.outfile = "%s.%s" % (self.outfile, self.outext)
        if not self.outdir == None:
            self.outfile = self.outdir + os.path.sep + self.outfile
        self.dprint("[_opt_update_output_fname] the output file is \"%s\"" % self.outfile)


    def _parse_opt_minimizer(self):
        """ minimizer parameters (forcefield, charges, minimizer steps)"""
        # get options 
        self.automini = self.args.automini
        self.nomini = self.args.nomini
        self.noextra = self.args.noextra
        self.norotamer = self.args.norotamer
        self.rotamer_conf = self.args.rotamer_conf
        self.rotamer_steps = self.args.rotamer_steps
        # minimizer
        self.sdsteps = self.args.sdsteps
        self.sdconv = self.args.sdconv
        self.cgsteps = self.args.cgsteps
        self.cgconv = self.args.cgconv
        self.forcefield = self.args.forcefield.lower()
        # minimizer extra
        self.sdsteps_extra = self.args.sdsteps_extra
        self.sdconv_extra = self.args.sdconv_extra
        self.cgsteps_extra = self.args.cgsteps_extra
        self.cgconv_extra = self.args.cgconv_extra
        self.forcefield_extra = self.args.forcefield_extra.lower()
        # validate force fields
        for ff in [self.forcefield, self.forcefield_extra]:
            if (not ff in self.ff_list):
                msg = ('Specified force field [%s] is not recognized. '
                   'Supported force fields  are: %s.') % (ff, self.ff_list_str)
                self.print_msg(msg)
                print(ff, self.ff_list)
                sys.exit(1)
        # charges model 
        self.chargemodel = self.args.chargemodel.lower()
        # charges (check them here, they're going to be used always)
        if not self.chargemodel in self.charge_models_list:
            msg = ('Specified charge model [%s] is not recognized. '
                   'Supported charge models are: %s.') % (self.chargemodel,
                   self.charge_models_list_str)
            self.print_msg(msg)
            sys.exit(1)
        # check parameters consistency
        self._check_mini_parms()
        self._check_mini_parms_extra()
        self._check_rotamer_parms()

    def _check_rotamer_parms(self):
        """ check rotamer search parameters"""
        if self.norotamer:
            #anyparm = (self.rotamer_conf > 0 or self.rotamer_steps > 0)
            anyparm = ((self.rotamer_conf is not None and self.rotamer_conf > 0) or (self.rotamer_steps is not None and self.rotamer_steps > 0))
            if anyparm:
                msg = ('Rotameric search disabled (--norotamer) but rotamer '
                      'parameters defined (--rotamer_conf/--rotamer_steps) requested')
                self.print_msg(msg, q=5)
            self.rotamer_conf = 0
            self.rotamer_steps = 0
        else:
            if self.rotamer_conf == None:
                self.rotamer_conf = self.default_rotamer_conf
            if self.rotamer_steps == None:
                self.rotamer_steps = self.default_rotamer_steps
            if (self.rotamer_conf == 0) and (self.rotamer_steps == 0):
                self.dprint("[parseRotamer] rotameric search disabled")
            if (self.rotamer_conf < 0) or (self.rotamer_steps < 0):
                msg = ('Number of rotamer conformers or steps negative!')
                self.print_msg(msg, q=5)

    def _check_mini_parms(self):
        # disable minimization
        if self.nomini:
            self.automini = None
            anysteps = ((self.sdsteps is not None and self.sdsteps > 0) or (self.cgsteps is not None and self.cgsteps > 0))
            if anysteps:
                msg = ('Minimization disabled (--nomini) but minimizer '
                      'steps defined (--sdsteps/--cgsteps) requested')
                self.print_msg(msg, q=5)
                # printmsg will exit by now
            self.sdsteps = 0
            self.cgsteps = 0
            self.dprint('[_parse_opt_minimizer] minimization disabled (--nomini)')
            return
        else:
            if self.sdsteps == None:
                self.sdsteps = self.default_sdsteps
            if self.sdconv == None:
                self.sdconv = self.default_sdconv
            if self.cgsteps == None:
                self.cgsteps = self.default_cgsteps
            if self.cgconv == None:
                self.cgconv = self.default_cgconv           
        # steepest descent
        if self.sdsteps == 0:
            self.dprint("[_parse_opt_minimizer] SD disabled")
        if self.sdsteps < 0:
            msg = ('SD steps are negative! (seriously? A negative number of steps '
                    'for the minimizer? Try again)')
            self.printmsg(msg)
        # conjugate gradient
        if self.cgsteps == 0:
            self.dprint("[_parse_opt_minimizer] CG disabled (0 steps)")
        if self.cgsteps < 0:
            msg = ('CG steps are negative! (seriously? A negative number of steps '
                    'for the minimizer? Try again)')
            self.printmsg(msg)
        if self.sdsteps == self.cgsteps == 0:
            msg = ('(SD and CG steps have been set to 0: '
                    'use --nomini to disable minimizer)')
            print(msg)
        # forcefield
        self.dprint("[_parse_opt_minimizer] forcefield set to [%s]"% self.forcefield.upper())
        if not self.forcefield in [x.lower() for x in list(self.ff_list.keys())]:
            msg = ('Specified forcefield [%s] is not recognized. '
                   'Supported forcefields are: %s.') % (self.forcefield,
                   self.ff_list_str)
            self.print_msg(msg)

    def _check_mini_parms_extra(self):
        # disable minimization
        if self.noextra:
            #anysteps = (self.sdsteps_extra > 0 or self.cgsteps_extra > 0)
            anysteps = ((self.sdsteps_extra is not None and self.sdsteps_extra > 0) or (self.cgsteps_extra is not None and self.cgsteps_extra > 0))
            if anysteps:
                msg = ('Extra minimization disabled (--noextra) but minimizer '
                      'steps defined (--sdsteps_extra/--cgsteps_extra) requested')
                self.print_msg(msg, q=5)
                # printmsg will exit by now
            self.sdsteps_extra = 0
            self.cgsteps_extra = 0
            self.dprint('[_parse_opt_minimizer] extra minimization disabled (--noextra)')
            return
        else:
            if self.sdsteps_extra == None:
                self.sdsteps_extra = self.default_sdsteps_extra
            if self.sdconv_extra == None:
                self.sdconv_extra = self.default_sdconv_extra
            if self.cgsteps_extra == None:
                self.cgsteps_extra = self.default_cgsteps_extra
            if self.cgconv_extra == None:
                self.cgconv_extra = self.default_cgconv_extra
        # steepest descent
        if self.sdsteps_extra == 0:
            self.dprint("[_parse_opt_minimizer] extra minimization SD disabled (0 steps)")
        if self.sdsteps_extra < 0:
            msg = ('SD steps are negative! (seriously? A negative number of steps '
                    'for the minimizer? Try again)')
            self.printmsg(msg)
        # conjugate gradient
        if self.cgsteps_extra == 0:
            self.dprint("[_parse_opt_minimizer] extra minimization CG disabled (0 steps)")
        if self.cgsteps_extra < 0:
            msg = ('CG steps are negative! (seriously? A negative number of steps '
                    'for the minimizer? Try again)')
            self.printmsg(msg)
        if self.sdsteps_extra == self.cgsteps_extra == 0:
            msg = ('(SD and CG steps have been set to 0: '
                    'use --nomini to disable minimizer)')
            print(msg)
        # forcefield
        self.dprint("[_parse_opt_minimizer] forcefield set to [%s]"% self.forcefield_extra.upper())
        if not self.forcefield_extra in [x.lower() for x in list(self.ff_list.keys())]:
            msg = ('Specified forcefield [%s] is not recognized. '
                   'Supported forcefields are: %s.') % (self.forcefield_extra,
                   self.ff_list_str)
            self.print_msg(msg)
    
    def _parse_opt_mol_processing(self):
        """ parse molecule processing options"""
        # mol processing 
        self.pH = self.args.pH
        self.nopH = self.args.nopH
        self.flipamide = not self.args.noflipamide
        self.stripsalts = not self.args.nostripsalts
        self.checkhydro = not self.args.nocheckhydro
        # self.fix_non_std_autodock = self.args.fix_non_std_autodock
        if self.args.noprocess:
            self.nopH = self.noflipamide = True
            self.stripsalts = self.checkhydro = self.flipamide = False
        # check that --nopH and --pH are not used at the same time
        if (not self.nopH == None) and (not self.pH == None):
            msg = ('--pH value for protonation specified together with '
                   '--nopH option: only one can be used at the time')
            self.print_msg(msg, q=6)
        if (self.pH == None) and (self.nopH == None):
            self.pH = self.default_pH
            self.nopH = self.default_nopH
            self.dprint ("[_parse_opt_mol_processing) pH and protonation settings set to default")
        if self.nopH:
            self.pH = None
            self.dprint ("[_parse_opt_mol_processing) protonation disabled")

        if not self.args.enumchiral.lower() in ['all', 'undefined', 'protomer', 'off']:
            msg = ('incorrect value for --enumchiral option '
                    '[%s], allowed:' % self.args.enumchiral.lower())
                      #  "|".join(['all', 'undefined', 'None']))
            self.print_msg(msg)

    def _parse_opt_exclude(self):
        """ parse and check that SMARTS patterns are valid"""
        # parse options
        self.exclude = self.args.exclude
        self.excludefromfile = self.args.excludefromfile
        self.blacklist = []
        # exclude
        if len(self.exclude):
            msg = "[_parse_opt_exclude] %d SMARTS pattern(s) specified with --exclude:\n" % len(self.exclude)
            patterns = ["   |%s|"%x for x in self.exclude]
            msg += "\n".join(patterns)
            self.dprint(msg)
            result = self.validate_smarts(self.exclude)
            if not result == True:
                msg = ('Pattern defined in --exclude (%d) is not a valid '
                        'SMARTS pattern |%s|') % (result[0], result[1])
                self.print_msg(msg)
            self.blacklist += self.exclude
        # exclude from file
        if self.excludefromfile:
            msg = None
            patterns = self.get_lines_from_file(self.excludefromfile)
            if len(patterns) == 0:
                msg = '[%s] is an empty file'
            result = self.validate_smarts(patterns)
            if not result == True:
                msg = ('Pattern at line %d in %s is not a valid SMARTS pattern '
                       '|%s|') % (result[0], self.excludefromfile, result[1])
            if not msg == None:
                self.print_msg(msg)
            self.blacklist += patterns
        self._init_mol_hound(self.blacklist)

    def _parse_opt_multiproc(self):
        """ parse options for multiprocessor/multicore and nice values """
        self.multiproc_max = self.args.multicore
        if self.multiproc_max > psutil.cpu_count():
            print("WARNING: the number of detected cores is smaller than the requested number of multiprocessing jobs.")
        self.nice = self.args.nice


    def _init_ob_params(self):
        """ initialize choices depending on OB features"""
        # forcefields list
        ff_list =  ob.vectorString()
        ob.OBPlugin.ListAsVector('forcefields', None, ff_list)
        self.ff_list = {}
        self.ff_list_str = []
        self.ff_list_names = []
        for k, v in [ x.split(None, 1) for x in ff_list ]:
            self.ff_list[k.lower()] = v[:-1]
            self.ff_list_str.append( "\"%s\" (%s)" % (k, v[:-1]))
            self.ff_list_names.append('"%s"' % k.lower())
        self.ff_list_str = ", ".join(self.ff_list_str)
        self.ff_list_names = ', '.join(self.ff_list_names)

        # charges list
        charge_models_list = ob.vectorString()
        ob.OBPlugin.ListAsVector('charges', None, charge_models_list)
        self.charge_models_list = {}
        self.charge_models_list_str = []
        self.charge_models_list_names = []
        for k, v in [x.split(None, 1) for x in charge_models_list]:
            self.charge_models_list[k] = v #[:-1]
            self.charge_models_list_str.append( "\"%s\" (%s)" % (k, v)) # [:-1]))
            self.charge_models_list_names.append('"%s"' % k)
        self.charge_models_list_str = ", ".join(self.charge_models_list_str)
        self.charge_models_list_names = ', '.join(self.charge_models_list_names)

    def is_multi(self, fname, firstonly=True):
        """ check if there is more than one structure
            NOTE (2014.9.15): this function is never used 
            for actual counting, that could turn out to be
            very expensing when processing large data files

            The function should be re-written to check 
            that empty files are recognized!
        """
        string_patterns = { 'smi': '\n', 'can': '\n', 
                'sdf' : '$$$$',
              'mol2': '@<TRIPOS>MOLECULE',
              'pdb': 'MODEL',
              }

        binary_patterns = {'cdx' : None,
                        }
        #name, ext = self.get_name_ext(fname)
        ext = self.informat
        if ext in string_patterns:
            patt = string_patterns.get(ext, None)
            if patt == None:
                print(('Warning! The [%s] format is a not '
                       'supported multi-structure file '
                       '(assuming multi-structure)' % ext))
                return True
            c = 0
            with open(fname, 'r') as f:
                for l in f:
                    if l.find(patt) >= 0:
                        c +=1
                        if c==2 and firstonly:
                            f.close()
                            return True
                # full count
                f.close()
                return c
            f.close()
        elif ext in binary_patterns:
            self.dprint("[is_multi] found CDX/binary")
            patt = binary_patterns.get(ext, None)
            # XXX this is an approximation
            # by definition, a CDX is considered a
            # multi-structure file
            return True
            
        return False

    def get_name_ext(self, fname):
        """ extract name and extension from the input file"""
        #name, ext = fname.rsplit('.', 1)
        name, ext = os.path.splitext(fname)
        return name, ext[1:].lower()

    def get_lines_from_file(self, fname):
        """ parse stripped lines from file"""
        try:
            fp = open(fname, 'r')
            lines = fp.readlines()
            fp.close()
            lines = [ x.strip() for x in lines if x.strip() ]
            return lines
        except IOError as e:
            msg = ('Error while reading [%s]: I/O error (%d) : %s') % (fname, e.errno, e.strerror)
            self.print_msg(msg)


    def validate_smarts(self, plist):
        """ check that every smarts pattern is """
        sys.stdout.flush()
        for i, p in enumerate(plist):
            self.dprint('[validate_smarts]: checking pattern[%d]:"%s"' % (i,p))
            #print "[validate_smarts]: checking %s\n" %p * 100
            validator = ob.OBSmartsPattern()
            validator.Init(p)
            if not validator.IsValid():
                self.dprint('[validate_smarts]: invalid pattern found')
                
                return i,p
        self.dprint('[validate_smarts]: all patterns are valid')
        return True

    def _init_mol_hound(self, blacklist):
        """ initialize mol_hound SMARTS matcher"""
        self.mol_hound = MolecularHound(blacklist)

    def filter_smarts(self, mol):
        """ basic SMARTS pattern matcher """
        if len(self.blacklist) == 0:
            return True
        result = self.mol_hound.sniff_mol(mol)
        if result == False:
            self.dprint("[filter_smarts] all filters passed")
            return True
        self.dprint("[filter_smarts] [%s] rejected by SMARTS pattern #%d" % (self.current_mol_name, result))
        self._rejected+=1
        return False

    def _init_docs(self):
        """ initialize documentation"""
        # args/opt initialization
        self.desc = ('Minimize molecular structure(s) in input file.\n\n'
                'By default, steepest descent minimization is performed, \n'
                'followed by conjugate gradients minimization.\nWhen necessary, 3D structure'
                ' is generated automatically using MMFF94.\n'
                'All (reasonable) OpenBabel input/output file formats are supported.\n\n'
                '      _(\-/)_\n'
                '     {(#b^d#)} \n'
                '     `-.(Y).-` \n') # JGS

        #usage = 'This. Is. Sparta!'
        #self.usage = 'This. Is. Sparta!'
        self.usage = None
        self.epilog = None

    def show_advanced_help(self):
        """ do what it says..."""
        cc = ''
        for k in sorted(self.charge_models_list.keys()):
            v = self.charge_models_list[k]
            cc += '        {0:9s} : {1:30s}\n'.format(k,v)
            #cc += '        %s\t: %s\n' % (k, v)

        ff = ''
        for k in sorted(self.ff_list.keys()):
            v = self.ff_list[k]
            ff += '        {0:9s} : {1:30s}\n'.format(k,v)
        self.parser.print_help()
        self._advanced_help_text = """
    
An example usage is the generation of ligand structures to be used for docking, 
with AutoDock, starting from SMILES files or ChemDraw schemes:

  %(progname)s --infile ligand_library.smi --outfile ligand_library.mol2
   
This will generate automatically 3D coordinates and save them in a multi-structure
Mol2 file.

Another usage could be to split a Mol2 file that already has 3D coordinates and 
generate single ligands suitable for docking with AutoDock (i.e., minimized 
conformation, protonation state, etc...):

  %(progname)s --infile nci_p0.0.mol2 --outfile dockings.pdbqt \
       --usemolname --end 10
   
This will convert the first 10 molecules found in the file Mol2, 
saving ZINC79106062.pdbqt, ZINC33512584.pdbqt, etc...

   

   %(progname)s is a general tool...
   PROTOCOL   When a molecule is read from a file, the following protocol is
   applied by default:

                                    [ START ]
                                        |
                                        |
                                  SMARTS filter
                                        |
                                        |
                               Strip salts/fragments
                                        |
                ________________________|_______________________
               |                        |                       |
               |                        |                       |
          1D formats               2D formats                3D formats
        (SMI, CDX,...)           (MOL2, SDF...)           (MOL2, SDF, PDB...)
               |                        |                       |
               |                        |                       |
         canonical SMI                  |                       |
               |                        |                       |
               |                        |                       |
           3D coordinates generation (mmff94)             check hydrogens
                         |                                      |
                         |______________________________________|
                                        |
                                        |
                      Check primary amide trans-conformation
                                        |
                                        |
                                 protonate for pH
                               (default: %(pH)1.1f)
                                        |
                                        |
                            Calculate partial charges
                         (default charge model: %(charge)s)
                                        |
                                        |
                           Steepest descent minimization
                          (default forcefield: %(forcefield)s)
                                        |
                                        |
                          Conjugate gradient minimization
                          (default forcefield: %(forcefield)s)
                                        |
                       Quick weighted conformational search 
                                        |
                          Post-optimization using UFF ???
                                        |
                                        |
                                 Save output file
                              (default format: %(outformat)s)
                                        |
                                        |
                                     [ END ]            


   In the default non-verbose mode, processing output is shown is presented as following:

                      .- Steepest-descent minimization
                      |
           molname:SD[x],CG[x]
                            |
                            '- Conjugate gradient minimization

   where 'x' can be 0 ( = converged) or  1 ( = max.iterations reached)
   Ideally, all minimizations should converge to guarantee that a true minimum
   conformation is reached. This is particularly important for docking.

   3D coordinates are generated automatically when 1D or 2D structures are found, then
   a series of structural clean-ups are performed.   

   SUPPORTED FORCEFIELDS
   The following forcefields are supported:
%(ff)s   
    SUPPORTED PARTIAL CHARGE MODELS
    The following partial charge models are available:
%(cc)s

EXAMPLES

   %(scriptname)s --infile ligand_library.smi --outfile ligand_library.mol2
       #  Generate 3D structures for all molecules found in the input file, adding hydrogens
       #  protonating groups at pH 7.4, minimizing coordinates and performing all default 
       #  structural cleanups, and saving the result in a single multi-structure file
        
   %(scriptname)s --infile nci_p0.0.mol2 --outfile docking.pdbqt --noprocess --nomini --end 100 --usemolname
       #  Convert first 100 molecules in the input file to separate PDBQT files ready 
       #  to be docked and using molecule names contained in the input file; no minimization is
       #  performed; all processing steps (salt stripping, hydrogens/protonation, trans-amide check)
       #  are disabled; ligands will be saved as: 
       #       docking_ZINC04783481.pdbqt
       #       docking_ZINC04786808.pdbqt
       #       docking_ZINC04786811.pdbqt
       #       ...
              
====================================================	 
%(progname)s (C)2021 ForliLab, Scripps Research 
        """ % { 'progname' :self.progname,
        'pH': self.default_pH, 'charge':self.default_charges, 'forcefield':self.default_forcefield,
        'outformat': self.default_out, 'cc': cc, 'ff':ff, 'scriptname': os.path.basename(sys.argv[0]),
        }  
        print(self._advanced_help_text)


    def _init_options(self):
        """ initialize command line options and information to properly
            pack the parser respecting proper order
        """
        print("\n\n ADD SHORT AND LONG VERSION! for INITIOPTIOHNS!\n\n")

        self._args_dict = { '--infile' : { 'help':'input file; '
                                  'format is guessed from file extension; ',
                                 #'this is the only required argument', 
                                 'action': 'store',
                                 'metavar' : 'INPUT_FILE[.EXT]', 'required' : True, 'type': str},

        '--single': {  'help':'process only Nth structure in the file' ,
                        'action' :  'store', 'metavar': 'STRUCTURE_NUMBER', 'type' : int},

        '--begin' : { 'help':'start processing from Nth molecule', # (default: first molecule)', 
                       'metavar' : 'FIRST_STRUCTURE_NUMBER', 'type': int, 'action':'store', 
                       'default': None},

        '--end' : { 'help': 'stop processing at Nth molecule ', # (default: last molecule)', 
                    'action' : 'store', 'metavar' : 'LAST_STRUCTURE_NUMBER', 'type': int, 
                    'default': None},

        '--byname' : { 'help' : 'process only molecules with specified name (note: cAse sEnSiTivE,'
                    ' spaces are not allowed; ...; )',
                    'metavar' : 'MOL_NAME', 'action' : 'store', 'type' : str },

        '--outfile' :  { 'help':('set output file; format is guessed from file extension; '
                                   'by default, output filename is \'INFILE_mini\' and same '
                                   'type as input; if input file format is not a 3D format '
                                   '(e.g. SMI) default %s is used; if an existing  directory '
                                   'is specified, output will be written in there') % self.default_out.upper(),
                        'metavar': 'OUTPUT_FILE[.EXT] ', 'action': 'store', 'type' :str, 'default':''},

        '--informat' : {'help':('force input file to be parsed as specified'),
                             'action':'store', 'metavar':'[mol2|sdf|smi|pdb|...]', 'type':str},

        '--outformat' : {'help':('force output file format to be written as specified'),
                             'action':'store', 'metavar' : '[mol2|sdf|pdb|...]', 'type':str, 'default':''},

        # Molecule naming schemes
        '--usemolname':{'help':('enable \'--split\' and use molecule names for saving splitted output files; '
                           'when names are not found or duplicated, progressive numbering is used; if --outfile\n'
                           'is used, names will have the outname string as prefix'
                           ),
                           'action':'store_true', 'default': self.default_usemolname},

        '--usemolnamesafe':{'help':('same as \'--usemolname\', but sanitize names to be valid/clean '
                           'filenames'),
                           'action':'store_true', 'default': self.default_usemolnamesafe},


        # '--fix_non_std_autodock':{'help':('fix atoms that are not in the standard AD4 force field with the closest existing atom type'),
        #                    'action':'store_true', 'default': self.default_fix_non_std_autodock},


        '--usefieldname' : {'help' : ('use specified field as molecule name (only for SDF, Mol2 '
                             'and PDB input files)')
                             #; NOTE: spaces must specified with HTML notation, '
                             #'i.e. \'Cpd ID\' -> \'Cpd%%20ID\' (thank you sys.argv!)')
                             , 'metavar': 'FIELD', 'action':'store', 
                             'default' : None},

        # Splitting
        '--split' : { 'help':('split each molecule as separate file, progressively numbered '
                                '(see \'--usemolname\' for using molecule names); '
                                'disabled by default if output format supports multi-structure'),
                                'default': self.default_split, 'action': 'store_true'},
        
        # minimizing 
        # XXX TODO XXX TODO XXX TODO
        # provide a minimum set of values good for everything
        # '--mini_quick' : { 'help':('steepest descent max. iterations; set to 0 to disable '
        #                 '(default %d).') % self.default_sdsteps, 'action': 'store', 'metavar':'[ %4g ]' % self.default_sdsteps,
        #                 'type' : int, 'default' : None},
        # '--mini_std' : {},
        # '--mini_aggressive' : {},
        # XXX TODO XXX TODO XXX TODO

        '--automini' : { 'help':('set minimization parameters automatically basing on each'
                            ' molecule complexity; if no argument is provided, by default \'%s\' '
                            'settings are used; each level has 10x more energy evaluations than the previous; '
                            'the formula used to calculate evaluations is : '
                            '(bonds/10) * (atoms/10) * rotatable_bonds * LEVEL') % self.default_heuristics,
                         'nargs':'?',  
                         'action': 'store', 
                         'const': "quick",
                         'metavar':'%s]' % "|".join(self.default_heuristics_list),
                        'type' : str, 'default' : None},

        '--sdsteps' : { 'help':('steepest descent max. iterations; set to 0 to disable '
                        '(default %d).') % self.default_sdsteps, 'action': 'store', 'metavar':'[ %4g ]' % self.default_sdsteps,
                        'type' : int, 'default' : None},

        '--sdconv' : {  'help':('set steepest descent convergence criterion. (default: %s)') % self.default_sdconv, 
                        'action': 'store', 'metavar':  '[ %4g ]' % self.default_sdconv, 'type': float,
                        'default': self.default_sdconv},

        '--cgsteps': {  'help': ('conjugate gradients max. iterations; set to 0 to disable '
                         '(default %d).' % self.default_cgsteps), 'action': 'store', 'metavar':'[ %3g ]' % self.default_cgsteps,
                         'type' :  int, 'default': None},

        '--cgconv' : { 'help':('set conjugate gradients descent convergece criterion. '
                             '(default: %s)' % self.default_cgconv), 'action': 'store',
                             'metavar': '[ %4g ]' % self.default_cgconv, 'type':float, 'default': None},

        '--sdsteps_extra' : { 'help':('steepest descent max. iterations (post-processing minimization with UFF?); set to 0 to disable '
                        '(default %d).') % self.default_sdsteps_extra, 'action': 'store', 'metavar':'[ %4g ]' % self.default_sdsteps_extra,
                        'type' : int, 'default' : None},

        '--sdconv_extra' : {  'help':('set steepest descent convergence criterion (post-processing minimization with UFF?). (default: %s)') % self.default_sdconv_extra, 
                        'action': 'store', 'metavar':  '[ %4g ]' % self.default_sdconv_extra, 'type': float,
                        'default': self.default_sdconv_extra},

        '--cgsteps_extra': {  'help': ('conjugate gradients max. iterations (post-processing minimization); set to 0 to disable '
                         '(default %d).' % self.default_cgsteps_extra), 'action': 'store', 'metavar':'[ %4g ]' % self.default_cgsteps_extra,
                         'type' :  int, 'default': None},

        '--cgconv_extra' : { 'help':('set conjugate gradients descent convergece criterion (post-processing minimization). '
                             '(default: %s)' % self.default_cgconv_extra), 'action': 'store',
                             'metavar': '[ %4g ]' % self.default_cgconv_extra, 'type':float, 'default': None},

        '--rotamer_conf' : { 'help':('number of random conformers to use in weighted conformational search '
                             '(default: %s)' % self.default_rotamer_conf), 'action': 'store',
                             'metavar': '[ %d ]' % self.default_rotamer_conf, 'type':int, 'default': None},
        
        '--rotamer_steps' : { 'help':('number of geometry minimization steps to per form during weighted '
                              'conformational search; this value should be small, otherwise calculation times '
                              'will increase considerably) (default: %s)' % self.default_rotamer_steps), 'action': 'store',
                             'metavar': '[ %s ]' % self.default_rotamer_steps, 'type':int, 'default': None},

        '--forcefield' : {'help': ('set forcefield for minimization (default: \'%s\'; see --help_advanced for '
                            'choices)') % (self.default_forcefield) , 
                         'default': self.default_forcefield, 'metavar':'[ %s ]' % self.default_forcefield,
                         'type': str}, 

        '--forcefield_extra' : {'help': ('set forcefield for post-minimization (default: \'%s\'; see --help_advanced for '
                            'choices)') % (self.default_forcefield_extra) , 
                         'default': self.default_forcefield_extra, 'metavar':'[ %s ]' % self.default_forcefield_extra,
                         'type': str}, 

        '--strict'     : {'help':'do not process molecules with atoms missing forcefield/charges '
                          'parameters (default: %s)' % self.default_strict, 'action':'store_true' },

        '--nomini' :{ 'help': 'disable minimization',
                    'default':self.default_nomini, 'action': 'store_true'},

        '--noextra' :{ 'help': 'disable extra minimization',
                    'default':self.default_noextra, 'action': 'store_true'},

        '--norotamer' :{ 'help': 'disable rotameric minimum search',
                    'default': self.default_norotamer, 'action': 'store_true'},

        '--log' :  { 'help': 'XXASADAD Asave minimization logs in specified file; by default the output '
                     'filename is used  (OUTFILE.log)', 'action': 'store', 'metavar':'LOGFILENAME'},

        '--verbose' : {'help': ('enable verbose mode'),
                                'default' : self.default_verbose, 'action': 'store_true'},

        '--noflipamide' :  { 'help': 'disable forcing cis conformation for non-cyclic primary/secondary amides', 
                             'action': 'store_true', 'default': self.default_noflipamide},

        '--nostripsalts' :  { 'help': 'disable removing salts/ions and small fragments (i.e. '
                            'keep largest fragment only)', 'action': 'store_true', 
                             'default': self.default_nostripsalts},

        '--nocheckhydro' :  { 'help': 'disable check for missing hydrogens', 'action': 'store_true', 
                            'default': self.default_nocheckhydro },

        '--noprocess' :  { 'help': 'disable all structure processing (equivalent to (--noflipamide) + '
                            '(--nocheckhydro) + (--nostripsalt) + (--nopH)', 'action': 'store_true', 
                            'default': self.default_noprocess },

        '--pH' :  { 'help': 'generate protonation state for requested pH (default: %1.1f)' % self.default_pH, 
                    'action': 'store', 'type': float, 'default': None, 'metavar':'[ 7.4 ]'}, # default set to None to perform check later

        '--nopH' :  { 'help': 'disable protonation state modifications', 
                    'action': 'store_true', 'default': None}, # default set to None to perform check later

        '--chargemodel' : { 'help': ('set charge model to compute atom partial '
                            'charges (default \'%s\'; see --help_advanced for options)')
                            % (self.default_charges), 'metavar': '[ %s ]' % self.default_charges,
                            'action': 'store', 'type': str, 'default': self.default_charges}, 
        '--enumchiral' : { 'help' : ('generate multiple stereoisomers for tetrahedral chiral'
                                'centers (note: input geometry will be ignored); \'undefined\': '
                                'enumerate only undefined chiral centers; '
                                '\'all\': enumerate all chiral centers; default: \'off\' (no enumeration)'),
                                'default': self.default_enumchiral, 'action' : 'store',
                                'metavar' : '[undefined|protomer|all|off]', 'type' : str, },

        '--maxenumchiral' : { 'help': ('maximum number of enantiomers to generate if --enumchiral'
                        ' is active; default: %d' % self.default_maxenumchiral),
                'default': self.default_maxenumchiral, 'action' : 'store', 
                'metavar': 'MAX_ENANTIOMERS', 'type':int},

        '--exclude' :  { 'help': ('skip molecules matching specified SMARTS pattern; '
                                 'multiple patterns can be repeated by using multiple \'--exclude\' '
                                'or using \'--excludefromfile\''), 
                    'action': 'append', 'metavar': 'SMARTS', 'type': str, 'default': []},

        '--excludefromfile' :  { 'help': ('skip molecules matching any of the SMARTS patterns'
                                        ' in each line of FILENAME'), 
                    'action': 'store', 'metavar': 'FILENAME', 'type': str, 'default': None},

        '--multicore' :  { 'help':('specify how many cores/cpu to use '
                             '(default: %d)' % self.default_multiproc_max), 'action': 'store',
                             'metavar': '[ %d ]' % self.default_multiproc_max, 'type':int, 
                            'default': self.default_multiproc_max},


        '--nice' : {'help': 'set nice level (not working: requires psutil)'},

        '--help_advanced' : {'help': 'show the advanced help with methods, etc...', 'action': 'store_true',
                    'default': False},

                    }

        # enforce the proper arguments order for the help message
        self._optparser_args_order = {  # in/out
                'INPUT/OUTPUT' : [ '--infile', '--outfile',  '--informat', 
                        '--outformat', '--usemolname', '--usemolnamesafe', '--usefieldname' ],
                # slicing
                'SLICING' : ['--single', '--begin', '--end', '--split', '--byname'],
                # mol processing 
                'CLEANUP': [ '--pH', '--nopH', '--noflipamide', '--nostripsalts', '--nocheckhydro', 
                            '--noprocess', '--enumchiral', '--maxenumchiral', '--exclude', '--excludefromfile',
                            # '--fix_non_std_autodock',
                            ],
                # minimizer 
                'ENERGY MINIMIZATION AND FORCEFIELD PARAMETERS' : [ '--automini', '--sdsteps', '--sdconv', '--cgsteps', '--cgconv', '--forcefield', 
                             '--sdsteps_extra', '--sdconv_extra', '--cgsteps_extra', '--cgconv_extra', '--forcefield_extra', 
                             '--rotamer_conf', '--rotamer_steps',

                '--nomini', '--noextra', '--norotamer', '--chargemodel', '--strict'],

                # multiprocessing
                'MULTIPROCESSING' : ['--multicore', '--nice'],

                # logging and info
                'LOGGING': ['--log', '--verbose', '--help_advanced'],
                     }

        self._optparser_groups_order = [
            ('INPUT/OUTPUT', ('Control input/output data files and types.')),
                             #'%s supports (in theory)  most  of the OpenBabel  file formats.\n'
                             #'3D structures are generated automatically for molecules and files\n'
                             #'that do not have 3D information (i.e., SMI, CDX).')),

            ('SLICING', ('Multi-structure files management')),
                        #'For files containing multiple molecules, it is possible to specify\n'
                        #'the range of molecules to be processed, split a multi-structure file in\n'
                        #'separate files.' )),
            ('CLEANUP', 'Pre- and post-processing cleanups performed on molecules'),

            ('ENERGY MINIMIZATION AND FORCEFIELD PARAMETERS',
                        'Options for the structure optimization'),

            ('MULTIPROCESSING', 'Manage multi-processing and parallel computation, nice, etc...'),

            ('LOGGING', 'Control logging, message reports and extra info'),
                        ]


    def _init_mol_parser(self):
        """ initialize the molecule parser"""
        self.ob_mol_parser = ob.OBConversion()
        """
        self.ob_mol_parser.SetInAndOutFormats(self.intype,self.outtype)
        self.dprint("[_init_mol_parser] initialized ob_mol_parser")
        if self.split or self.single:
            # files will be read and written separately
            return
        self.dprint("[_init_mol_parser] Open combined in/output parser")
        self.ob_mol_parser.OpenInAndOutFiles(self.infile, self.outfile)
        """
        self.ob_mol_parser.SetInAndOutFormats(self.intype,self.intype) #outtype)
        self.dprint("[_init_mol_parser] initialized ob_mol_parser: %s,%s"% (self.intype, self.intype))
        #if self.split or self.single:
        #    # files will be read and written separately
        #    return
        #self.dprint("[_init_mol_parser] Open combined in/output parser")
        #self.ob_mol_parser.OpenInAndOutFiles(self.infile, self.outfile)

    def read_mol(self, first=False):
        """ initialze the molecule parser and get the first molecule"""
        """ read the next molecule from the molecular parser
            or report when the end is reached
        """
        mol = ob.OBMol()
        # XXX This part is not very clear...
        if first:
            #if (not self.split) and (not self.single):
            #    more = self.ob_mol_parser.Read(mol)
            #else:
            #    more = self.ob_mol_parser.ReadFile(mol, self.infile)
            if self.split or (not self.single == None):
                more = self.ob_mol_parser.ReadFile(mol, self.infile)
            else:
                more = self.ob_mol_parser.ReadFile(mol, self.infile)
        else:
            more = self.ob_mol_parser.Read(mol)
        if not more:
            return False
        self._counter += 1
        return mol

    def get_mol_string(self, mol):
        """ generate string to submit to the queue"""
        return self.ob_mol_parser.WriteString(mol) #, self.intype)
    
    def get_mol_title(self, mol):
        """ get the molecule name stored in the OBMol obj
            and sanitize it if required: spaces will be converted to
            underscores, then unacceptable characters will be
            discarded
        """
        title = mol.GetTitle()
        if not self.safenames:
            return title
        title = title.replace(' ', '_')
        validname = ''.join(c for c in title if c in self._validFnameChars)
        return validname

    def get_mol_name(self, mol, data={}): #, field=None):
        """ find the appropriate unique name for the molecule
            accordingly to the naming policy selected;
            if the molecule lacks a proper name or the
            progressive numbering is asked, a new name
            is genearated and accounting performed

            perform also the filter by name, returning None
            if no match with byname value is found
        """
        # get the molecule name
        mname = ''
        #print "FFFFIEL", self.args.usefieldname
        if not self.args.usefieldname == None:
            self.args.usefieldname = self.args.usefieldname.replace('%20', ' ')
            self.dprint('[get_mol_name] Search field name :"%s"' % self.args.usefieldname)
            # by field 
            try:
                mname = data[self.args.usefieldname].strip()
                self.dprint('[get_mol_name] Found name :"%s"' % mname)
                mol.SetTitle(mname)
            except KeyError:
                mname = ''
                self.dprint('[get_mol_name] WARNING: field name not found')
        elif self.usemolname:
            # by mol name
            mname = self.get_mol_title(mol)
            if mname == '':
                mname = 'MOL'
        else:
            # by default
            mname = "MOL_%d" % self._counter
        if mname.strip() == '':
            mname = 'MOL'
        # filter and accept only molecule by name
        if self.args.byname:
            if not self.get_mol_title(mol).strip() == self.args.byname:
                msg = ("[get_mol_name] Skipping molecule by name filter (%s) [%s]")
                self.dprint(msg % (self.get_mol_title(mol), self.args.byname) )
                return None
        if not mname in list(self.mol_names_dict.keys()):
            self.mol_names_dict[mname] = -1
        self.mol_names_dict[mname] += 1
        if self.mol_names_dict[mname] == 0:
            name = mname
        else:
            name = "%s_%d" % (mname, self.mol_names_dict[mname])
        self.dprint("[get_mol_name] Name is now [%s]" % name)
        return name

    def get_mol_data(self, mol):
        """ parse generic data (e.g., SDF, Mol2 fields...)
            [ this code is shamelessly borrowed from pybel ]
        """
        data = mol.GetData()
        answer = [ x for x in data if
            x.GetDataType() == ob.PairData or
            x.GetDataType() == ob.CommentData]    
        raw = [ ( x.GetAttribute(), x.GetValue() ) for x in answer ]
        return dict(raw)        
        
    def write_output(self, mol):
        """ write processed molecule accordingly to the file
            output policy: update the multi-structure file
            or create a new file 
        """
        if self.split or self.single:
            fname = self.outfile % (self.current_mol_name)
            self.dprint( "[write_output] Saving new file: %s"% fname)
            self.ob_mol_parser.WriteFile(mol, fname)
        else:
            self.dprint( "[write_output] Updating file: %s"% self.outfile)
            self.ob_mol_parser.Write(mol)

    def _close_mol_parser(self):
        """ close output file used for multi-structure savings"""
        if not self.split:
            self.ob_mol_parser.CloseOutFile()

    def heuristic_mini_parms(self, mol, level='quick'):
        """ try to guess """
        level_settings = {'quick':1, 'accurate':10, 'long':100, 'extreme':1000} 
        rot = mol.NumRotors()
        atoms = mol.NumHvyAtoms()
        bonds = mol.NumBonds()
        print("[PARM HEURISTIC]  %s R:%d  A:%d  B:%d"% (mol.GetTitle(), rot, atoms, bonds), end=' ')
        sdsteps = max( ( (bonds/10.) * (atoms/10.) * rot), 10) * level_settings[level]
        cgsteps = sdsteps / 2
        sdsteps_extra = sdsteps * 0.75
        cgsteps_extra = cgsteps * 0.75
        print("==> SD steps [ %d ] CG steps [ %d ]   SD_extra steps [ %d ] CG_extra steps [ %d ]" % ( sdsteps, cgsteps, sdsteps_extra, cgsteps_extra))
        # return

        self.scrub_options['sdsteps'] = int(sdsteps)
        self.scrub_options['cgsteps'] = int(cgsteps)
        self.scrub_options['sdsteps_extra'] = int(sdsteps_extra)
        self.scrub_options['cgsteps_extra'] = int(cgsteps_extra)
        self.scrub_options['rotamerConf'] = int(rot)


    def _init_threads(self):
        """ """
        # print("CALLED", self.forcefield)
        self.scrub_options = { 'ff':self.forcefield, 'sdsteps':self.sdsteps,
            'sdconv':self.sdconv, 'cgsteps':self.cgsteps, 
            'cgconv':self.cgconv, 'ff_extra':self.forcefield_extra, 
                'sdsteps_extra':self.sdsteps_extra, 
                'sdconv_extra':self.sdconv_extra, 
                'cgsteps_extra':self.cgsteps_extra,
                'cgconv_extra':self.cgconv_extra, 
                'rotamer_conf' : self.rotamer_conf,
                'rotamer_geom_steps' : self.rotamer_steps,
                'pH':self.pH, 'chargemodel':self.chargemodel,
                'flip_cis_amide': self.flipamide,
                'stripsalts':self.stripsalts, 
                'checkhydro':self.checkhydro,   'name': None, #self.current_mol_name, 
                'verbose':self.verbose, 'auto':True
            }
        # feeding queue
        self.queue_in = mp.JoinableQueue(maxsize=self.multiproc_max)
        # writer queue
        self.queue_out = mp.Queue(maxsize=-1)
        # create processes and start them
        for i in range(self.multiproc_max):
           s = ScrubMultiprocess(self.queue_in, self.queue_out, nice=None, 
                in_type=self.intype, out_type=self.outtype,  scrub_opts = self.scrub_options)
           s.start()
        # XXX this has to be fixed for split, etc...
        if self.split or self.single:
            fname = None
        else:
            fname = self.outfile
        writer = MolWriterMultiProcess(queue=self.queue_out, threads=self.multiproc_max,
                    fname=fname, ext=self.outtype)
        writer.start()
        print("%d threads initialized" % self.multiproc_max)


    def _close_threads(self):
        """ """
        for i in range(self.multiproc_max):
            self.queue_in.put((None,None))
        #self.queue_out.put((None, None))
        
    def _init_loop(self):
        """ initialize variables, molecule input file parser and threads"""
        self._counter = 1 # general counter
        self._rejected = 0
        self._processed = 0
        self._init_mol_parser()
        #self._init_threads()
        
    
    def start(self):
        """ main loop does the actual job"""
        self.dprint('====================[ START ]====================')
        # start!
        self._init_loop()

        scrub_options = { 'ff':self.forcefield, 'sdsteps':self.sdsteps,
            'sdconv':self.sdconv, 'cgsteps':self.cgsteps, 
            'cgconv':self.cgconv, 'ff_extra':self.forcefield_extra, 
                'sdsteps_extra':self.sdsteps_extra, 
                'sdconv_extra':self.sdconv_extra, 
                'cgsteps_extra':self.cgsteps_extra,
                'cgconv_extra':self.cgconv_extra, 
                'rotamer_conf' : self.rotamer_conf,
                'rotamer_geom_steps' : self.rotamer_steps,
                'pH':self.pH, 'chargemodel':self.chargemodel,
                'flip_cis_amide': self.flipamide,
                'stripsalts':self.stripsalts, 
                'checkhydro':self.checkhydro,   'name': None, #self.current_mol_name, 
                'verbose':self.verbose, 'auto':True
            }

        output_fp = open(self.outfile, "w")
        scrub = Scrub(**scrub_options)

        raw_mol = self.read_mol(first=True)
        if not self.begin == 1:
            print("Seeking begin molecule %d... " % (self.begin))
        while raw_mol:
            # getting the data is more efficient (before enumerating chirals)
            # check that there are atoms (useful for CDX files)
            if not raw_mol.NumAtoms():
                self.dprint("No atoms, skipping molecule")
                raw_mol = self.read_mol()
                continue
            # skip by counter
            if self._counter < self.begin:
                raw_mol = self.read_mol()
                self.dprint("[start] molecule[%s] skipped by counting" % self._counter)
                continue
            molData = self.get_mol_data(raw_mol) 
            # enumerate chiral structures
            for mol in ChiralEnumerator(raw_mol, self.args.enumchiral,
                                self.args.maxenumchiral, debug=self.verbose):
                # TODO consider if grouping these functions in an "EvaluateMol" wrapper?
                # naming
                self.current_mol_name = self.get_mol_name(mol, molData)
                if self.current_mol_name == None:
                    # skip by name
                    continue
                msg = "-------------[processing mol.%d: %s]---------" 
                if not self.automini is None:
                    print("\n\nGUESSING PARMS!!!!\n\n")
                    self.heuristic_mini_parms(mol, self.automini)
                args = (self._counter, self.current_mol_name)
                self.dprint(msg % args)
                # skip by SMARTS
                if not self.filter_smarts(mol):
                    molRaw = self.read_mol()
                    self.dprint("[start] molecule[%s] skipped by SMARTS" % self.current_mol_name)
                    continue
                mol_string = self.get_mol_string(mol)
                #print("start scrub: %s" % mol.GetTitle())
                output_string = scrub.process(mol_string, self.intype, self.outtype)
                output_fp.write(output_string)
                self._processed += 1
                #print("scrubbed successfully %d" % self._processed)
                #self.add_mol_to_queue(mol, self.current_mol_name)
            # terminate by counter
            if self.end and (self._counter == self.end):
                break
            raw_mol = self.read_mol()
        output_fp.close()
        # poison pill
        # self._close_threads()
        print("[ DONE ]Total structures processed: %d" % (self._processed))

    def add_mol_to_queue(self, mol, name):
        """ """
        string = self.get_mol_string(mol)
        self.queue_in.put((string,name), block=True)

    # TODO not used?
    #def setStopCriterion(self):
    #    for i in range(self.multiproc_max): 
    #        self.queue_in.put((None,None))
    #    self.queue_out.put(None)
    #    #self.queue_in.close()
    #    #self.queue_in.join_thread()
         
###### XXX XXX HIC SUNT LEONES
if __name__ == '__main__':
    MiniMee()
    """
    print "- COUNT THE MOLECULES IN THE FILE? TO HAVE A PROPER ZERO-PADDING?"
    print "- ADD EXAMPLES"
    print "- ADD INCLUDE (opposite than EXCLUDE"
    print "PH RANGES..."
    print "TAUTOMERS"
    print "SULFONAMIDES FLIP?"
    print "NONPRIMARY AMIDE FLIPS?"
    print "- ADD AUTOFILTER FOR REACTIVE GROUPS"
    print "preserve properties  (GETDATA)"
    print "RENAME atoms to have unique names"
    print "\n\n\n"


    print "GET IDEAS:http://openbabel.org/docs/dev/Command-line_tools/babel.html#append-option " 
    print "\n Implement UNIQUE ATOM NAMES"
    #print "-- enumerate chirality"
    """
