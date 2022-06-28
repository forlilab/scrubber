from operator import itemgetter

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem


from ..common import ScrubberClass, UniqueMoleculeContainer, mol2smi
# from scrubber.core.base import ScrubberClass, UniqueMoleculeContainer, mol2smi
from .base import (
    MoleculeTransformations,
    MaxResultsException,
    MaxIterException,
    MolecularReactionsLogger,
)

# this file contains the definition of the named default reactions used
REACTIONS_FILE = "reactions.txt"


class Reactor(ScrubberClass, MoleculeTransformations):
    """perform chemical reactions to functionalize molecules
    # the two modes should be:
    #   - "all", (exhaustive, all sites saturated)
    #   - "single", ( enumeration of all mono-modified sites )

    TODO protect warhead atoms to prevent polymerization and endless reactions!
    https://www.rdkit.org/docs/GettingStartedInPython.html#protecting-atoms
    """

    default_init = {
        "mode": "all",
        "reaction_file": None,
        "verbose": False,
        # "suppress_rdkit_warnings": True,
    }
    # default_options = {"reactions_list": None, "keep_only_reacted": False}
    # PROTEC
    def __init__(
        self,
        mode: str = "all",
        reactions_file: str = None,
        keep_only_reacted: bool = False,
        verbose: bool = False,
        # suppress_rdkit_warnings: bool = False,
        _stop_at_defaults: bool = False,
    ):
        self.mode = mode
        self.reactions_file = reactions_file
        self.keep_only_reacted = keep_only_reacted
        self.verbose = verbose
        # self.suppress_rdkit_warnings = suppress_rdkit_warnings
        if _stop_at_defaults:
            return
        MoleculeTransformations().__init__(
            self,
            self.verbose,
            self.suppress_rdkit_warnings,
        )
        self.__init_reactions(self.reactions_file)

    def __init_reactions(self, reaction_file=None, only_custom=False):
        """initialize the default reactions for the reactor
        by default, standard reactions from the data file are parsed, then
        custom reactions can be added if requested. If 'custom_only' flag is
        used, then default reaction file is not parsed"""
        self.reactions = {}
        if not custom_only:
            default_reactions = self.get_datafile(REACTIONS_FILE)
            reactions, errors = self.__parse_reaction_file(default_reactions)
            for rxn_obj, _, _, tag in reactions:
                self.reactions[tag] = rxn_obj
        if not reactions_file is None:
            reactions, errors = self.__parse_reaction_file(reaction_file)
            for rxn_obj, _, _, tag in reactions:
                self.reactions[tag] = rxn_obj

    def add_reaction(self, reaction_string, name=None):
        """add a custom reaction from a valid SMIRKS string, optionally"""
        rxn_obj, rxn_left, rxn_right, tag = self.__parse_reaction_line(line)
        if not name is None:
            tag = name
        if tag == "":
            i = 0
            while (tag in self.reactions) or (i == 0):
                tag = "rxn_%d" % (len(self.reactions) + i)
                i += 1
        self.reactions[tag] = rxn_obj
        return tag

    def call_Reaction_here(self, reaction_name):
        """this is where the reaction occurs"""
        # this is a class function with 0 indentation
        pass
