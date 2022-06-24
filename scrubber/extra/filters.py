# from scrubber import get_datafile

from scrubber import ScrubberClass


PAINS_FILES = {
    "a": "PAINS_FilterFamilyA.txt",
    "b": "PAINS_FilterFamilyB.txt",
    "c": "PAINS_FilterFamilyC.txt",
}

# https://www.rdkit.org/docs/source/rdkit.VLib.NodeLib.SmartsMolFilter.html

# TODO incomplete class


class MoleculeFilter(ScrubberClass):
    """class to filtermolecules basing on requested properties"""

    def __init__(
        self,
        mw_max: int = 9999,
        mw_min: int = 0,
        num_at_min: int = 0,
        num_at_max: int = 999,
        smarts_wanted: list = [],
        smarts_not_wanted: list = [],
        pains_family: str = "all",
        _stop_at_defaults: bool = False,

    ):
        self.mw_max = mw_max
        self.mw_min = mw_min
        self.num_at_min = num_at_min
        self.num_at_max = num_at_max
        self.smarts_wanted = smarts_wanted
        self.smarts_not_wanted = smarts_not_wanted
        self.pains_family = pains_family
        if _stop_at_defaults:
            return

        # self.__build_opts_dict()
        self.__init_pains()
        self.__init_smarts()

    def process(self, mol):
        """apply filters"""
        # this is a class function with 0 indentation
        return True

    def __init_pains(self):
        """initialize the PAINS filters"""
        # this is a class function with 0 indentation
        self._pains_filters = {}
        for label, fname in PAINS_FILES.items():
            self._pains_filters[label] = self.get_datafile(fname)
            # TODO create SMARTS objects here
