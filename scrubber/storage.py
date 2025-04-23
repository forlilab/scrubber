import multiprocessing

from rdkit import Chem

class SMIMolSupplierWrapper:
    """RDKit SMI molecule supplier wrapper."""

    def __init__(
        self,
        filename: str,
        sanitize: bool = True,
        titleLine: bool = False,
        queue_err: multiprocessing.Queue = None,
        discarded_input_fname: str = None,
        is_enamine_cxsmiles: bool = False,
        _stop_at_defaults: bool = False,
    ):
        self.filename = filename
        self.sanitize = sanitize
        self.titleLine = titleLine
        self.queue_err = queue_err
        self.discarded_input_fname = discarded_input_fname
        self.is_enamine_cxsmiles = is_enamine_cxsmiles
        if _stop_at_defaults:
            return
        self.fp_input = open(filename, "r")
        self.fp_errors = None
        self._buff = []
        if self.titleLine:
            self.fp_input.readline() # ditch first line
        # print("INITIALIZED", self.titleLine)
        # print("INITIALIZED", self.queue_err)

    def _close_fp(self):
        """close all potential file pointers"""
        self.fp_input.close()
        if not self.fp_errors is None:
            self.fp_errors.close()

    def __iter__(self):
        self.fp_input.seek(0)
        if self.titleLine:
            self.fp_input.readline() # ditch first line
        return self

    def reset(self): # same interface as rdkit.Chem.SDMolSupplier
        self.fp_input.seek(0)

    def __next__(self):
        """iterator step"""
        while True:
            line = self.fp_input.readline()
            # readline() returns "" at the end of file, "\n" for empty lines in the file
            if not line:
                self.fp_input.close()
                if not self.fp_errors is None:
                    self.fp_errors.close()
                raise StopIteration
            # skip empty lines
            if not line.strip():
                continue
            try:
                if self.is_enamine_cxsmiles:
                    smiles, name, _ = line.split("\t", maxsplit=2)
                    mol = Chem.MolFromSmiles(smiles, sanitize=self.sanitize)
                    mol.SetProp("_Name", name) 
                else:
                    mol = Chem.MolFromSmiles(line, sanitize=self.sanitize)
                if mol is None:
                    if not self.queue_err is None:
                        self.queue_err.put(("input", line), block=True)
                    if not self.discarded_input_fname is None:
                        if self.fp_errors is None:
                            self.fp_errors = open(self.discarded_input_fname)
                        self.fp_errors.write(line + "\n")
                return mol
            except:
                if not self.queue_err is None:
                    self.queue_err.put(("input", line), block=True)
                if not self.discarded_input_fname is None:
                    if self.fp_errors is None:
                        self.fp_errors = open(self.discarded_input_fname)
                    self.fp_errors.write(line + "\n")
                return None
