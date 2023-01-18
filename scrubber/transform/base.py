from string import ascii_uppercase, ascii_lowercase
from collections import defaultdict

from rdkit import RDLogger
from rdkit import Chem #, AtomValenceException, KekulizeException
from rdkit.Chem.PropertyMol import PropertyMol

from ..common import UniqueMoleculeContainer, mol2smi, copy_mol_properties


class MolecularReactionsLogger(object):
    """graph to store molecular reactions"""

    def __init__(self, verbose=False, fname=None):
        """"""
        self.verbose = verbose
        self.logfile = None
        if not fname is None:
            self.logfile = open(fname, "wb")
        self.graph = defaultdict(list)
        self.inverse_graph = defaultdict(list)
        self.trajectory = {}
        self.data = {}  # defaultdict(lambda:0)
        self.labels_pool = ascii_uppercase + ascii_lowercase
        self.label_counter = -1
        self._largest_string = 0
        self._first = None
        self._counter = 0

    def add(self, reagent_mol, product_mol, transformation, iteration=0, kept=False):
        """register a molecule in the molecular graph"""
        self._counter += 1
        reagent, reagent_name, reagent_count = self._create_name(reagent_mol)
        self._largest_string = max(len(reagent), self._largest_string)
        if self._first is None:
            self._first = reagent
        if not product_mol is None:
            product, product_name, product_count = self._create_name(product_mol)
            self._largest_string = max(len(product), self._largest_string)
        else:
            product = None
            product_name = None
            product_count = None
        # store the direct reaction
        self.trajectory[iteration] = (
            reagent,
            reagent_name,
            reagent_count,
            product,
            product_name,
            product_count,
            transformation,
            kept,
        )
        direct = (product, transformation, iteration, kept)
        self.graph[reagent].append(direct)
        # store the inverse reaction
        if not product_mol is None:
            inverse = (reagent, transformation, iteration, kept)
            if inverse in self.inverse_graph[product]:
                if self.verbose:
                    print("WARNING: duplicate reaction:", inverse)
            self.inverse_graph[product].append(inverse)
        if not self.logfile is None:
            # TODO implement a file logging system
            pass
        return reagent_name, product_name

    def _create_name(self, mol):
        """ """
        smi = mol2smi(mol)
        if not smi in self.data:
            self.label_counter += 1
            label = self.labels_pool[self.label_counter % len(self.labels_pool)]
            self.data[smi] = [0, label]
        self.data[smi][0] += 1
        name = self.data[smi][1]
        count = self.data[smi][0]
        return smi, name, count

    def get_info(self, mol):
        """retrieve SMILES and molecule name (e.g.: A, B, ...), if registered,
        None otherwise"""
        smi = mol2smi(mol)
        if not smi in self.data:
            return None
        return smi, self.data[smi][1]

    def _get_current_label(self):
        """DEBUG FUNCTION TO RETRIVE CURRENT LABEL"""
        return self.labels_pool[self.label_counter % len(self.labels_pool)]

    def get_reagent(self, product):
        """given a product, return reagent and transformation to obtain it"""
        smi = mol2smi(product)
        return self.inverse_graph[smi]

    def get_product(self, reagent):
        """given a reagent return the list of products and transformations"""
        smi = mol2smi(reagent)
        return self.graph[smi]

    def print_history(self):
        """print a text version of the history"""
        print("\n\n===============================================")
        print("Reactions trajectory")
        # print("steps registered", self._counter)
        print("===============================================")
        format_line_full = (
            "%%5d - %%%ds ( %%2s ) [ %%3d ] "
            "  ---->     %%%ds ( %%2s ) [ %%3d ]  \t  |  %%s "
            % (
                self._largest_string,
                self._largest_string,
            )
        )
        format_line_partial = (
            "%%5d - %%%ds ( %%2s ) [ %%3d ] "
            "  --X->     %s                 \t  |  %%s "
            % (
                self._largest_string,
                (self._largest_string) * " ",
            )
        )

        idx = 0
        for idx, (iter_step, result) in enumerate(self.trajectory.items()):
            if not result[3] is None:
                print(
                    format_line_full
                    % (
                        iter_step,
                        result[0],
                        result[1],
                        result[2],
                        result[3],
                        result[4],
                        result[5],
                        result[6],
                    ),
                )
            else:
                print(
                    format_line_partial
                    % (
                        iter_step,
                        result[0],
                        result[1],
                        result[2],
                        # result[3],
                        # result[4],
                        # result[5],
                        result[6],
                    ),
                )

        print("steps printed", idx + 1)


class MoleculeTransformations(object):
    """Generic class to handle transformations

    Children classes need to implement a process() method that uses
    internally the __exhaustive_reaction() method to perform the actual
    transformations.
    """

    def __init__(self, verbose: bool = False, suppress_rdkit_warnings: bool = True):
        self.verbose = verbose
        self._set_rdkit_log(enable=not suppress_rdkit_warnings)
        self.reaction_log = None

    def process(self, mol, *args):
        """prototype of the process method used to perform the specific tranformations
        The reaction_log object is an instance of the MolecularReactionsLogger used to track
        """
        self.reaction_log = None
        raise NotImplemented

    def _set_rdkit_log(self, enable=False):
        """convenience function to enable/disable the RDKit logging messages during the process"""
        if enable:
            RDLogger.EnableLog("rdApp.*")
            if self.verbose:
                print("[VERBOSE] RDKit logging *ENABLED*")
        else:
            RDLogger.DisableLog("rdApp.*")
            if self.verbose:
                print("[VERBOSE] RDKit logging *DISABLED*")

    def _exhaustive_reaction(
        self,
        reaction_list: list,
        keep_all: bool = True,
        max_results: int = 1000,
        max_iter: int = 10000,
        extra_reagent=None,
        reaction_log=None,
        preserve_mol_properties:bool=True,
    ) -> tuple:
        """perform the RDKit reaction. The format of the reaction is [(rxn_id,
        reaction)] This method processes the molecular pool (self.mol_pool) but
        does not add results to it, letting the calling method to take care of
        that (allowing for filtering etc.)
        https://github.com/tentrillion/rdkit-tutorials/blob/master/notebooks/003_SMARTS_ReactionsExamples.ipynb

        TODO this function should/could be used to attach covalent residues or warheads
        extra_reagent = covalent residue/warhead...

        reaction_list               : list containing the RDKit reaction objects to be applied
        keep_all                    :
        preserve_mol_properties     : copy all molecular properties of the reagents into the generated molecules
        """
        # counter for the max_iter counting
        _iter = 0
        # results = UniqueMoleculeContainer()
        reaction_products = UniqueMoleculeContainer()
        consumed_reagents = UniqueMoleculeContainer()
        visited = UniqueMoleculeContainer(ignore_chirality=True)
        visited_pairs = set()
        try:
            for rxn_idx, (rxn_pattern, rxn_obj) in enumerate(reaction_list):
                self._iterations += 1
                if not keep_all:
                    if len(reaction_products):
                        reagents_pool = reaction_products.copy()
                        reaction_products.clear()
                    else:
                        reagents_pool = self.mol_pool.copy()
                    reagents_pool -= consumed_reagents
                else:
                    reagents_pool = self.mol_pool.copy()
                    reagents_pool += reaction_products
                    reaction_products.clear()
                    reagents_pool += consumed_reagents
                    # print("WE START FROM HERE", len(reagents_pool), reagents_pool[0], "\n",mol2smi(reagents_pool[0]))
                if len(reagents_pool) > max_results:
                    raise MaxResultsException
                # loop this reaction until either no reagents are left, or no products are generated
                while True:
                    # reaction_products.sealed = True
                    for reagent in reagents_pool:
                        ###### MULTIPLICITY REDUCTION 1 (TO TEST)
                        if not (mol2smi(reagent, False),rxn_pattern) in visited_pairs:
                            visited_pairs.add((mol2smi(reagent,False),rxn_pattern))
                        else:
                            continue
                        self._iterations += 1
                        _iter += 1
                        if _iter > max_iter:
                            raise MaxIterException
                        products = rxn_obj.RunReactants((reagent,))
                        if self.verbose:
                            print("ReactionId[%s]: %d products" %(rxn_idx, len(products)))
                        if len(products):
                            consumed_reagents.add(reagent)
                        for p, *_ in products:
                            self._iterations += 1
                            # TODO add cache check here to check if discarded
                            # molecules are going to be processed here
                            try:
                                Chem.SanitizeMol(p)
                                # skip molecules that have been discarded already
                                # if mol2smi(p,False) in self._discarded:
                                #     continue
                                ###### MULTIPLICITY REDUCTION 2 (TO TEST)
                                if p in visited:
                                    continue
                                visited.add(p)
                                ###### MULTIPLICITY REDUCTION 3 (TO TEST)
                                if not (p,rxn_pattern) in visited_pairs:
                                    visited_pairs.add((mol2smi(p,False),rxn_pattern))
                                else:
                                    continue
                                visited_pairs.add((p, rxn_pattern))
                                ##############################################
                                if preserve_mol_properties:
                                    copy_mol_properties(reagent, p, strict=False, include_name=True)
                                reaction_products.add(PropertyMol(p))
                                if not reaction_log is None:
                                    reaction_log.add(
                                        reagent,
                                        p,
                                        rxn_pattern,
                                        iteration=self._iterations,
                                    )
                            except Chem.AtomValenceException:
                                if self.verbose:
                                    print("[exhaustive enumeration] Product rejected: valence violation")
                            except Chem.KekulizeException:
                                if self.verbose:
                                    print("[exhaustive enumeration] Product rejected: kekule violation")
                            # except Exception as error:
                            #     print("\n\n", error)
                            #     success = False
                    if len(reaction_products):
                        reagents_pool = UniqueMoleculeContainer(reaction_products)
                        reaction_products = UniqueMoleculeContainer()
                    else:
                        reaction_products = reagents_pool
                        break
        except MaxResultsException:
            success = False
            print("WARNING: maximum number of results reached (%d)" % max_results)
        except MaxIterException:
            success = False
            print("WARNING: maximum number of iterations reached (%d)" % _iter)
        if keep_all:
            reagents_pool += consumed_reagents
            reagents_pool += reaction_products
        return reagents_pool, True

    def _parse_reaction_file(self, datafile: str) -> tuple:
        """parse a datafile by stripping comment and empty lines
        that are passed to the function __parse_reaction_line method.
        the line format is the following:
                  *[ ]>>[ ]* [tag]

        the first part must be a valid RDKit SMIRKS reaction or SMARTS
        transformation, with or without spacing between the reaction pattern
        (">>"), followed by a tag.

        The tag can be a float value (the pKa of the protomer transformation), or the name of the transformation (tautomers)

        The reaction is tested before being accepted
        If a failed reaction is found, the code exits

        """
        reactions = []
        # for a future GUI...
        errors = []
        with open(datafile, "r") as fp:
            for idx, line in enumerate(fp.readlines()):
                if line[0] == "#" or not line.strip():
                    continue
                try:
                    rxn_obj, rxn_left, rxn_right, tag = self._parse_reaction_line(line)
                    reactions.append((rxn_obj, rxn_left, rxn_right, tag))
                # TODO convert to reaction failure exception
                except Exception as exc:
                    msg = "ERROR: invalid reaction definition at line [%d]:\n%s" % (
                        idx,
                        line,
                    )
                    errors.append(msg)
        return reactions, errors

    def _parse_reaction_line(self, line):
        """parse a reaction string line
        the line format is the following:
                  *[ ]>>[ ]* [tag]

        the first part must be a valid RDKit SMIRKS reaction or SMARTS
        transformation, with or without spacing between the reaction pattern
        (">>"), followed by a tag.

        The tag can be a float value (the pKa of the protomer transformation),
        or the name of the transformation (tautomers)

        The reaction is tested before being accepted
        If a failed reaction is found, the code exits
        """
        rxn_left, rxn_right = line.split(">>")
        rxn_right, tag = rxn_right.split(None, 1)
        tag = tag.strip()
        rxn_string = "%s >> %s" % (rxn_left, rxn_right)
        rxn_obj = Chem.AllChem.ReactionFromSmarts(rxn_string)
        return rxn_obj, rxn_left, rxn_right, tag


class MaxResultsException(Exception):
    """custom class to manage abrupt interruption of iterations due to maximum
    number of protomers generated"""

    pass


class MaxIterException(Exception):
    """custom class to manage abrupt interruption of iterations due to maximum
    number of protomers generated"""

    pass


def exhaustive_reaction(
    mol,
    reaction_list: list,
    keep_all: bool = True,
    max_results: int = 1000,
    max_iter: int = 10000,
    reaction_log=None,
    preserve_mol_properties:bool=True,
    verbose=False,
) -> tuple:
    """perform the RDKit reaction. The format of the reaction is [(rxn_id,
    reaction)] This method processes the molecular pool (self.mol_pool) but
    does not add results to it, letting the calling method to take care of
    that (allowing for filtering etc.)
    https://github.com/tentrillion/rdkit-tutorials/blob/master/notebooks/003_SMARTS_ReactionsExamples.ipynb

    reaction_list               : list containing the RDKit reaction objects to be applied
    keep_all                    :
    preserve_mol_properties     : copy all molecular properties of the reagents into the generated molecules
    """
    # counter for the max_iter counting
    _iter = 0
    # results = UniqueMoleculeContainer()
    reaction_products = UniqueMoleculeContainer()
    consumed_reagents = UniqueMoleculeContainer()
    visited = UniqueMoleculeContainer(ignore_chirality=True)
    visited_pairs = set()
    for rxn_idx, (rxn_pattern, rxn_obj) in enumerate(reaction_list):
        print("outer loop in standalone exhaustive_reaction: %d %s" % (rxn_idx, rxn_pattern))
        if not keep_all:
            if len(reaction_products):
                reagents_pool = reaction_products.copy()
                reaction_products.clear()
            else:
                reagents_pool = UniqueMoleculeContainer([mol])
            reagents_pool -= consumed_reagents
        else:
            reagents_pool = UniqueMoleculeContainer([mol])
            reagents_pool += reaction_products
            reaction_products.clear()
            reagents_pool += consumed_reagents
            # print("WE START FROM HERE", len(reagents_pool), reagents_pool[0], "\n",mol2smi(reagents_pool[0]))
        if len(reagents_pool) > max_results:
            print("WARNING: maximum number of results reached (%d)" % max_results)
            break
        # loop this reaction until either no reagents are left, or no products are generated
        while True:
            # reaction_products.sealed = True
            for reagent in reagents_pool:
                ###### MULTIPLICITY REDUCTION 1 (TO TEST)
                if not (mol2smi(reagent, False),rxn_pattern) in visited_pairs:
                    visited_pairs.add((mol2smi(reagent,False),rxn_pattern))
                else:
                    continue
                _iter += 1
                if _iter > max_iter:
                    print("WARNING: maximum number of iterations reached (%d)" % _iter)
                    break
                products = rxn_obj.RunReactants((reagent,))
                if verbose:
                    print("ReactionId[%s]: %d products" %(rxn_idx, len(products)))
                if len(products):
                    consumed_reagents.add(reagent)
                for p, *_ in products:
                    # TODO add cache check here to check if discarded
                    # molecules are going to be processed here
                    try:
                        Chem.SanitizeMol(p)
                        # skip molecules that have been discarded already
                        # if mol2smi(p,False) in self._discarded:
                        #     continue
                        ###### MULTIPLICITY REDUCTION 2 (TO TEST)
                        if p in visited:
                            continue
                        visited.add(p)
                        ###### MULTIPLICITY REDUCTION 3 (TO TEST)
                        if not (p,rxn_pattern) in visited_pairs:
                            visited_pairs.add((mol2smi(p,False),rxn_pattern))
                        else:
                            continue
                        visited_pairs.add((p, rxn_pattern))
                        ##############################################
                        if preserve_mol_properties:
                            copy_mol_properties(reagent, p, strict=False, include_name=True)
                        reaction_products.add(PropertyMol(p))
                        if not reaction_log is None:
                            reaction_log.add(
                                reagent,
                                p,
                                rxn_pattern,
                                iteration=self._iterations,
                            )
                    except Chem.AtomValenceException:
                        if verbose:
                            print("[exhaustive enumeration] Product rejected: valence violation")
                    except Chem.KekulizeException:
                        if verbose:
                            print("[exhaustive enumeration] Product rejected: kekule violation")
                    # except Exception as error:
                    #     print("\n\n", error)
                    #     success = False
            if len(reaction_products):
                reagents_pool = UniqueMoleculeContainer(reaction_products)
                reaction_products = UniqueMoleculeContainer()
            else:
                reaction_products = reagents_pool
                break
    if keep_all:
        reagents_pool += consumed_reagents
        reagents_pool += reaction_products
    return reagents_pool

