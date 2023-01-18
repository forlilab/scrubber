from operator import itemgetter

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.EnumerateStereoisomers import (
    EnumerateStereoisomers,
    StereoEnumerationOptions,
    GetStereoisomerCount,
)

from ..common import ScrubberBase, UniqueMoleculeContainer, mol2smi

from .base import MoleculeTransformations
from .base import MaxResultsException
from .base import MaxIterException
from .base import MolecularReactionsLogger
from .base import exhaustive_reaction


# TODO add pro-chiral patterns?

PH_DATAFILE = "phmodel.txt"
TAUTOMERS_DATAFILE = "tautomers.txt"


class MoleculeIsomers(ScrubberBase, MoleculeTransformations):
    """Class to apply chemical transformations:
    - stereoisomers (chirality enumeration)
    - tautomers
    - protomers"""

    def __init__(
        self,
        ## STEREOISOMERS
        stereo_enum: bool = False,
        stereo_max_results: int = 50,
        stereo_gen3d: bool = False,  # MOVE TO STEREO?
        ## PROTOMERS
        proto_enum: bool = True,
        # TODO have two separate options for single pH and pH range?
        proto_pH: float = 7.4,
        proto_max_results: int = 50,
        proto_keep_all: bool = False,
        proto_max_net_charge: int = 5,
        proto_neutralize_only: bool = False,
        ## TAUTOMERS
        tauto_enum: bool = True,
        tauto_max_results: int = 50000,
        # tauto_keep_all: bool = True,
        tauto_protect_aromatic: bool = True,
        tauto_protect_amide: bool = True,
        ## GENERIC
        add_hydrogens: bool = True,
        max_cycles: int = 10,
        # max_iter:50,
        verbose: bool = False,
        ph_datafile: str = None,
        tauto_datafile: str = None,
        suppress_rdkit_warnings: bool = True,
        _stop_at_defaults: bool = False,
    ):
        self.stereo_enum = stereo_enum
        self.stereo_max_results = stereo_max_results
        self.stereo_gen3d = stereo_gen3d
        self.proto_enum = proto_enum
        self.proto_pH = proto_pH
        self.proto_max_results = proto_max_results
        self.proto_keep_all = proto_keep_all
        self.proto_max_net_charge = proto_max_net_charge
        self.proto_neutralize_only = proto_neutralize_only
        self.tauto_enum = tauto_enum
        self.tauto_max_results = tauto_max_results
        # self.tauto_keep_all = tauto_keep_all
        self.tauto_protect_aromatic = tauto_protect_aromatic
        self.tauto_protect_amide = tauto_protect_amide
        self.add_hydrogens = add_hydrogens
        self.max_cycles = max_cycles
        self.verbose = verbose
        self.ph_datafile = ph_datafile
        self.tauto_datafile = tauto_datafile
        self.suppress_rdkit_warnings = suppress_rdkit_warnings
        if _stop_at_defaults:
            return
        valid_isomer = [False, "undefined", "all"]
        if not self.stereo_enum in valid_isomer:
            raise ValueError(
                "steroisomer enumeration can be any of these [%s], not %s"
                % (
                    ", ".join([str(x) for x in valid_isomer]),
                    str(self.stereo_enum),
                )
            )
        MoleculeTransformations.__init__(
            self,
            verbose=self.verbose,
            suppress_rdkit_warnings=self.suppress_rdkit_warnings,
        )
        self.__init_data(
            ph_datafile=self.ph_datafile,
            tauto_datafile=self.tauto_datafile,
        )
        # TODO disabled; figure out how to make it work across neutrality threshold...
        # if isinstance(self.proto_pH, int):
        #     self.proto_pH = [self.proto_pH, self.proto_pH]
        # elif isinstance(self.proto_pH, list):
        #     self.proto_pH = [self.proto_pH[0], self.proto_pH[1]]

    def __init_data(self, ph_datafile: str = None, tauto_datafile: str = None) -> None:
        """initialize data files for transformations"""
        self.__init_tautomers(tauto_datafile)
        self.__init_protomers(ph_datafile)

    def __init_protomers(self, fname: str = None) -> None:
        """initialize the transformations for protomers"""
        self.__ph_model = {}
        if fname is None:
            fname = self.get_datafile(PH_DATAFILE)
        reactions, errors = self._parse_reaction_file(fname)
        for rxn, rxn_left, _, pKa in reactions:
            pKa = float(pKa)
            self.__ph_model[rxn_left] = [rxn, pKa]
            if self.verbose:
                print(
                    "[VERBOSE] PROTOMERS: added model [ %s | pKa %2.2f ]"
                    % (rxn_left, pKa)
                )

    def __init_tautomers(self, fname: str = None) -> None:
        """initialize the transformations for tautomers"""
        self.__tauto_model = {}
        if fname is None:
            fname = self.get_datafile(TAUTOMERS_DATAFILE)
        reactions, errors = self._parse_reaction_file(fname)
        for rxn, _, _, name in reactions:
            self.__tauto_model[name] = [rxn]
            if self.verbose:
                print("[VERBOSE] TAUTOMERS: added model [%s] [ %s ]" % (name, rxn))

    def process(
        self,
        mol,
    ) -> list:
        """apply all transformations"""
        success_record = []
        self.reaction_log = MolecularReactionsLogger(self.verbose)
        self.mol_pool = UniqueMoleculeContainer([mol])
        # disposable set of molecules that have been rejected
        # TODO: do we need to keep track of this? Use a UniqueMolContainer?
        # TODO: use this to save molecules?
        self._discarded = set()
        ##### protomers/tautomers loop
        self.__process_iter = 0
        self._iterations = 0
        while True:
            if (not self._iterations == 0) and (self.mol_pool.sealed == True):
                break
            if self.proto_enum:
                if self.verbose:
                    print(
                        "[VERBOSE] Calling protomer generation, round:",
                        self.__process_iter,
                    )
                success = self.enumerate_protomers(
                    self.proto_pH,
                    self.proto_max_results,
                    self.proto_keep_all,
                    self.proto_max_net_charge,
                    self.proto_neutralize_only,
                )
                # store the status of the last operation
                success_record.append((success, self._iterations, "protomer"))
                if self.verbose:
                    print(
                        "protomer enum [iter:%d] success: %s (%d results)"
                        % (self.__process_iter, success, len(self.mol_pool))
                    )
            # seal the container to track for modifications
            self.mol_pool.sealed = True
            if self.tauto_enum:
                if self.verbose:
                    print(
                        "[VERBOSE] Calling tautomer generation round:",
                        self.__process_iter,
                    )
                success = self.enumerate_tautomers(
                    self.tauto_max_results,
                    self.tauto_protect_aromatic,
                    self.tauto_protect_amide,
                )
                success_record.append((success, self._iterations, "tautomer"))
                if self.verbose:
                    print(
                        "[VERBOSE] Tautomer enumeration [iter:%d] success: %s (%d results)"
                        % (self.__process_iter, success, len(self.mol_pool))
                    )
            self.__process_iter += 1
            if self.__process_iter > self.max_cycles:
                if self.verbose:
                    print("[VERBOSE] Seal not broken, no more things to do")
                # TODO is this a failure? (arguably, yes)
                success_record.append(
                    (False, self._iterations, "proto/tauto loop interrupted")
                )
                break
            if (not self.proto_enum) or (not self.tauto_enum):
                if self.verbose:
                    print("[VERBOSE] either tautomers or protomers are not requested")
                break
            elif self.mol_pool.sealed:
                if self.verbose:
                    print("[VERBOSE] MolPool container is still sealed (done)")
                break
            else:
                self.mol_pool.sealed = True
        ##### isomers
        if self.stereo_enum is not False:
            success = self.enumerate_stereoisormers(
                mode=self.stereo_enum,
                max_results=self.stereo_max_results,
                gen3d=self.stereo_gen3d,
            )
            success_record.append((success, self._iterations, "stereoisomer"))
            if self.verbose:
                print(
                    "[VERBOSE] Enantiomers enumeration success",
                    success,
                    len(self.mol_pool),
                )
        if self.add_hydrogens:
            # add back hydrogens here
            for old_mol in self.mol_pool:
                new_mol = Chem.AddHs(old_mol)
                self.mol_pool.add(new_mol, replace=True)
            # post = set(list(self.mol_pool._UniqueMoleculeContainer__data.keys()))
        # TODO return a single boolean
        if self.verbose:
            print("SUCCESSS", success_record)
        del self._discarded
        return success_record

    def enumerate_stereoisormers(self, mode=None, max_results=100, gen3d=True) -> bool:
        """
        source https://www.rdkit.org/docs/source/rdkit.Chem.EnumerateStereoisomers.html
        """
        # TODO add random enumeration
        unassigned_only = True
        if mode == "all":
            unassigned_only = False
        opts = StereoEnumerationOptions(
            unique=True,
            tryEmbedding=gen3d,
            onlyUnassigned=unassigned_only,
            maxIsomers=max_results,
        )
        results = UniqueMoleculeContainer()
        # process results and register the information of this transformation
        for m in self.mol_pool:
            if self.verbose:
                print("[VERBOSE] processing molecule", mol2smi(m))
            for ent in EnumerateStereoisomers(m, options=opts):
                self._iterations += 1
                results.add(ent)
                self.reaction_log.add(
                    reagent_mol=m,
                    product_mol=ent,
                    transformation="stereoisomer-gen",
                    iteration=self._iterations,
                    kept=True,
                )
        if len(results):
            # updating counters for naming
            for mol in results:
                self._set_scrubber_property(mol, "Scrubber_stereo_count")
            self.mol_pool = results
            if self.verbose:
                print("[VERBOSE] %d enantiomers" % len(results))
            return True
        return False

    def enumerate_protomers(
        self,
        pH: float = 7.4,  # range, too?
        max_results: int = 50,
        keep_all: bool = False,
        max_net_charge: int = None,
        neutralize: bool = False,
    ) -> bool:
        """enumerate protomer species for the molecule(s) present in the
        molecules pool (self.mol_pool) at the input pH.  If the process does
        not create any new molecules, then the input will be retained. By
        default, only the final products after exhaustive transformations will
        be retained (i.e.: NH2-X-CCOH -> (NH3+)-X-(COO-) ), but if the
        'keep_all' option is used, then all intermediate states will be kept (
        i.e.: NH2-X-COOH -> (NH3+)-X-COOH, NH2-X-(COO-), (NH3+)-X-(COO-) )

        pH              :   upper bound (inclusive) of the pH at which the tranformation
                            will be applied

        max_results     :   maximum number of results to be generated; if this
                            number is exceeded, the enumeration is stopped

        unique          :   remove duplicate results (for symmetric molecules)

        keep_all        :   by default, only the final results of exhaustive
                            transformations are returned; if keep_all is used, also
                            intermediate protonation states are returned

        max_net_charge  :   maximum absolute net charge accepted; transformed
                            results that generate species with absolute formal charges
                            higher than this value will be rejected. NOTE this filter
                            is applied only after all results have been enumerated, to
                            prevent the generation of zwitterions with a lower formal
                            charge than the intermediate protonation states

        neutralize      :   try neutralizing groups that can be transformed in
                            their neutral form. NOTE: when this option is used all other
                            options are ignored
        """
        if neutralize:
            container = UniqueMoleculeContainer()
            for mol in self._mol_pool:
                container.add(self._neutralize_atoms(mol))
            # if no molecule is truly modified, just return
            if container == self.mol_pool:
                return
            self.mol_pool = container
            return
        # find which transformations apply
        proto_reactions = []
        for pattern, (rxn, rxn_pH) in self.__ph_model.items():
            if rxn_pH >= pH:
                continue
            # TODO fix after figuring out how to deal with neutrality threshold
            # if pH[0] == pH[1]:
            #     if rxn_pH >= pH[1]:
            #         continue
            # else:
            #     if pH[0] <= rxn_pH >= pH[1]:
            #         continue
            proto_reactions.append((pattern, rxn, rxn_pH))
        # sort transformations by decreasing pKa (TODO: maybe useful for the future)
        proto_reactions = sorted(proto_reactions, key=itemgetter(2))
        proto_reactions = [(x[0], x[1]) for x in proto_reactions]
        # enumerate all transformations
        results, success = self._exhaustive_reaction(
            proto_reactions,
            keep_all,
            max_results,
            reaction_log=self.reaction_log,
        )
        if len(results) > 0:
            # filter by net charge
            if not max_net_charge is None:
                pre_count = len(results)
                remove = []
                for m in results:
                    if abs(Chem.GetFormalCharge(m)) > max_net_charge:
                        remove.append(m)
                self._discarded.update(remove)
                # for r in remove:
                #     print(">>DISCARDED:", mol2smi(r))
                results.remove_mols(remove)
            if self.verbose:
                print(
                    "[VERBOSE] %d results generated (success:%s)"
                    % (len(results), success)
                )
                if not pre_count == len(results):
                    print(
                        "[VERBOSE] %d results filtered out by charge criterion"
                        % (pre_count - len(results))
                    )
                print("[VERBOSE]", results)
                for mol in results:
                    self._set_scrubber_property(mol, "Scrubber_proto_count")
            if keep_all:
                self.mol_pool += results
            else:
                self.mol_pool = results
        return success

    def enumerate_tautomers(
        self,
        max_results: int = 50,
        # keep_all: int = True, # TODO all tautomers must be kept?
        protect_aromatic: bool = True,
        protect_amide: bool = True,
    ) -> bool:
        """enumerate tautomers"""
        # print("CALLED WITH", protect_aromatic, protect_amide )

        def __check_property_violation(_type, pool1, pool2):
            """helper function to check for violations"""
            if _type == "aromatic":
                check_func = self.__count_aromatics
                msg = ("*** REJECTED *** aromatic decrease",)
            elif _type == "amide":
                check_func = self.__count_amides
                msg = ("*** REJECTED *** amide decrease",)
            count_list1 = [check_func(x) for x in pool1]
            count_list2 = [check_func(x) for x in pool2]
            joined = count_list1 + count_list2
            if len(joined) == 0:
                print("JOINED EMPTY?", joined)
                # return pool
            max_type_value = max(count_list1 + count_list2)
            for pool, count_list in ((pool1, count_list1), (pool2, count_list2)):
                delete = []
                for idx, count in enumerate(count_list):
                    if count < max_type_value:
                        delete.append(idx)
                        self._iterations += 1
                        self.reaction_log.add(
                            pool[idx],
                            None,
                            "*** REJECTED *** aromatic decrease",
                            self._iterations,
                            kept=False,
                        )
                        if self.verbose:
                            print(
                                "[VERBOSE] molecule [%s] rejected by %s count (current: %d | max %d)"
                                % (_type, mol2smi(results[idx]), count, max_type_value)
                            )
                self._discarded.update(pool.remove_from_indices(delete))
            # TODO check if this is even necessary, we modify objects in place...
            return pool1, pool2

        reactions = [(x[0], x[1][0]) for x in self.__tauto_model.items()]
        results, success = self._exhaustive_reaction(
            reaction_list=reactions,
            keep_all=True,
            max_results=max_results,
            reaction_log=self.reaction_log,
            preserve_mol_properties=True,
        )
        if self.verbose:
            print(
                "[VERBOSE] %d tautomers generated (success: %s)"
                % (len(results), success)
            )
        if len(results):
            if protect_aromatic:
                results, self.mol_pool = __check_property_violation(
                    "aromatic", results, self.mol_pool
                )
            else:
                print("WARNING: aromatic protection disabled.")
            if protect_amide:
                results, self.mol_pool = __check_property_violation(
                    "amide", results, self.mol_pool
                )
            else:
                print("WARNING: amide protection disabled.")
        for mol in results:
            self._set_scrubber_property(mol, "Scrubber_tauto_count")
        if len(results):
            self.mol_pool += results
        return success

    def __count_aromatics(self, mol):
        """count the number of aromatic atoms in a molecule"""
        return sum([atom.GetIsAromatic() for atom in mol.GetAtoms()])

    def __count_amides(self, mol):
        """count the number of amides in a molecule"""
        # does NOT match 2-Pyridone (intentionally)
        pattern = Chem.MolFromSmarts("[OX1,SX1]=[CX3][NX3]")
        return len(mol.GetSubstructMatches(pattern))

    def _neutralize_atoms(self, mol):
        """Neutralie charged molecules by atom
        source: http://www.rdkit.org/docs/Cookbook.html#neutralizing-molecules
        """
        pattern = Chem.MolFromSmarts(
            "[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]"
        )
        at_matches = mol.GetSubstructMatches(pattern)
        at_matches_list = [y[0] for y in at_matches]
        if len(at_matches_list) > 0:
            for at_idx in at_matches_list:
                atom = mol.GetAtomWithIdx(at_idx)
                chg = atom.GetFormalCharge()
                hcount = atom.GetTotalNumHs()
                atom.SetFormalCharge(0)
                atom.SetNumExplicitHs(hcount - chg)
                atom.UpdatePropertyCache()
        return mol

    def _set_scrubber_property(self, mol, property_name):
        """update scrubber incremental value properties, e.g.:
        stereoisomer/tautomer count, etc."""
        if not mol.HasProp(property_name):
            mol.SetIntProp(property_name, 0)
            current = 0
        else:
            current = mol.GetIntProp(property_name)
        current += 1
        # print("UPDATING CRRENT!", current, property_name)
        mol.SetIntProp(property_name, current)


if __name__ == "__main__":
    import sys

    verbose = False
    suppress_rdkit_warnings = True
    isomer_enum = False
    proto_enum = False
    tauto_enum = False
    proto_keep_all = False
    if "-v" in sys.argv:
        verbose = True
        suppress_rdkit_warnings = False
    if "-i" in sys.argv:
        isomer_enum = "undefined"
    elif "-I" in sys.argv:
        isomer_enum = "all"
    if "-p" in sys.argv:
        proto_enum = True
        proto_keep_all = False
    if "-P" in sys.argv:
        proto_enum = True
        proto_keep_all = True
    if "-t" in sys.argv:
        tauto_enum = True

    # try:
    #     infile = sys.argv[1]
    # except:
    #     infile = "tests/noreact.sdf"
    infile = sys.argv[1]
    try:
        ph_datafile = sys.argv[2]
    except:
        ph_datafile = "data/test_model.txt"
    try:
        tauto_datafile = sys.argv[3]
    except:
        tauto_datafile = "data/tautomers.txt"

    mt = MoleculeIsomers(
        verbose=verbose,
        ph_datafile=ph_datafile,
        tauto_datafile=tauto_datafile,
        suppress_rdkit_warnings=suppress_rdkit_warnings,
    )
    # parse molecule
    with Chem.SDMolSupplier(infile, removeHs=True) as supp:
        for mol in supp:
            print("\n===============================================")
            print("Processing [ %s ]" % mol.GetProp("_Name"))
            print(">SMI:", mol2smi(mol))
            print("-----------------------------------------------")
            outcome = mt.process(
                mol,
                stereo_enum=isomer_enum,
                stereo_max_results=50,
                proto_enum=proto_enum,
                proto_pH=7.4,
                proto_max_results=50,
                proto_keep_all=proto_keep_all,
                proto_max_net_charge=5,
                tauto_enum=tauto_enum,
                tauto_max_results=500,
                # tauto_keep_all=True,
                tauto_protect_aromatic=True,
                tauto_protect_amide=True,
                ## GENERIC
                add_hydrogens=True,
                gen3d=False,
            )
            #
            #
            # DEBUG/REPORT
            #
            #
            graph = mt.reaction_log
            traj = graph.trajectory
            print("\n===============================================")
            print(
                "Accepted products: %d ( %d total steps )"
                % (len(mt.mol_pool), mt._iterations)
            )
            print("-----+--------------------------------------------")
            print("name | SMILES")
            print("-----+--------------------------------------------")
            for mol in mt.mol_pool:
                smi, label = graph.get_info(mol)
                print("  %s  | %s " % (label, smi))
            print("-----------------------------------------------")
        if verbose:
            print(" Success |   operation   | iteration")
            print("---------------------------------------")
            for s, i, o in outcome:
                print("  %5s  | %13s | %d" % (str(s), o, i))
            ######################
            print("steps registered", graph._counter)
            graph.print_history()

            # print("\n\n===============================================")
            # print("Reactions trajectory")
            # print("===============================================")
            # format_line = ("%%5d - %%%ds ( %%2s ) [ %%3d ]  \t   "
            #         "--> %%%ds ( %%2s ) [ %%3d ]  \t  |  %%s " % (
            #     graph._largest_string,
            #     graph._largest_string,
            # ))
            # for idx, (iter_step, result) in enumerate(traj.items()):
            #     print(
            #         format_line
            #         % (
            #             iter_step,
            #             result[0],
            #             result[1],
            #             result[2],
            #             result[3],
            #             result[4],
            #             result[5],
            #             result[6],
            #         ),
            #     )
            # print("steps printed", idx + 1)


def enumerate_tautomers(
    mol,
    reactions,
    max_results: int = 50,
    # keep_all: int = True, # TODO all tautomers must be kept?
    protect_aromatic: bool = True,
    protect_amide: bool = True,
    verbose=False,
) -> list:
    """enumerate tautomers"""
    # print("CALLED WITH", protect_aromatic, protect_amide )

    def check_property_violation(pool, querymol):
        """helper function to check for violations"""
        count_list = [len(mol.GetSubstructMatches(querymol)) for mol in pool]
        max_count = max(count_list)
        delete = []
        for idx, count in enumerate(count_list):
            if count < max_count:
                delete.append(idx)
                #self.reaction_log.add(
                #    pool[idx],
                #    None,
                #    "*** REJECTED *** aromatic decrease",
                #    self._iterations,
                #    kept=False,
                #)
                if verbose:
                    print(
                        "[VERBOSE] molecule [%s] rejected by %s count (current: %d | max %d)"
                        % (mol2smi(results[idx]), Chem.MolToSmiles(querymol), count, max_count)
                    )
        pool.remove_from_indices(delete)

    results = exhaustive_reaction(
        mol,
        reaction_list=reactions,
        keep_all=True,
        max_results=max_results,
        reaction_log=None,
        preserve_mol_properties=True,
        verbose=verbose,
    )
    if verbose:
        print("%d tautomers generated" % len(results))
    check_property_violation(results, Chem.MolFromSmarts("[a]"))
    check_property_violation(results, Chem.MolFromSmarts("[OX1,SX1]=[CX3][NX3]"))
    return results
