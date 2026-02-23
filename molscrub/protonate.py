import pathlib
from importlib.resources import files
from .common import UniqueMoleculeContainer
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdChemReactions
from rdkit.Chem import rdFMCS
import numpy as np
import joblib
import pandas as pd
from rdkit.Chem import Descriptors

datapath = files("molscrub") / "data"
default_tautomers_fn = datapath / "tautomers.txt"
default_pka_reactions_fn = datapath / "pka_reactions.txt"

etr1_path = datapath  / "ETR_rules_rdkit_desc.joblib"

class AcidBaseConjugator:
    def __init__(self, pka_reactions, pka_model="rules"):
        self.pka_reactions = pka_reactions
        self.pka_model = pka_model

    def __call__(self, input_mol, ph_range_low, ph_range_high):
        if ph_range_low > ph_range_high:
            raise ValueError("ph_range_low must be lesser than or equal to ph_range_high")

        # protonate according to selected model
        if self.pka_model == "rules":
            return self.protonate_with_rules(input_mol, ph_range_low, ph_range_high)
        elif self.pka_model == "etr1":
            return self.protonate_with_model(input_mol, ph_range_low, ph_range_high)
        else:
            raise ValueError("Unknown pKa model given. I don't know what to do.")


    def protonate_with_model(self, input_mol, ph_range_low, ph_range_high):
        """
        return appropriate protonated mols with the ML model. 
        Currently only one model is supported, but this may be expanded
        in the future. 
        """

        rxn_info = self.get_rxn_info(input_mol)

        if len(rxn_info) == 0:
            # no rxn happened based on the rules
            return [input_mol]
        
        mol_list = []

        for rxn in rxn_info:
            ml_pka = self.calculate_pka(input_mol, rxn, model="etr1")[0]

            if ph_range_high < ml_pka:
                if rxn["direction"] == "gain_h":
                    mol_list.append(rxn["product"])
                else: 
                    mol_list.append(input_mol)

            elif ph_range_low > ml_pka:
                if rxn["direction"] == "lose_h":
                    mol_list.append(rxn["product"])
                else: 
                    mol_list.append(input_mol)
            else: 
                # keep both
                mol_list.append(input_mol)
                mol_list.append(rxn["product"])

        return mol_list


    def protonate_with_rules(self,input_mol, ph_range_low, ph_range_high):
        """
        return appropriate protonated mols based on the smarts rxn given in the
        default pka_reactions file or a custom file provided by the user. 
        """

        mol_list = [input_mol]
        for r in self.pka_reactions:
            if ph_range_high < r["pka"]:
                mol_list = [convert_exhaustive(mol, r["rxn_gain_h"]) for mol in mol_list]
            elif ph_range_low > r["pka"]:
                mol_list = [convert_exhaustive(mol, r["rxn_lose_h"]) for mol in mol_list]
            else: # keep both states for each transformation
                tmp = UniqueMoleculeContainer()
                for mol in mol_list:
                    tmp.add(mol)
                    convert_recursive(mol, r["rxn_gain_h"], tmp)
                    convert_recursive(mol, r["rxn_lose_h"], tmp)
                mol_list = [mol for mol in tmp]
        for mol in mol_list:
            copy_mol_props(input_mol, mol)
        return mol_list

    def mol_comparisons(self, mol1, mol2):
        """
        Check if 2 mols are the same, using inchi keys
        """

        inchi1 = Chem.MolToInchiKey(mol1)
        inchi2 = Chem.MolToInchiKey(mol2)

        return inchi1 == inchi2


    def _one_hot(self, size, index):
        return np.eye(size)[index]

    def find_protonation_site_with_mcs(self, original: Chem.Mol, reacted: Chem.Mol):
        """
        Finds the atom index in `original` that changed protonation state
        by computing maximal structural overlap (MCS) and building
        an optimal atom mapping.
        
        Returns:
            original atom index that changed
            or None if no unique site found
        """

        # Compute MCS
        mcs = rdFMCS.FindMCS(
            [original, reacted],
            bondCompare=rdFMCS.BondCompare.CompareOrder,
            atomCompare=rdFMCS.AtomCompare.CompareElements,
            ringMatchesRingOnly=True,
            completeRingsOnly=True,
        )

        if mcs.canceled or mcs.numAtoms == 0:
            return None

        mcs_mol = Chem.MolFromSmarts(mcs.smartsString)

        # Get atom matches
        match_orig = original.GetSubstructMatch(mcs_mol)
        match_react = reacted.GetSubstructMatch(mcs_mol)

        if not match_orig or not match_react:
            return None

        # Build mapping: reacted_idx to original_idx
        react_to_orig = dict(zip(match_react, match_orig))

        # Detect changed atom
        changed_atoms = []

        for r_idx, o_idx in react_to_orig.items():

            atom_orig = original.GetAtomWithIdx(o_idx)
            atom_react = reacted.GetAtomWithIdx(r_idx)

            if (
                atom_orig.GetFormalCharge() != atom_react.GetFormalCharge()
                or atom_orig.GetTotalNumHs() != atom_react.GetTotalNumHs()
            ):
                changed_atoms.append(o_idx)

        if len(changed_atoms) == 1:
            return changed_atoms[0]

        # If multiple changed atoms, return None 
        return None

    def get_rxn_info(self, mol) -> list:
        """
        this function runs a mol object through all possible reactions and
        returns information in the form of a dictionary on the reactions that succeed. 

        dictionary keys:
        ----------------
        product: the product

        rule_pka: the pka of the reaction as given by the pka_reactions.txt file

        rxn_name: name of the reaction as given in the pka_reactions.txt file

        rxn: the actual reaction in rdkit format

        direction: whether the reaction that proceeds is the forward or backward reaction, as given in the file
        
        protonated_atom: the index of the atom that is protonated by the pka rxn. 
        
        rxn_1hot_encoding: array of 1-hot encodings indicating which reaction succeeded. 
        """


        size = len(self.pka_reactions)
        reacted_mols = []


        for i,r in enumerate(self.pka_reactions):
            temp_forward = convert_exhaustive(mol, r["rxn_gain_h"])
            temp_forward = Chem.MolFromSmiles(Chem.MolToSmiles(temp_forward))
            if not self.mol_comparisons(mol, temp_forward):
                changed_atom = self.find_protonation_site_with_mcs(mol, temp_forward)

                reacted_mols.append({"product": temp_forward, 
                                     "rule_pka": r["pka"], 
                                     "rxn_name": r["name"], 
                                     "rxn": r["rxn_gain_h"], 
                                     "direction": "gain_h", 
                                     "protonated_atom": changed_atom, 
                                     "rxn_1hot_encoding": self._one_hot(size, i)})

            temp_backward = convert_exhaustive(mol, r["rxn_lose_h"])
            temp_backward = Chem.MolFromSmiles(Chem.MolToSmiles(temp_backward))
            if not self.mol_comparisons(mol, temp_backward):
                changed_atom = self.find_protonation_site_with_mcs(mol, temp_backward)
                reacted_mols.append({"product":temp_backward, 
                                     "rule_pka": r["pka"], 
                                     "rxn_name":r["name"], 
                                     "smarts_rxn": r["rxn_lose_h"], 
                                     "direction": "lose_h", 
                                     "protonated_atom": changed_atom, 
                                     "rxn_1hot_encoding": self._one_hot(size, i)})
        
        return reacted_mols

    def calculate_pka(self, mol, rxn_info, model="etr1"):
        """calculates pkas from the given model

        for now only one model is supported, but this may be updated in the future
        """

        if model == "etr1":
            model = joblib.load(etr1_path)
            x = self._prepare_data_for_model(mol, rxn_info)
        else:
            raise ValueError("invalid model chosen for pka prediction.")

        pka = model.predict(x)
        return pka
        

    def _prepare_data_for_model(self, mol, rxn_info):
        """
        prepares data for input into the ml_model, using pandas
        """

        # rules
        df = pd.DataFrame([rxn_info])
        df = df[["rule_pka", "rxn_1hot_encoding"]]


        # expand 1hot encoding
        expanded_cols = df['rxn_1hot_encoding'].apply(pd.Series)
        expanded_cols.columns = [f'encoding{i+1}' for i in range(expanded_cols.shape[1])] 

        # rdkit descriptors
        descriptors = self.getMolDescriptors(mol)
        desc_df = pd.DataFrame([descriptors])

        x = pd.concat((desc_df.reset_index(drop=True), expanded_cols.reset_index(drop=True), df["rule_pka"].reset_index(drop=True)), axis=1)


        # this is needed because that's how the model names it. 
        x.rename({"rule_pka":"base_pka"}, axis=1, inplace=True)

        return x


    def getMolDescriptors(self, mol, missingVal=None):
        ''' calculate the full list of descriptors for a molecule
        
            missingVal is used if the descriptor cannot be calculated
        '''
        res = {}
        for nm,fn in Descriptors._descList:
            # some of the descriptor fucntions can throw errors if they fail, catch those here:
            try:
                val = fn(mol)
            except:
                # print the error message:
                import traceback
                traceback.print_exc()
                # and set the descriptor value to whatever missingVal is
                val = missingVal
            res[nm] = val
        return res

    @classmethod 
    def from_default_data_files(cls, model="rules"):
        pka_reactions = cls.parse_reaction_file(default_pka_reactions_fn)
        return cls(pka_reactions, pka_model = model)

    @staticmethod
    def parse_reaction_file(datafile: str) -> list:
        """the line format is the following:
               SMARTS <<>> SMARTS NAME PKA_VALUE
           the space between SMARTS and <<>> is optional
        """
        reactions = []
        name_set = set()
        with open(datafile, "r") as fp:
            for line in fp:
                line = line.strip()
                if len(line) == 0 or line[0] == "#":
                    continue
                rxn_left, rxn_right = line.split("<<>>")
                rxn_right, name, pka = rxn_right.split()
                if name in name_set:
                    raise ValueError("reaction name must be unique") 
                name_set.add(name)
                r = {}
                r["name"] = name
                r["pka"] = float(pka)
                r["rxn_lose_h"] = AllChem.ReactionFromSmarts("%s >> %s" % (rxn_left, rxn_right))
                r["rxn_gain_h"] = AllChem.ReactionFromSmarts("%s >> %s" % (rxn_right, rxn_left))
                reactions.append(r)
        return reactions


class Tautomerizer:
    def __init__(self, reactions, keepmax_smarts, nr_rounds=2):
        self.reactions = reactions
        self.keepmax_smarts = keepmax_smarts
        self.nr_rounds = nr_rounds

    @classmethod
    def from_default_data_files(cls):
        reactions, keepmax_smarts = cls.parse_tautomers_config_file(default_tautomers_fn)
        return cls(reactions, keepmax_smarts)

    @classmethod
    def from_reactions_filename(cls, filename, nr_rounds=2):
        reactions, keepmax_smarts = cls.parse_tautomers_config_file(filename)
        return cls(reactions, keepmax_smarts, nr_rounds)

    def __call__(self, input_mol):
        tautomers = UniqueMoleculeContainer([input_mol])
        for roundid in range(self.nr_rounds):
            tmp = UniqueMoleculeContainer()
            for mol in tautomers:
                for r in self.reactions:
                    uniq = UniqueMoleculeContainer()
                    rxn = r['rxn']
                    products = react_and_sanitize(mol, rxn)
                    for product in products:
                        tmp.add(product)
            for mol in tmp:
                tautomers.add(mol)

        # count occurences of each SMARTS
        smarts_count = [[0]*len(tautomers) for _ in self.keepmax_smarts]
        for i, smarts in enumerate(self.keepmax_smarts):
            smarts_mol = Chem.MolFromSmarts(smarts["smarts"])
            for j, mol in enumerate(tautomers):
                smarts_count[i][j] = len(mol.GetSubstructMatches(smarts_mol))

        # select tautomers that have the max count of each SMARTS
        best_of_all_counts = False # fewer tautomers if set to True
        is_selected = [True] * len(tautomers)
        for index in range(len(self.keepmax_smarts)):
            fn = self.keepmax_smarts[index]["fn"]
            count = smarts_count[index]
            current_best = fn([count[i] for i in range(len(tautomers)) if is_selected[i] or best_of_all_counts])
            for j in range(len(tautomers)):
                is_selected[j] = is_selected[j] and (count[j] == current_best)
        
        output = [tautomers[j] for j in range(len(tautomers)) if is_selected[j]]
        for mol in output:
            copy_mol_props(input_mol, mol)
        return output

    @staticmethod
    def parse_tautomers_config_file(fname):
        reactions = []
        keepmax_smarts = []
        with open(fname) as f:
            for line in f:
                line = line.strip()
                if len(line) == 0 or line[0] == "#":
                    continue
                if line.startswith("KEEPMAX_SMARTS") or line.startswith("KEEPMIN_SMARTS"):
                    _, smarts, name = line.split()
                    fn = max if line[4:7] == "MAX" else min
                    keepmax_smarts.append({"smarts": smarts, "name": name, "fn": fn})
                else:
                    smirks, name = line.split()
                    reactions.append({"rxn": rdChemReactions.ReactionFromSmarts(smirks), "name": name, "smarts": smirks})
        return reactions, keepmax_smarts
            

def react_and_sanitize(mol, rxn):
    nr_react = rxn.GetNumReactantTemplates()
    nr_prod = rxn.GetNumProductTemplates()
    #regenerate 
    mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol)) # fixes second round tautomerization bugs
    if nr_react != 1 or nr_prod != 1:
        raise RuntimeError("reaction must be single reactant -> single product")
    output_products = []
    products = rxn.RunReactants((mol,))
    for product in products:
        product = product[0] # nr products == NumProductTemplates == 1
        try:
            s = Chem.SanitizeMol(product)
            #product.UpdatePropertyCache()
            # loading a fresh molecule detects errors that updating the property cache doesn't
            product = Chem.MolFromSmiles(Chem.MolToSmiles(product))
            if product is None:
                continue
        except Chem.AtomValenceException as e:
            continue
        except Chem.KekulizeException as e:
            continue
        except Exception as e:
            print("uncought exception", e, type(e))
            continue
        output_products.append(product)
    return output_products


def convert_recursive(mol, rxn, container):
    for product in react_and_sanitize(mol, rxn):
        container.add(product)
        convert_recursive(product, rxn, container)


def convert_exhaustive(mol, rxn):
    """
        Returns exactly one molecule:
            - the product when the reaction occurs and sanitization succeeds,
            - the input mol otherwise.
        The returned product will have all substructures reacted.
    """
    nr_react = rxn.GetNumReactantTemplates()
    nr_prod = rxn.GetNumProductTemplates()
    if nr_react != 1 or nr_prod != 1:
        raise RuntimeError("reaction must be single reactant -> single product" )
    # in case of multiple reactive substructures, maxProducts=1 returns a product with
    # only one reacted substructure. If maxProducts was not set to 1, we would get a tuple
    # of products, and each would have one (different) reacted substructure.
    # We want to take any of these products and make it react again, recursively, to
    # produce a single product will all matching substructures reacted. 
    products_list = rxn.RunReactants((mol,), maxProducts=1)
    if len(products_list) == 0:
        return mol
    elif len(products_list) == 1:
        products = products_list[0]
        product = products[0] # nr products == NumProductTemplates == 1
        try:
            Chem.SanitizeMol(product)
        # the following exceptions arise often with Chem.SanitizeMol
        except Chem.AtomValenceException:
            return mol
        except Chem.KekulizeException:
            return mol
        # run reaction again with the product to return 
        return convert_exhaustive(product, rxn)
    else:
        raise RuntimeError("RunReactants(maxProducts=1) returned %d products" % len(products_list))


def enumerate_pka_fewer_combos(mol, pka_reactions, ph_range_low, ph_range_high):
    """ if the pka of a reaction is between ph_range_low and pk_range_high, and the
        molecule has multiple substructures that are affected by the reaction,
        no states will be returned with a subset of the substructures reacted.
        There will be states in which all the substructures reacted, and states
        in which zero substructures reacted. This limitation is for combinations of
        states within each reaction, not across reactions. 
    """
    if ph_range_low > ph_range_high:
        raise ValueError("ph_range_low must be lesser than or equal to ph_range_high")
    mol_set = [mol]
    for r in pka_reactions:
        # tmp is set() and not list() because when both gain_H and lose_H are
        # applied but don't react, the input molecule will be returned twice.
        # No need to use smiles to check if the molecule is the same because
        # the molecule object identity is the same.
        tmp = set()
        if ph_range_low <= r["pka"]:
            for mol in mol_set:
                tmp.add(convert_exhaustive(mol, r["rxn_gain_h"]))
        if ph_range_high >= r["pka"]:
            for mol in mol_set:
                tmp.add(convert_exhaustive(mol, r["rxn_lose_h"]))
        mol_set = [mol for mol in tmp]
    return mol_set

def copy_mol_props(original_mol, target_mol):
    if original_mol.HasProp("_Name"):
        target_mol.SetProp("_Name", original_mol.GetProp("_Name"))
    for prop_name in original_mol.GetPropNames():
        target_mol.SetProp(prop_name, original_mol.GetProp(prop_name))
