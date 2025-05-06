import pathlib
from importlib.resources import files
from .common import UniqueMoleculeContainer
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdChemReactions

datapath = files("molscrub") / "data"
default_tautomers_fn = datapath / "tautomers.txt"
default_pka_reactions_fn = datapath / "pka_reactions.txt"

class AcidBaseConjugator:
    def __init__(self, pka_reactions):
        self.pka_reactions = pka_reactions

    def __call__(self, input_mol, ph_range_low, ph_range_high):
        if ph_range_low > ph_range_high:
            raise ValueError("ph_range_low must be lesser than or equal to ph_range_high")
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

    @classmethod 
    def from_default_data_files(cls):
        pka_reactions = cls.parse_reaction_file(default_pka_reactions_fn)
        return cls(pka_reactions)

    @classmethod
    def from_reactions_filename(cls, fname):
        pka_reactions = cls.parse_reaction_file(fname)
        return cls(pka_reactions) 

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
                    products = react_and_sanitize(mol, r["rxn"])
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
                    reactions.append({"rxn": rdChemReactions.ReactionFromSmarts(smirks), "name": name})
        return reactions, keepmax_smarts
            

def react_and_sanitize(mol, rxn):
    nr_react = rxn.GetNumReactantTemplates()
    nr_prod = rxn.GetNumProductTemplates()
    if nr_react != 1 or nr_prod != 1:
        raise RuntimeError("reaction %s must be single reactant -> single product" % name)
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
        raise RuntimeError("reaction %s must be single reactant -> single product" % name)
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
