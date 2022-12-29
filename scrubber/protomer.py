from .common import UniqueMoleculeContainer
from rdkit import Chem
from rdkit.Chem import AllChem

def parse_reaction_file(datafile: str) -> list:
    """parse a datafile by stripping comment and empty lines
    that are passed to the function __parse_reaction_line method.
    the line format is the following:
              *[ ]>>[ ]* [tag]

    the first part must be a valid RDKit SMIRKS reaction or SMARTS
    transformation, with or without spacing between the reaction pattern
    (">>"), followed by a tag.

    The tag can be a float value (the pKa of the protomer transformation), or the name of the transformation (tautomers)

    """
    reactions = []
    # for a future GUI...
    with open(datafile, "r") as fp:
        for idx, line in enumerate(fp.readlines()):
            if line[0] == "#" or not line.strip():
                continue
            rxn_left, rxn_right = line.split(">>")
            rxn_right, tag = rxn_right.split(maxsplit=1)
            tag = tag.strip()
            #rxn_string = "%s >> %s" % (rxn_left, rxn_right)
            #rxn_obj = Chem.AllChem.ReactionFromSmarts(rxn_string)
            #reactions.append((rxn_obj, rxn_left, rxn_right, tag))
            reactions.append((rxn_left, rxn_right, tag))
    return reactions

def build_pka_reactions(reactions):
    pka_reactions = []
    name_set = set()
    for (rxn_left, rxn_right, tag) in reactions:
        name, pka = tag.split()
        if name in name_set:
            raise ValueError("reaction name must be unique") 
        name_set.add(name)
        r = {}
        r["name"] = name
        r["pka"] = float(pka)
        r["rxn_lose_h"] = AllChem.ReactionFromSmarts("%s >> %s" % (rxn_left, rxn_right))
        r["rxn_gain_h"] = AllChem.ReactionFromSmarts("%s >> %s" % (rxn_right, rxn_left))
        pka_reactions.append(r)
    return pka_reactions

def convert_recursive(mol, rxn, container):
    nr_react = rxn.GetNumReactantTemplates()
    nr_prod = rxn.GetNumProductTemplates()
    if nr_react != 1 or nr_prod != 1:
        raise RuntimeError("reaction %s must be single reactant -> single product" % name)
    for products in rxn.RunReactants((mol,)):
        product = products[0] # nr products == NumProductTemplates == 1
        try:
            Chem.SanitizeMol(product)
        except Chem.AtomValenceException:
            continue
        except Chem.KekulizeException:
            continue
        container.add(product)
        convert_recursive(product, rxn, container)

def convert_exhaustive(mol, rxn):
    """
        - If the reaction occurs, return the product.
        - If the reaction does not occur or sanitization fails, return the input mol.
        - If multiple substructures react, the returned product will have all groups converted
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


def enumerate_pka(mol, pka_reactions, ph_range_low, ph_range_high):
    if ph_range_low > ph_range_high:
        raise ValueError("ph_range_low must be lesser than or equal to ph_range_high")
    mol_list = [mol]
    for r in pka_reactions:
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
    return mol_list


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
