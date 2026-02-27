# The dataset used to train the model comes from:
#  https://github.com/czodrowskilab/Machine-learning-meets-pKa
#
# The dataset is distributed under the CC-BY-4.0 license (also included). 
# This license applies ONLY to the dataset, and not to the rest of the molscrub or training
# code. That is still under GPL3.0. 


from rdkit import Chem
from rdkit.Chem import PandasTools
import pandas as pd
from molscrub import AcidBaseConjugator
from molscrub.protonate import convert_exhaustive, convert_all_single_sites
import itertools

import argparse
import pickle


def getPkaFromMol(mol: Chem.Mol, target_atom_index: int) -> list:
    conjugator = AcidBaseConjugator.from_default_data_files()

    size = len(conjugator.pka_reactions)
    reacted_mols = []

    for i,r in enumerate(conjugator.pka_reactions):
        temp_forward = convert_exhaustive(mol, r["rxn_gain_h"])
        temp_forward = Chem.MolFromSmiles(Chem.MolToSmiles(temp_forward))
        if not conjugator.mol_comparisons(mol, temp_forward):
            changed_atom = conjugator.find_protonation_site_with_mcs(mol, temp_forward)
            if changed_atom == target_atom_index:
                reacted_mols.append((temp_forward, r["pka"], r["name"], "forward", changed_atom, conjugator._one_hot(size, i)))

        temp_backward = convert_exhaustive(mol, r["rxn_lose_h"])
        temp_backward = Chem.MolFromSmiles(Chem.MolToSmiles(temp_backward))
        if not conjugator.mol_comparisons(mol, temp_backward):
            changed_atom = conjugator.find_protonation_site_with_mcs(mol, temp_backward)
            if changed_atom == target_atom_index:
                reacted_mols.append((temp_backward, r["pka"], r["name"], "backward", changed_atom, conjugator._one_hot(size, i)))
    
    return reacted_mols

def convert_single_sites_allrxns(mol):
    conjugator = AcidBaseConjugator.from_default_data_files()

    reacted_mols = []

    seen_smiles = set()

    size = len(conjugator.pka_reactions)

    for i,r in enumerate(conjugator.pka_reactions):
        temp_forward = convert_all_single_sites(mol, r["rxn_gain_h"])

        for m in temp_forward:    
            smi = Chem.MolToSmiles(m, canonical=True)
            if smi not in seen_smiles:
                seen_smiles.add(smi)
                changed_atom = conjugator.find_protonation_site_with_mcs(mol, m)
                reacted_mols.append({"original": mol, 
                                     "product": m, 
                                     "rule_pka": r["pka"], 
                                     "rxn_name": r["name"], 
                                     "direction": "gain_h", 
                                     "protonated_atom": changed_atom, 
                                     "rxn_1hot_encoding": conjugator._one_hot(size, i)})


        temp_backward = convert_all_single_sites(mol, r["rxn_lose_h"])
        for m in temp_backward:    
            smi = Chem.MolToSmiles(m, canonical=True)
            if smi not in seen_smiles:
                seen_smiles.add(smi)
                changed_atom = conjugator.find_protonation_site_with_mcs(mol, m)
                reacted_mols.append({"original": mol, 
                                     "product": m, 
                                     "rule_pka": r["pka"], 
                                     "rxn_name": r["name"], 
                                     "direction": "lose_h", 
                                     "protonated_atom": changed_atom, 
                                     "rxn_1hot_encoding": conjugator._one_hot(size, i)})
    
    return reacted_mols


def generate_biased_data():
    conjugator = AcidBaseConjugator.from_default_data_files()

    new_data = [{"smiles": "CC1CN(C)CCN1C", "marvin_atom": 3, "pKa": 9}, 
                {"smiles": "CC1C[NH+](C)CCN1C", "marvin_atom": 3, "pKa": 9},
                {"smiles": "CC1CN(C)CC[NH+]1C", "marvin_atom": 3, "pKa": 5},
                {"smiles": "CC1C[NH+](C)CC[NH+]1C", "marvin_atom": 3, "pKa": 5},
                {"smiles": "CCCNCCNC", "marvin_atom": 3, "pKa": 9}, 
                {"smiles": "CCC[NH2+]CCNC", "marvin_atom": 3, "pKa": 9},
                {"smiles": "CCCNCC[NH2+]C", "marvin_atom": 3, "pKa": 5},
                {"smiles": "CCC[NH2+]CC[NH2+]C", "marvin_atom": 3, "pKa": 5},
                {"smiles": "CCCNCCNCC=C", "marvin_atom": 3, "pKa": 9}, 
                {"smiles": "CCC[NH2+]CCNCC=C", "marvin_atom": 3, "pKa": 9},
                {"smiles": "CCCNCC[NH2+]CC=C", "marvin_atom": 3, "pKa": 5},
                {"smiles": "CCC[NH2+]CC[NH2+]CC=C", "marvin_atom": 3, "pKa": 5},
                {"smiles": "CCCNCCNCCC", "marvin_atom": 3, "pKa": 9}, 
                {"smiles": "CCC[NH2+]CCNCCC", "marvin_atom": 3, "pKa": 9},
                {"smiles": "CCCNCC[NH2+]CCC", "marvin_atom": 3, "pKa": 5},
                {"smiles": "CCC[NH2+]CC[NH2+]CCC", "marvin_atom": 3, "pKa": 5}
                ]


    tempDF = pd.DataFrame(new_data)
    tempDF["ROMol"] = tempDF.smiles.apply(Chem.MolFromSmiles)
    pkas2 = tempDF.ROMol.apply(convert_single_sites_allrxns).to_list()
    new_df = pd.DataFrame(list(itertools.chain.from_iterable(pkas2))).rename(columns={"original": "ROMol", 
                                                                                    "product": "products", 
                                                                                    "rule_pka": "base_pkas", 
                                                                                    "rxn_1hot_encoding": "rxn_encoding" })
    # columns=["ROMol", "products", "base_pkas", "rxn_name", "direction", "protonated_atom", "rxn_encoding"]
    # for x in pkas2:
    #     print(x[0].values(), x[1].values() if len(x) == 2 else x[0].values())
    new_df = new_df[new_df["protonated_atom"] == 3]


    new_df["charge_diff"] = new_df.apply(lambda x: conjugator._charge_diff(x.ROMol, x.protonated_atom), axis=1)
    new_df["pKa"] = tempDF["pKa"].to_list()
    return new_df

if __name__ == "__main__":

    #get the output filename to save model
    parser = argparse.ArgumentParser(description="train model with given dataset")
    parser.add_argument("--dataset", type=str, default="dataset.sdf", help="input sdf dataset to train on.")
    parser.add_argument("--model_out", type=str, help="output filename to save the model (pickle)")
    args = parser.parse_args()

    # Load data
    sdf_path = args.dataset
    all_df = PandasTools.LoadSDF(sdf_path)
    all_df["smiles"] = all_df.ROMol.apply(Chem.MolToSmiles)
    all_df[["pKa", "marvin_pKa"]] = all_df[["pKa", "marvin_pKa"]].astype(float)
    all_df[["marvin_atom"]] = all_df[["marvin_atom"]].astype(int)

    pkas = all_df.apply(lambda x: getPkaFromMol(x.ROMol, x.marvin_atom), axis=1).to_list()


    conjugator = AcidBaseConjugator.from_default_data_files()
    
    new_df = pd.DataFrame([zip(*x) for x in pkas], columns=["products", "base_pkas", "rxn_name", "direction", "protonated_atom", "rxn_encoding"])

    data1 = pd.concat((all_df, new_df), axis=1).dropna()

    # select single reactions and expand. 
    data2 = data1[data1.direction.apply(len) == 1]
    data2[new_df.columns.values] = data2[new_df.columns.values].map(lambda x: x[0])
    data3 = data2[['pKa', "marvin_atom", 'ID', 'ROMol', 'smiles', 'products', 'base_pkas', 'rxn_name', 'direction', 'protonated_atom', 'rxn_encoding']]
    data3["pka_diff"] = data3.pKa - data3.base_pkas
    data3["charge_diff"] = data3.apply(lambda x: conjugator._charge_diff(x.ROMol, x.protonated_atom), axis=1)

    biased_df = generate_biased_data()

    descriptors = data3.ROMol.apply(conjugator.getMolDescriptors)
    biased_descriptors = biased_df.ROMol.apply(conjugator.getMolDescriptors)


    desc_data = pd.DataFrame(descriptors.to_list())
    # desc_data = desc_data[['fr_COO','fr_COO2','fr_C_S','NumHeteroatoms','NumAromaticHeterocycles', 'MaxPartialCharge', 'MinAbsPartialCharge', 'fr_sulfonamd', 'fr_Ar_N', 'PEOE_VSA9', 'PEOE_VSA8', 'fr_NH1', 'PEOE_VSA4', 'fr_aniline', 'fr_morpholine', 'PEOE_VSA7', 'PEOE_VSA2', 'qed', 'BCUT2D_MRHI', 'fr_nitrile', 'MinEStateIndex', 'MaxAbsPartialCharge', 'fr_halogen', 'BCUT2D_LOGPLOW', 'SMR_VSA1']]
    data4 = pd.concat((data3.reset_index(drop=True), desc_data), axis=1)

    biased_desc_data = pd.DataFrame(biased_descriptors.to_list())
    # new_desc_data = new_desc_data[['fr_COO','fr_COO2','fr_C_S','NumHeteroatoms','NumAromaticHeterocycles', 'MaxPartialCharge', 'MinAbsPartialCharge', 'fr_sulfonamd', 'fr_Ar_N', 'PEOE_VSA9', 'PEOE_VSA8', 'fr_NH1', 'PEOE_VSA4', 'fr_aniline', 'fr_morpholine', 'PEOE_VSA7', 'PEOE_VSA2', 'qed', 'BCUT2D_MRHI', 'fr_nitrile', 'MinEStateIndex', 'MaxAbsPartialCharge', 'fr_halogen', 'BCUT2D_LOGPLOW', 'SMR_VSA1']]
    biased_data4 = pd.concat((biased_df.reset_index(drop=True), biased_desc_data), axis=1)


    # generate x and y
    x = data4[data4.columns[13:]]
    expanded_cols = data4['rxn_encoding'].apply(pd.Series)
    expanded_cols.columns = [f'encoding{i+1}' for i in range(expanded_cols.shape[1])] # Rename columns if needed

    x = pd.concat([x, expanded_cols], axis=1)

    x["base_pka"] = data4.base_pkas
    x["charge_diff"] = data4.charge_diff

    x_biased = biased_data4[biased_data4.columns[9:]]
    expanded_cols = biased_data4['rxn_encoding'].apply(pd.Series)
    expanded_cols.columns = [f'encoding{i+1}' for i in range(expanded_cols.shape[1])] # Rename columns if needed
    x_biased = pd.concat([x_biased, expanded_cols], axis=1)

    x_biased["base_pka"] = biased_data4.base_pkas
    x_biased["charge_diff"] = biased_data4.charge_diff

    y = data4.pKa
    y_biased = biased_data4.pKa

    # train using the biased dataset :

    from sklearn.ensemble import ExtraTreesRegressor
    from sklearn.model_selection import train_test_split
    from sklearn.metrics import root_mean_squared_error

    x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.20)

    # biased dataset
    x_train = pd.concat((x_train.reset_index(drop=True), x_biased.reset_index(drop=True)), axis=0)
    y_train = pd.concat((y_train.reset_index(drop=True), y_biased.reset_index(drop=True)), axis=0)

    etr = ExtraTreesRegressor(n_estimators=600,  n_jobs=-1)
    etr.fit(x_train, y_train)
    # l_model.score(x_test, )
    y_predict = etr.predict(x_test)
    mse = root_mean_squared_error(y_test, y_predict)
    print(f"RMSE on test data: {mse:.4f}")
    print(f"Score: {etr.score(x_test, y_test)}")

    with open(args.model_out, 'wb') as f:
        pickle.dump(etr, f )