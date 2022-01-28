from rdkit import Chem 
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

import pandas as pd
import os, sys

def get_columns():
    rdkit2D = ["BalabanJ", "BertzCT", "Ipc", "HallKierAlpha", "Kappa1", "Kappa2", "Kappa3", "Chi0", "Chi1", "Chi0n", "Chi1n", "Chi2n", "Chi3n", "Chi4n", "Chi0v", "Chi1v", "Chi2v", "Chi3v", "Chi4v", "MolLogP", "MolMR", "MolWt", "ExactMolWt", "HeavyAtomCount", "HeavyAtomMolWt", "NHOHCount", "NOCount", "NumHAcceptors", "NumHDonors", "NumHeteroatoms", "NumRotatableBonds", "NumValenceElectrons", "NumAmideBonds", "NumAromaticRings", "NumSaturatedRings", "NumAliphaticRings", "NumAromaticHeterocycles", "NumSaturatedHeterocycles", "NumAliphaticHeterocycles", "NumAromaticCarbocycles", "NumSaturatedCarbocycles", "NumAliphaticCarbocycles", "RingCount", "FractionCSP3", "NumSpiroAtoms", "NumBridgeheadAtoms", "TPSA", "LabuteASA", "PEOE_VSA1", "PEOE_VSA2", "PEOE_VSA3", "PEOE_VSA4", "PEOE_VSA5", "PEOE_VSA6", "PEOE_VSA7", "PEOE_VSA8", "PEOE_VSA9", "PEOE_VSA10", "PEOE_VSA11", "PEOE_VSA12", "PEOE_VSA13", "PEOE_VSA14", "SMR_VSA1", "SMR_VSA2", "SMR_VSA3", "SMR_VSA4", "SMR_VSA5", "SMR_VSA6", "SMR_VSA7", "SMR_VSA8", "SMR_VSA9", "SMR_VSA10", "SlogP_VSA1", "SlogP_VSA2", "SlogP_VSA3", "SlogP_VSA4", "SlogP_VSA5", "SlogP_VSA6", "SlogP_VSA7", "SlogP_VSA8", "SlogP_VSA9", "SlogP_VSA10", "SlogP_VSA11", "SlogP_VSA12", "EState_VSA1", "EState_VSA2", "EState_VSA3", "EState_VSA4", "EState_VSA5", "EState_VSA6", "EState_VSA7", "EState_VSA8", "EState_VSA9", "EState_VSA10", "EState_VSA11", "VSA_EState1", "VSA_EState2", "VSA_EState3", "VSA_EState4", "VSA_EState5", "VSA_EState6", "VSA_EState7", "VSA_EState8", "VSA_EState9", "VSA_EState10", "MQNs_atom_counts_c", "MQNs_atom_counts_f", "MQNs_atom_counts_cl", "MQNs_atom_counts_br", "MQNs_atom_counts_i", "MQNs_atom_counts_s", "MQNs_atom_counts_p", "MQNs_atom_counts_an", "MQNs_atom_counts_cn", "MQNs_atom_counts_ao", "MQNs_atom_counts_co", "MQNs_atom_counts_hac", "MQNs_bond_counts_asb", "MQNs_bond_counts_adb", "MQNs_bond_counts_atb", "MQNs_bond_counts_csb", "MQNs_bond_counts_cdb", "MQNs_bond_counts_ctb", "MQNs_bond_counts_rbc", "MQNs_polarity_counts_hbam", "MQNs_polarity_counts_hba", "MQNs_polarity_counts_hbdm", "MQNs_polarity_counts_hbd", "MQNs_polarity_counts_negc", "MQNs_polarity_counts_posc", "MQNs_topology_counts_asv", "MQNs_topology_counts_adv", "MQNs_topology_counts_atv", "MQNs_topology_counts_aqv", "MQNs_topology_counts_cdv", "MQNs_topology_counts_ctv", "MQNs_topology_counts_cqv", "MQNs_topology_counts_r3", "MQNs_topology_counts_r4", "MQNs_topology_counts_r5", "MQNs_topology_counts_r6", "MQNs_topology_counts_r7", "MQNs_topology_counts_r8", "MQNs_topology_counts_rg", "MQNs_topology_counts_rgIO", "MQNs_topology_counts_afrc","MQNs_topology_counts_bfrc"]
    return rdkit2D

def mol2_to_mol(ligand, tmp_folder_name, name):
    ligand_name = os.path.basename(os.path.normpath(ligand)).split(".", 1)[0].lower()
    mol_path = os.path.join(tmp_folder_name, ligand_name + ".mol")

    command = "babel -imol2 " + ligand + " -omol " + mol_path
    print(command)
    os.system(command)

    return mol_path

def get_descriptors(mol_path, name):
    
    indx = name

    df = pd.DataFrame(columns = get_columns())

    m = Chem.MolFromMolFile(mol_path, sanitize=False)

    sanitize = False
    try:    
        Chem.SanitizeMol(m)
        sm = Chem.MolToSmiles(m)
        print("Sanitize Mol: True")
        sanitize = True
    except ValueError as e:
        sm = str(e)
        print("Sanitize Mol: False")

    df.loc[str(indx), "BalabanJ"] = Descriptors.BalabanJ(m)
    df.loc[str(indx), "BertzCT"] = Descriptors.BertzCT(m)
    df.loc[str(indx), "Ipc"] = Descriptors.Ipc(m)
    df.loc[str(indx), "HallKierAlpha"] = Descriptors.HallKierAlpha(m)
    df.loc[str(indx), "Kappa1"] = Descriptors.Kappa1(m)
    df.loc[str(indx), "Kappa2"] = Descriptors.Kappa2(m)
    df.loc[str(indx), "Kappa3"] = Descriptors.Kappa3(m)
    df.loc[str(indx), "Chi0"] = Descriptors.Chi0(m)
    df.loc[str(indx), "Chi1"] = Descriptors.Chi1(m)

    if sanitize:
        df.loc[str(indx), "Chi0n"] = Descriptors.Chi0n(m)
    else:
        df.loc[str(indx), "Chi0n"] = 0

    if sanitize:
        df.loc[str(indx), "Chi1n"] = Descriptors.Chi1n(m)
    else:
        df.loc[str(indx), "Chi1n"] = 0

    if sanitize:
        df.loc[str(indx), "Chi2n"] = Descriptors.Chi2n(m)
    else:
        df.loc[str(indx), "Chi2n"] = 0

    if sanitize:
        df.loc[str(indx), "Chi3n"] = Descriptors.Chi3n(m)
    else:
        df.loc[str(indx), "Chi3n"] = 0

    if sanitize:
        df.loc[str(indx), "Chi4n"] = Descriptors.Chi4n(m)
    else:
        df.loc[str(indx), "Chi4n"] = 0

    if sanitize:
        df.loc[str(indx), "Chi0v"] = Descriptors.Chi0v(m)
    else:
        df.loc[str(indx), "Chi0v"] = 0

    if sanitize:
        df.loc[str(indx), "Chi1v"] = Descriptors.Chi1v(m)
    else:
        df.loc[str(indx), "Chi1v"] = 0

    if sanitize:
        df.loc[str(indx), "Chi2v"] = Descriptors.Chi2v(m)
    else:
        df.loc[str(indx), "Chi2v"] = 0

    if sanitize:
        df.loc[str(indx), "Chi3v"] = Descriptors.Chi3v(m)
    else:
        df.loc[str(indx), "Chi3v"] = 0

    if sanitize:
        df.loc[str(indx), "Chi4v"] = Descriptors.Chi4v(m)
    else:
        df.loc[str(indx), "Chi4v"] = 0

    if sanitize:
        df.loc[str(indx), "MolLogP"] = Descriptors.MolLogP(m)
    else:
        df.loc[str(indx), "MolLogP"] = 0

    if sanitize:
        df.loc[str(indx), "MolMR"] = Descriptors.MolMR(m)
    else:
        df.loc[str(indx), "MolMR"] = 0

    if sanitize:
        df.loc[str(indx), "MolWt"] = Descriptors.MolWt(m)
    else:
        df.loc[str(indx), "MolWt"] = 0

    if sanitize:
        df.loc[str(indx), "ExactMolWt"] = Descriptors.ExactMolWt(m)
    else:
        df.loc[str(indx), "ExactMolWt"] = 0

    df.loc[str(indx), "HeavyAtomCount"] = Descriptors.HeavyAtomCount(m)
    df.loc[str(indx), "HeavyAtomMolWt"] = Descriptors.HeavyAtomMolWt(m)

    if sanitize:
        df.loc[str(indx), "NHOHCount"] = Descriptors.NHOHCount(m)
    else:
        df.loc[str(indx), "NHOHCount"] = 0

    df.loc[str(indx), "NOCount"] = Descriptors.NOCount(m)

    if sanitize:
        df.loc[str(indx), "NumHAcceptors"] = Descriptors.NumHAcceptors(m)
    else:
        df.loc[str(indx), "NumHAcceptors"] = 0

    if sanitize:
        df.loc[str(indx), "NumHDonors"] = Descriptors.NumHDonors(m)
    else:
        df.loc[str(indx), "NumHDonors"] = 0

    df.loc[str(indx), "NumHeteroatoms"] = Descriptors.NumHeteroatoms(m)

    if sanitize:
        df.loc[str(indx), "NumRotatableBonds"] = Descriptors.NumRotatableBonds(m)
    else:
        df.loc[str(indx), "NumRotatableBonds"] = 0

    if sanitize:
        df.loc[str(indx), "NumValenceElectrons"] = Descriptors.NumValenceElectrons(m)
    else:
        df.loc[str(indx), "NumValenceElectrons"] = 0

    if sanitize:
        df.loc[str(indx), "NumAmideBonds"] = rdMolDescriptors.CalcNumAmideBonds(m)
    else:
        df.loc[str(indx), "NumAmideBonds"] = 0

    df.loc[str(indx), "NumAromaticRings"] = Descriptors.NumAromaticRings(m)
    df.loc[str(indx), "NumSaturatedRings"] = Descriptors.NumSaturatedRings(m)
    df.loc[str(indx), "NumAliphaticRings"] = Descriptors.NumAliphaticRings(m)    

    df.loc[str(indx), "NumAromaticHeterocycles"] = Descriptors.NumAromaticHeterocycles(m)
    df.loc[str(indx), "NumSaturatedHeterocycles"] = Descriptors.NumSaturatedHeterocycles(m)
    df.loc[str(indx), "NumAliphaticHeterocycles"] = Descriptors.NumAliphaticHeterocycles(m)

    df.loc[str(indx), "NumAromaticCarbocycles"] = Descriptors.NumAromaticCarbocycles(m)
    df.loc[str(indx), "NumSaturatedCarbocycles"] = Descriptors.NumSaturatedCarbocycles(m)
    df.loc[str(indx), "NumAliphaticCarbocycles"] = Descriptors.NumAliphaticCarbocycles(m)

    if sanitize:
        df.loc[str(indx), "RingCount"] = Descriptors.RingCount(m)
    else:
        df.loc[str(indx), "RingCount"] = 0

    if sanitize:
        df.loc[str(indx), "FractionCSP3"] = Descriptors.FractionCSP3(m)
    else:
        df.loc[str(indx), "FractionCSP3"] = 0

    df.loc[str(indx), "NumSpiroAtoms"] = rdMolDescriptors.CalcNumSpiroAtoms(m)
    df.loc[str(indx), "NumBridgeheadAtoms"] = rdMolDescriptors.CalcNumBridgeheadAtoms(m)

    if sanitize:
        df.loc[str(indx), "TPSA"] = Descriptors.TPSA(m)
    else:
        df.loc[str(indx), "TPSA"] = 0

    df.loc[str(indx), "LabuteASA"] = Descriptors.LabuteASA(m)

    if sanitize:
        df.loc[str(indx), "PEOE_VSA1"] = Descriptors.PEOE_VSA1(m)
    else:
        df.loc[str(indx), "PEOE_VSA1"] = 0

    if sanitize:
        df.loc[str(indx), "PEOE_VSA2"] = Descriptors.PEOE_VSA2(m)
    else:
        df.loc[str(indx), "PEOE_VSA2"] = 0

    if sanitize:
        df.loc[str(indx), "PEOE_VSA3"] = Descriptors.PEOE_VSA3(m)
    else:
        df.loc[str(indx), "PEOE_VSA3"] = 0

    if sanitize:
        df.loc[str(indx), "PEOE_VSA4"] = Descriptors.PEOE_VSA4(m)
    else:
        df.loc[str(indx), "PEOE_VSA4"] = 0

    if sanitize:
        df.loc[str(indx), "PEOE_VSA5"] = Descriptors.PEOE_VSA5(m)
    else:
        df.loc[str(indx), "PEOE_VSA5"] = 0

    if sanitize:
        df.loc[str(indx), "PEOE_VSA6"] = Descriptors.PEOE_VSA6(m)
    else:
        df.loc[str(indx), "PEOE_VSA6"] = 0

    if sanitize:
        df.loc[str(indx), "PEOE_VSA7"] = Descriptors.PEOE_VSA7(m)
    else:
        df.loc[str(indx), "PEOE_VSA7"] = 0

    if sanitize:
        df.loc[str(indx), "PEOE_VSA8"] = Descriptors.PEOE_VSA8(m)
    else:
        df.loc[str(indx), "PEOE_VSA8"] = 0

    if sanitize:
        df.loc[str(indx), "PEOE_VSA9"] = Descriptors.PEOE_VSA9(m)
    else:
        df.loc[str(indx), "PEOE_VSA9"] = 0

    if sanitize:
        df.loc[str(indx), "PEOE_VSA10"] = Descriptors.PEOE_VSA10(m)
    else:
        df.loc[str(indx), "PEOE_VSA10"] = 0

    if sanitize:
        df.loc[str(indx), "PEOE_VSA11"] = Descriptors.PEOE_VSA11(m)
    else:
        df.loc[str(indx), "PEOE_VSA11"] = 0

    if sanitize:
        df.loc[str(indx), "PEOE_VSA12"] = Descriptors.PEOE_VSA12(m)
    else:
        df.loc[str(indx), "PEOE_VSA12"] = 0

    if sanitize:
        df.loc[str(indx), "PEOE_VSA13"] = Descriptors.PEOE_VSA13(m)
    else:
        df.loc[str(indx), "PEOE_VSA13"] = 0

    if sanitize:
        df.loc[str(indx), "PEOE_VSA14"] = Descriptors.PEOE_VSA14(m)
    else:
        df.loc[str(indx), "PEOE_VSA14"] = 0

    if sanitize:
        df.loc[str(indx), "SMR_VSA1"] = Descriptors.SMR_VSA1(m)
    else:
        df.loc[str(indx), "SMR_VSA1"] = 0

    if sanitize:
        df.loc[str(indx), "SMR_VSA2"] = Descriptors.SMR_VSA2(m)
    else:
        df.loc[str(indx), "SMR_VSA2"] = 0

    if sanitize:
        df.loc[str(indx), "SMR_VSA3"] = Descriptors.SMR_VSA3(m)
    else:
        df.loc[str(indx), "SMR_VSA3"] = 0

    if sanitize:
        df.loc[str(indx), "SMR_VSA4"] = Descriptors.SMR_VSA4(m)
    else:
        df.loc[str(indx), "SMR_VSA4"] = 0

    if sanitize:
        df.loc[str(indx), "SMR_VSA5"] = Descriptors.SMR_VSA5(m)
    else:
        df.loc[str(indx), "SMR_VSA5"] = 0

    if sanitize:
        df.loc[str(indx), "SMR_VSA6"] = Descriptors.SMR_VSA6(m)
    else:
        df.loc[str(indx), "SMR_VSA6"] = 0

    if sanitize:
        df.loc[str(indx), "SMR_VSA7"] = Descriptors.SMR_VSA7(m)
    else:
        df.loc[str(indx), "SMR_VSA7"] = 0

    if sanitize:
        df.loc[str(indx), "SMR_VSA8"] = Descriptors.SMR_VSA8(m)
    else:
        df.loc[str(indx), "SMR_VSA8"] = 0

    if sanitize:
        df.loc[str(indx), "SMR_VSA9"] = Descriptors.SMR_VSA9(m)
    else:
        df.loc[str(indx), "SMR_VSA9"] = 0

    if sanitize:
        df.loc[str(indx), "SMR_VSA10"] = Descriptors.SMR_VSA10(m)
    else:
        df.loc[str(indx), "SMR_VSA10"] = 0

    if sanitize:
        df.loc[str(indx), "SlogP_VSA1"] = Descriptors.SlogP_VSA1(m)
    else:
        df.loc[str(indx), "SlogP_VSA1"] = 0

    if sanitize:
        df.loc[str(indx), "SlogP_VSA2"] = Descriptors.SlogP_VSA2(m)
    else:
        df.loc[str(indx), "SlogP_VSA2"] = 0

    if sanitize:
        df.loc[str(indx), "SlogP_VSA3"] = Descriptors.SlogP_VSA3(m)
    else:
        df.loc[str(indx), "SlogP_VSA3"] = 0

    if sanitize:
        df.loc[str(indx), "SlogP_VSA4"] = Descriptors.SlogP_VSA4(m)
    else:
        df.loc[str(indx), "SlogP_VSA4"] = 0

    if sanitize:
        df.loc[str(indx), "SlogP_VSA5"] = Descriptors.SlogP_VSA5(m)
    else:
        df.loc[str(indx), "SlogP_VSA5"] = 0

    if sanitize:
        df.loc[str(indx), "SlogP_VSA6"] = Descriptors.SlogP_VSA6(m)
    else:
        df.loc[str(indx), "SlogP_VSA6"] = 0

    if sanitize:
        df.loc[str(indx), "SlogP_VSA7"] = Descriptors.SlogP_VSA7(m)
    else:
        df.loc[str(indx), "SlogP_VSA7"] = 0

    if sanitize:
        df.loc[str(indx), "SlogP_VSA8"] = Descriptors.SlogP_VSA8(m)
    else:
        df.loc[str(indx), "SlogP_VSA8"] = 0

    if sanitize:
        df.loc[str(indx), "SlogP_VSA9"] = Descriptors.SlogP_VSA9(m)
    else:
        df.loc[str(indx), "SlogP_VSA9"] = 0

    if sanitize:
        df.loc[str(indx), "SlogP_VSA10"] = Descriptors.SlogP_VSA10(m)
    else:
        df.loc[str(indx), "SlogP_VSA10"] = 0

    if sanitize:
        df.loc[str(indx), "SlogP_VSA11"] = Descriptors.SlogP_VSA11(m)
    else:
        df.loc[str(indx), "SlogP_VSA11"] = 0

    if sanitize:
        df.loc[str(indx), "SlogP_VSA12"] = Descriptors.SlogP_VSA12(m)
    else:
        df.loc[str(indx), "SlogP_VSA12"] = 0

    if sanitize:
        df.loc[str(indx), "EState_VSA1"] = Descriptors.EState_VSA1(m)
    else:
        df.loc[str(indx), "EState_VSA1"] = 0

    if sanitize:
        df.loc[str(indx), "EState_VSA2"] = Descriptors.EState_VSA2(m)
    else:
        df.loc[str(indx), "EState_VSA2"] = 0

    if sanitize:
        df.loc[str(indx), "EState_VSA3"] = Descriptors.EState_VSA3(m)
    else:
        df.loc[str(indx), "EState_VSA3"] = 0

    if sanitize:
        df.loc[str(indx), "EState_VSA4"] = Descriptors.EState_VSA4(m)
    else:
        df.loc[str(indx), "EState_VSA4"] = 0

    if sanitize:
        df.loc[str(indx), "EState_VSA5"] = Descriptors.EState_VSA5(m)
    else:
        df.loc[str(indx), "EState_VSA5"] = 0

    if sanitize:
        df.loc[str(indx), "EState_VSA6"] = Descriptors.EState_VSA6(m)
    else:
        df.loc[str(indx), "EState_VSA6"] = 0

    if sanitize:
        df.loc[str(indx), "EState_VSA7"] = Descriptors.EState_VSA7(m)
    else:
        df.loc[str(indx), "EState_VSA7"] = 0

    if sanitize:
        df.loc[str(indx), "EState_VSA8"] = Descriptors.EState_VSA8(m)
    else:
        df.loc[str(indx), "EState_VSA8"] = 0

    if sanitize:
        df.loc[str(indx), "EState_VSA9"] = Descriptors.EState_VSA9(m)
    else:
        df.loc[str(indx), "EState_VSA9"] = 0

    if sanitize:
        df.loc[str(indx), "EState_VSA10"] = Descriptors.EState_VSA10(m)
    else:
        df.loc[str(indx), "EState_VSA10"] = 0

    if sanitize:
        df.loc[str(indx), "EState_VSA11"] = Descriptors.EState_VSA11(m)
    else:
        df.loc[str(indx), "EState_VSA11"] = 0

    if sanitize:
        df.loc[str(indx), "VSA_EState1"] = Descriptors.VSA_EState1(m)
    else:
        df.loc[str(indx), "VSA_EState1"] = 0

    if sanitize:
        df.loc[str(indx), "VSA_EState2"] = Descriptors.VSA_EState2(m)
    else:
        df.loc[str(indx), "VSA_EState2"] = 0

    if sanitize:
        df.loc[str(indx), "VSA_EState3"] = Descriptors.VSA_EState3(m)
    else:
        df.loc[str(indx), "VSA_EState3"] = 0

    if sanitize:
        df.loc[str(indx), "VSA_EState4"] = Descriptors.VSA_EState4(m)
    else:
        df.loc[str(indx), "VSA_EState4"] = 0

    if sanitize:
        df.loc[str(indx), "VSA_EState5"] = Descriptors.VSA_EState5(m)
    else:
        df.loc[str(indx), "VSA_EState5"] = 0

    if sanitize:
        df.loc[str(indx), "VSA_EState6"] = Descriptors.VSA_EState6(m)
    else:
        df.loc[str(indx), "VSA_EState6"] = 0

    if sanitize:
        df.loc[str(indx), "VSA_EState7"] = Descriptors.VSA_EState7(m)
    else:
        df.loc[str(indx), "VSA_EState7"] = 0

    if sanitize:
        df.loc[str(indx), "VSA_EState8"] = Descriptors.VSA_EState8(m)
    else:
        df.loc[str(indx), "VSA_EState8"] = 0

    if sanitize:
        df.loc[str(indx), "VSA_EState9"] = Descriptors.VSA_EState9(m)
    else:
        df.loc[str(indx), "VSA_EState9"] = 0

    if sanitize:
        df.loc[str(indx), "VSA_EState10"] = Descriptors.VSA_EState10(m)
    else:
        df.loc[str(indx), "VSA_EState10"] = 0

    if sanitize:
        df.loc[str(indx), "MQNs_atom_counts_c"] = rdMolDescriptors.MQNs_(m)[0]
    else:
        df.loc[str(indx), "MQNs_atom_counts_c"] = 0

    if sanitize:
        df.loc[str(indx), "MQNs_atom_counts_f"] = rdMolDescriptors.MQNs_(m)[1]
    else:
        df.loc[str(indx), "MQNs_atom_counts_f"] = 0

    if sanitize:
        df.loc[str(indx), "MQNs_atom_counts_cl"] = rdMolDescriptors.MQNs_(m)[2]
    else:
        df.loc[str(indx), "MQNs_atom_counts_cl"] = 0

    if sanitize:
        df.loc[str(indx), "MQNs_atom_counts_br"] = rdMolDescriptors.MQNs_(m)[3]
    else:
        df.loc[str(indx), "MQNs_atom_counts_br"] = 0

    if sanitize:
        df.loc[str(indx), "MQNs_atom_counts_i"] = rdMolDescriptors.MQNs_(m)[4]
    else:
        df.loc[str(indx), "MQNs_atom_counts_i"] = 0

    if sanitize:
        df.loc[str(indx), "MQNs_atom_counts_s"] = rdMolDescriptors.MQNs_(m)[5]
    else:
        df.loc[str(indx), "MQNs_atom_counts_s"] = 0

    if sanitize:
        df.loc[str(indx), "MQNs_atom_counts_p"] = rdMolDescriptors.MQNs_(m)[6]
    else:
        df.loc[str(indx), "MQNs_atom_counts_p"] = 0

    if sanitize:
        df.loc[str(indx), "MQNs_atom_counts_an"] = rdMolDescriptors.MQNs_(m)[7]
    else:
        df.loc[str(indx), "MQNs_atom_counts_an"] = 0

    if sanitize:
        df.loc[str(indx), "MQNs_atom_counts_cn"] = rdMolDescriptors.MQNs_(m)[8]
    else:
        df.loc[str(indx), "MQNs_atom_counts_cn"] = 0

    if sanitize:
        df.loc[str(indx), "MQNs_atom_counts_ao"] = rdMolDescriptors.MQNs_(m)[9]
    else:
        df.loc[str(indx), "MQNs_atom_counts_ao"] = 0

    if sanitize:
        df.loc[str(indx), "MQNs_atom_counts_co"] = rdMolDescriptors.MQNs_(m)[10]
    else:
        df.loc[str(indx), "MQNs_atom_counts_co"] = 0

    if sanitize:
        df.loc[str(indx), "MQNs_atom_counts_hac"] = rdMolDescriptors.MQNs_(m)[11]
    else:
        df.loc[str(indx), "MQNs_atom_counts_hac"] = 0

    if sanitize:
        df.loc[str(indx), "MQNs_bond_counts_asb"] = rdMolDescriptors.MQNs_(m)[12]
    else:
        df.loc[str(indx), "MQNs_bond_counts_asb"] = 0

    if sanitize:
        df.loc[str(indx), "MQNs_bond_counts_adb"] = rdMolDescriptors.MQNs_(m)[13]
    else:
        df.loc[str(indx), "MQNs_bond_counts_adb"] = 0

    if sanitize:
        df.loc[str(indx), "MQNs_bond_counts_atb"] = rdMolDescriptors.MQNs_(m)[14]
    else:
        df.loc[str(indx), "MQNs_bond_counts_atb"] = 0

    if sanitize:
        df.loc[str(indx), "MQNs_bond_counts_csb"] = rdMolDescriptors.MQNs_(m)[15]
    else:
        df.loc[str(indx), "MQNs_bond_counts_csb"] = 0

    if sanitize:
        df.loc[str(indx), "MQNs_bond_counts_cdb"] = rdMolDescriptors.MQNs_(m)[16]
    else:
        df.loc[str(indx), "MQNs_bond_counts_cdb"] = 0

    if sanitize:
        df.loc[str(indx), "MQNs_bond_counts_ctb"] = rdMolDescriptors.MQNs_(m)[17]
    else:
        df.loc[str(indx), "MQNs_bond_counts_ctb"] = 0

    if sanitize:
        df.loc[str(indx), "MQNs_bond_counts_rbc"] = rdMolDescriptors.MQNs_(m)[18]
    else:
        df.loc[str(indx), "MQNs_bond_counts_rbc"] = 0

    if sanitize:
        df.loc[str(indx), "MQNs_polarity_counts_hbam"] = rdMolDescriptors.MQNs_(m)[19]
    else:
        df.loc[str(indx), "MQNs_polarity_counts_hbam"] = 0

    if sanitize:
        df.loc[str(indx), "MQNs_polarity_counts_hba"] = rdMolDescriptors.MQNs_(m)[20]
    else:
        df.loc[str(indx), "MQNs_polarity_counts_hba"] = 0

    if sanitize:
        df.loc[str(indx), "MQNs_polarity_counts_hbdm"] = rdMolDescriptors.MQNs_(m)[21]
    else:
        df.loc[str(indx), "MQNs_polarity_counts_hbdm"] = 0

    if sanitize:
        df.loc[str(indx), "MQNs_polarity_counts_hbd"] = rdMolDescriptors.MQNs_(m)[22]
    else:
        df.loc[str(indx), "MQNs_polarity_counts_hbd"] = 0

    if sanitize:
        df.loc[str(indx), "MQNs_polarity_counts_negc"] = rdMolDescriptors.MQNs_(m)[23]
    else:
        df.loc[str(indx), "MQNs_polarity_counts_negc"] = 0

    if sanitize:
        df.loc[str(indx), "MQNs_polarity_counts_posc"] = rdMolDescriptors.MQNs_(m)[24]
    else:
        df.loc[str(indx), "MQNs_polarity_counts_posc"] = 0

    if sanitize:
        df.loc[str(indx), "MQNs_topology_counts_asv"] = rdMolDescriptors.MQNs_(m)[25]
    else:
        df.loc[str(indx), "MQNs_topology_counts_asv"] = 0

    if sanitize:
        df.loc[str(indx), "MQNs_topology_counts_adv"] = rdMolDescriptors.MQNs_(m)[26]
    else:
        df.loc[str(indx), "MQNs_topology_counts_adv"] = 0

    if sanitize:
        df.loc[str(indx), "MQNs_topology_counts_atv"] = rdMolDescriptors.MQNs_(m)[27]
    else:
        df.loc[str(indx), "MQNs_topology_counts_atv"] = 0

    if sanitize:
        df.loc[str(indx), "MQNs_topology_counts_aqv"] = rdMolDescriptors.MQNs_(m)[28]
    else:
        df.loc[str(indx), "MQNs_topology_counts_aqv"] = 0

    if sanitize:
        df.loc[str(indx), "MQNs_topology_counts_cdv"] = rdMolDescriptors.MQNs_(m)[29]
    else:
        df.loc[str(indx), "MQNs_topology_counts_cdv"] = 0

    if sanitize:
        df.loc[str(indx), "MQNs_topology_counts_ctv"] = rdMolDescriptors.MQNs_(m)[30]
    else:
        df.loc[str(indx), "MQNs_topology_counts_ctv"] = 0

    if sanitize:
        df.loc[str(indx), "MQNs_topology_counts_cqv"] = rdMolDescriptors.MQNs_(m)[31]
    else:
        df.loc[str(indx), "MQNs_topology_counts_cqv"] = 0

    if sanitize:
        df.loc[str(indx), "MQNs_topology_counts_r3"] = rdMolDescriptors.MQNs_(m)[32]
    else:
        df.loc[str(indx), "MQNs_topology_counts_r3"] = 0

    if sanitize:
        df.loc[str(indx), "MQNs_topology_counts_r4"] = rdMolDescriptors.MQNs_(m)[33]
    else:
        df.loc[str(indx), "MQNs_topology_counts_r4"] = 0

    if sanitize:
        df.loc[str(indx), "MQNs_topology_counts_r5"] = rdMolDescriptors.MQNs_(m)[34]
    else:
        df.loc[str(indx), "MQNs_topology_counts_r5"] = 0

    if sanitize:
        df.loc[str(indx), "MQNs_topology_counts_r6"] = rdMolDescriptors.MQNs_(m)[35]
    else:
        df.loc[str(indx), "MQNs_topology_counts_r6"] = 0

    if sanitize:
        df.loc[str(indx), "MQNs_topology_counts_r7"] = rdMolDescriptors.MQNs_(m)[36]
    else:
        df.loc[str(indx), "MQNs_topology_counts_r7"] = 0

    if sanitize:
        df.loc[str(indx), "MQNs_topology_counts_r8"] = rdMolDescriptors.MQNs_(m)[37]
    else:
        df.loc[str(indx), "MQNs_topology_counts_r8"] = 0

    if sanitize:
        df.loc[str(indx), "MQNs_topology_counts_rg"] = rdMolDescriptors.MQNs_(m)[38]
    else:
        df.loc[str(indx), "MQNs_topology_counts_rg"] = 0

    if sanitize:
        df.loc[str(indx), "MQNs_topology_counts_rgIO"] = rdMolDescriptors.MQNs_(m)[39]
    else:
        df.loc[str(indx), "MQNs_topology_counts_rgIO"] = 0

    if sanitize:
        df.loc[str(indx), "MQNs_topology_counts_afrc"] = rdMolDescriptors.MQNs_(m)[40]
    else:
        df.loc[str(indx), "MQNs_topology_counts_afrc"] = 0

    if sanitize:
        df.loc[str(indx), "MQNs_topology_counts_bfrc"] = rdMolDescriptors.MQNs_(m)[41]
    else:
        df.loc[str(indx), "MQNs_topology_counts_bfrc"] = 0

    #df.loc[str(indx), "sanitize"] = sanitize
    #df.loc[str(indx), "protein"] = protein

    df["pdb"] = name
    df.index = df.pop("pdb")

    return df


def get_rdkit2d147_descriptors(ligand, tmp_folder_name, name):
    
    mol_path = mol2_to_mol(ligand = ligand, tmp_folder_name = tmp_folder_name, name = name)
    dt = get_descriptors(mol_path = mol_path, name = name)

    return dt