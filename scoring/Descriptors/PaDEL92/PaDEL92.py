
import os, sys
import pandas as pd

def get_columns():
    padel92 = ["nAtom", "nHeavyAtom", "nH", "nB", "nC", "nN", "nO", "nS", "nP", "nF", "nCl", "nBr", "nI", "nX", "nBonds", "nBonds2", "nBondsS", "nBondsS2", "nBondsS3", "nBondsD", "nBondsD2", "nBondsT", "nBondsQ", "nBondsM", "nRing", "n3Ring", "n4Ring", "n5Ring", "n6Ring", "n7Ring", "n8Ring", "n9Ring", "n10Ring", "n11Ring", "n12Ring", "nG12Ring", "nFRing", "nF4Ring", "nF5Ring", "nF6Ring", "nF7Ring", "nF8Ring", "nF9Ring", "nF10Ring", "nF11Ring", "nF12Ring", "nFG12Ring", "nTRing", "nT4Ring", "nT5Ring", "nT6Ring", "nT7Ring", "nT8Ring", "nT9Ring", "nT10Ring", "nT11Ring", "nT12Ring", "nTG12Ring", "nHeteroRing", "n3HeteroRing", "n4HeteroRing", "n5HeteroRing", "n6HeteroRing", "n7HeteroRing", "n8HeteroRing", "n9HeteroRing", "n10HeteroRing", "n11HeteroRing", "n12HeteroRing", "nG12HeteroRing", "nFHeteroRing", "nF4HeteroRing", "nF5HeteroRing", "nF6HeteroRing", "nF7HeteroRing", "nF8HeteroRing", "nF9HeteroRing", "nF10HeteroRing", "nF11HeteroRing", "nF12HeteroRing", "nFG12HeteroRing", "nTHeteroRing", "nT4HeteroRing", "nT5HeteroRing", "nT6HeteroRing", "nT7HeteroRing", "nT8HeteroRing", "nT9HeteroRing", "nT10HeteroRing", "nT11HeteroRing", "nT12HeteroRing", "nTG12HeteroRing"]
    return padel92

def mol2_to_smi(ligand, tmp_folder_name, name):
    ligand_name = os.path.basename(os.path.normpath(ligand)).split(".", 1)[0].lower()
    smi_path = os.path.join(tmp_folder_name, ligand_name + ".smi")

    command = "babel -imol2 " + ligand + " -osmi " + smi_path
    print(command)
    os.system(command)

    return smi_path

def execute_padel_descriptors(smi_path, tmp_folder_name, name):
    base_path = os.path.join(__file__.rsplit("/", 1)[0], "PaDEL-Descriptor")

    jar_path = os.path.join(base_path, "PaDEL-Descriptor.jar")
    config_path = os.path.join(os.path.join(base_path, "PaDEL_Configuration"), "Default_Configuration")
    descriptors_path = os.path.join(os.path.join(base_path, "PaDEL_Configuration"), "Descriptors_Default.xml")
    result_path = os.path.join(tmp_folder_name, name + ".padel")

    command = "java -jar -Xmx20G " + jar_path + " -2d -descriptortypes " + descriptors_path + " -dir " + smi_path + " -file " + result_path
    
    print(command)
    os.system(command)

    return result_path

def get_descriptors(result_path, name):
    dt = pd.read_csv(result_path)
    dt = dt.iloc[:, 1:len(dt.columns) + 1]

    dt["pdb"] = name
    dt.index = dt.pop("pdb")

    return dt

def get_pd92_descriptors(ligand, tmp_folder_name, name):
    
    smi_path = mol2_to_smi(ligand = ligand, tmp_folder_name = tmp_folder_name, name = name)    
    result_path = execute_padel_descriptors(smi_path = ligand, tmp_folder_name = tmp_folder_name, name = name)
    
    dt = get_descriptors(result_path = result_path, name = name)
    
    return dt
