import scoring.Descriptors.NNscore350.oddt
from scoring.Descriptors.NNscore350.oddt.scoring.descriptors.binana import binana_descriptor
import six
import os, sys

import pandas as pd

def get_binana_titles():
    titles = binana_descriptor().titles
    return titles

def mol2_to_sdf(ligand, tmp_folder_name, name):
    ligand_name = os.path.basename(os.path.normpath(ligand)).split(".", 1)[0].lower()
    sdf_path = os.path.join(tmp_folder_name, ligand_name + ".sdf")

    command = "babel -imol2 " + ligand + " -osdf " + sdf_path
    print(command)
    os.system(command)

    return sdf_path

def mol2_to_smi(ligand, tmp_folder_name, name):
    ligand_name = os.path.basename(os.path.normpath(ligand)).split(".", 1)[0].lower()
    smi_path = os.path.join(tmp_folder_name, ligand_name + ".smi")

    command = "babel -imol2 " + ligand + " -osmi " + smi_path
    print(command)
    os.system(command)

    return smi_path

def get_columns():
    return ['NN' + str(i+1) for i in range(350)]

def get_nn350_descriptors(receptor, ligand, tmp_folder_name, name, num = 58):

    # Load data
    receptor_extension = receptor.rsplit(".", 1)[1]
    rcptr = six.next(scoring.Descriptors.NNscore350.oddt.toolkit.readfile(receptor_extension, receptor, lazy = True))        
    #smi_path = mol2_to_smi(ligand = ligand, tmp_folder_name = tmp_folder_name, name = name)    
    sdf_path = mol2_to_sdf(ligand = ligand, tmp_folder_name = tmp_folder_name, name = name)    
    
    #ligand_extension = ligand.rsplit(".", 1)[1]
    #ligand_extension = smi_path.rsplit(".", 1)[1]
    ligand_extension = sdf_path.rsplit(".", 1)[1]
    
    #lgnd = six.next(scoring.Descriptors.NNscore350.oddt.toolkit.readfile(ligand_extension, ligand))
    #lgnd = six.next(scoring.Descriptors.NNscore350.oddt.toolkit.readfile(ligand_extension, smi_path))
    lgnd = six.next(scoring.Descriptors.NNscore350.oddt.toolkit.readfile(ligand_extension, sdf_path, lazy = True))    
    
    nnscore_descriptors = binana_descriptor(rcptr)
    features = nnscore_descriptors.build([lgnd])[0].tolist()

    dt = pd.DataFrame(columns = get_columns())
    dt.loc[name,:] = features

    dt["pdb"] = name
    dt.index = dt.pop("pdb")

    return dt
