import os, sys
import pandas as pd
import time
import numpy as np

from scoring.Descriptors.AminoAcid20 import AminoAcid20
from scoring.Descriptors.DSSP34 import DSSP34
from scoring.Descriptors.RDKit2D147 import RDKit2D147
from scoring.Descriptors.RDKit3D11 import RDKit3D11
from scoring.Descriptors.PaDEL92 import PaDEL92
from scoring.Descriptors.Sasa10 import Sasa10
from scoring.Descriptors.Vina58 import Vina58
from scoring.Descriptors.NNscore350 import NNscore350

def get_features_names():
    names = ["amino20", "dssp34", "nn350", "pd92", "rdkt2d147", "rdkt3d11", "sasa10", "vina58", "deltavina20"]
    return names

def get_descriptor(descriptor_name, tmp_folder_name, receptor, ligand, fasta, receptor_mol2, idx):
    if descriptor_name == "amino20":
        return AminoAcid20.get_amino20_descriptors(fasta = fasta, name = idx)
    elif descriptor_name == "dssp34":
        return DSSP34.get_dssp34_descriptors(receptor = receptor, tmp_folder_name = tmp_folder_name, name = idx)
    elif descriptor_name == "rdkt2d147":  
        return RDKit2D147.get_rdkit2d147_descriptors(ligand = ligand, tmp_folder_name = tmp_folder_name, name = idx)
    elif descriptor_name == "rdkt3d11":  
        return RDKit3D11.get_rdkit3d11_descriptors(ligand = ligand, tmp_folder_name = tmp_folder_name, name = idx)
    elif descriptor_name == "pd92":  
        return PaDEL92.get_pd92_descriptors(ligand = ligand, tmp_folder_name = tmp_folder_name, name = idx)
    elif descriptor_name == "sasa10":  
        return Sasa10.get_sasa10_descriptors(receptor = receptor, ligand = ligand, tmp_folder_name = tmp_folder_name, name = idx)
    elif descriptor_name == "vina58":  
        dt_features, dt_score = Vina58.get_vina58_descriptors(receptor = receptor, ligand = ligand, tmp_folder_name = tmp_folder_name, name = idx)
        return dt_features
    elif descriptor_name == "deltavina20":  
        dt_features, dt_score = Vina58.get_vina58_descriptors(receptor = receptor, ligand = ligand, tmp_folder_name = tmp_folder_name, name = idx)
        return dt_score
    elif descriptor_name == "nn350":  
        if not receptor_mol2 is None:
            receptor = receptor_mol2                        
        return NNscore350.get_nn350_descriptors(receptor = receptor, ligand = ligand, tmp_folder_name = tmp_folder_name, name = idx)
    else:
        print(descriptor_name + " not found")

#def get_descriptors(descriptor_names, tmp_folder_name, receptor, ligand, fasta, idx, features_selected = []):
def get_descriptors(descriptor_names, tmp_folder_name, pdblist, features_selected = []):    
    
    final_descriptors_set = pd.DataFrame()
    descriptors_set = pd.DataFrame()
    
    final_times = pd.DataFrame()
    
    # Dataframe Times columns
    times_columns = get_features_names()
    times_columns.append("total")
    times = pd.DataFrame(columns = times_columns)    

    for ix in pdblist.index.tolist():
        time_total = 0
        receptor = pdblist.loc[ix, "receptor"]
        ligand = pdblist.loc[ix, "ligand"]
        fasta = pdblist.loc[ix, "fasta"]
        receptor_mol2 = pdblist.loc[ix, "receptor_mol2"]
        idx = os.path.basename(os.path.normpath(receptor)).split(".", 1)[0].lower() + "+" + os.path.basename(os.path.normpath(ligand)).split(".", 1)[0].lower()
        
        n = 0
        #for name in names in order:
        for i in range(0, len(descriptor_names)):
            # Start Time
            inicia_exec = time.time()

            n+=1        
            descriptor_name = descriptor_names[i]
            
            if n == 1:
                descriptors_set = get_descriptor(descriptor_name = descriptor_name, tmp_folder_name = tmp_folder_name, receptor = receptor, ligand = ligand, fasta = fasta, receptor_mol2 = receptor_mol2, idx = idx)
            else:
                descriptors_set = descriptors_set.merge(get_descriptor(descriptor_name = descriptor_name, tmp_folder_name = tmp_folder_name, receptor = receptor, ligand = ligand, fasta = fasta, receptor_mol2 = receptor_mol2, idx = idx), on="pdb", how="inner")

            # Finish Time
            fim_exec = time.time()
            time_descriptor = fim_exec - inicia_exec

            print("\n Descriptor: " + descriptor_name + " Total Time: " + str(time_descriptor))
            
            times.loc[idx, descriptor_name] = time_descriptor
            time_total += time_descriptor

        # Total Time
        times.loc[idx, "total"] = time_total
        
        # Features Selected
        if not len(features_selected) == 0:        
            descriptors_set = descriptors_set.loc[:, features_selected]

        # Add into final set
        final_descriptors_set = final_descriptors_set.append(descriptors_set)
        final_times = final_times.append(times)
        times = times.iloc[0:0]
    
    print(final_descriptors_set)
    print(final_times)

    return final_descriptors_set, final_times


def save_descriptors(descriptor_names, output, name, tmp_folder_name, pdblist):
    final_descriptors_set = pd.DataFrame()

    n = 0
    #for name in order:
    for i in range(0, len(descriptor_names)):        
        
        n += 1

        # Start Time
        inicia_exec = time.time()
               
        descriptor_name = descriptor_names[i]
        
        print("Descriptor: " + descriptor_name)
        
        descriptors_set = pd.DataFrame()

        #n = 0
        for ix in pdblist.index.tolist():
            
            receptor = pdblist.loc[ix, "receptor"]
            ligand = pdblist.loc[ix, "ligand"]            
            receptor_mol2 = pdblist.loc[ix, "receptor_mol2"]
            
            fasta = None
            if "amino20" in descriptor_names:
                fasta = pdblist.loc[ix, "fasta"]            

            idx = os.path.basename(os.path.normpath(receptor)).split(".", 1)[0].lower()
            descriptors_set = descriptors_set.append(get_descriptor(descriptor_name = descriptor_name, tmp_folder_name = tmp_folder_name, receptor = receptor, ligand = ligand, fasta = fasta, receptor_mol2 = receptor_mol2, idx = idx))

        print(descriptors_set)
        
        # Save descriptors
        descriptors_set.to_csv(os.path.join(output, name + "_" + descriptor_name + ".csv"))

        # Finish Time
        fim_exec = time.time()
        time_total = fim_exec - inicia_exec
        print("\nDescriptor: " + descriptor_name + " Total Time: " + str(time_total))
    
        if n == 1:          
            final_descriptors_set = descriptors_set
        else:                    
            final_descriptors_set = final_descriptors_set.merge(descriptors_set, on="pdb", how="inner")
        
    # Save final descriptors set
    if len(descriptor_names) > 1:
        final_descriptors_set.to_csv(os.path.join(output, name + "_all.csv"))
        