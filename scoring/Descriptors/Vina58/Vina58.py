import os, sys
import numpy as np
import pandas as pd

def get_columns():
    return ['F' + str(i+1) for i in range(58)]

def get_vina_terms():
    vina_terms = ['ad4_solvation(charge=T)','ad4_solvation(charge=F)','electrostatic(x=1)','electrostatic(x=2)','gauss(0,0.3)','gauss(0.5,0.3)','gauss(1,0.3)', 'gauss(1.5,0.3)', 'gauss(2,0.3)', 'gauss(2.5,0.3)', 'gauss(0,0.5)', 'gauss(1,0.5)', 'gauss(2,0.5)', 'gauss(0,0.7)', 'gauss(1,0.7)', 'gauss(2,0.7)', 'gauss(0,0.9)', 'gauss(1,0.9)', 'gauss(2,0.9)', 'gauss(3,0.9)', 'gauss(0,1.5)', 'gauss(1,1.5)', 'gauss(2,1.5)', 'gauss(3,1.5)', 'gauss(4,1.5)', 'gauss(0,2)', 'gauss(1,2)', 'gauss(2,2)', 'gauss(3,2)', 'gauss(4,2)', 'gauss(0,3)', 'gauss(1,3)', 'gauss(2,3)', 'gauss(3,3)', 'gauss(4,3)', 'repulsion(0.4)', 'repulsion(0.2)', 'repulsion(0.0)', 'repulsion(-0.2)', 'repulsion(-0.4)','repulsion(-0.6)', 'repulsion(-0.8)', 'repulsion(-1.0)', 'hydrophobic(0.5,1)', 'hydrophobic(0.5,1.5)', 'hydrophobic(0.5,2)', 'hydrophobic(0.5,3)', 'non_hydrophobic(0.5,1.5)', 'vdw(4,8)', 'non_dir_h_bond(-0.7,0)', 'non_dir_h_bond(-0.7,0.2)', 'non_dir_h_bond(-0.7,0.4)', 'num_tors', 'num_rotors', 'num_heavy_atoms', 'num_hydrophobic_atoms', 'ligand_max_num_h_bonds', 'ligand_length']
    return vina_terms
    
def runVina(receptor_pdbqt, ligand_pdbqt, tmp_folder_name):
    
    score_path = os.path.join(tmp_folder_name, "score_v1.tmp")
    cmd = "$VINADIR/vina --receptor " + receptor_pdbqt + " --ligand " + ligand_pdbqt + \
          " --score_only > " + score_path
    os.system(cmd)
    
    vinalist = []
    with open(score_path, "r") as f:
        for lines in f:
            if lines[0:4] in ["Affi", "Term"]:
                vinalist.append(float(lines.split()[1]))
    
    return vinalist
    
def prepare_receptor(receptor, tmp_folder_name):
    #protpdbqt = fprot + ".pdbqt"
    fprot, __ = os.path.splitext(os.path.basename(receptor))
    receptor_pdbqt = os.path.join(tmp_folder_name, fprot + ".pdbqt")

    out1_path = os.path.join(tmp_folder_name, "out1.tmp")

    cmd = "$MGLPY $MGLUTIL/prepare_receptor4.py -r "  + receptor + \
          " -o " + receptor_pdbqt + " -U 'nphs' > " + out1_path
    os.system(cmd)

    return receptor_pdbqt

def prepare_ligand(ligand, tmp_folder_name):
    #ligpdbqt = flig + ".pdbqt"
    flig, __ = os.path.splitext(os.path.basename(ligand))
    ligand_pdbqt = os.path.join(tmp_folder_name, flig + ".pdbqt")
    
    out2_path = os.path.join(tmp_folder_name, "out2.tmp")

    cmd = "$MGLPY $MGLUTIL/prepare_ligand4.py -l " + ligand  + \
          " -o " + ligand_pdbqt +  " -U 'nphs' > " + out2_path
    os.system(cmd)

    return ligand_pdbqt

def featureVina(receptor, ligand, tmp_folder_name, name):    
    #protpdbqt = fprot + ".pdbqt"
    #protpdbqt = fprot + ".pdbqt"
    
    receptor_pdbqt = prepare_receptor(receptor = receptor, tmp_folder_name = tmp_folder_name)
    ligand_pdbqt = prepare_ligand(ligand = ligand, tmp_folder_name = tmp_folder_name)
    vinalist = runVina(receptor_pdbqt = receptor_pdbqt, ligand_pdbqt = ligand_pdbqt, tmp_folder_name = tmp_folder_name)
    
    # convert vina score and interaction term to be pKd unit
    c = -0.73349
    vinalist[0:53] = [c * i for i in vinalist[0:53]]
            
    return vinalist

def get_vina58_descriptors(receptor, ligand, tmp_folder_name, name, num = 58):
    vinalist = featureVina(receptor = receptor, ligand = ligand, tmp_folder_name = tmp_folder_name, name = name)
    
    vinaScore = vinalist[0]
    vinaFeatures = vinalist[1:]

    # Vina Features
    dt_features = pd.DataFrame(columns = get_columns())
    dt_features.loc[name,:] = vinaFeatures

    dt_features["pdb"] = name
    dt_features.index = dt_features.pop("pdb")

    # Vina Score
    dt_score = pd.DataFrame(columns = ["vina"])
    dt_score.loc[name,"vina"] = vinaScore

    dt_score["pdb"] = name
    dt_score.index = dt_score.pop("pdb")

    return dt_features, dt_score

