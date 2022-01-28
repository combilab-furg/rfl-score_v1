from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors3D
from rdkit.Chem import rdMolDescriptors

import pandas as pd
import os, sys

randomSeed = 123456

def get_columns():
    rdkit3D = ["PBF", "PMI1", "PMI2", "PMI3", "NPR1", "NPR2", "RadiusOfGyration", "InertialShapeFactor", "Eccentricity", "Asphericity", "SpherocityIndex"]
    return rdkit3D

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

    m = AllChem.MolFromMolFile(mol_path, sanitize=False)

    sanitize = False
    molecule_conformers = True
    
    try:        
        AllChem.SanitizeMol(m)
        sm = AllChem.MolToSmiles(m)
        print("Sanitize Mol: True")
        sanitize = True
    except ValueError as e:
        sm = str(e)
        print("Sanitize Mol: False")
    
    if sanitize:

        m = AllChem.AddHs(m)

        try:    
            AllChem.EmbedMolecule(m, randomSeed=randomSeed)
        except:
            # Bad bond order
            molecule_conformers = False

        if molecule_conformers:
            try:
                df.loc[str(indx), "PBF"] = rdMolDescriptors.CalcPBF(m)
            except:
                df.loc[str(indx), "PBF"] = 0
                molecule_conformers = False

            try:
                df.loc[str(indx), "PMI1"] = Descriptors3D.PMI1(m)
            except:
                df.loc[str(indx), "PMI1"] = 0
                molecule_conformers = False

            try:
                df.loc[str(indx), "PMI2"] = Descriptors3D.PMI2(m)
            except:
                df.loc[str(indx), "PMI2"] = 0
                molecule_conformers = False

            try:			
                df.loc[str(indx), "PMI3"] = Descriptors3D.PMI3(m)
            except:
                df.loc[str(indx), "PMI3"] = 0
                molecule_conformers = False

            try:			
                df.loc[str(indx), "NPR1"] = Descriptors3D.NPR1(m)
            except:
                df.loc[str(indx), "NPR1"] = 0
                molecule_conformers = False

            try:			
                df.loc[str(indx), "NPR2"] = Descriptors3D.NPR2(m)
            except:
                df.loc[str(indx), "NPR2"] = 0
                molecule_conformers = False

            try:			
                df.loc[str(indx), "RadiusOfGyration"] = Descriptors3D.RadiusOfGyration(m)
            except:
                df.loc[str(indx), "RadiusOfGyration"] = 0
                molecule_conformers = False

            try:			
                df.loc[str(indx), "InertialShapeFactor"] = Descriptors3D.InertialShapeFactor(m)
            except:
                df.loc[str(indx), "InertialShapeFactor"] = 0
                molecule_conformers = False

            try:		
                df.loc[str(indx), "Eccentricity"] = Descriptors3D.Eccentricity(m)
            except:
                df.loc[str(indx), "Eccentricity"] = 0
                molecule_conformers = False

            try:		
                df.loc[str(indx), "Asphericity"] = Descriptors3D.Asphericity(m)
            except:
                df.loc[str(indx), "Asphericity"] = 0
                molecule_conformers = False

            try:		
                df.loc[str(indx), "SpherocityIndex"] = Descriptors3D.SpherocityIndex(m)
            except:
                df.loc[str(indx), "SpherocityIndex"] = 0
                molecule_conformers = False
        else:
            df.loc[str(indx), "PBF"] = 0
            df.loc[str(indx), "PMI1"] = 0
            df.loc[str(indx), "PMI2"] = 0
            df.loc[str(indx), "PMI3"] = 0
            df.loc[str(indx), "NPR1"] = 0
            df.loc[str(indx), "NPR2"] = 0
            df.loc[str(indx), "RadiusOfGyration"] = 0
            df.loc[str(indx), "InertialShapeFactor"] = 0
            df.loc[str(indx), "Eccentricity"] = 0
            df.loc[str(indx), "Asphericity"] = 0
            df.loc[str(indx), "SpherocityIndex"] = 0
    else:
        #AllChem.EmbedMolecule(m, randomSeed=randomSeed)

        df.loc[str(indx), "PBF"] = 0

        df.loc[str(indx), "PMI1"] = 0
        df.loc[str(indx), "PMI2"] = 0
        df.loc[str(indx), "PMI3"] = 0

        df.loc[str(indx), "NPR1"] = 0
        df.loc[str(indx), "NPR2"] = 0

        df.loc[str(indx), "RadiusOfGyration"] = 0
        df.loc[str(indx), "InertialShapeFactor"] = 0
        df.loc[str(indx), "Eccentricity"] = 0
        df.loc[str(indx), "Asphericity"] = 0
        df.loc[str(indx), "SpherocityIndex"] = 0

    #df.loc[str(indx), "random_seed"] = randomSeed
    #df.loc[str(indx), "sanitize"] = sanitize
    #df.loc[str(indx), "molecule_conformers"] = molecule_conformers
    #df.loc[str(indx), "protein"] = protein    

    df["pdb"] = name
    df.index = df.pop("pdb")

    return df

def get_rdkit3d11_descriptors(ligand, tmp_folder_name, name):
    
    mol_path = mol2_to_mol(ligand = ligand, tmp_folder_name = tmp_folder_name, name = name)
    dt = get_descriptors(mol_path = mol_path, name = name)

    return dt