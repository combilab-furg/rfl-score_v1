import os, sys
import pandas as pd
import re

# https://github.com/PDB-REDO/dssp/issues/1
# CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1

cryst1 = 'CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1'

def get_columns():
    dssp = ['D' + str(i+1) for i in range(34)]
    return dssp

def add_cryst1_line(receptor, tmp_folder_name, name):
    
    pdb_tmp_path = os.path.join(tmp_folder_name, name + '_temp.pdb');
    
    flag = 0;
    with open(pdb_tmp_path, 'a') as tmp:
        with open(receptor, 'r') as f:
            for line in f:
                if ('ATOM' in line) and flag == 0:                    
                    line = cryst1 + '\n' + line
                    flag = flag + 1
                tmp.write(line)
    
    return pdb_tmp_path

def execute_dssp(receptor, tmp_folder_name, name):

    pdb_tmp_path = add_cryst1_line(receptor, tmp_folder_name, name);
    
    command = "dssp -i " + pdb_tmp_path + " -o " + os.path.join(tmp_folder_name, name + ".dssp")
    print(command)
    os.system(command)    

def get_descriptors(tmp_folder_name, name):
    file_path = os.path.join(tmp_folder_name, name + ".dssp")

    df = pd.DataFrame(columns = get_columns())

    # Read file
    f = open(file_path, "r")
    lines = f.readlines()

    indx = name

    # Line 6
    line = re.split(r'\s+', lines[6])

    # TOTAL NUMBER OF RESIDUES
    df.loc[indx, "D1"] = line[1]

    # NUMBER OF CHAINS
    df.loc[indx, "D2"] = line[2]

    # NUMBER OF SS-BRIDGES_TOTAL
    df.loc[indx, "D3"] = line[3]

    # NUMBER OF SS-BRIDGES_INTRACHAIN
    df.loc[indx, "D4"] = line[4]

    # NUMBER OF SS-BRIDGES_INTERCHAIN
    df.loc[indx, "D5"] = line[5]


    # Line 7
    line = re.split(r'\s+', lines[7])

    # ACCESSIBLE SURFACE OF PROTEIN (ANGSTROM**2)
    if not line[1] == "ACCESSIBLE":
        df.loc[indx, "D6"] = line[1]
    else:
        df.loc[indx, "D6"] = 0

    # Line 8
    line = re.split(r'\s+', lines[8])

    # TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(J)
    df.loc[indx, "D7"] = line[1]

    # TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(J)  , SAME NUMBER PER 100 RESIDUES
    df.loc[indx, "D8"] = line[2]


    # Line 9
    line = re.split(r'\s+', lines[9])

    # TOTAL NUMBER OF HYDROGEN BONDS IN     PARALLEL BRIDGES
    df.loc[indx, "D9"] = line[1]

    # TOTAL NUMBER OF HYDROGEN BONDS IN     PARALLEL BRIDGES, SAME NUMBER PER 100 RESIDUES
    df.loc[indx, "D10"] = line[2]


    # Line 10
    line = re.split(r'\s+', lines[10])

    # TOTAL NUMBER OF HYDROGEN BONDS IN ANTIPARALLEL BRIDGES
    df.loc[indx, "D11"] = line[1]

    # TOTAL NUMBER OF HYDROGEN BONDS IN ANTIPARALLEL BRIDGES, SAME NUMBER PER 100 RESIDUES 
    df.loc[indx, "D12"] = line[2]


    # Line 11
    line = re.split(r'\s+', lines[11])

    # TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I-5)
    df.loc[indx, "D13"] = line[1]

    # TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I-5), SAME NUMBER PER 100 RESIDUES                              .
    df.loc[indx, "D14"] = line[2]


    # Line 12
    line = re.split(r'\s+', lines[12])

    # TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I-4)
    df.loc[indx, "D15"] = line[1]

    # TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I-4), SAME NUMBER PER 100 RESIDUES                              .
    df.loc[indx, "D16"] = line[2]

    # Line 13
    line = re.split(r'\s+', lines[13])

    # TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I-3)
    df.loc[indx, "D17"] = line[1]

    # TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I-3), SAME NUMBER PER 100 RESIDUES                              .
    df.loc[indx, "D18"] = line[2]


    # Line 14
    line = re.split(r'\s+', lines[14])

    # TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I-2)
    df.loc[indx, "D19"] = line[1]

    # TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I-2), SAME NUMBER PER 100 RESIDUES                              .
    df.loc[indx, "D20"] = line[2]


    # Line 15
    line = re.split(r'\s+', lines[15])

    # TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I-1)
    df.loc[indx, "D21"] = line[1]

    # TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I-1), SAME NUMBER PER 100 RESIDUES                              .
    df.loc[indx, "D22"] = line[2]


    # Line 16
    line = re.split(r'\s+', lines[16])

    # TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I+0)
    df.loc[indx, "D23"] = line[1]

    # TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I+0), SAME NUMBER PER 100 RESIDUES                              .
    df.loc[indx, "D24"] = line[2]


    # Line 17
    line = re.split(r'\s+', lines[17])

    # TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I+1)
    df.loc[indx, "D25"] = line[1]

    # TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I+1), SAME NUMBER PER 100 RESIDUES                              .
    df.loc[indx, "D26"] = line[2]


    # Line 18
    line = re.split(r'\s+', lines[18])

    # TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I+2)
    df.loc[indx, "D27"] = line[1]

    # TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I+2), SAME NUMBER PER 100 RESIDUES                              .
    df.loc[indx, "D28"] = line[2]


    # Line 19
    line = re.split(r'\s+', lines[19])

    # TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I+3)
    df.loc[indx, "D29"] = line[1]

    # TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I+3), SAME NUMBER PER 100 RESIDUES                              .
    df.loc[indx, "D30"] = line[2]


    # Line 20
    line = re.split(r'\s+', lines[20])

    # TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I+4)
    df.loc[indx, "D31"] = line[1]

    # TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I+4), SAME NUMBER PER 100 RESIDUES                              .
    df.loc[indx, "D32"] = line[2]


    # Line 21
    line = re.split(r'\s+', lines[21])

    # TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I+5)
    df.loc[indx, "D33"] = line[1]

    # TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I+5), SAME NUMBER PER 100 RESIDUES                              .
    df.loc[indx, "D34"] = line[2]

    f.close()
    
    df["pdb"] = name
    df.index = df.pop("pdb")
    
    return df

def get_dssp34_descriptors(receptor, tmp_folder_name, name):
    
    execute_dssp(receptor = receptor, tmp_folder_name = tmp_folder_name, name = name)
    dt = get_descriptors(tmp_folder_name = tmp_folder_name, name = name)
    
    return dt
