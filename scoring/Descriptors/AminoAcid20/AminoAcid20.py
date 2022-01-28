from Bio.SeqUtils.ProtParam import ProteinAnalysis 
from scoring.Descriptors.AminoAcid20 import tools as tls
import pandas as pd

def get_columns():
    amino_acid = ["A", "G", "M", "S", "C", "H", "N", "T", "D", "I", "P", "V", "E", "K", "Q", "W", "F", "L", "R", "Y"]
    return amino_acid

def get_amino20_descriptors(fasta, name):

    columns = get_columns().copy()
    dataset = pd.DataFrame(columns = columns)

    fasta_string = tls.return_fasta_string(fasta)
    print(name)
    print("FASTA String: " + name)
    print(fasta_string)
    X = ProteinAnalysis(fasta_string)

    index = name
    #dataset.loc[name, "protein"] = protein

    for aa in columns:
        print(aa + ": " + str(X.get_amino_acids_percent()[aa]))
        dataset.loc[index, aa] = X.get_amino_acids_percent()[aa]
    
    dataset["pdb"] = name
    dataset.index = dataset.pop("pdb")

    return dataset


