#!/usr/bin/env python

import os, sys, time
import pandas as pd
import shutil

from scoring.Descriptors import Descriptors as descriptors
from optparse import OptionParser

def check_folder(tmp_folder_name, create = True):
    if os.path.exists(tmp_folder_name):
        shutil.rmtree(tmp_folder_name)
    
    if create:
        os.makedirs(tmp_folder_name)

# options for read files
parser = OptionParser()
parser.add_option("-r", "--receptor",  help="receptor file path")
parser.add_option("-f", "--fasta",  help="fasta file path")
parser.add_option("-l", "--ligand", help="ligand file path")
parser.add_option("-n", "--name", help="name file")
parser.add_option("-o", "--output", help="output path")
parser.add_option("-g", "--group",  help="a list of receptor-ligand-fasta files")
parser.add_option("-x", "--features_list", help="features list")
parser.add_option("-m", "--receptor_mol2",  help="receptor mol2 file path")

(options, args) = parser.parse_args()
options =  options.__dict__

# Simple
if options['name'] is None:
    print("Please specify name file")
    sys.exit()
elif options['output'] is None:
    print("Please specify output file")
    sys.exit()
elif options['features_list'] is None:
    print("Please specify features list")
    sys.exit()

name = options['name']
output = options['output']

# Check features
features = options['features_list'].split(" ")
features = list(set(features))

features_names = descriptors.get_features_names()

x = list(set(features) - set(features_names))

if len(x) > 0:
    print("\"" + x[0] + "\" features not found")
    print("Features:")
    print(features_names)
    sys.exit()

# assign args to pdblist
#pdblist = []
pdblist = pd.DataFrame(columns = ["receptor", "receptor_mol2", "ligand", "fasta"])

if options['group'] is None:
    idx = len(pdblist)

    if options['receptor'] is None:
        print("Please specify receptor file")
        sys.exit()
    else:
        pdblist.loc[idx, "receptor"] = options['receptor']

    if options['ligand'] is None:
        print("Please specify ligand file")
        sys.exit()
    else:
        pdblist.loc[idx, "ligand"] = options['ligand']
    
    # Check fasta file
    pdblist.loc[idx, "fasta"] = None
    
    if "amino20" in features and options['fasta'] is None:
        print("Please specify fasta file")
        sys.exit()
    else:
        if not options['fasta'] is None:        
            pdblist.loc[idx, "fasta"] = options['fasta']
    
    # Check receptor_mol2 file
    if not options['receptor_mol2'] is None:        
        pdblist.loc[idx, "receptor_mol2"] = options["receptor_mol2"]
    else:
        pdblist.loc[idx, "receptor_mol2"] = None
else:    
    # one or more protein-ligand complex
    with open(options['group'], 'r') as f:        
        line = 0
        for lines in f:
            idx = len(pdblist)

            line += 1
            n_params = lines.split(",")

            if n_params[0] is None or n_params[0].strip() == "":
                print("Please specify receptor file. Line: " + str(line))
                sys.exit()
            else:
                pdblist.loc[idx, "receptor"] = n_params[0].strip()

            if n_params[1] is None or n_params[1].strip() == "":
                print("Please specify ligand file. Line: " + str(line))
                sys.exit()
            else:
                pdblist.loc[idx, "ligand"] = n_params[1].strip()

            if "amino20" in features:
                if len(n_params) < 3 or n_params[2] is None or n_params[2].strip() == "":
                    print("Please specify fasta file. Line: " + str(line))
                    sys.exit()
                else:
                    pdblist.loc[idx, "fasta"] = n_params[2].strip()
            else:
                pdblist.loc[idx, "fasta"] = None

            if len(n_params) == 4 and not n_params[3] is None and not n_params[3].strip() == "":
                pdblist.loc[idx, "receptor_mol2"] = n_params[3].strip()
            else:
                pdblist.loc[idx, "receptor_mol2"] = None

            print(n_params)
            #if n == 3 and len(n_params) < 3:
            #    print("Please specify fasta file. Line: " + str(line))
            #    sys.exit()
            
                #pdblist.append(lines.split()[0:n])
    
print(pdblist)

tmp_folder_name = "tmp"
check_folder(tmp_folder_name = tmp_folder_name)

descriptors.save_descriptors(descriptor_names = features, output = output, name = name, tmp_folder_name = tmp_folder_name, pdblist = pdblist)

# Delete foldet tmp
check_folder(tmp_folder_name = tmp_folder_name, create = False)