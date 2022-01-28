#!/usr/bin/env python

import os, sys, time
import pandas as pd
import shutil

from scoring.Descriptors import Descriptors as features
from scoring.FeatureSelection import FeatureSelection 
from scoring.RScript import RScript

from optparse import OptionParser

def check_folder(tmp_folder_name, create = True):
    if os.path.exists(tmp_folder_name):
        shutil.rmtree(tmp_folder_name)
    
    if create:
        os.makedirs(tmp_folder_name)

# Start Time 
inicia_exec = time.time()

# options for read files
parser = OptionParser()
parser.add_option("-r", "--receptor",  help="receptor file path")
parser.add_option("-f", "--fasta",  help="fasta file path")
parser.add_option("-l", "--ligand", help="ligand file path")
parser.add_option("-n", "--name", help="name file")
parser.add_option("-o", "--output", help="output path")
parser.add_option("-g", "--group",  help="a list of receptor-ligand-fasta files")
parser.add_option("-t", "--tmp_path",  help="tmp folder path")
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

name = options['name']
output = options['output']

# assign args to pdblist
#pdblist = []
pdblist = pd.DataFrame(columns = ["receptor", "receptor_mol2", "ligand", "fasta"])
if options['group'] is None:
    idx = len(pdblist)

    if options['receptor'] is None:
        print("Please specify receptor file")
        sys.exit()
    elif options['ligand'] is None:
        print("Please specify ligand file")
        sys.exit()
    elif options['fasta'] is None:
        print("Please specify fasta file")
        sys.exit()
    else:
        
        pdblist.loc[idx, "receptor"] = options['receptor']
        pdblist.loc[idx, "ligand"] = options['ligand']
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
            
            if n_params[2] is None or n_params[2].strip() == "":
                print("Please specify fasta file. Line: " + str(line))
                sys.exit()
            else:
                pdblist.loc[idx, "fasta"] = n_params[2].strip()
        
            if len(n_params) == 4 and not n_params[3] is None and not n_params[3].strip() == "":
                pdblist.loc[idx, "receptor_mol2"] = n_params[3].strip()
            else:
                pdblist.loc[idx, "receptor_mol2"] = None

print(pdblist)

tmp_folder_name = "tmp"

if not options['tmp_path'] is None:
    tmp_folder_name = os.path.join(options['tmp_path'], "tmp")

check_folder(tmp_folder_name = tmp_folder_name)

# Load Features Selected
features_selected = FeatureSelection.load_features_selected_160()
print(len(features_selected))

#dt = features.get_descriptors(descriptor_names = ["amino20"], tmp_folder_name = tmp_folder_name, receptor = receptor, ligand = ligand, fasta = fasta, idx = name, features_selected = features_selected)
#dt = features.get_descriptors(descriptor_names = ["amino20", "dssp34", "rdkt2d147", "rdkt3d11", "pd92", "sasa10", "vina58", "deltavina20", "nn350"], tmp_folder_name = tmp_folder_name, receptor = receptor, ligand = ligand, fasta = fasta, idx = name, features_selected = features_selected)
#dt, times = features.get_descriptors(descriptor_names = features.get_features_names(), tmp_folder_name = tmp_folder_name, receptor = receptor, ligand = ligand, fasta = fasta, idx = name, features_selected = features_selected)
dt, times = features.get_descriptors(descriptor_names = features.get_features_names(), tmp_folder_name = tmp_folder_name, pdblist = pdblist, features_selected = features_selected)

# Times 
times_path = os.path.join(output, name + "_descriptors_times.csv")
times.to_csv(times_path)

# Input
input_path = os.path.join(output, name + "_input.csv")

# Save
dt.to_csv(input_path)

# Output
output_path = os.path.join(output, name + "_output.csv")

# Script
script_path = RScript.make_model_script(input_path = input_path, output_path = output_path, n_features = len(features_selected), tmp_folder_name = tmp_folder_name)

# run R script
cmd = "R CMD BATCH " + script_path
print(cmd)
os.system(cmd)

check_folder(tmp_folder_name = tmp_folder_name, create = False)

# Finish Time
fim_exec = time.time()
time_total = fim_exec - inicia_exec

print("\nTotal Time: " + str(time_total))

# Save
time_total_dt = pd.DataFrame(columns = ["process", "time"])

time_total_dt.loc[0, "process"] = "sf250"
time_total_dt.loc[0, "time"] = time_total

time_total_path = os.path.join(output, name + "_process_time.csv")
time_total_dt.to_csv(time_total_path)