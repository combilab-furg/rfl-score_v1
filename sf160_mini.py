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
parser.add_option("-n", "--name", help="name file")
parser.add_option("-i", "--input", help="input path")
parser.add_option("-o", "--output", help="output path")
parser.add_option("-t", "--tmp_path",  help="tmp folder path")

(options, args) = parser.parse_args()
options =  options.__dict__

# Simple
if options['name'] is None:
    print("Please specify name file")
    sys.exit()
elif options['input'] is None:
    print("Please specify input file")
    sys.exit()
elif options['output'] is None:
    print("Please specify output file")
    sys.exit()

name = options['name']
input_path = options['input']
output = options['output']

tmp_folder_name = "tmp"

if not options['tmp_path'] is None:
    tmp_folder_name = os.path.join(options['tmp_path'], "tmp")

check_folder(tmp_folder_name = tmp_folder_name)

# Load Features Selected
features_selected = FeatureSelection.load_features_selected_160()
print(len(features_selected))

# Output
output_path = os.path.join(output, name + "_output.csv")

# Script
script_path = RScript.make_model_script(input_path = input_path, output_path = output_path, n_features = len(features_selected), tmp_folder_name = tmp_folder_name)

# run R script
cmd = "R CMD BATCH " + script_path
print(cmd)
os.system(cmd)

check_folder(tmp_folder_name = tmp_folder_name, create = False)