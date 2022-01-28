import os, sys
import pandas as pd

file_name_160 = os.path.join(__file__.rsplit("/", 1)[0], "features_selected_160.csv")
file_name_250 = os.path.join(__file__.rsplit("/", 1)[0], "features_selected_250.csv")


def load_features_selected_160():
    
    features = pd.read_csv(file_name_160)
    features_selected_160 = features.iloc[:,0].tolist()

    return features_selected_160

def load_features_selected_250():
    
    features = pd.read_csv(file_name_250)
    features_selected_250 = features.iloc[:,0].tolist()

    return features_selected_250