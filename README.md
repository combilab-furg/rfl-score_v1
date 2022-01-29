# RFL-Score v1

This function defined as RFL-Score was trained using the PDBBind 2018 (general)  and CSAR (decoys  and NRC HiQ). From these complexes we extract molecular descriptors using BiNANA, RDKit 2D/3D, SASA, PaDEL, Vina, totaling 723 descriptors. Then, we selected the most promising attributes using LassoCV and parametrized a Random Forest (RF) algorithm using the GridSearchCV. This RF model was validated using CASF 2013 and 2016.

## Requirements
---------

  * Python 3.7+
  * MGLTools 1.5.6
  * vina4dv
  * DSSP 3.0.0+
  * OpenJDK 11.0.10+
  * Pandas 1.1.3+
  * Biopython 1.78+
  * RDKit 2020.09.1+
  * OpenBabel 2.4
  * Joblib 0.17.0+
  * Scikit-learn 0.23.2+
  * ODDT 0.7+
  * R 4.0.3
  * Random Forest Package 4.6.14+

## Install
---------

### Conda Environment and Packages

### Extras

Corresponds to software not found as an Anaconda package.

**MGLTools** was developed in the Sanner lab at the Center for Computational Structural Biology (CCRB) and visualization and analysis of molecular structures. The installer can be downloaded from https://ccsb.scripps.edu/mgltools/downloads/

The software **vina4dv** is a fork of AutoDock Vina modified to output features for **DeltaVinaRF20 score**. The repository can be downloaded from https://github.com/chengwang88/vina4dv

For both **MGLTools** and **vina4dv** are necessary to create environment variables in **.bashrc** (Linux) or **.bash_profile** (macOS) file in the home directory:

````
# Set MGLTools-1.5.6
export PATH=$PATH:/home/myuser/MGLTools-1.5.6/bin
export MGL=/home/myuser/MGLTools-1.5.6/
export MGLPY=$MGL/bin/python
export MGLUTIL=$MGL/MGLToolsPckgs/AutoDockTools/Utilities24/

# Set vina4v
export VINADIR=/home/myuser/vina4dv/build/linux/release/ 
````

## Scoring Function
---------

### Score

### Descriptors

## Reference
---------

