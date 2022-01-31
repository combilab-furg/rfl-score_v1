# RFL-Score v1

This function defined as RFL-Score was trained using the PDBBind 2018 (general)  and CSAR (decoys  and NRC HiQ). From these complexes we extract molecular descriptors using BiNANA, RDKit 2D/3D, SASA, PaDEL-Descriptor, Vina, totaling 723 descriptors. Then, we selected the most promising attributes using LassoCV and parametrized a Random Forest (RF) algorithm using the GridSearchCV. This RF model was validated using CASF 2013 and 2016.

## Requirements

  * Python 3.7+
  * MGLTools 1.5.6
  * vina4dv
  * DSSP 3.0.0+
  * JRE 11.0.10+
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

### Conda Environment
The current scoring function is executed with **Python 3.7**. The most simple way to install the SF is by creating an **Anaconda Environment**. To download the last version of **Anaconda Individal Edition** visit his official page: https://www.anaconda.com/produacts/individual

A proper environment in **Anaconda** is useful to install the SF. To create a new environment with **Python 3.7** execute the following statement in **Anaconda's prompt**:
````
conda create -n 3_7 python=3.7
````
To activate the new environment execute the following statement:
````
source activate 3_7
````
The new environment can be activated automatically adding the previous statement in the **.bashrc** file.

### Conda Packages

**Pandas**
It is a package for data analysis and manipulation tool. To install in **Conda** execute the following statement:
````
conda install -c anaconda pandas
````
**Biopython**
It is a set of freely available tools for biological computation. To install in **Conda** execute the following statement:
````
conda install -c anaconda biopython
````
**RDKit**
It is a cheminformatics software to create 2D and 3D molecular features for machine learning applications. To install in **Conda** execute the following statement:
````
conda install -c rdkit rdkit
````
**OpenBabel 2.4**
It is a chemical software designed to speak the many languages of chemical data. The **version 2.4** is mandatory to install the SF. To install in **Conda** execute the following statement:
````
conda install -c openbabel openbabel 
````
**ODDT (Open Drug Discovery Toolkit)**
It is a modular and comprehensive toolkit for use in cheminformatics, molecular modeling etc. To install in **Conda** execute the following statement:
````
conda install -c oddt oddt
````
**Scikit-learn**
It is a toolkit for predictive data analysis. To install in **Conda** execute the following statement:
````
conda install -c anaconda scikit-learn
````
**Joblib**
It is a suite to provide lightweight pipelining in Python. To install in **Conda** execute the following statement:
````
conda install -c anaconda joblib
````
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
The **DSSP (hydrogen bond estimation algorithm)** program was designed by Wolfgang Kabsch and Chris Sander and calculates the secondary structure assignment given the 3D structure of a protein. To install the **DSSP** run the following statement:
````
sudo apt-get install -y dssp
````
The **Java Runtime Environment** (JRE) is a software layer to run Java programs. **PaDEL-Descriptor** is a program written in Java, it is inside the implementation of the scoring function and it is used to calculate ligand features. To install the **JRE** run the following statement:
````
sudo apt install default-jre
````
**R** is an integrated suite of software facilities for data science. The scoring function developed for this work was trained with **R**. To install the **R** run the following statement:
````
sudo apt install r-base
````
The  machine learning **Random Forest** algorithm was used to train **RFL-Score**, therefore it is necessary to install its corresponding package within **R**. In the **R** console execute the following statement:
````
install.packages("randomForest")
````
## Scoring Function

The current repository includes the SF and features generation scripts.

### Score
To calculate the score of a protein-ligand complex.

Params:
  * **-r:** Protein's file path. **PDB** format is mandatory.
  * **-l:** Ligand's file path. **MOL2** format is mandatory.
  * **-o:** Output path.
  * **-n:** Output file name.
  * **-g:** Protein-ligand list file path. 
  * **-t:** Temporary folder path.

A score example:
````
python sf160.py -r examples/1a30/1a30_protein.pdb -f examples/1a30/1A30.fasta.txt -l examples/1a30/1a30_ligand.mol2 -n 1a30 -o /home/oarrua/Desktop/Test -t /home/oarrua/Desktop/Test
````
A protein-ligand list example:
````
python sf160.py -g examples/list_example -n all -o /home/oarrua/Desktop/Test
````

### Descriptors
To calculate a set of features of a protein-ligand complex.

**Feature Set**
* amino20: Related features of amino acid percentages.
* dssp34:  Related features of the secondary  structure  of  proteins.
* nn350: NNscore 2.0 features
* pd92: A set ligand features of PaDEL-Descriptors.
* rdkt2d147: A set of 2D ligand descriptors of RDKit.
* rdkt3d11: A set of 3D ligand descriptors of RDKit.
* sasa10: Related features of solvent-accessible  surface  area.
* vina58: Vina  implemented  terms.
* deltavina20: Autodock Vina score.

Params:
  * **-r:** Protein's file path. **PDB** format is mandatory.
  * **-f:** Protein Sequence's file path. **FASTA** format is mandatory.
  * **-l:** Ligand's file path. **MOL2** format is mandatory.
  * **-o:** Output path.
  * **-n:** Output file name.
  * **-g:** Protein-ligand list file path. 
  * **-t:** Temporary folder path.
  
A simple example:
````
python sf723_descriptors.py -r examples/1a30/1a30_protein.pdb -f examples/1a30/1A30.fasta.txt -l examples/1a30/1a30_ligand.mol2 -n 1a30 -o /home/oarrua/Desktop/Test -x "amino20 dssp34 nn350 pd92 rdkt2d147 rdkt3d11 sasa10 vina58 deltavina20"
````
A protein-ligand list example:
````
python sf723_descriptors.py -g examples/list_example_with_amino20 -n all -o /home/oarrua/Desktop/Test -x "amino20 dssp34 nn350 pd92 rdkt2d147 rdkt3d11 sasa10 vina58 deltavina20"
````


### Videos
## Reference
