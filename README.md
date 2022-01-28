# RFL-Score v1

This function defined as RFL-Score was trained using the PDBBind 2018 (general)  and CSAR (decoys  and NRC HiQ). From these complexes we extract molecular descriptors using BiNANA, RDKit 2D/3D, SASA, PaDEL, Vina, totaling 723 descriptors. Then, we selected the most promising attributes using LassoCV and parametrized a Random Forest (RF) algorithm using the GridSearchCV. This RF model was validated using CASF 2013 and 2016.

## Requirements
<<<<<<< HEAD


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

### Essential

### Conda Environment and Packages

### Extras

## Scoring Function

### Score

### Descriptors

## Reference

=======
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

### Essential

### Conda Environment and Packages

### Extras

## Scoring Function
---------

### Score

### Descriptors

## Reference
---------
>>>>>>> 88e3bf040c8598a749e66ce77eb2c7a13eb6a5d3
