# Whole Cell Network Reconstruction for CHO cells

## Requirements
setuptools>=65.4

## Solver Install
Set up the qMINOS solver. You will need the qminos file, which can be obtained for academic use from Prof. Michael Saunders at Stanford University.


&emsp;i) download the qminos file into a specified solver directory, which we refer to here as "solver_parent_directory".

&emsp;ii) the solver can be installed using the conda environment as follows:

```console
conda activate <your_environment_name>
 ```
 
```console
python install_solver_conda.py <path/to/solver/solver_parent_directory> qminos
```
