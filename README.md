# Whole Cell Network Reconstruction for CHO cells

## Solver Install
Set up the qMINOS solver. You will need the qminos file, which can be obtained for academic use from Prof. Michael Saunders at Stanford University.

&emsp;i) download the qminos file into a specified solver directory, which we refer to here as "solver_parent_directory".

&emsp;ii) the solver can be installed using the human_me Makefile as follows:
 
```console
make -C <path/to/human_me/> install-qminos SOLVER_PATH=<path/to/solver/solver_parent_directory>
```
