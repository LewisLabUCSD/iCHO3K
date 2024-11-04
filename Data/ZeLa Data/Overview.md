# Overview
This folder contains all of the information required to do the labelling and consumption/production rate analysis.

## Script Overview

`Labelling.ipynb`
The script that generates the enrichment analysis of the TCA cycle intermediates normalised to 
the extracellular glucose fraction. If the measured value reaches 1 that means that the metabolite
is both a) at a steady-state labelling, and b) completely derived from glucose. 

`Rate_calculations.ipynb`
This script determines uptake/production and growth rates for the CHO-S wt and CHO-ZeLa cell lines. 
We conducted a Bayesian regression to improve the confidence intervals of the parameters of interest. 
The unit of measurement being mol/mL_cell/h was used to adjust for differences in cell diameter between cell lines. Additionally, the bayesian regression struggled to fit the model when using cell count. However, both models are included.

There is a summary of the sampling intervals in the jupyter notebook.

## Data Overview

`extracellular_labelling_mvd`
GC-MS measured extracellular glucose concentrations directly after addtion of the label. This was 
used in order to normalise the intracellular metabolite pools.

`Results_17_11_23.csv`
LC-MS measurements of the different isotopes of intracellular metabolites.

`Timcouse.xlsx`
Cell characterisation, extracellular metabolite concentrations, and times of sampling.