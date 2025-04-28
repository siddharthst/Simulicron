# Simulations

This repository contains four directories:

* **figures**, which contains the R code to generate three figures (the 4 panels of Figure 1, the 2 panels of Figure 2, and the 2 panels of Figure 3). The figures (in pdf format) will be generated in this directory. 
* **results**, an empty directory in which the simulation results (pickle files) will be generated
* **run**, which contains the python code to run the simulations. Three scripts, corresponding to the data for the 3 figures.
* **src**, the python code for the simulation software.

## Dependencies and requirements
``
Simulations and figures were designed to be run on a standard Linux distribution (24.04.2 LTS, desktop or server version). 

### For simulations

Requirments:
* Python >3.6
* numpy >1.18.1
* pandas >1.0.3

**warning**: simulations were run on a decent server with more CPU and RAM than a desk computer. Change the "maxProcceses" (sic) parameter in Default.parameters to prevent crashes if necessary. 

**note**: simulations are quite computer intensive, then can take several hours/days to run depending on the number of CPUs involved. 

### For figures

Figures were drawn with a recent version of R (v4.1.2). The only dependency is package reticulate (version 1.39.0) to read the .pickle files. 


## Running simulations

Reproducing the figures from the paper requires two independent steps: 
1. Run the simulations:
```
  cd run
  python3 Figure1.py
  python3 Figure2.py
  python3 Figure3.py
```

2. Generate the figures
```
  cd ../figures
  Rscript Figure1-draw.R
  Rscript Figure2-draw.R
  Rscript Figure3-draw.R
```
