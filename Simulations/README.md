# Simulations

This repository contains four directories:

* **figures**, which contains the R code to generate three figures (the 4 panels of Figure 1, the 2 panels of Figure 2, and the 2 panels of Figure 3). The figures (in pdf format) will be generated in this directory. 
* **results**, an empty directory in which the simulation results (pickle files) will be generated
* **run**, which contains the python code to run the simulations. Three scripts, corresponding to the data for the 3 figures.
* **src**, the python code for the simulation software. This directory contains its own README file with the requirements and dependencies for the simulations.

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
