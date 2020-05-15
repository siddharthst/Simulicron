  

# Readme

## Setting up the development enviroment

Requirments:

- Python >3.6

- numpy >1.18.1

- pandas >1.0.3

- click >7.1.2

  

Conda is used throughout to provde package consitency and enviroment tracking. Ideally a conda installation will preserve the system native python installation and will be installed in the local directory. I would suggest using MiniForge instead of MiniConda and traditional Anaconda distribution.

- To begin, get the latest release

`$ wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge-pypy3-Linux-x86_64.sh`

- Install it with defaults

`$ ./Miniforge-pypy3-Linux-x86_64.sh`

- Clone the current repository

`$ git clone git@github.com:siddharthst/Simulicron.git`

- Use the `environment_working.yaml` in working directory to initialize a new enviroment.

`$ conda env create --name simulicron --file=environment_working.yaml`

- Activate the simulicron enviroment

`$ conda activate simulicron`

  

## Using simulicron

- The tool still lacks valid binary packaging or setup scheme, so working directly with source code is required. However, there is a file `/src/simulicronalpha/frontEnd.py` which acts as a wrapper around the main program.

- Testing it is trivial, with `frontEnd.py --help` displaying options. Currently, only sane options are displayed, and more options can be explored in `popsim.py`. If `frontEnd.py` is executed without any argument, it will simply display help.

- The output is stored in a [python pickle file](https://docs.python.org/3/library/pickle.html).

 ## Parsing output
 - The tool outputs a pickle file in form of a list of tuples. This list can be expanded into a DataFrame or other desired data structure. Each element of list can be further divided into following components (in order - acting as a pseudo row of a DataFrame) :
 
 **Simulation state**|**Fixed TEs**|**Transient TEs**|**Lost TEs**|**Simulation exit generation**|**Total number of transposition events in simulation**|**Average copy number**|**Copy number variance**
:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:
String – Fixed, Lost or in Flux|List of fixed transposon insertions (datatype int - containing transposon ID)|List of transposon insertions in flux (datatype int - containing transposon ID)|List of transposon insertions lost in simulation (datatype int - containing transposon ID)|Integer|Integer – total number of transposition events|List of type float (average copy number at each generation)|List of type float (copy number variance at each generation) 

### Work in progress

- Interactive parser for the output file.

