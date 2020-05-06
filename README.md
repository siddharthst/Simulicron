
# Readme
## Setting up the development enviroment
Requirments:
- Python >3.6
- numpy  >1.18.1
- pandas >1.0.3
- click  >7.1.2

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

### Work in progress
- Interactive parser for the output file.

