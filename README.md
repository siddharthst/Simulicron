  

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

