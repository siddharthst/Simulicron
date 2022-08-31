  

# Readme

## Setting up the development enviroment

Requirments:

- Python >3.6

- numpy >1.18.1

- pandas >1.0.3

Conda is used throughout to provde package consitency and enviroment tracking. Ideally a conda installation will preserve the system native python installation and will be installed in the local directory. I would suggest using MiniForge instead of MiniConda or traditional Anaconda distribution.

- To begin, get the latest release

`$ wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge-pypy3-Linux-x86_64.sh`

- Install it with defaults

`$ ./Miniforge-pypy3-Linux-x86_64.sh`

- Clone the current repository

`$ git clone git@github.com:siddharthst/Simulicron.git`

- Use the `environment_working.yaml` in working directory to initialize a new enviroment.

`$ conda env create --name simulicron --file=Simulicron_env.yaml`

- Activate the simulicron enviroment

`$ conda activate simulicron`

  

## Using simulicron

- The best way to use Simulicron is by invoking the functions directly in a python script. Code used for generating figures in the folder `./Figures` can be adapted as a boilerplate.


- The output is stored in a [python pickle file](https://docs.python.org/3/library/pickle.html).

 ## Parsing output
 - The tool outputs a pickle file in form of a dictionary with following keys:

     simulationState,
    i,
    numberOfTranspositionEvents,
    averageCopyNumber,
    varianceCopyNumber,
    TEfamilyCountArrRes,
    TEfamilyVarArrRes,
    TEregulationArrRes,
    avgFitness,
    HMTgen,
    eta,
    NumberOfTransposonInsertions,
    FrequencyOfInsertions,
    ExcisionRates,
    tau,
    selPen,
    piRNAindices,
    overlap,
    fitnessFunction,
    epistasisCoefficient,
    TEset,
    TECoreOverlap,


**Key**|** Value**
:-----:|:-----:
`State`| FLUX or LOSS (STRING)
`Generatrion`| Simulation exit generation (INT)
`NTE`| Total number of transposition events (INT)
`AvgCopyNum`| Average copy number per generation [LIST]
`CopyNumVar`| Average copy number variation per generation [LIST]
`TEfamilyCN`| Average copy number per generation per TE {DICT},
`TEfamilyVR`| Average copy number variation per generation per TE {DICT},
`TEfamilyRg`| Average regulation per TE {DICT},

