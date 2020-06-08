  

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

- The tool still lacks valid binary packaging or setup scheme, so working directly with source code is required. However, there is a file `/src/simulicronalpha/frontEnd.py` which acts as a wrapper around the main program. (This method is not recommended)

- The best way to use Simulicron is by invoking the functions directly in a python script. This can be done using the boilerplate code below:

```
# Init
import os
import sys
import pickle

# Path to source
module_path = os.path.abspath(os.path.join("../src/simulicronalpha/"))
if module_path not in sys.path:
    sys.path.append(module_path)


# Import the main function
from popSim import runBatch

# Output file name and location
output = "output.pickle"

# numberOfSimulations = Total number of simulations (int)
# baseSelection = Selection pen. at each site (Float)
# numberOfInsertionSites = Number of insertion sites (int)
# numberOfChromosomes = Number of chromosomes (int)
# baseRecombinationRate = Recombination rate at each insertion site (float)
# baseTau = Regulatory coef. (float)
# numberOfPiRNA = Number of piRNA (int) 
# piPercentage = Total piRNA length w.r.t to insertion sites (int)
# NumberOfIndividual = Number of individuals (int)
# NumberOfTransposonTypes = Number of transposon types (int)
# NumberOfInsertionsPerType = Number of insetions per type (List of type int)
# FrequencyOfInsertions = Number of insetions per type (List of type int)
# ExcisionRates = Excision rate per type (List of type float)
# RepairRates = Repair rate per type (List of type float)
# InsertionRates = Insertion rate per type (List of type float) 
# NumberOfGenerations = Number of generations (int)
# numberOfThreads = Number of threads (int)


sim   = runBatch(    
    numberOfSimulations=50,
    baseSelection=0,
    baseInsertionProb=1,
    numberOfInsertionSites=10000,
    numberOfChromosomes=6,
    baseRecombinationRate=0.01,
    baseTau=1,
    numberOfPiRNA=6,
    piPercentage=3,
    enablePiRecombination=False,
    NumberOfIndividual=1000,
    NumberOfTransposonTypes=1,
    NumberOfInsertionsPerType=[1],
    FrequencyOfInsertions=[1.0],
    ExcisionRates=[1.0],
    RepairRates=[1],
    InsertionRates=[1],
    NumberOfGenerations=10000,
    numberOfThreads=8,
)

with open(output, "wb") as f:
    pickle.dump((sim), f)
```


- The output is stored in a [python pickle file](https://docs.python.org/3/library/pickle.html).

 ## Parsing output
 - The tool outputs a pickle file in form of a dictionary with following keys:

**Key**|** Value**
:-----:|:-----:
`State`| FLUX or LOSS (STRING)
`Generatrion`| Simulation exit generation (INT)
`NTE`| Total number of transposition events (INT)
`AvgCopyNum`| Average copy number per generation [LIST]
`CopyNumVar`| Average copy number variation per generation [LIST]
### Work in progress

- Interactive parser for the output file.

