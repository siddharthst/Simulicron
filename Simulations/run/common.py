from generateSim import generatePopulation, generateGenome, initHGT
from popSim import runSim
import random
import pickle


# Wrapper function for multiprocessing
def worker(parameters, outputdir):
	    
    # Generate genome and population
    genome, piset, piIndice, rates = generateGenome(
        numberOfInsertionSites = parameters["NumberOfInsertionSites"],
        numberOfChromosomes    = parameters["NumberOfChromosomes"],
        baseRecombinationRate  = parameters["RecombinationRate"],
        baseSelection          = parameters["selectionPenalty"],
        baseTau                = parameters["tau"],
        disablePiRecombination = True,
    )
        
    population, transposons, TEset = generatePopulation(
        genome,
        piIndice,
        NumberOfIndividual     = parameters["Individuals"],
        NumberOfTransposonTypes = 2,
        NumberOfInsertionsPerType = [1, 0],
        FrequencyOfInsertions  =[
            parameters["FrequencyOfInsertionMain"],
            parameters["FrequencyOfInsertionHGT"],
        ],
        ExcisionRates          = [
            parameters["ExcisionRateMain"], 
            parameters["ExcisionRateHGT"],
        ],
        RepairRates            = [1, 1],
        InsertionRates         = [1, 1],
    )
    
    # run the simulation
    result = runSim(
        genomeMatrix           = genome,
        populationMatrix       = population,
        transposonMatrix       = transposons,
        TEset                  = TEset,
        NumberOfTransposonInsertions=2,
        generations            = parameters["Generations"],
        genMap                 = rates,
        piRNAindices           = piIndice,
        simHGT                 = None,
        HMTgen                 = parameters["HGTgeneration"],
        NumberOfInsertionsPerType = None,
        FrequencyOfInsertions  =[
            parameters["FrequencyOfInsertionMain"],
            parameters["FrequencyOfInsertionHGT"],
        ],
        ExcisionRates          = [
            parameters["ExcisionRateMain"], 
            parameters["ExcisionRateHGT"],
        ],
        RepairRates            = [1, 1],
        InsertionRates         = [1, 1],
        eta                    = parameters["eta"],
        tau                    = parameters["tau"],
        selPen                 = parameters["selectionPenalty"],
        epistasisCoefficient   = parameters["epistasisCoefficient"],
        fitnessFunction        = 2
    )
    
    ff = outputdir + '%030x' % random.randrange(16**30) + "-" + parameters["saveSuffix"] + ".pickle"
    with open(ff, "wb") as f:
        pickle.dump((result), f)
        
    return ff


