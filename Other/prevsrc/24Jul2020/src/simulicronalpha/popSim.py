import numpy as np
from numpy import cumsum
from numpy import concatenate as c
import pandas as pd
import random
import time

from generateSim import (
    generateGenome,
    generatePopulation,
)
from recombination import recombination
from fitness import calculateFitness
from transposition import transposition
from checkCopyNumber import checkCopyNumber

# Current multiprocessing implementation
import multiprocessing

# For future cluster based implementation
# import ray
# ray.init()

# For future cluster based implementation
# @ray.remote()
def runSim(
    genomeMatrix,
    populationMatrix,
    transposonMatrix,
    TEset,
    NumberOfTransposonInsertions,
    generations,
    genMap,
    piSet,
    simHGT=None,
    HGTgen=None,
):
    # ------------------#
    # lambda/macros
    flatten = lambda *n: (
        e for a in n for e in (flatten(*a) if isinstance(a, (tuple, list)) else (a,))
    )
    # ------------------#
    # ------------------#
    # for storing transposons which are fixed
    fixedTE = []
    # for storing transposons which are not fixed
    unfixedTE = []
    # for storing transposons which are lost
    lostTE = []
    # For storing the average copy number per generation
    averageCopyNumber = []
    # For storing the copy number variance per generation
    varianceCopyNumber = []
    # Per family statistics
    TEfamilyCountArr = []
    TEfamilyVarArr = []
    TEregulationArr = []
    # Create a insertionSiteSepecific array
    insertionSiteFrequencyArray = np.zeros(
        (len(genomeMatrix), NumberOfTransposonInsertions,)
    )
    transposonMatrixCopy = transposonMatrix
    populationMatrixCopy = populationMatrix

    # Create a fixed size array for storing transposon information
    # Also add a "counter" - with a large fixed array, it is not
    # possible to use len(transposons) to find the last filled index
    transposonMatrixCopy = np.append(
        transposonMatrix, np.zeros((10000000, 6), dtype=object), axis=0,
    )
    numberOfTranspositionEvents = len(transposonMatrix)

    # Calculate piRNA coordinates
    piCoord = []
    for i in piSet.values():
        piCoord = piCoord + (list(range(i[0], i[1])))

    # Calculate the CN and CNV for generation 0
    (copyNumber, varianceNumber, TEfamilyCount, TEfamilyVar,) = checkCopyNumber(
        populationMatrixCopy, TEset, transposonMatrixCopy
    )
    averageCopyNumber.append(copyNumber)
    varianceCopyNumber.append(varianceNumber)
    TEfamilyCountArr.append(TEfamilyCount)
    TEfamilyVarArr.append(TEfamilyVar)

    # Driver loop
    for i in range(generations):
        # Start the clock
        # timeStart = time.time()
        # print(i)
        # Create arrays to store information
        populationV1 = []
        populationV2 = []
        populationFit = []
        populationRegulation = {k: [] for k in range(1, len(TEset.keys()) + 1)}
        for k in list(range(populationMatrixCopy.shape[0])):
            fitness = list(populationMatrixCopy[0:, 2])

            p1, p2 = random.choices(
                list(range(populationMatrixCopy.shape[0])), weights=fitness, k=2,
            )

            # Since recombination function only accepts arrays,
            # checking and forcing type conversion as needed
            # for the respective alleles
            if populationMatrixCopy[p1, 0] == 0 and populationMatrixCopy[p1, 1] == 0:
                v1 = 0
            else:
                cP1V1 = populationMatrixCopy[p1, 0]
                cP1V2 = populationMatrixCopy[p1, 1]
                if isinstance(cP1V1, list):
                    cP1V1 = np.asarray(cP1V1)
                else:
                    cP1V1 = np.asarray([cP1V1])
                if isinstance(cP1V2, list):
                    cP1V2 = np.asarray(cP1V2)
                else:
                    cP1V2 = np.asarray([cP1V2])

                v1 = recombination(genMap, transposonMatrixCopy, v1=cP1V1, v2=cP1V2,)

            if populationMatrixCopy[p2, 0] == 0 and populationMatrixCopy[p2, 1] == 0:
                v2 = 0
            else:
                cP2V1 = populationMatrixCopy[p2, 0]
                cP2V2 = populationMatrixCopy[p2, 1]
                if isinstance(cP2V1, list):
                    cP2V1 = np.asarray(cP2V1)
                else:
                    cP2V1 = np.asarray([cP2V1])
                if isinstance(cP2V2, list):
                    cP2V2 = np.asarray(cP2V2)
                else:
                    cP2V2 = np.asarray([cP2V2])

                v2 = recombination(genMap, transposonMatrixCopy, v1=cP2V1, v2=cP2V2,)

            if v1 == 0 and v2 == 0:
                indFitness = 1
                for key in TEset.keys():
                    populationRegulation[key].append(0)

            else:
                (
                    v1,
                    v2,
                    transposonMatrixCopy,
                    TEset,
                    numberOfTranspositionEvents,
                    RegulationStrength
                ) = transposition(
                    transposonMatrix=transposonMatrixCopy,
                    genomeMatrix=genomeMatrix,
                    NumberOfTransposonInsertions=NumberOfTransposonInsertions,
                    TEset=TEset,
                    piCoord=piCoord,
                    numberOfTranspositionEvents=numberOfTranspositionEvents,
                    v1=v1,
                    v2=v2,
                )
                indFitness = calculateFitness(transposonMatrixCopy, v1, v2)
                for key in TEset.keys():
                    populationRegulation[key].append(RegulationStrength[key])

            populationV1.append(v1)
            populationV2.append(v2)
            populationFit.append(indFitness)
        # Return i+2 since i start at 0 = generation 1
        # and the condition check happens at generation
        # n-1, hence i + 1 + 1

        # Check if there are no transposons left
        if all(np.array_equal(v, [0, 0]) for v in np.c_[populationV1, populationV2]):
            # Make the family information "flat"
            # https://stackoverflow.com/questions/5946236/how-to-merge-multiple-dicts-with-same-key
            dict1 = {}
            dict2 = {}
            dict3 = {}
            for k in TEset.keys():
                dict1[k] = tuple(dict1[k] for dict1 in TEfamilyCountArr)
                dict2[k] = tuple(dict2[k] for dict2 in TEfamilyVarArr)
                dict3[k] = tuple(dict3[k] for dict3 in TEregulationArr)
            TEfamilyCountArrRes = dict1
            TEfamilyVarArrRes = dict2
            TEregulationArrRes = dict3
            return {
                "State": "LOSS",
                "Generatrion": i + 2,
                "NTE": numberOfTranspositionEvents,
                "AvgCopyNum": averageCopyNumber,
                "CopyNumVar": varianceCopyNumber,
                "TEfamilyCN": TEfamilyCountArrRes,
                "TEfamilyVR": TEfamilyVarArrRes,
                "TEfamilyRg": TEregulationArrRes,
            }
        else:
            pass
        # Major bug in numpy - forced to use pandas
        # Refer to the question
        # https://stackoverflow.com/questions/60210897
        # populationMatrixCopy = pd.DataFrame(
        #    [populationV1, populationV2, populationFit]
        # ).T.to_numpy()
        populationMatrixCopy = np.array(
            [populationV1, populationV2, populationFit,], dtype="object"
        ).T
        # Regulation strength for each family
        for key in TEset.keys():
            populationRegulation[key] = sum(populationRegulation[key]) / len(
                populationMatrixCopy
            )
        TEregulationArr.append(populationRegulation)
        # (
        #     compyNumber,
        #     varianceNumber,
        #     insertionSiteFrequencyArray,
        # ) = statistics(
        #     populationMatrixCopy,
        #     transposonMatrixCopy,
        #     TEset,
        #     insertionSiteFrequencyArray,
        # )
        (copyNumber, varianceNumber, TEfamilyCount, TEfamilyVar,) = checkCopyNumber(
            populationMatrixCopy, TEset, transposonMatrixCopy
        )
        averageCopyNumber.append(copyNumber)
        varianceCopyNumber.append(varianceNumber)
        TEfamilyCountArr.append(TEfamilyCount)
        TEfamilyVarArr.append(TEfamilyVar)
        # Stop the timer
        # timeStop = time.time()
        # Terminate the loop if it runs for more than
        # 8 seconds!
        # if (timeStop-timeStart) > 10.0:
        #    return {
        #        "State": "TIMEXC",
        #        "Generatrion": i + 2,
        #        "NTE": numberOfTranspositionEvents,
        #        "AvgCopyNum": averageCopyNumber,
        #        "CopyNumVar": varianceCopyNumber,
        #        "TEfamilyCN": TEfamilyCountArr,
        #        "TEfamilyVR": TEfamilyVarArr,
        #    }

    # Quit simulation if there in a transient state
    # i.e. no loss
    dict1 = {}
    dict2 = {}
    dict3 = {}
    for k in TEset.keys():
        dict1[k] = tuple(dict1[k] for dict1 in TEfamilyCountArr)
        dict2[k] = tuple(dict2[k] for dict2 in TEfamilyVarArr)
        dict3[k] = tuple(dict3[k] for dict3 in TEregulationArr)
    TEfamilyCountArrRes = dict1
    TEfamilyVarArrRes = dict2
    TEregulationArrRes = dict3
    return {
        "State": "FLUX",
        "Generatrion": i + 2,
        "NTE": numberOfTranspositionEvents,
        "AvgCopyNum": averageCopyNumber,
        "CopyNumVar": varianceCopyNumber,
        "TEfamilyCN": TEfamilyCountArrRes,
        "TEfamilyVR": TEfamilyVarArrRes,
        "TEfamilyRg": TEregulationArrRes,
    }


# Creating a generator for the dataset
def createData(
    numberOfSimulations=1000,
    baseSelection=1,
    baseInsertionProb=1,
    numberOfInsertionSites=1000,
    numberOfChromosomes=6,
    baseRecombinationRate=0.01,
    baseTau=1,
    numberOfPiRNA=6,
    piPercentage=3,
    enablePiRecombination=False,
    NumberOfIndividual=1000,
    NumberOfTransposonTypes=2,
    NumberOfInsertionsPerType=[2, 2],
    FrequencyOfInsertions=[0.5, 0.5],
    ExcisionRates=[0, 0],
    RepairRates=[1, 1],
    InsertionRates=[1, 1],
    HardyWeinberg=False,
    NumberOfGenerations=10000,
):
    for i in range(numberOfSimulations):
        gen, piset, piIndice, rate2Map = generateGenome(
            baseSelection=baseSelection,
            baseInsertionProb=1,
            numberOfInsertionSites=numberOfInsertionSites,
            numberOfChromosomes=numberOfChromosomes,
            baseRecombinationRate=baseRecombinationRate,
            baseTau=baseTau,
            numberOfPiRNA=numberOfPiRNA,
            piPercentage=piPercentage,
            enablePiRecombination=False,
        )
        pop, tr, TEset = generatePopulation(
            genomeMatrix=gen,
            piRNAindices=piIndice,
            NumberOfIndividual=NumberOfIndividual,
            NumberOfTransposonTypes=NumberOfTransposonTypes,
            NumberOfInsertionsPerType=NumberOfInsertionsPerType,
            FrequencyOfInsertions=FrequencyOfInsertions,
            ExcisionRates=ExcisionRates,
            RepairRates=RepairRates,
            InsertionRates=InsertionRates,
            HardyWeinberg=False,
            numberOfPiRNA=numberOfPiRNA,
        )

        yield (
            (
                gen,
                pop,
                tr,
                TEset,
                NumberOfTransposonTypes,
                NumberOfGenerations,
                rate2Map,
                piset,
            )
        )


def runBatch(
    numberOfSimulations=1,
    baseSelection=1,
    baseInsertionProb=1,
    numberOfInsertionSites=1000,
    numberOfChromosomes=6,
    baseRecombinationRate=0.01,
    baseTau=1,
    numberOfPiRNA=6,
    piPercentage=3,
    enablePiRecombination=False,
    NumberOfIndividual=1000,
    NumberOfTransposonTypes=2,
    NumberOfInsertionsPerType=[2, 2],
    FrequencyOfInsertions=[0.5, 0.5],
    ExcisionRates=[0, 0],
    RepairRates=[1, 1],
    InsertionRates=[1, 1],
    HardyWeinberg=False,
    NumberOfGenerations=100,
    numberOfThreads=1,
):
    dataSet = createData(
        numberOfSimulations=numberOfSimulations,
        baseSelection=baseSelection,
        baseInsertionProb=baseInsertionProb,
        numberOfInsertionSites=numberOfInsertionSites,
        numberOfChromosomes=numberOfChromosomes,
        baseRecombinationRate=baseRecombinationRate,
        baseTau=baseTau,
        numberOfPiRNA=numberOfPiRNA,
        piPercentage=piPercentage,
        enablePiRecombination=enablePiRecombination,
        NumberOfIndividual=NumberOfIndividual,
        NumberOfTransposonTypes=NumberOfTransposonTypes,
        NumberOfInsertionsPerType=NumberOfInsertionsPerType,
        FrequencyOfInsertions=FrequencyOfInsertions,
        ExcisionRates=ExcisionRates,
        RepairRates=RepairRates,
        InsertionRates=InsertionRates,
        HardyWeinberg=False,
        NumberOfGenerations=NumberOfGenerations,
    )
    inputSet = []
    for i in dataSet:
        inputSet.append(i)

    with multiprocessing.Pool(processes=numberOfThreads) as pool:
        results = pool.starmap(runSim, inputSet)

    return results