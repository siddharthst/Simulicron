import numpy as np
from numpy import cumsum
from numpy import concatenate as c
import pandas as pd
import random
import time

from generateSim import generateGenome, generatePopulation, initHGT
from recombination import recombination
from fitness import calculateFitness
from transposition import transposition
from checkCopyNumber import checkCopyNumber
from TEpiOverlap import TEpiOverlap

# Current multiprocessing implementation
import multiprocessing

# For future cluster based implementation
# import ray
# ray.init()

# Aux return function


def returnHelper(TEset, TEfamilyCountArr, TEfamilyVarArr, TEregulationArr):
    # Make the family information "flat"
    # https://stackoverflow.com/questions/5946236/how-to-merge-multiple-dicts-with-same-key
    dict1 = {}
    dict2 = {}
    dict3 = {}
    for k in TEset.keys():
        dict1[k] = tuple(dict1[k] for dict1 in TEfamilyCountArr)
        dict2[k] = tuple(dict2[k] for dict2 in TEfamilyVarArr)
        dict3[k] = tuple(dict3[k] for dict3 in TEregulationArr)
    return (dict1, dict2, dict3)


def coreReturn(
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
    populationArray,
    transposonMatrix,
):
    return {
        "State": simulationState,
        "Generatrion": i + 2,
        "NTE": numberOfTranspositionEvents,
        "AvgCopyNum": averageCopyNumber,
        "CopyNumVar": varianceCopyNumber,
        "TEfamilyCN": TEfamilyCountArrRes,
        "TEfamilyVR": TEfamilyVarArrRes,
        "TEfamilyRg": TEregulationArrRes,
        "AvgFit": avgFitness,
        "HGTGen": HMTgen,
        "ETA": eta,
        "NTI": NumberOfTransposonInsertions,
        "Freq": FrequencyOfInsertions,
        "TRate": ExcisionRates,
        "Tau": tau,
        "selPen": selPen,
        "piRNA": piRNAindices,
        "TEpi": overlap,
        "FitnessFunction": fitnessFunction,
        "epistasisCoefficient": epistasisCoefficient,
        "TEset": TEset,
        "populationArray": populationArray,
        "transposonMatrix": transposonMatrix,
    }


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
    piRNAindices,
    simHGT=None,
    HMTgen=None,
    NumberOfInsertionsPerType=None,
    FrequencyOfInsertions=None,
    ExcisionRates=None,
    RepairRates=None,
    InsertionRates=None,
    eta=0.0,
    tau=None,
    selPen=None,
    maxAvgTE=None,
    epistasisCoefficient=0.0,
    fitnessFunction=1,
):
    # ------------------#
    # ------------------#
    # ------------------#
    # For storing population state
    populationArray = []
    # for storing transposons which are fixed
    fixedTE = []
    # for storing transposons which are not fixed
    unfixedTE = []
    # for storing transposons which are lost
    lostTE = []
    # For storing the average fitness per generation
    avgFitness = []
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
        (
            len(genomeMatrix),
            NumberOfTransposonInsertions,
        )
    )
    transposonMatrixCopy = transposonMatrix
    populationMatrixCopy = populationMatrix

    # Create a fixed size array for storing transposon information
    # Also add a "counter" - with a large fixed array, it is not
    # possible to use len(transposons) to find the last filled index
    transposonMatrixCopy = np.append(
        transposonMatrix,
        np.zeros((10000000, 6), dtype=object),
        axis=0,
    )
    numberOfTranspositionEvents = len(transposonMatrix)

    # Calculate the CN and CNV for generation 0
    (
        copyNumber,
        varianceNumber,
        TEfamilyCount,
        TEfamilyVar,
    ) = checkCopyNumber(populationMatrixCopy, TEset, transposonMatrixCopy)
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
        if HMTgen == i:
            (
                populationMatrixCopy,
                transposonMatrixCopy,
                TEset,
                numberOfTranspositionEvents,
            ) = initHGT(
                populationMatrixCopy,
                transposonMatrixCopy,
                genomeMatrix,
                TEset,
                piRNAindices,
                numberOfTranspositionEvents,
                FrequencyOfInsertions[-1],
                ExcisionRates[-1],
                RepairRates[-1],
                InsertionRates[-1],
            )
        for k in list(range(populationMatrixCopy.shape[0])):
            fitness = list(populationMatrixCopy[0:, 2])
            p1, p2 = random.choices(
                list(range(populationMatrixCopy.shape[0])),
                weights=fitness,
                k=2,
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

                v1 = recombination(
                    genMap,
                    transposonMatrixCopy,
                    v1=cP1V1,
                    v2=cP1V2,
                )

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

                v2 = recombination(
                    genMap,
                    transposonMatrixCopy,
                    v1=cP2V1,
                    v2=cP2V2,
                )

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
                    RegulationStrength,
                ) = transposition(
                    transposonMatrix=transposonMatrixCopy,
                    genomeMatrix=genomeMatrix,
                    NumberOfTransposonInsertions=NumberOfTransposonInsertions,
                    TEset=TEset,
                    piCoord=piRNAindices,
                    numberOfTranspositionEvents=numberOfTranspositionEvents,
                    v1=v1,
                    v2=v2,
                    eta=eta,
                )
                indFitness = calculateFitness(
                    transposonMatrixCopy,
                    v1,
                    v2,
                    fitnessFunction=fitnessFunction,
                    epistasisCoefficient=epistasisCoefficient,
                )
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
            TEfamilyCountArrRes, TEfamilyVarArrRes, TEregulationArrRes = returnHelper(
                TEset, TEfamilyCountArr, TEfamilyVarArr, TEregulationArr
            )
            simulationState = "LOSS"
            overlap = "NA"
            return coreReturn(
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
                populationArray,
                transposonMatrix,
            )
        else:
            pass
        # Generate population matrix for next iteration
        populationMatrixCopy = np.array(
            [
                populationV1,
                populationV2,
                populationFit,
            ],
            dtype="object",
        ).T
        # Append population state into matrix
        populationArray.append(populationMatrixCopy)
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
        (
            copyNumber,
            varianceNumber,
            TEfamilyCount,
            TEfamilyVar,
        ) = checkCopyNumber(populationMatrixCopy, TEset, transposonMatrixCopy)
        averageCopyNumber.append(copyNumber)
        varianceCopyNumber.append(varianceNumber)
        TEfamilyCountArr.append(TEfamilyCount)
        TEfamilyVarArr.append(TEfamilyVar)
        avgFitness.append(sum((populationMatrixCopy[:, 2]) / len(populationMatrixCopy)))

        # Exit the simulation of TE copy number exceeds the threshold
        if maxAvgTE != None:
            if maxAvgTE > copyNumber:
                (
                    TEfamilyCountArrRes,
                    TEfamilyVarArrRes,
                    TEregulationArrRes,
                ) = returnHelper(
                    TEset, TEfamilyCountArr, TEfamilyVarArr, TEregulationArr
                )
                overlap = TEpiOverlap(
                    populationMatrixCopy,
                    transposonMatrixCopy,
                    TEset,
                    piRNAindices,
                )
                simulationState = "ATMAX"
                return coreReturn(
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
                    populationArray,
                    transposonMatrix,
                )
    # Quit simulation if there in a transient state
    # i.e. no loss
    TEfamilyCountArrRes, TEfamilyVarArrRes, TEregulationArrRes = returnHelper(
        TEset, TEfamilyCountArr, TEfamilyVarArr, TEregulationArr
    )
    # Return the TE locations
    overlap = TEpiOverlap(
        populationMatrixCopy,
        transposonMatrixCopy,
        TEset,
        piRNAindices,
    )
    simulationState = "FLUX"
    return coreReturn(
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
        populationArray,
        transposonMatrix,
    )
