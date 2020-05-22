import numpy as np
from numpy import cumsum
from numpy import concatenate as c
import pandas as pd
import random

from generateSim import (
    generateGenome,
    generatePopulation,
    generateTransposon,
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
    generations=200,
    ignoreFixation=True,
):
    # ------------------#
    # lambda/macros
    flatten = lambda *n: (
        e
        for a in n
        for e in (
            flatten(*a) if isinstance(a, (tuple, list)) else (a,)
        )
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
    # Create a insertionSiteSepecific array
    insertionSiteFrequencyArray = np.zeros(
        (len(genomeMatrix), NumberOfTransposonInsertions,)
    )

    transposonMatrixCopy = transposonMatrix
    populationMatrixCopy = populationMatrix
    for i in range(generations):
        populationV1 = []
        populationV2 = []
        populationFit = []
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
            if (
                populationMatrixCopy[p1, 0] == 0
                and populationMatrixCopy[p1, 1] == 0
            ):
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
                    genomeMatrix[0:, 2],
                    transposonMatrixCopy,
                    v1=cP1V1,
                    v2=cP1V2,
                )

            if (
                populationMatrixCopy[p2, 0] == 0
                and populationMatrixCopy[p2, 1] == 0
            ):
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
                    genomeMatrix[0:, 2],
                    transposonMatrixCopy,
                    v1=cP2V1,
                    v2=cP2V2,
                )

            if v1 == 0 and v2 == 0:
                indFitness = 1

            else:
                v1, v2, transposonMatrixCopy, TEset = transposition(
                    transposonMatrix=transposonMatrixCopy,
                    genomeMatrix=genomeMatrix,
                    NumberOfTransposonInsertions=NumberOfTransposonInsertions,
                    TEset=TEset,
                    v1=v1,
                    v2=v2,
                )
                indFitness = calculateFitness(
                    transposonMatrixCopy, v1, v2
                )

            populationV1.append(v1)
            populationV2.append(v2)
            populationFit.append(indFitness)
        # Return i+2 since i start at 0 = generation 1
        # and the condition check happens at generation
        # n-1, hence i + 1 + 1

        # Check if there are no transposons left
        if all(
            np.array_equal(v, [0, 0])
            for v in np.c_[populationV1, populationV2]
        ):
            return (
                "LOSS",
                fixedTE,
                unfixedTE,
                lostTE,
                i + 2,
                transposonMatrixCopy.size / 4 - 1,
                averageCopyNumber,
                varianceCopyNumber,
            )
        else:
            continue

        # Major bug in numpy - forced to use pandas
        # Refer to the question
        # https://stackoverflow.com/questions/60210897
        populationMatrixCopy = pd.DataFrame(
            [populationV1, populationV2, populationFit]
        ).T.to_numpy()
        (
            compyNumber,
            varianceNumber,
            insertionSiteFrequencyArray,
        ) = statistics(
            populationMatrixCopy,
            transposonMatrixCopy,
            TEset,
            insertionSiteFrequencyArray,
        )
        copyNumber, varianceNumber = checkCopyNumber(
            populationMatrixCopy
        )
        averageCopyNumber.append(copyNumber)
        varianceCopyNumber.append(varianceNumber)
    # Quit simulation if there in a transient state
    # i.e. no loss
    return (
        "FLUX",
        fixedTE,
        unfixedTE,
        lostTE,
        i + 2,
        transposonMatrixCopy.size / 4 - 1,
        averageCopyNumber,
        varianceCopyNumber,
    )


# Creating a generator for the dataset
def createData(
    numberOfSimulations=1000,
    numberOfChromosomes=4,
    numberOfInsertionSites=1000,
    baseRecombinationRate=0.1,
    NumberOfIndividual=1000,
    InsertIntoOne=False,
    InsertIntoAll=False,
    NumberOfTransposonInsertions=2,
    NumberOfGenerations=100000,
    baseSelection=1,
    baseExcision=1,
    baseRepair=1,
    baseInsertion=1,
    consecutiveTransposons=False,
    changeRecombination=False,
    baseTrRecombination=0.1,
    insertionFrequency=False,
    HardyWeinberg=False,
):
    for i in range(numberOfSimulations):
        gen = generateGenome(
            numberOfInsertionSites=numberOfInsertionSites,
            numberOfChromosomes=numberOfChromosomes,
            baseRecombinationRate=baseRecombinationRate,
            baseSelection=baseSelection,
        )
        tr, gen, TEset = generateTransposon(
            genomeArray=gen,
            NumberOfTransposonInsertions=NumberOfTransposonInsertions,
            baseExcision=baseExcision,
            baseRepair=baseRepair,
            baseInsertion=baseInsertion,
            consecutiveTransposons=consecutiveTransposons,
            changeRecombination=changeRecombination,
            RecombinationRate=baseTrRecombination,
            insertionFrequency=insertionFrequency,
            NumberOfIndividual=NumberOfIndividual,
        )
        pop = generatePopulation(
            transposonMatrix=tr,
            NumberOfIndividual=NumberOfIndividual,
            NumberOfTransposonInsertions=NumberOfTransposonInsertions,
            InsertIntoOne=InsertIntoOne,
            InsertIntoAll=InsertIntoAll,
            insertionFrequency=insertionFrequency,
            HardyWeinberg=HardyWeinberg,
        )

        yield (
            (
                gen,
                pop,
                tr,
                TEset,
                NumberOfTransposonInsertions,
                NumberOfGenerations,
            )
        )


def runBatch(
    numberOfSimulations=1000,
    numberOfChromosomes=4,
    numberOfInsertionSites=1000,
    baseRecombinationRate=0.1,
    NumberOfIndividual=1000,
    InsertIntoOne=False,
    InsertIntoAll=False,
    NumberOfTransposonInsertions=2,
    NumberOfGenerations=100000,
    baseSelection=1,
    baseExcision=1,
    baseRepair=1,
    baseInsertion=1,
    consecutiveTransposons=False,
    changeRecombination=False,
    baseTrRecombination=0.1,
    insertionFrequency=False,
    HardyWeinberg=False,
    numberOfThreads=1,
):
    dataSet = createData(
        numberOfSimulations,
        numberOfChromosomes,
        numberOfInsertionSites,
        baseRecombinationRate,
        NumberOfIndividual,
        InsertIntoOne,
        InsertIntoAll,
        NumberOfTransposonInsertions,
        NumberOfGenerations,
        baseSelection,
        baseExcision,
        baseRepair,
        baseInsertion,
        consecutiveTransposons,
        changeRecombination,
        baseTrRecombination,
        insertionFrequency,
        HardyWeinberg,
    )
    inputSet = []
    for i in dataSet:
        inputSet.append(i)

    with multiprocessing.Pool(processes=numberOfThreads) as pool:
        results = pool.starmap(runSim, inputSet)

    return results
