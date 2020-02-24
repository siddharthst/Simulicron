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
    generations=100000,
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
            )

        # Check if all members of population contain transposon
        if not any(
            any(k == 0 for k in z)
            for z in [
                list(flatten(v.flatten().tolist()))
                for v in np.c_[populationV1, populationV2]
            ]
        ):
            counter = 0
            for TE in list(TEset):
                if all(
                    bool(set(k).intersection(TEset[TE]))
                    for k in populationV1 + populationV2
                ):
                    # Check if the TE appears in
                    # unfixed array and remove it
                    # from that array
                    unfixedTE = list(filter((TE).__ne__, unfixedTE))
                    fixedTE.append(TE)
                    counter += 1
                elif any(
                    bool(set(k).intersection(TEset[TE]))
                    for k in populationV1 + populationV2
                ):
                    unfixedTE.append(TE)
                else:
                    # Check if the TE appears in
                    # unfixed array and remove it
                    # from that array
                    unfixedTE = list(filter((TE).__ne__, unfixedTE))
                    lostTE.append(TE)
                    del TEset[TE]

                # If all transposons are fixed
                if counter == len(TEset):
                    # Cleanup the output
                    fixedTE = list(set(fixedTE))
                    return (
                        "FIXED",
                        fixedTE,
                        unfixedTE,
                        lostTE,
                        i + 2,
                        transposonMatrixCopy.size / 4 - 1,
                    )
        # Major bug in numpy - forced to use pandas
        # Refer to the question
        # https://stackoverflow.com/questions/60210897
        populationMatrixCopy = pd.DataFrame(
            [populationV1, populationV2, populationFit]
        ).T.to_numpy()

    # Quit simulation if there in a transient state
    # i.e. no fixation or loss
    return (
        "FLUX",
        fixedTE,
        unfixedTE,
        lostTE,
        i + 2,
        transposonMatrixCopy.size / 4 - 1,
    )


# Creating a generator for the dataset
def createData(
    numberOfSimulations=1000,
    numberOfChromosomes=4,
    numberOfInsertionSites=1000,
    baseRecombinationRate=0.1,
    NumberOfIndividual=1000,
    InsertIntoOne=True,
    NumberOfTransposonInsertions=2,
    NumberOfGenerations=100000,
    baseSelection=1,
    baseTransposition=1,
    consecutiveTransposons=False,
    changeRecombination=False,
    baseTrRecombination=0.1,
    insertionFrequency=False,
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
            baseTransposition=baseTransposition,
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
            insertionFrequency=insertionFrequency,
        )

        yield ((gen, pop, tr, TEset, NumberOfTransposonInsertions,))


def runBatch(
    numberOfSimulations=1000,
    numberOfChromosomes=4,
    numberOfInsertionSites=1000,
    baseRecombinationRate=0.1,
    NumberOfIndividual=1000,
    InsertIntoOne=True,
    NumberOfTransposonInsertions=2,
    NumberOfGenerations=100000,
    baseSelection=1,
    baseTransposition=1,
    consecutiveTransposons=False,
    changeRecombination=False,
    baseTrRecombination=0.1,
    insertionFrequency=False,
):
    dataSet = createData(
        numberOfSimulations,
        numberOfChromosomes,
        numberOfInsertionSites,
        baseRecombinationRate,
        NumberOfIndividual,
        InsertIntoOne,
        NumberOfTransposonInsertions,
        NumberOfGenerations,
        baseSelection,
        baseTransposition,
        consecutiveTransposons,
        changeRecombination,
        baseTrRecombination,
        insertionFrequency,
    )
    inputSet = []
    for i in dataSet:
        inputSet.append(i)

    with multiprocessing.Pool(processes=8) as pool:
        results = pool.starmap(runSim, inputSet)

    return results

