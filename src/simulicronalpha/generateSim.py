import numpy as np
from numpy import cumsum
from numpy import concatenate as c
import random

def generateGenome(
    baseSelection=1,
    numberOfInsertionSites=1000,
    numberOfChromosomes=10,
    baseRecombinationRate=0.01,
):
    if baseSelection == 0.0:
        SelectionCoef = np.zeros(numberOfInsertionSites)
    else:
        SelectionCoef = np.random.normal(
            -0.02, 0.01, numberOfInsertionSites
        )
    insertionProbability = np.random.uniform(
        0.01, 0.99, numberOfInsertionSites
    )
    RecombinationRates = np.full(
        numberOfInsertionSites, baseRecombinationRate, dtype=float
    )
    chromosomeLocation = np.random.choice(
        np.arange(RecombinationRates.size),
        replace=False,
        size=numberOfChromosomes,
    )
    # Only insert more chromosomes if needed
    if numberOfChromosomes > 1:
        RecombinationRates[chromosomeLocation] = 0.5
    genome = np.vstack(
        (SelectionCoef, insertionProbability, RecombinationRates)
    ).T
    return genome


def generatePopulation(
    NumberOfIndividual=1000, NumberOfTransposonInsertions=2
):
    population = np.zeros((NumberOfIndividual, 3), dtype=np.ndarray)
    infectedIndividuals = np.random.choice(
        range(1, NumberOfIndividual),
        NumberOfTransposonInsertions,
        replace=False,
    )
    counter = 1
    for i in list(range(NumberOfTransposonInsertions)):
        allele = random.choice([0, 1])
        population[infectedIndividuals[i]][allele] = [counter]
        counter += 1
    population[0:, 2] = 1
    return population


# def generateTransposon(
#     genomeArray, baseTransposition=1, NumberOfTransposonInsertions=2
# ):
#     transposons = np.zeros(
#         (NumberOfTransposonInsertions + 1, 4), dtype=np.ndarray
#     )
#     insertionSites = np.random.choice(
#         np.arange(genomeArray.shape[0]),
#         replace=False,
#         size=NumberOfTransposonInsertions,
#     )
#     counter = 1
#     for i in insertionSites:
#         transposons[counter][3] = (
#             0
#             if baseTransposition == 0
#             else np.random.uniform(0.02, 0.03)
#         )
#         transposons[counter][2] = genomeArray[i][0]
#         transposons[counter][1] = i
#         transposons[counter][0] = "%030x" % random.randrange(16 ** 30)
#         counter += 1
#     return transposons


def generateTransposon(
    genomeArray,
    baseTransposition=1,
    NumberOfTransposonInsertions=2,
    consecutiveTransposons=False,
    changeRecombination=False,
    baseRecombination=0.01,
):
    transposons = np.zeros(
        (NumberOfTransposonInsertions + 1, 4), dtype=np.ndarray
    )
    if consecutiveTransposons == True:
        # The starting position is padded to prevent
        # index overflow
        start = random.choice(
            range(
                1 + NumberOfTransposonInsertions,
                genomeArray.shape[0] - NumberOfTransposonInsertions,
            )
        )
        insertionSites = np.arange(
            start, start + NumberOfTransposonInsertions
        )
    else:
        insertionSites = np.random.choice(
            np.arange(genomeArray.shape[0]),
            replace=False,
            size=NumberOfTransposonInsertions,
        )

    # If there is a requirment to change the recombination rate
    # at transposon insertion site
    if changeRecombination == True:
        genomeArray[insertionSites+1, 2] = baseRecombination

    counter = 1
    for i in insertionSites:
        transposons[counter][3] = (
            0
            if baseTransposition == 0
            else np.random.uniform(0.02, 0.03)
        )
        transposons[counter][2] = genomeArray[i][0]
        transposons[counter][1] = i
        transposons[counter][0] = "%030x" % random.randrange(16 ** 30)
        counter += 1

    # Create sets to store transposon propogation
    TEset = {}
    for i in NumberOfTransposonInsertions:
        tempSet = set()
        TEset[i] = tempSet.add(i)
    return transposons, genomeArray, TEset
