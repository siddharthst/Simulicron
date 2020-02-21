import numpy as np
from numpy import cumsum
from numpy import concatenate as c
import random


def generateGenome(
    baseSelection=None,
    baseInsertionProb=None,
    numberOfInsertionSites=1000,
    numberOfChromosomes=10,
    baseRecombinationRate=0.01,
):
    # Define selection penalty for the insertion site
    if baseSelection == None:
        SelectionCoef = np.zeros(numberOfInsertionSites)
    elif baseSelection == "Random":
        SelectionCoef = np.random.normal(
            -0.02, 0.01, numberOfInsertionSites
        )
    else:
        SelectionCoef = np.full(numberOfInsertionSites, baseSelection)

    # Define insertion probability for the insertion site
    if baseInsertionProb == None:
        insertionProbability = np.zeros(numberOfInsertionSites)
    elif baseInsertionProb == "Random":
        insertionProbability = np.random.uniform(
            0.01, 0.99, numberOfInsertionSites
        )
    else:
        insertionProbability = np.full(
            numberOfInsertionSites, baseInsertionProb
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
    transposonMatrix,
    NumberOfIndividual=1000,
    NumberOfTransposonInsertions=2,
    InsertIntoOne=True,
):
    population = np.zeros((NumberOfIndividual, 3), dtype=np.ndarray)
    # Set base fitness
    population[0:, 2] = 1

    if InsertIntoOne == True:
        infectedIndividual = random.choice(
            range(1, NumberOfIndividual)
        )
        # Choose the maternal or paternal chromosome
        allele = random.choice([0, 1])
        population[infectedIndividual][allele] = list(
            range(1, NumberOfTransposonInsertions + 1)
        )
        population[infectedIndividual][2] = 1 + sum(
            [
                transposonMatrix[i][2]
                for i in range(1, len(transposonMatrix))
            ]
        )

    else:
        infectedIndividuals = np.random.choice(
            range(1, NumberOfIndividual),
            NumberOfTransposonInsertions,
            replace=False,
        )

        # Insert transposons and change fitness
        counter = 1
        for i in list(range(NumberOfTransposonInsertions)):
            allele = random.choice([0, 1])
            population[infectedIndividuals[i]][allele] = [counter]
            population[infectedIndividuals[i]][2] = 1 + (
                transposonMatrix[i + 1][2]
            )
            counter += 1

    return population


def generateTransposon(
    genomeArray,
    baseTransposition=1,
    NumberOfTransposonInsertions=2,
    consecutiveTransposons=False,
    changeRecombination=False,
    RecombinationRate=0.01,
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
    # This changes the recombination on index position + 1
    if changeRecombination == True:
        genomeArray[insertionSites[-2], 2] = RecombinationRate

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
    for i in list(range(NumberOfTransposonInsertions)):
        TEset[i + 1] = set([i + 1])
    return transposons, genomeArray, TEset
