import numpy as np
from numpy import cumsum
from numpy import concatenate as c
import random
from random import shuffle
from math import exp


def generateGenome(
    baseSelection=None,
    baseInsertionProb=1,
    numberOfInsertionSites=1000,
    numberOfChromosomes=10,
    baseRecombinationRate=0.01,
):
    # Define selection penalty for the insertion site
    if baseSelection == None:
        SelectionCoef = np.zeros(numberOfInsertionSites)
    elif baseSelection == "Random":
        SelectionCoef = np.random.normal(-0.02, 0.01, numberOfInsertionSites)
    else:
        SelectionCoef = np.full(numberOfInsertionSites, baseSelection)

    # Define insertion probability for the insertion site
    if baseInsertionProb == None:
        insertionProbability = np.zeros(numberOfInsertionSites)
    elif baseInsertionProb == "Random":
        insertionProbability = np.random.uniform(0.01, 0.99, numberOfInsertionSites)
    else:
        insertionProbability = np.full(numberOfInsertionSites, baseInsertionProb)

    RecombinationRates = np.full(
        numberOfInsertionSites, baseRecombinationRate, dtype=float
    )
    chromosomeLocation = np.random.choice(
        np.arange(RecombinationRates.size), replace=False, size=numberOfChromosomes,
    )
    # Only insert more chromosomes if needed
    if numberOfChromosomes > 1:
        RecombinationRates[chromosomeLocation] = 0.5
    genome = np.vstack((SelectionCoef, insertionProbability, RecombinationRates)).T
    return genome


def generatePopulation(
    transposonMatrix,
    NumberOfIndividual=1000,
    NumberOfTransposonInsertions=2,
    InsertIntoOne=False,
    InsertIntoAll=False,
    insertionFrequency=False,
    HardyWeinberg=False,
):
    population = np.zeros((NumberOfIndividual, 3), dtype=np.ndarray)
    # Set base fitness
    population[0:, 2] = 1

    if HardyWeinberg != False and insertionFrequency != False:
        # Calculate the number of insertions
        NumhomozygousInsertion = int(NumberOfIndividual * (insertionFrequency ** 2))
        NumheterozygousInsertion = int(
            NumberOfIndividual * (2 * insertionFrequency * (1 - insertionFrequency))
        )
        indices = list(range(1, NumberOfIndividual))
        shuffle(indices)
        homozygousInsertionSites = indices[0:NumhomozygousInsertion]
        heterozygousInsertionSites = indices[
            NumhomozygousInsertion : NumhomozygousInsertion + NumheterozygousInsertion
        ]
        # Insert transposons and change fitness
        counter = 1
        for i in homozygousInsertionSites:
            population[i][0] = [counter]
            population[i][1] = [counter]
            population[i][2] = exp(2 * transposonMatrix[1][2])

        for i in heterozygousInsertionSites:
            allele = random.choice([0, 1])
            population[i][allele] = [counter]
            population[i][2] = exp(transposonMatrix[1][2])

    elif insertionFrequency != False:
        # Calculate the number of insertions
        numberOfInsertions = int(NumberOfIndividual * insertionFrequency)
        infectedIndividuals = np.random.choice(
            range(1, numberOfInsertions), NumberOfTransposonInsertions, replace=False,
        )
        # Insert transposons and change fitness
        counter = 1
        for i in infectedIndividuals:
            allele = random.choice([0, 1])
            population[i][allele] = [counter]
            population[i][2] = exp(transposonMatrix[1][2])

    elif InsertIntoAll != False:
        # Choose the maternal or paternal chromosome
        for i in range(NumberOfIndividual):
            allele = random.choice([0, 1])
            transposon = random.choice(range(1, NumberOfTransposonInsertions + 1))
            population[i][allele] = [transposon]
            population[i][2] = exp(transposonMatrix[transposon][2])

    elif InsertIntoOne == True:
        infectedIndividual = random.choice(range(1, NumberOfIndividual))
        # Choose the maternal or paternal chromosome
        allele = random.choice([0, 1])
        population[infectedIndividual][allele] = list(
            range(1, NumberOfTransposonInsertions + 1)
        )
        population[infectedIndividual][2] = 1 + sum(
            [transposonMatrix[i][2] for i in range(1, len(transposonMatrix))]
        )

    else:
        infectedIndividuals = np.random.choice(
            range(1, NumberOfIndividual), NumberOfTransposonInsertions, replace=False,
        )

        # Insert transposons and change fitness
        counter = 1
        for i in list(range(NumberOfTransposonInsertions)):
            allele = random.choice([0, 1])
            population[infectedIndividuals[i]][allele] = [counter]
            population[infectedIndividuals[i]][2] = 1 + (transposonMatrix[i + 1][2])
            counter += 1

    return population


def generateTransposon(
    genomeArray,
    baseExcision=0,
    baseRepair=1,
    baseInsertion=1,
    NumberOfTransposonInsertions=2,
    consecutiveTransposons=False,
    changeRecombination=False,
    RecombinationRate=0.01,
    insertionFrequency=False,
    NumberOfIndividual=None,
):
    # if insertionFrequency != False:
    #    NumberOfTransposonInsertions = int(
    #        NumberOfIndividual * insertionFrequency
    #    )

    transposons = np.zeros((NumberOfTransposonInsertions + 1, 6), dtype=np.ndarray)
    if consecutiveTransposons == True:
        # The starting position is padded to prevent
        # index overflow
        start = random.choice(
            range(
                1 + NumberOfTransposonInsertions,
                genomeArray.shape[0] - NumberOfTransposonInsertions,
            )
        )
        insertionSites = np.arange(start, start + NumberOfTransposonInsertions)
    else:
        insertionSites = np.random.choice(
            np.arange(genomeArray.shape[0]),
            replace=False,
            size=NumberOfTransposonInsertions,
        )

    # If there is a requirment to change the recombination rate
    # at transposon insertion site
    # This changes the recombination on index position
    if changeRecombination == True:
        genomeArray[insertionSites[-2], 2] = RecombinationRate

    counter = 1
    for i in insertionSites:
        transposons[counter][5] = baseInsertion
        transposons[counter][4] = baseRepair
        transposons[counter][3] = (
            0 if baseExcision == 0 else baseExcision
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
