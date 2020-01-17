import numpy as np
from numpy import cumsum
from numpy import concatenate as c
import random


def generateGenome(
    numberOfInsertionSites=1000,
    numberOfChromosomes=10,
    baseRecombinationRate=0.01,
):
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
    population[0:, 2] = np.exp(
        np.random.uniform(0.6, 1.0, NumberOfIndividual)
    )
    return population


def generateTransposon(genomeArray, NumberOfTransposonInsertions=2):
    transposons = np.zeros(
        (NumberOfTransposonInsertions + 1, 4), dtype=np.ndarray
    )
    insertionSites = np.random.choice(
        np.arange(genomeArray.shape[0]),
        replace=False,
        size=NumberOfTransposonInsertions,
    )
    counter = 1
    for i in insertionSites:
        transposons[counter][3] = np.random.uniform(0.01, 0.02)
        transposons[counter][2] = genomeArray[i][0]
        transposons[counter][1] = i
        transposons[counter][0] = "%030x" % random.randrange(16 ** 30)
        counter += 1
    return transposons


def generateFitness(
    populationMatrix, transposonMatrix, fitnessFunction=1
):
    allele1Index = np.nonzero(populationMatrix[0:, 0])[0]
    allele2Index = np.nonzero(populationMatrix[0:, 1])[0]
    for i in np.concatenate((allele1Index, allele2Index), axis=0):
        transposonContent = np.trim_zeros(
            np.hstack(
                [populationMatrix[i][0]] + [populationMatrix[i][1]]
            )
        )
        SelectionCoef = transposonMatrix[transposonContent, 2]
        if fitnessFunction == 1:
            populationMatrix[i][2] = populationMatrix[i][2] - np.exp(
                -sum(SelectionCoef)
            )


def recombination2(rates, transposonMatrix, v1=None, v2=None):
    rec = np.random.uniform(size=len(rates)) < rates
    start = [0 if (np.random.uniform() < 0.5) else 1]
    whichcol = 1 + cumsum(c((start, rec))) % 2
    allele1 = np.where(whichcol == 1)[0]
    allele2 = np.where(whichcol == 2)[0]
    r1 = v1[(np.isin(transposonMatrix[[v1], 1], allele1))[0]]
    r2 = v2[(np.isin(transposonMatrix[[v2], 1], allele2))[0]]
    # print (allele1)
    # print (allele2)
    # print (np.isin(transposonMatrix[[v1], 1], allele1))
    # print (np.isin(transposonMatrix[[v2], 1], allele2))
    # print (type(r1))
    # print (type(r2))
    # print (r1)
    # print (r2)
    if r1.size == 0 and r2.size == 0:
        return 0
    elif r1.size == 0 and r2.size != 0:
        return r2
    elif r1.size != 0 and r2.size == 0:
        return r1
    else:
        return c([r1, r2])


def recombination(rates, transposonMatrix, v1=None, v2=None):
    rec = np.random.uniform(size=len(rates)) < rates
    start = [0 if (np.random.uniform() < 0.5) else 1]
    whichcol = 1 + cumsum(c((start, rec))) % 2
    allele1 = set(np.where(whichcol == 1)[0])
    allele2 = set(np.where(whichcol == 2)[0])
    transposon1 = transposonMatrix[[v1], 1].tolist()[0]
    transposon2 = transposonMatrix[[v2], 1].tolist()[0]
    haplotype1 = [x in allele1 for x in transposon1]
    haplotype2 = [x in allele2 for x in transposon2]
    r1 = v1[haplotype1]
    r2 = v2[haplotype2]
    if r1.size == 0 and r2.size == 0:
        return 0
    elif r1.size == 0 and r2.size != 0:
        return r2
    elif r1.size != 0 and r2.size == 0:
        return r1
    else:
        return c([r1, r2])


def transposition(transposonMatrix, genomeMatrix, v1=None, v2=None):
    if v1 == 0:
        allele1 = []
        allele1Sites = []
        allele1Index = []
    else:
        allele1 = transposonMatrix[[v1], 3].tolist()[0]
        allele1Sites = transposonMatrix[[v1], 1].tolist()[0]
        allele1Index = v1
    if v2 == 0:
        allele2 = []
        allele2Sites = []
        allele2Index = []
    else:
        allele2 = transposonMatrix[[v2], 3].tolist()[0]
        allele2Sites = transposonMatrix[[v2], 1].tolist()[0]
        allele2Index = v2

    transposonIndices =  np.array(allele1Index + allele2Index)
    transposonRates = np.array(allele1 + allele2)
    filledSites = allele1Sites + allele2Sites
    Transoposecheck = transposonRates > np.random.uniform(
        0, 0.1, len(transposonRates)
    )

    if not any(Transoposecheck):
        return (v1, v2)
    else:
        transposonsToTranspose = transposonIndices[Transoposecheck]
        emptySiteIndices = [
            x
            for x in list(range(len(genomeMatrix)))
            if x not in filledSites
        ]
        emptySitesProb = genomeMatrix[emptySiteIndices, 1].tolist()
        probSum = sum(emptySitesProb)
        InsertionProb = [float(i) / probSum for i in emptySitesProb]
        sites = np.random.choice(
            emptySiteIndices,
            size=sum(Transoposecheck),
            replace=False,
            p=InsertionProb,
        )
        progenyAllele = random.choices(
            ["v1", "v2"], k=sum(Transoposecheck)
        )
        for i in list(range(len(transposonsToTranspose))):
            transposonToAdd = [
                "%030x" % random.randrange(16 ** 30),
                sites[i],
                genomeMatrix[sites[i]][0],
                transposonMatrix[transposonsToTranspose[i], 3],
            ]
            transposonMatrix = np.vstack(
                [transposonMatrix, transposonToAdd,]
            )
            if progenyAllele[i] == "v1":
                allele1Index.append(len(transposonMatrix)-1)
            if progenyAllele[i] == "v2":
                allele2Index.append(len(transposonMatrix)-1)
    if (allele1Index == []):
        allele1Index = 0
    if (allele2Index == []):
        allele2Index = 0            
    return (allele1Index, allele2Index)
