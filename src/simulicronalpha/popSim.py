import numpy as np
from numpy import cumsum
from numpy import concatenate as c
import random

# Current multiprocessing implementation
import multiprocessing

# For future cluster based implementation
# import ray
# ray.init()

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


def calculateFitness(
    transposonMatrix, currentFitness, v1, v2, fitnessFunction=1,
):
    cV1 = v1
    cV2 = v2
    if cV1 == 0:
        cV1 = np.asarray([])
    else:
        cV1 = np.asarray(v1)
    if cV2 == 0:
        cV2 = np.asarray([])
    else:
        cV2 = np.asarray(v2)
    teContent = c([cV1, cV2]).astype(int)
    penalties = transposonMatrix[teContent, 2]
    if fitnessFunction == 1:
        return currentFitness - np.exp(sum(penalties))


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
    # print (haplotype1)
    # print (r1)
    # print (allele1)
    # print (haplotype2)
    # print (r2)
    # print (allele2)
    if r1.size == 1:
        if r1 == [0]:
            r1 = np.asarray([])
    if r2.size == 1:
        if r2 == [0]:
            r2 = np.asarray([])
    if r1.size == 0 and r2.size == 0:
        return 0
    elif r1.size == 0 and r2.size != 0:
        return list(r2)
    elif r1.size != 0 and r2.size == 0:
        return list(r1)
    else:
        return list(c([r1, r2]))


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

    transposonIndices = np.array(allele1Index + allele2Index)
    transposonRates = np.array(allele1 + allele2)
    filledSites = allele1Sites + allele2Sites
    Transoposecheck = transposonRates > np.random.uniform(
        0, 1, len(transposonRates)
    )
    if not any(Transoposecheck):
        return (v1, v2, transposonMatrix)
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
                [
                    transposonMatrix,
                    np.asarray(transposonToAdd, object),
                ]
            )
            if progenyAllele[i] == "v1":
                allele1Index.append(len(transposonMatrix) - 1)
            if progenyAllele[i] == "v2":
                allele2Index.append(len(transposonMatrix) - 1)
    if allele1Index == []:
        allele1Index = 0
    if allele2Index == []:
        allele2Index = 0

    return (allele1Index, allele2Index, transposonMatrix)

# For future cluster based implementation
# @ray.remote()
def runSim(
    genomeMatrix, populationMatrix, transposonMatrix, generations=100,
):
    transposonMatrixCopy = transposonMatrix
    populationMatrixCopy = populationMatrix
    for i in list(range(generations)):
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
            baseFitness = random.choice(
                [
                    populationMatrixCopy[p1, 2],
                    populationMatrixCopy[p2, 2],
                ]
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
                indFitness = baseFitness

            else:
                v1, v2, transposonMatrixCopy = transposition(
                    transposonMatrixCopy, genomeMatrix, v1=v1, v2=v2
                )
                indFitness = calculateFitness(
                    transposonMatrixCopy, baseFitness, v1, v2
                )

            populationV1.append(v1)
            populationV2.append(v2)
            populationFit.append(indFitness)
        populationMatrixCopy = np.vstack(
            (populationV1, populationV2, populationFit)
        ).T
        if all(v == 0 for v in populationV1) and all(
            v == 0 for v in populationV2
        ):
            return (0, i)
        if all(v != 0 for v in populationV1) or all(
            v != 0 for v in populationV2
        ):
            return (1, i)
    # print(transposonMatrixCopy.size / 4)
    # print ('p2v1', populationMatrixCopy[p2, 0])
    # print ('c2v1', cP2V1)
    # print ('c2v1Type', type(populationMatrixCopy[p2, 0]))
    # print ('p2v2', populationMatrixCopy[p2, 1])
    # print ('c2v2', cP2V2)
    # print ('c2v2Type', type(populationMatrixCopy[p2, 1]))

def runBatch(numberOfGenerations=10000):
    results = []
    for i in range(numberOfGenerations):
        gen = generateGenome(numberOfInsertionSites=1000)
        pop = generatePopulation()
        tr = generateTransposon(gen)
        





