import numpy as np
from numpy import cumsum
from numpy import concatenate as c
import random


def transposition(
    transposonMatrix,
    genomeMatrix,
    NumberOfTransposonInsertions,
    TEset,
    v1,
    v2,
):
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

        # Choose the allele for tranposition
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
            # Code to track the genealogy of TE
            for k in range(NumberOfTransposonInsertions):
                if (transposonsToTranspose[i] in TEset[k+1]):
                    TEset[k+1].add(len(transposonMatrix) - 1)
                    pass
            
            # Assign TE to the choosed allele
            if progenyAllele[i] == "v1":
                allele1Index.append(len(transposonMatrix) - 1)
            if progenyAllele[i] == "v2":
                allele2Index.append(len(transposonMatrix) - 1)
    if allele1Index == []:
        allele1Index = 0
    if allele2Index == []:
        allele2Index = 0

    return (allele1Index, allele2Index, transposonMatrix, TEset)
