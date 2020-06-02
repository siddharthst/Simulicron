import numpy as np
from numpy import cumsum
from numpy import concatenate as c
from regulation import regulation
import random


def transposition(
    transposonMatrix, genomeMatrix, NumberOfTransposonInsertions, TEset, v1, v2,
):
    if v1 == 0:
        allele1 = []
        allele1Sites = []
        allele1Index = []
        allele1ExicsionRate = []
        allele1RepairRate = []
        allele1InsertionRate = []
    else:
        allele1ExicsionRate = transposonMatrix[[v1], 3].tolist()[0]
        allele1RepairRate = transposonMatrix[[v1], 4].tolist()[0]
        allele1InsertionRate = transposonMatrix[[v1], 5].tolist()[0]
        allele1Sites = transposonMatrix[[v1], 1].astype(int).tolist()[0]
        allele1Index = v1

    if v2 == 0:
        allele2 = []
        allele2Sites = []
        allele2Index = []
        allele2ExicsionRate = []
        allele2RepairRate = []
        allele2InsertionRate = []

    else:
        allele2ExicsionRate = transposonMatrix[[v2], 3].tolist()[0]
        allele2RepairRate = transposonMatrix[[v2], 4].tolist()[0]
        allele2InsertionRate = transposonMatrix[[v2], 5].tolist()[0]
        allele2Sites = transposonMatrix[[v2], 1].astype(int).tolist()[0]
        allele2Index = v2

    transposonIndices = np.array(allele1Index + allele2Index)
    # Find effective ExcisionRates
    transposonExcisionRates = regulation(
        transposons=transposonIndices,
        TEset=TEset,
        transposonMatrix=transposonMatrix,
        genomeMatrix=genomeMatrix,
    )
    transposonRepairRates = np.array(allele1RepairRate + allele2RepairRate)
    transposonInsertionRates = np.array(allele1InsertionRate + allele2InsertionRate)

    exicsionCheck = transposonExcisionRates > np.random.uniform(
        0, 1, len(transposonExcisionRates)
    )
    repairCheck = transposonRepairRates > np.random.uniform(
        0, 1, len(transposonRepairRates)
    )
    insertionCheck = transposonInsertionRates > np.random.uniform(
        0, 1, len(transposonInsertionRates)
    )
    # Find the intersection
    # Transoposecheck = [
    #     x for x in exicsionCheck if x in repairCheck and x in insertionCheck
    # ]
    Transoposecheck = exicsionCheck & repairCheck & insertionCheck
    # Return alleles as they are if no transposition happens
    if not any(Transoposecheck):
        return (v1, v2, transposonMatrix, TEset)
    else:
        transposonsToTranspose = transposonIndices[Transoposecheck]
        # If transposons share the same site - replace the old transposon
        # with new transposon
        # Choose the allele for tranposition
        progenyAllele = random.choices(["v1", "v2"], k=sum(Transoposecheck))
        # probSum = sum(emptySitesProb)
        # InsertionProb = [float(i) / probSum for i in emptySitesProb]
        sites = np.random.choice(list(range(len(genomeMatrix))), size=sum(Transoposecheck), replace=False)
        for i in list(range(len(transposonsToTranspose))):
            transposonToAdd = [
                0,
                sites[i],
                genomeMatrix[sites[i]][0],
                transposonMatrix[transposonsToTranspose[i], 3],
                transposonMatrix[transposonsToTranspose[i], 4],
                transposonMatrix[transposonsToTranspose[i], 5],
            ]

            transposonMatrix = np.vstack(
                [transposonMatrix, np.asarray(transposonToAdd, object),]
            )
            # Code to track the genealogy of TE
            for k in range(NumberOfTransposonInsertions):
                if transposonsToTranspose[i] in TEset[k + 1]:
                    TEset[k + 1].add(len(transposonMatrix) - 1)
                    pass

            # Assign TE to the choosen allele
            if progenyAllele[i] == "v1":
                # First check if the transposon is replacing another transposon on same allele
                # If yes, we need to remove the replaced transposon from the said allele
                if sites[i] in allele1Sites:
                    vIndex = allele1Sites.index(sites[i])
                    del allele1Index[vIndex]
                allele1Index.append(len(transposonMatrix) - 1)
            if progenyAllele[i] == "v2":
                # Same as above
                if sites[i] in allele2Sites:
                    vIndex = allele2Sites.index(sites[i])
                    del allele2Index[vIndex]
                allele2Index.append(len(transposonMatrix) - 1)

    if allele1Index == []:
        allele1Index = 0
    if allele2Index == []:
        allele2Index = 0

    return (allele1Index, allele2Index, transposonMatrix, TEset)
