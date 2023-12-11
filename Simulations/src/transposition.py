# ~ Copyright or © or Copr. CNRS, contributor(s) : Siddharth Sing Tomar (2020-2023)

# ~ Contact: siddharth.egce@gmail.com , arnaud.le-rouzic@universite-paris-saclay.fr

# ~ This software is a computer program whose purpose is to [describe
# ~ functionalities and technical features of your software].

# ~ This software is governed by the CeCILL license under French law and
# ~ abiding by the rules of distribution of free software.  You can  use, 
# ~ modify and/ or redistribute the software under the terms of the CeCILL
# ~ license as circulated by CEA, CNRS and INRIA at the following URL
# ~ "http://www.cecill.info". 

# ~ As a counterpart to the access to the source code and  rights to copy,
# ~ modify and redistribute granted by the license, users are provided only
# ~ with a limited warranty  and the software's author,  the holder of the
# ~ economic rights,  and the successive licensors  have only  limited
# ~ liability. 

# ~ In this respect, the user's attention is drawn to the risks associated
# ~ with loading,  using,  modifying and/or developing or reproducing the
# ~ software by the user in light of its specific status of free software,
# ~ that may mean  that it is complicated to manipulate,  and  that  also
# ~ therefore means  that it is reserved for developers  and  experienced
# ~ professionals having in-depth computer knowledge. Users are therefore
# ~ encouraged to load and test the software's suitability as regards their
# ~ requirements in conditions enabling the security of their systems and/or 
# ~ data to be ensured and,  more generally, to use and operate it in the 
# ~ same conditions as regards security. 

# ~ The fact that you are presently reading this means that you have had
# ~ knowledge of the CeCILL license and that you accept its terms.



import numpy as np
from numpy import cumsum
from numpy import concatenate as c
from regulation import regulation
import random


def transposition(
    transposonMatrix,
    genomeMatrix,
    NumberOfTransposonInsertions,
    TEset,
    piCoord,
    numberOfTranspositionEvents,
    v1,
    v2,
    eta,
):
    # If any of the haplotyle is empty, convert the appr.
    # parameters to empty lists
    if v1 == 0:
        allele1Index = []
        allele1RepairRate = []
        allele1InsertionRate = []
        allele1Sites = []
    else:
        allele1RepairRate = transposonMatrix[[v1], 4].tolist()[0]
        allele1InsertionRate = transposonMatrix[[v1], 5].tolist()[0]
        allele1Sites = transposonMatrix[[v1], 1].tolist()[0]
        allele1Index = v1

    if v2 == 0:
        allele2Index = []
        allele2RepairRate = []
        allele2InsertionRate = []
        allele2Sites = []

    else:
        allele2RepairRate = transposonMatrix[[v2], 4].tolist()[0]
        allele2InsertionRate = transposonMatrix[[v2], 5].tolist()[0]
        allele2Sites = transposonMatrix[[v2], 1].tolist()[0]
        allele2Index = v2

    transposonIndices = np.array(allele1Index + allele2Index)
    # Find effective ExcisionRates
    transposonExcisionRates, RegulationStrengthSet = regulation(
        transposons=transposonIndices,
        TEset=TEset,
        transposonMatrix=transposonMatrix,
        genomeMatrix=genomeMatrix,
        piRNAindices=piCoord,
        eta=eta
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

    Transoposecheck = exicsionCheck & repairCheck & insertionCheck
    # Return alleles as they are if no transposition happens
    if not any(Transoposecheck):
        return (
            v1,
            v2,
            transposonMatrix,
            TEset,
            numberOfTranspositionEvents,
            RegulationStrengthSet,
        )
    else:
        # Create a vector to choose sites from
        genomeSites = list(range(len(genomeMatrix)))
        transposonsToTranspose = transposonIndices[Transoposecheck]
        # If transposons share the same site - replace the old transposon
        # with new transposon
        # Choose the allele for tranposition
        progenyAllele = random.choices(["v1", "v2"], k=sum(Transoposecheck))
        # probSum = sum(emptySitesProb)
        # InsertionProb = [float(i) / probSum for i in emptySitesProb]
        sites = random.sample(genomeSites, sum(Transoposecheck))
        for i in list(range(len(transposonsToTranspose))):
            numberOfTranspositionEvents += 1
            transposonMatrix[numberOfTranspositionEvents, 0] = transposonMatrix[
                transposonsToTranspose[i], 0
            ]
            transposonMatrix[numberOfTranspositionEvents, 1] = sites[i]
            transposonMatrix[numberOfTranspositionEvents, 2] = genomeMatrix[sites[i]][0]
            transposonMatrix[numberOfTranspositionEvents, 3] = transposonMatrix[
                transposonsToTranspose[i], 3
            ]
            transposonMatrix[numberOfTranspositionEvents, 4] = transposonMatrix[
                transposonsToTranspose[i], 4
            ]
            transposonMatrix[numberOfTranspositionEvents, 5] = transposonMatrix[
                transposonsToTranspose[i], 5
            ]

            # transposonMatrix = np.vstack(
            #    [transposonMatrix, np.asarray(transposonToAdd, object),]
            # )
            # transposonMatrix = np.append(
            #    transposonMatrix, np.array([transposonToAdd]), axis=0
            # )
            # Code to track the genealogy of TE
            for k in range(NumberOfTransposonInsertions):
                if transposonsToTranspose[i] in TEset[k + 1]:
                    TEset[k + 1].add(numberOfTranspositionEvents)
                    pass

            # Assign TE to the choosen allele
            if progenyAllele[i] == "v1":
                # First check if the transposon is replacing another transposon on same allele
                # If yes, we need to remove the replaced transposon from the said allele
                if sites[i] in allele1Sites:
                    vIndex = allele1Sites.index(sites[i])
                    del allele1Index[vIndex]
                allele1Index.append(numberOfTranspositionEvents)
            if progenyAllele[i] == "v2":
                # Same as above
                if sites[i] in allele2Sites:
                    vIndex = allele2Sites.index(sites[i])
                    del allele2Index[vIndex]
                allele2Index.append(numberOfTranspositionEvents)

    if allele1Index == []:
        allele1Index = 0
    if allele2Index == []:
        allele2Index = 0

    return (
        allele1Index,
        allele2Index,
        transposonMatrix,
        TEset,
        numberOfTranspositionEvents,
        RegulationStrengthSet,
    )
