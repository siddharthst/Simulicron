import numpy as np


def stats(
    populationMatrix,
    transposonMatrix,
    TEset,
    insertionSiteFrequencyArray,
):
    # Soft reset the array - there should be a better way to do this
    # procedure
    insertionSiteFrequencyArray = np.zeros(insertionSiteFrequencyArray.shape)

    # Count the transposons at each position
    # With respect to each individual
    for i in range(populationMatrix.shape[0]):
        allele1 = populationMatrix[i][0]
        allele2 = populationMatrix[i][1]
        if allele1 == 0 and allele2 == 0:
            continue

        # Create a new dictionary to store the transposon membership
        TEcontent = {
            k: [] for k in range(1, len(TEset.keys()) + 1)
        }

        if allele1 == 0:
            allele1 = []
        else:
            for l in TEset.keys():
                TEcontent[l] = [z for z in allele1 if z in TEset[l]]

        if allele2 == 0:
            allele2 = []
        else:
            for l in TEset.keys():
                TEcontent[l] = TEcontent[l] + [z for z in allele2 if z in TEset[l]]

        for m in TEcontent.keys():
            for n in TEcontent[m]:
                insertionSiteFrequencyArray[int(transposonMatrix[n,1])][m-1] += 1

    return insertionSiteFrequencyArray


#bool(set(allele1) & TEset)
