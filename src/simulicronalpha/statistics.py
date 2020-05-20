import numpy as np


def statistics(
    populationMatrix,
    transposonMatrix,
    TEset,
    insertionSiteFrequencyArray,
):
    # Create a empty matrix where rows =  individuals
    # and columns = transposons
    for i in len(populationMatrix):
        allele1 = populationMatrix[i][0]
        allele2 = populationMatrix[i][1]
        if allele1 == 0 and allele2 == 0:
            continue
        else:
            # Check the progenitor of TE
            if allele1 == 0:
                pass
            else:
                bool(set(allele1) & TEset)
                
