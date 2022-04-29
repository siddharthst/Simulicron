import numpy as np


def regulation(
    transposons,
    TEset,
    transposonMatrix,
    genomeMatrix,
    piRNAindices,
    fast=False,
    eta=0.0,
    regulationStr=0.0,
):
    # eta = co-regulation coefficient
    # Create an empty array to store transposon locations
    TElocations = []
    # Create an empty array to store their actual excision rates
    TEexcision = []
    # Create an empty array to store the family information
    TEfamily = []
    # Create an empty array to store the tau
    GenomeTau = []
    # Create an empty array to store effective excision rates
    TEexcEffective = np.zeros(len(transposons))
    # Covert transposon IDs to location
    TElocations = transposonMatrix[transposons, 1].astype(int).tolist()
    # Find their actual exision rates
    TEexcision = transposonMatrix[transposons, 3]
    # Create empty set to store regulation values for each family
    TERegulationSet = {}
    for i in range(len(TEset)):
        TERegulationSet[i + 1] = 0.0
    # Implementing a "faster" version of regulation when only
    # classical piRNA effect is being studied (Kofler, 2019)
    if fast == True:
        pass

    else:
        # print ("-------")
        TEfamily = transposonMatrix[transposons, 0]
        # TEfamily = np.array(TEfamily)
        # Fill the tau array
        GenomeTau = genomeMatrix[TElocations, 3]
        # Identify the families present and create a set
        # Also create a list to store effect regulation factor
        tauList = []
        TEfamilySet = set(TEfamily)
        # Calculate effective excision rates per family
        for key in TERegulationSet.keys():
            if key in TEfamilySet:
                indices = (TEfamily == key).nonzero()[0]
                netTau = sum(GenomeTau[indices])
                if netTau > 1:
                    netTau = 1
                TERegulationSet[key] = netTau
                tauList.append(netTau)
            else:
                tauList.append(0)
                TERegulationSet[key] = 0
        # Adjust transposition rates 
        for keys in TERegulationSet.keys():
            TERegulationSet[keys] = TERegulationSet[keys] + (
                sum(
                    [
                        value
                        for key, value in TERegulationSet.items()
                        if key != keys
                    ]
                )
                * eta 
            )
            if TERegulationSet[keys] > 1:
                TERegulationSet[keys] = 1.0
            indices = (TEfamily == keys).nonzero()[0]
            TEexcEffective[indices] = (TEexcision[indices] - (
                TEexcision[indices] * TERegulationSet[keys]
                )) / (1. + len(transposons) * regulationStr)
                
        return (TEexcEffective.tolist(), TERegulationSet)
