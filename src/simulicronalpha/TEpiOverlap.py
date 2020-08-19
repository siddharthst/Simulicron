import numpy as np

def TEpiOverlap(populationMatrix, transposonMatrix, TEset, piIndices):
    # Convert the list to flat list and remove 0
    # if present.
    # Also create a dictionary counting 
    # each family type
    TEfamilyPiCount = {k: [] for k in range(1, len(TEset.keys()) + 1)}
    
    for i in range(len(populationMatrix)):
        v1 = populationMatrix[i, 0] if type(populationMatrix[i, 0]) is list else []
        v2 = populationMatrix[i, 1] if type(populationMatrix[i, 1]) is list else []
        if (v1 + v2 == []):
            for k in TEset.keys():
                TEfamilyPiCount[k].append(0)
        else:
            # Check the occupancy in piRNA for each family
            transposons = v1 + v2
            TElocations = transposonMatrix[transposons, 1]
            TEfamily = transposonMatrix[transposons, 0]         
            for k in TEset.keys():
                indices = (TEfamily == k).nonzero()[0]
                familyLocations = TElocations[indices]
                piOverlap = sum(el in piIndices for el in familyLocations)
                TEfamilyPiCount[k].append(piOverlap)
    for k in TEset.keys():
        TEfamilyPiCount[k] = sum(TEfamilyPiCount[k])/len(populationMatrix)
    return(TEfamilyPiCount)
