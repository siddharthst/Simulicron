import numpy as np


def checkCopyNumber(populationMatrix, TEset, transposonMatrix):
    # Convert the list to flat list and remove 0
    # if present.
    # Also create a dictionary counting 
    # each family type
    TEfamilyCount = {k: [] for k in range(1, len(TEset.keys()) + 1)}
    TEfamilyVar = {k: [] for k in range(1, len(TEset.keys()) + 1)}
    TETotal = []
    for i in range(len(populationMatrix)):
        v1 = populationMatrix[i, 0] if type(populationMatrix[i, 0]) is list else []
        v2 = populationMatrix[i, 1] if type(populationMatrix[i, 1]) is list else []
        if (v1 + v2 == []):
            TETotal.append(0)
            for k in TEset.keys():
                TEfamilyCount[k].append(0)
        else:
            # Check the occupancy for each family
            TEfamilies = transposonMatrix[v1+v2, 0].tolist()
            TETotal.append(len(TEfamilies))
            for k in TEset.keys():
                TEfamilyCount[k].append(TEfamilies.count(k))
    # Calculate variance per family and sum the count
    for k in TEset.keys():
        TEfamilyVar[k].append(np.var(TEfamilyCount[k]))
        TEfamilyCount[k] = sum(TEfamilyCount[k])/len(populationMatrix)
    return sum(TETotal)/len(populationMatrix), np.var(TETotal), TEfamilyCount, TEfamilyVar