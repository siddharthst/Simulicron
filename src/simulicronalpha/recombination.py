import numpy as np
from numpy import cumsum
from numpy import concatenate as c
import random

def recombination(rates, transposonMatrix, v1, v2):
    # Empty vectors to store result
    r1 = []
    r2 = []
    # Creating lambda (macro)
    # "Match" does not exist in python
    match = lambda a, b: [ b.index(x) if x in b else 0 for x in a ]
    # Get the postion of transposons 
    positionV1 = transposonMatrix[v1,1].astype(int).tolist()
    positionV2 = transposonMatrix[v2,1].astype(int).tolist()
    # Create a 2D array storing information of transposons 
    # and their location
    # This will keep track of transposons and location during sorting
    posMatrix = np.transpose([v1 + v2, positionV1 + positionV2])
    posMatrix = posMatrix[posMatrix[:,1].argsort()]
    # This step sorts the locations and adds another location, 
    # 0 if not already present
    unqiquePos = list(set(positionV1+positionV2) | set([0]))
    unqiquePos.sort()
    # Calculate the effective rate from genome map
    effectiveRates = 0.5*(1-np.exp(-2*np.diff(rates[unqiquePos])))
    print (effectiveRates)
    # Performing "Recombination"
    rec = (np.random.uniform(size=len(effectiveRates)) < effectiveRates)
    # Select the direction to start from
    start = [(np.random.uniform() < 0.5)]
    # Concat. start and recombination
    # Also remove the added 0 and start from
    # whichhaplo and uniqpos
    whichhaplo = 1 + cumsum(c((start, rec))) % 2
    whichhaplo = np.delete(whichhaplo, 0)
    del unqiquePos[0]
    print (whichhaplo)
    # Generating the haplotype
    # Also checking if there is no transposon left
    # In the case above, return a (int) 0
    # Else return the array containing transposons
    if (positionV1 == [0]):
        r1 = [0]
    else:
        checkRec = (whichhaplo[(match(positionV1, unqiquePos))] == 1)
        #checkRec = 
        print (checkRec)
        #r1 = transposonMatrix[whichhaplo[match(positionV1, unqiquePos)] == 1]
    if (positionV2 == [0]):
        r2 = [0]
    else:
        pass
        #r2 = transposonMatrix[whichhaplo[match(positionV2, unqiquePos)] == 2]
    #print (r1)
    #print (r2)





    '''
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

'''