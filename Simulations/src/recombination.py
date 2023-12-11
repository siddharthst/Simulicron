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
import random


def recombination(rates, transposonMatrix, v1, v2):
    # Empty vectors to store result
    r1 = []
    r2 = []
    # Creating lambda (macro)
    # "Match" does not exist in python
    # match = lambda a, b: [ b.index(x) if x in b else 0 for x in a ]
    # Get the postion of transposons
    positionV1 = transposonMatrix[v1, 1].astype(int).tolist()
    positionV2 = transposonMatrix[v2, 1].astype(int).tolist()
    # This step sorts the locations and adds another location,
    # 0 if not already present
    unqiquePos = list(set(positionV1 + positionV2) | set([0]))
    unqiquePos.sort()
    # Calculate the effective rate from genome map
    effectiveRates = 0.5 * (1 - np.exp(-2 * np.diff(rates[unqiquePos])))
    # print (effectiveRates)
    # Performing "Recombination"
    rec = np.random.uniform(size=len(effectiveRates)) < effectiveRates
    # Select the direction to start from
    start = [(np.random.uniform() < 0.5)]
    # Concat. start and recombination
    # Also remove the added 0 and start from
    # whichhaplo and uniqpos
    whichhaplo = 1 + cumsum(c((start, rec))) % 2
    whichhaplo = np.delete(whichhaplo, 0)
    del unqiquePos[0]
    unqiquePos = np.asarray(unqiquePos)
    # Generating the haplotype
    # Also checking if there is no transposon left
    # In the case above, return a (int) 0
    # Else return the array containing transposons
    if positionV1 == [0]:
        r1 = []
    else:
        if not any(whichhaplo == 1):
            pass
        else:
            pos = set(unqiquePos[whichhaplo == 1])
            for i in v1:
                if (transposonMatrix[i, 1]) in pos:
                    r1.append(i)
    if positionV2 == [0]:
        r2 = []
    else:
        if not any(whichhaplo == 2):
            pass
        else:
            pos = set(unqiquePos[whichhaplo == 2])
            for i in v2:
                if (transposonMatrix[i, 1]) in pos:
                    r2.append(i)
    # Merge to create gamate
    r = r1 + r2
    # Return 0 if no transposon remains
    if not r:
        return 0
    else:
        return r
