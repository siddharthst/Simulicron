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


def TEpiOverlap(
    populationMatrix, transposonMatrix, TEset, piIndices, SingleFamily=False
):
    if SingleFamily == False:
        # Convert the list to flat list and remove 0
        # if present.
        # Also create a dictionary counting
        # each family type
        TEfamilyPiCount = {k: [] for k in range(1, len(TEset.keys()) + 1)}

        for i in range(len(populationMatrix)):
            v1 = populationMatrix[i, 0] if type(populationMatrix[i, 0]) is list else []
            v2 = populationMatrix[i, 1] if type(populationMatrix[i, 1]) is list else []
            if v1 + v2 == []:
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
            TEfamilyPiCount[k] = sum(TEfamilyPiCount[k]) / len(populationMatrix)
        return TEfamilyPiCount

    if SingleFamily == True:
        # Create an empty array to store the piRNA insertions for
        # each individual
        TECount = []
        for i in range(len(populationMatrix)):
            v1 = populationMatrix[i, 0] if type(populationMatrix[i, 0]) is list else []
            v2 = populationMatrix[i, 1] if type(populationMatrix[i, 1]) is list else []
            if v1 + v2 == []:
                TECount.append(0)
            else:
                transposons = v1 + v2
                TElocations = transposonMatrix[transposons, 1]
                piOverlap = sum(el in piIndices for el in TElocations)
                TECount.append(piOverlap)
        return sum(TECount) / len(populationMatrix)
