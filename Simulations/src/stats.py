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
