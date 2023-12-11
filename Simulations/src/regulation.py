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


def regulation(
    transposons,
    TEset,
    transposonMatrix,
    genomeMatrix,
    piRNAindices,
    fast=False,
    eta=0.0,
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
            TEexcEffective[indices] = TEexcision[indices] - (
                TEexcision[indices] * TERegulationSet[keys]
                )
                
        return (TEexcEffective.tolist(), TERegulationSet)
