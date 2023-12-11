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

def calculateFitness(
    transposonMatrix, v1, v2, fitnessFunction=1, epistasisCoefficient=1
):
    cV1 = v1
    cV2 = v2
    if cV1 == 0:
        cV1 = np.asarray([])
    else:
        cV1 = np.asarray(v1)
    if cV2 == 0:
        cV2 = np.asarray([])
    else:
        cV2 = np.asarray(v2)
    teContent = c([cV1, cV2]).astype(int)
    penalties = transposonMatrix[teContent, 2]
    if fitnessFunction == 1:
        return np.exp(sum(penalties))
    if fitnessFunction == 2:
        # This follows the assumption that all sites share same
        # selection pressure
        w = np.exp(sum(penalties) + (0.5 * epistasisCoefficient * (penalties[0])**2 * len(penalties)**2))
        return (w)
