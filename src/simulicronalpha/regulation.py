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
    TElocations = (
        transposonMatrix[transposons, 1].astype(int).tolist()
    )
    # Find their actual exision rates
    TEexcision = transposonMatrix[transposons, 3]

    # Implementing a "faster" version of regulation when only
    # classical piRNA effect is being studied (Kofler, 2019)
    if fast == True:
        # Check if any transposon is present in the piRNA cluster
        # if any(e in TElocations for e in piRNAindices):
        #   # Check which family has piRNA insertion
        #   inPiRNA = np.in1d(TElocations,piRNAindices)
        #   transposonsInPiRNA = transposons[inPiRNA]
        #   # Above ^
        #   # Check if at least one member from each family is present
        #   # in piRNA
        #   count = 0
        #   for i in TEset.keys():
        #     if bool(set(transposonsInPiRNA) & TEset[i]):
        #       count += 1
        #   if count == len(TEset.keys()):
        #     # All excision rates in this condition are
        #     # effectivly zero
        #     return ((TEexcision * 0).tolist())
        #   else:
        #     # Find which TE family is in piRNA and change the
        #     # effective excision rates to 0
        #     for i in TEset.keys():
        #       if bool(set(transposonsInPiRNA) & TEset[i]):
        #         counter = 0
        #         for k in transposons:
        #           if k in TEset[i]:
        #             TEexcision[i] = 0
        #           counter += 1
        #     return (TEexcision.tolist())
        # else:
        #   # Return the normal excision rates
        #   return(transposonMatrix[transposons,3].tolist())
        pass

    else:
        # print ("-------")
        # Identify the family
        # for i in transposons:
        #  TEfamily.append([a for a, b in TEset.items() if i in b][0])
        TEfamily = transposonMatrix[transposons, 0]
        # TEfamily = np.array(TEfamily)
        # Fill the tau array
        GenomeTau = genomeMatrix[TElocations, 3]
        # Identify the families present and create a set
        # Also create a list to store effect regulation factor
        tauList = []
        TEfamilySet = set(TEfamily)

        # Calculate effective excision rates per family
        for i in list(TEfamilySet):
            indices = (TEfamily == i).nonzero()[0]
            netTau = sum(GenomeTau[indices])
            if netTau > 1:
                netTau = 1
            tauList.append(netTau)
        # Implement co-regulation if only one out of n
        # families is present in piRNA/KRAB-ZfP
        # And if only one family is present, no need
        # to go through co-regulation
        if (len(TEfamilySet)) == 1:
            # print ('True1')
            TEexcEffective = TEexcision - (TEexcision * tauList[0])
        elif tauList == [1.0] * len(tauList):
            # print ('True2')
            # If there are multiple families - but they
            # are completely regulated by piRNA
            TEexcEffective = TEexcision - (TEexcision * tauList[0])
        elif tauList == [0.0] * len(tauList):
            # print ('True3')
            # If there are multiple families - but they
            # are completely NOT regulated by piRNA
            TEexcEffective = TEexcision
        else:
            # print ('True4')
            # Only one family is regulated, either partially or
            # completely - coregulation is only applicable in
            # this condition
            # The effective transposition rate in this case would be
            # u = 1- (tau(n) + (sum(tau(every other family) * eta)))
            for i in range(len(tauList)):
                netTau = tauList[i] - (
                    (sum([x for k, x in enumerate(tauList) if k != i]))
                    * eta
                )
                indices = (TEfamily == list(TEfamilySet)[i]).nonzero()[0]
                # TEexcEffective[indices] = netTau
                TEexcEffective[indices] = TEexcision[indices] - (TEexcision[indices] * netTau)
        # print (TEfamily)
        # print (TEfamilySet)
        # print (tauList)
        # print (TEexcision)
        # print (TEexcEffective)
        # print ("-------")
        return TEexcEffective.tolist()

