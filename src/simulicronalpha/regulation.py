import numpy as np

def regulation(transposons, TEset, transposonMatrix, genomeMatrix, piRNAindices, fast=False):
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
  TElocations = transposonMatrix[transposons,1].astype(int).tolist()
  # Find their actual exision rates
  TEexcision = transposonMatrix[transposons, 3]

  # Implementing a "faster" version of regulation when only
  # classical piRNA effect is being studied (Kofler, 2019)
  if fast == True:
    # Check if any transposon is present in the piRNA cluster
    if any(e in TElocations for e in piRNAindices):
      # Check which family has piRNA insertion
      inPiRNA = np.in1d(TElocations,piRNAindices)
      transposonsInPiRNA = transposons[inPiRNA]
      # Above ^
      # Check if at least one member from each family is present
      # in piRNA
      count = 0
      for i in TEset.keys():
        if bool(set(transposonsInPiRNA) & TEset[i]):
          count += 1 
      if count == len(TEset.keys()):
        # All excision rates in this condition are 
        # effectivly zero 
        return ((TEexcision * 0).tolist())
      else:
        # Find which TE family is in piRNA and change the
        # effective excision rates to 0
        for i in TEset.keys():
          if bool(set(transposonsInPiRNA) & TEset[i]):
            counter = 0
            for k in transposons:
              if k in TEset[i]:
                TEexcision[i] = 0
              counter += 1 
        return (TEexcision.tolist())
    else:
      # Return the normal excision rates
      return(transposonMatrix[transposons,3].tolist())

  else:
    # Identify the family
    #for i in transposons:
    #  TEfamily.append([a for a, b in TEset.items() if i in b][0]) 
    TEfamily = transposonMatrix[transposons, 0]
    #TEfamily = np.array(TEfamily)
    # Fill of the tau array
    GenomeTau = genomeMatrix[TElocations,3]
    # Calculate effective excision rates
    for i in list(set(TEfamily)):
      indices = (TEfamily == i).nonzero()[0]
      netTau = sum(GenomeTau[indices])
      if netTau > 1:
        netTau = 1
      TEexcEffective[indices] = TEexcision[indices] - (TEexcision[indices] * netTau)
    
    return(TEexcEffective.tolist())

    

    