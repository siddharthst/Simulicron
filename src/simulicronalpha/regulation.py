import numpy as np

def regulation(transposons, TEset, transposonMatrix, genomeMatrix):
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
  # Identify the family
  for i in transposons:
    TEfamily.append([a for a, b in TEset.items() if i in b][0]) 
  TEfamily = np.array(TEfamily)
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

    

    