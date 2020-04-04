import numpy as np

def checkCopyNumber(populationMatrix):
    # The local flatten lambda
    flatten = lambda *n: (e for a in n for e in (flatten(*a) if isinstance(a, (tuple, list)) else (a,)))
    # Convert the list to flat list and remove 0
    v1 = list(filter((0).__ne__,list(flatten(populationMatrix[:,0].tolist()))))
    v2 = list(filter((0).__ne__,list(flatten(populationMatrix[:,1].tolist()))))
    # Return the average number of transposon copies per individual
    return(len(v1+v2)/len(populationMatrix))