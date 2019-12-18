import numpy as np

def fitness(genomeFrame, populationFrame, function=1):
    TEcontent = (
        populationFrame["Insertion_Father"] + populationFrame["Insertion_Mother"]
    )
    selectionCoef = []
    for i in TEcontent:
        selectionCoef.append(
            genomeFrame[genomeFrame["InsertionSiteID"] == i]["SelectionCoef"].values[0]
        )
    if function == 1:
        fitnessValue = np.exp(sum(selectionCoef))
    if function == 2:
        fitnessValue = 1 + sum(selectionCoef)
    if function == 3:
        selectionCoef = np.array(selectionCoef)
        fitnessValue = np.prod(1 + selectionCoef)
    return fitnessValue
