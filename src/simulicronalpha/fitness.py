import numpy as np


def fitness(
    self, genomeFrame=None, populationFrame=None, transposonFrame=None, function=1
):
    TEcontent = (
        transposonFrame[
            transposonFrame["TID"].isin(populationFrame["TEfather"].tolist())
        ]["InsertionSite"].tolist()
        + transposonFrame[
            transposonFrame["TID"].isin(populationFrame["TEmother"].tolist())
        ]["InsertionSite"].tolist()
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
