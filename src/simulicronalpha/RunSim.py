import ClassDeco as cd
import pandas as pd
import random
import numpy as np
import uuid
from coolname import generate_slug
from numpy import random, random_intel
import logging
import coloredlogs
from InitSim import initSim

# Logger
logger = logging.getLogger(__name__)
coloredlogs.install(level="DEBUG")

###################################################################################################
# /////|   Class   |///////////////////////////////////////////////////////////////////////////////
###################################################################################################
class RunSim(initSim):
    def resetFitness(self, Function=1):
        if Function == 1:
            self.PopFrame["NetFitness"] = np.exp(self.PopFrame["NetFitness"])
        if Function == 2 or Function == 3:
            self.PopFrame["NetFitness"] = 1 + (self.PopFrame["NetFitness"])

    # /////|   Recombination   |///////////////////////////////////////////////////////////////////////

    def recombination(self, genomeOrg=None):
        # Create insertion list for the progeny
        TEid = []
        # Create a copy of genomeFrame
        genomeCopy = self.GenFrame.copy(deep=True)
        # Added lookup to transposon frame
        # Instead of looking into population database the insertion
        # sites are now accessed directly from the transposon data-
        # base
        TEid_Father = genomeOrg["TEfather"]
        TEid_Mother = genomeOrg["TEmother"]

        if TEid_Father[0] == 0 and TEid_Mother[0] == 0:
            return [0]
        else:
            initParent = random.choice(["M", "F"])
            surrogateParent = "M" if ("M" != initParent) else "F"
            RandArray = np.random_intel.uniform(0, 1.0, self.GenFrame.shape[0])
            switch = genomeCopy["RecombinationRate"] > RandArray
            # counter = collections.Counter(switch)
            # print(counter)
            genomeCopy["Progeny"] = np.where(
                switch.cumsum() % 2 == 0, initParent, surrogateParent
            )
            if all(v == 0 for v in TEid_Father):
                pass
            else:
                for i in TEid_Father:
                    insertionSite = self.TranspFrame.loc[self.TranspFrame["TID"] == i][
                        "InsertionSite"
                    ].values[0]
                    se = genomeCopy[genomeCopy["InsertionSite"] == insertionSite][
                        "Progeny"
                    ].values[0]
                    if se == "M":
                        TEid.append(i)

            if all(v == 0 for v in TEid_Mother):
                pass
            else:
                for i in TEid_Mother:
                    insertionSite = self.TranspFrame.loc[self.TranspFrame["TID"] == i][
                        "InsertionSite"
                    ].values[0]
                    se = genomeCopy[genomeCopy["InsertionSite"] == insertionSite][
                        "Progeny"
                    ].values[0]
                    if se == "F":
                        TEid.append(i)

            # print(genomeCopy['Progeny'].value_counts())
        if not TEid:
            TEid.append(0)
        return TEid

    # /////|   Fitness   |/////////////////////////////////////////////////////////////////////////////

    def fitness(
        self, TIDm, TIDf, populationFrame=None, transposonFrame=None, function=1
    ):
        TEcontent = (
            self.TranspFrame[self.TranspFrame["TID"].isin(TIDm)][
                "InsertionSite"
            ].tolist()
            + self.TranspFrame[self.TranspFrame["TID"].isin(TIDf)][
                "InsertionSite"
            ].tolist()
        )
        selectionCoef = []
        for i in TEcontent:
            selectionCoef.append(
                self.GenFrame[self.GenFrame["InsertionSiteID"] == i][
                    "SelectionCoef"
                ].values[0]
            )
        if function == 1:
            fitnessValue = np.exp(sum(selectionCoef))
        if function == 2:
            fitnessValue = 1 + sum(selectionCoef)
        if function == 3:
            selectionCoef = np.array(selectionCoef)
            fitnessValue = np.prod(1 + selectionCoef)
        return fitnessValue

    # /////|   Transposition   |///////////////////////////////////////////////////////////////////////

    def transposition(self, TIDm, TIDf):
        # Skip the check for empty transposon content. Do it in the parent
        TEfather = TIDm
        TEmother = TIDf
        filledSites = list(
            filter(
                (0).__ne__,
                self.TranspFrame[self.TranspFrame["TID"].isin([TIDm])][
                    "InsertionSite"
                ].tolist()
                + self.TranspFrame[self.TranspFrame["TID"].isin([TIDf])][
                    "InsertionSite"
                ].tolist(),
            )
        )
        print (filledSites)
        GenomicSites = pd.Series(
            self.GenFrame.InsertionProbability.values,
            index=self.GenFrame.InsertionSiteID,
        ).copy(deep=True)
        GenomicSites.drop(filledSites, inplace=True)
        for i in filledSites:
            transpositionRate = self.TranspFrame.loc[
                self.TranspFrame["TID"] == i, "TraRate"
            ]
            if transpositionRate > np.random_intel.uniform(0, 1.0):
                insertionSite = np.random_intel.choice(
                    GenomicSites.index.values,
                    1,
                    p=GenomicSites.InsertionProbability.values,
                )
                newParent = random.choice(["M", "F"])
                GenomicSites.drop(insertionSite, inplace=True)
                TID = uuid.uuid4().hex
                row = pd.Series(
                    {
                        "TID": TID,
                        "InsertionSite": insertionSite,
                        "TraRate": transpositionRate,
                        "Name": generate_slug(),
                        "Class": self.TranspFrame.loc[
                            self.TranspFrame["TID"] == i, "Class"
                        ],
                        "Traceback": [i],
                        "Generation": 1,
                        "Parent": newParent,
                    }
                )
                if newParent == "M":
                    TEfather.append(TID)
                else:
                    TEmother.append(TID)
                self.TranspFrame = self.TranspFrame.append(row, ignore_index=True)
        return (TEfather, TEmother)

    # /////|   SimRunner   |//////////////////////////////////////////////////////////////////////////
    # /////|   & selection   |////////////////////////////////////////////////////////////////////////
    def runSimulation(self, genMax=10000):
        SimFrame = self.PopFrame.copy(deep=True)
        for i in list(range(genMax)):
            print(i)
            currentPop = pd.DataFrame()
            for k in list(range(SimFrame.shape[0])):
                parentFrame = SimFrame.sample(n=2, weights="NetFitness")
                TIDm = self.recombination(parentFrame.iloc[0])
                TIDf = self.recombination(parentFrame.iloc[1])
                TIDm, TIDf = self.transposition(TIDm, TIDf)
                rowPop = pd.Series(
                    {
                        "PID": uuid.uuid4().hex,
                        "Fitness": 0,
                        "NetFitness": self.fitness(TIDm, TIDf),
                        "Name": generate_slug(),
                        "Sex": "H",
                        "Lineage": ["0"],
                        "Generation": 1,
                        "TEmother": TIDm,
                        "TEfather": TIDf,
                    }
                )
                currentPop = currentPop.append(rowPop, ignore_index=True)
        return 0
