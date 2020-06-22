import ClassDeco as cd
import pandas as pd
import random
import numpy as np
import uuid
from coolname import generate_slug


import logging
import coloredlogs

# Logger
logger = logging.getLogger(__name__)
coloredlogs.install(level="DEBUG")

########################################################################################################################
# /////|   Class   |////////////////////////////////////////////////////////////////////////////////////////////////////
########################################################################################################################
class initSim:
    """Class to initialize the progenitor (F0) population.

    initSim primarily creates an object which stores information about the 
    initial population. Moreover the inbuilt methods will return pandas
    dataframe for direct modification as desired.
    
    :return: Object
    :rtype: object type initSim
    """

    @cd.auto_assign_arguments
    def __init__(
        self,
        tcount=1,
        ttype=[1],
        popsize=100,
        insize=2000,
        trate=[0.02],
        tparent=["Mother"],
        nChr=50,
    ):
        """Class constructor
        :param tcount: Number of transposons to be present in initial population
        :param ttype: Type of transposons (class 1 or 2)
        :param popsize: Population size
        :param insize: Number of insertion sites
        :param trate: Transposition rates
        :param tparent: Parent carrying the transposon insertion (Mother/Father)
        :param nChr: Number of chromosomes
        :type tcount: int
        :type ttype: list [int]
        :type popsize: int
        :type insize: int
        :type trate: list [float]
        :type tparent: list [string]
        :type nChr: int
        """
        # Dataframe definations
        self.TranspFrame = pd.DataFrame(
            columns=[
                "TID",
                "Name",
                "Class",
                "Traceback",
                "Generation",
                "Parent",
                "TraRate",
            ]
        )
        self.PopFrame = pd.DataFrame(
            columns=[
                "PID",
                "Fitness",
                "Name",
                "Sex",
                "Lineage",
                "Generation",
                "TEfather",
                "TEmother",
            ]
        )
        self.GenFrame = pd.DataFrame(
            columns=[
                "InsertionSiteID",
                "InsertionProbability",
                "RecombinationRate",
                "SelectionCoef",
                "Filled",
            ]
        )

    def generateGenome(self):
        NumberInsertionSites = self.insize
        NumberChromosomes = self.nChr
        ###-----------Function main-----------###
        genome = list(range(1, NumberInsertionSites))
        InsertionSiteColumn = list(range(1, NumberInsertionSites))
        RecombinationRate = [0.01] * len(genome)
        chrLocation = np.random.choice(
            list(range(1, NumberInsertionSites - 10)), NumberChromosomes, replace=False,
        )
        SelectionCoef = np.random.normal(-0.02, 0.01, len(genome))
        insertionProbability = np.random.uniform(0.01, 0.99, len(genome))
        insertionProbability /= insertionProbability.sum()
        for i in chrLocation:
            RecombinationRate[i] = 0.5

        genomeDict = {
            "InsertionSiteID": InsertionSiteColumn,
            "InsertionProbability": insertionProbability,
            "RecombinationRate": RecombinationRate,
            "SelectionCoef": SelectionCoef,
        }
        genome = pd.DataFrame(genomeDict)

        # ###-----------Function plot-----------###
        # ### Plot the simulated variables
        # sns.set(style="ticks", palette="muted", color_codes=True)
        # # Set up the matplotlib figure
        # f, axes = plt.subplots(1, 3, figsize=(15, 7))
        # chartRecombinationRates = sns.distplot(
        #     RecombinationRate,
        #     ax=axes[0],
        #     kde=False,
        #     axlabel="Histogram - Recombination rate",
        # )
        # chartInsertionProbability = sns.distplot(
        #     insertionProbability,
        #     ax=axes[1],
        #     kde=False,
        #     axlabel="Histogram - Insertion probability",
        # )
        # chartSelectionCoef = sns.distplot(
        #     SelectionCoef,
        #     ax=axes[2],
        #     kde=False,
        #     axlabel="Histogram - Selection coefficient",
        # )

        ###---------------Return--------------###
        self.GenFrame = genome
        return genome

    # Init transposons
    def initT(self):
        """Method to dataframe containing initial transposon population
        
        :return: Dataframe containing Transposon information
        :rtype: Dataframe
        """
        if self.tcount > len(self.trate):
            logger.info(
                "Mismatch between transposon count and transposition rates. Using default for each transposon count!"
            )
            self.trate = [0.02] * self.tcount
        if self.tcount > len(self.ttype):
            logger.info(
                "Mismatch between transposon count and transposon types. Using default for each transposon count!"
            )
            self.ttype = [1] * self.tcount
        if self.tcount > len(self.tparent):
            logger.info(
                "Mismatch between transposon count and transposon parent. Using default for each transposon count!"
            )
            self.tparent = ["Mother"] * self.tcount

        # Create random filled insertion sites
        inSiteArray = random.sample(range(1, self.insize - 10), self.tcount)

        for i in range(0, self.tcount):
            row = pd.Series(
                {
                    "TID": uuid.uuid4().hex,
                    "InsertionSite": inSiteArray[i],
                    "TraRate": self.trate[i],
                    "Name": generate_slug(),
                    "Class": self.ttype[i],
                    "Traceback": ["0"],
                    "Generation": 1,
                    "Parent": self.tparent[i],
                }
            )
            self.TranspFrame = self.TranspFrame.append(row, ignore_index=True)
            self.TranspFrame["InsertionSite"] = self.TranspFrame[
                "InsertionSite"
            ].astype(int)
        return self.TranspFrame

    # Init population
    def initPG(self):
        """Method to create initial population and their respective genomes
        
        :return: tuple(population,genome)
            WHERE
            population is population dataframe
            genome is genome dataframe
        :rtype: DataFrame
        """

        # Create transposon insertions in randomly selected individuals
        IndividualToInsert = random.sample(list(range(1, self.popsize)), self.tcount)
        TIDlist = self.TranspFrame.TID.tolist()
        TIDcounter = 0
        Parent = "0"
        insertion_Father = 0
        insertion_Mother = 0
        FitnessPen = 0
        TEfather = "0"
        TEmother = "0"
        for i in range(self.popsize):
            # In case this (un)lucky individual has transposon insertion
            if i in IndividualToInsert:
                TE = TIDlist[TIDcounter]
                TIDcounter += 1
                insertionSiteID = self.TranspFrame[self.TranspFrame["TID"] == TE][
                    "InsertionSite"
                ].values[0]
                Parent = self.TranspFrame[self.TranspFrame["TID"] == TE][
                    "Parent"
                ].values[0]
                FitnessPen = self.GenFrame[
                    self.GenFrame["InsertionSiteID"] == insertionSiteID
                ]["SelectionCoef"].values[0]

                if Parent == "Mother":
                    insertion_Mother = self.TranspFrame[self.TranspFrame["TID"] == TE][
                        "InsertionSite"
                    ].values[0]
                    TEmother = TE

                if Parent == "Father":
                    insertion_Father = self.TranspFrame[self.TranspFrame["TID"] == TE][
                        "InsertionSite"
                    ].values[0]
                    TEfather = TE

            else:
                TE = "0"
                Parent = "0"
                insertion_Father = 0
                insertion_Mother = 0
                FitnessPen = 0
                TEmother = 0
                TEfather = 0

            # Populate the population!
            # Define intial fitness
            fitness = random.uniform(0.6, 1.0)
            rowPop = pd.Series(
                {
                    "PID": uuid.uuid4().hex,
                    "Fitness": fitness,
                    "NetFitness": fitness + FitnessPen,
                    "Name": generate_slug(),
                    "Sex": "H",
                    "Lineage": ["0"],
                    "Generation": 1,
                    # "Insertion_Father": [insertion_Father],
                    # "Insertion_Mother": [insertion_Mother],
                    "TEmother": [TEmother],
                    "TEfather": [TEfather],
                }
            )
            self.PopFrame = self.PopFrame.append(rowPop, ignore_index=True)

        self.PopFrame["Lineage"] = self.PopFrame["Lineage"].astype("object")
        # self.PopFrame["Insertion_Father"] = self.PopFrame["Insertion_Father"].astype(
        #     "object"
        # )
        # self.PopFrame["Insertion_Mother"] = self.PopFrame["Insertion_Mother"].astype(
        #     "object"
        # )
        return self.PopFrame

    def createSim(self):
        """Method to generate the initial simulation dataset
        
        :return: List
            WHERE 
            index 0 is transposon dataframe
            index 1 is population dataframe
            index 2 is genome dataframe
        :rtype: list
        """

        genome = self.generateGenome()
        transposon = self.initT()
        population = self.initPG()
        return [transposon, population, genome]


########################################################################################################################
