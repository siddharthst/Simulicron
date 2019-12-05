# Imports
import inspect
import random
import tempfile
import uuid
from functools import wraps

import coloredlogs
import dask as dd
import logging
import numpy
import pandas as pd
from coolname import generate_slug

# Logger
logger = logging.getLogger(__name__)
coloredlogs.install(level='DEBUG')

# StackOverflow snippet #1
########################################################################################################################
# /////|   Decorator   |////////////////////////////////////////////////////////////////////////////////////////////////
########################################################################################################################


def auto_assign_arguments(function):

    @wraps(function)
    def wrapped(self, *args, **kwargs):
        _assign_args(self, list(args), kwargs, function)
        function(self, *args, **kwargs)

    return wrapped


########################################################################################################################
# /////|   Utils   |////////////////////////////////////////////////////////////////////////////////////////////////////
########################################################################################################################


def _assign_args(instance, args, kwargs, function):

    def set_attribute(instance, parameter, default_arg):
        if not (parameter.startswith("_")):
            setattr(instance, parameter, default_arg)

    def assign_keyword_defaults(parameters, defaults):
        for parameter, default_arg in zip(reversed(parameters),
                                          reversed(defaults)):
            set_attribute(instance, parameter, default_arg)

    def assign_positional_args(parameters, args):
        for parameter, arg in zip(parameters, args.copy()):
            set_attribute(instance, parameter, arg)
            args.remove(arg)

    def assign_keyword_args(kwargs):
        for parameter, arg in kwargs.items():
            set_attribute(instance, parameter, arg)

    def assign_keyword_only_defaults(defaults):
        return assign_keyword_args(defaults)

    def assign_variable_args(parameter, args):
        set_attribute(instance, parameter, args)

    POSITIONAL_PARAMS, VARIABLE_PARAM, _, KEYWORD_DEFAULTS, _, KEYWORD_ONLY_DEFAULTS, _ = inspect.getfullargspec(
        function)
    POSITIONAL_PARAMS = POSITIONAL_PARAMS[1:]  # remove 'self'

    if (KEYWORD_DEFAULTS):
        assign_keyword_defaults(parameters=POSITIONAL_PARAMS,
                                defaults=KEYWORD_DEFAULTS)
    if (KEYWORD_ONLY_DEFAULTS):
        assign_keyword_only_defaults(defaults=KEYWORD_ONLY_DEFAULTS)
    if (args): assign_positional_args(parameters=POSITIONAL_PARAMS, args=args)
    if (kwargs): assign_keyword_args(kwargs=kwargs)
    if (VARIABLE_PARAM):
        assign_variable_args(parameter=VARIABLE_PARAM, args=args)


########################################################################################################################


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

    @auto_assign_arguments
    def __init__(self,
                 tcount=1,
                 ttype=[1],
                 popsize=100,
                 insize=2000,
                 gensize=1000,
                 trate=[0.02],
                 tpenalty=[-0.02],
                 tparent=["Mother"]):
        """Class constructor
        :param tcount: Number of transposons to be present in initial population
        :param ttype: Type of transposons (class 1 or 2)
        :param popsize: Population size
        :param insize: Number of insertion sites
        :param gensize: Size of genome on cM
        :param trate: Transposition rates
        :param tpenalty: Selection penalty for each transposon
        :param tparent: Parent carrying the transposon insertion (Mother/Father)
        :type tcount: int
        :type ttype: list [int]
        :type popsize: int
        :type insize: int
        :type gensize: int
        :type trate: list [float]
        :type tpenalty: list [float]
        :type arg1: tparent [string]
        """
        # Dataframe definations
        self.TranspFrame = pd.DataFrame(columns=[
            'TID', 'Position', 'SelCo', 'Name', 'Class', 'Traceback',
            'Generation', 'Parent', "TraRate"
        ])
        self.PopFrame = pd.DataFrame(columns=[
            'PID', 'Fitness', 'Name', 'Sex', 'Lineage', 'Generation',
            'TEfather', 'TEmother', 'Insertion_Father', 'Insertion_Mother'
        ])

    # Init transposons
    def initT(self):
        """Method to dataframe containing initial transposon population
        
        :return: Dataframe containing Transposon information
        :rtype: Dataframe
        """

        if (self.tcount > len(self.tpenalty)):
            logger.info(
                "Mismatch between transposon count and selection penalties. Using default for each transposon count!"
            )
            self.tpenalty = [-0.02] * self.tcount
        if (self.tcount > len(self.trate)):
            logger.info(
                "Mismatch between transposon count and transposition rates. Using default for each transposon count!"
            )
            self.trate = [0.02] * self.tcount
        if (self.tcount > len(self.ttype)):
            logger.info(
                "Mismatch between transposon count and transposon types. Using default for each transposon count!"
            )
            self.ttype = [1] * self.tcount
        if (self.tcount > len(self.tparent)):
            logger.info(
                "Mismatch between transposon count and transposon parent. Using default for each transposon count!"
            )
            self.tparent = ["Mother"] * self.tcount

        # Create random filled insertion sites
        inSiteArray = random.sample(range(1, self.insize), self.tcount)

        for i in range(0, self.tcount):
            row = pd.Series({
                'TID': uuid.uuid4().hex,
                'Position': inSiteArray[i],
                'TraRate': self.trate[i],
                'SelCo': self.tpenalty[i],
                'Name': generate_slug(),
                'Class': self.ttype[i],
                'Traceback': ['0'],
                'Generation': 1,
                'Parent': self.tparent[i]
            })
            self.TranspFrame = self.TranspFrame.append(row, ignore_index=True)
        return (self.TranspFrame)

    # Init population and genome
    def initPG(self):
        """Method to create initial population and their respective genomes
        
        :return: tuple(population,genome)
            WHERE
            population is population dataframe
            genome is genome dataframe
        :rtype: DataFrame
        """

        # Create transposon insertions in randomly selected individuals
        IndividualToInsert = random.sample(list(range(1, self.popsize)),
                                           self.tcount)
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
            if (i in IndividualToInsert):
                TE = TIDlist[TIDcounter]
                TIDcounter += 1
                Parent = self.TranspFrame[self.TranspFrame['TID'] ==
                                          TE]['Parent'].values[0]
                FitnessPen = self.TranspFrame[self.TranspFrame['TID'] ==
                                              TE]['SelCo'].values[0]
                if (Parent == "Mother"):
                    insertion_Mother = self.TranspFrame[
                        self.TranspFrame['TID'] == TE]['Position'].values[0]
                    TEmother = TE

                if (Parent == "Father"):
                    insertion_Father = self.TranspFrame[
                        self.TranspFrame['TID'] == TE]['Position'].values[0]
                    TEfather = TE

            else:
                TE = '0'
                Parent = "0"
                insertion_Father = 0
                insertion_Mother = 0
                FitnessPen = 0
                TEmother = "0"
                TEfather = "0"

            # Populate the population!
            # Define intial fitness
            fitness = random.uniform(0.6, 1.0)
            rowPop = pd.Series({
                'PID': uuid.uuid4().hex,
                'Fitness': fitness,
                'NetFitness': fitness + FitnessPen,
                'Name': generate_slug(),
                'Sex': 'H',
                'Lineage': ['0'],
                'Generation': 1,
                'Insertion_Father': [insertion_Father],
                'Insertion_Mother': [insertion_Mother],
                'TEmother': [TEmother],
                'TEfather': [TEfather]
            })
            self.PopFrame = self.PopFrame.append(rowPop, ignore_index=True)

        self.PopFrame['Lineage'] = self.PopFrame['Lineage'].astype('object')
        self.PopFrame['Insertion_Father'] = self.PopFrame[
            'Insertion_Father'].astype('object')
        self.PopFrame['Insertion_Mother'] = self.PopFrame[
            'Insertion_Mother'].astype('object')
        return (self.PopFrame)

    def createSim(self):
        """Method to generate the initial simulation dataset
        
        :return: List
            WHERE 
            index 0 is transposon dataframe
            index 1 is population dataframe
            index 2 is genome dataframe
        :rtype: list
        """

        transposon = self.initT()
        genome = self.initPG()
        return ([transposon, genome])

########################################################################################################################
