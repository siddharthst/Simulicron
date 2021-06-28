import numpy as np
from numpy import cumsum
from numpy import concatenate as c
import random
from random import shuffle
from math import exp
from fitness import calculateFitness


def generateGenome(
    baseSelection=None,
    baseInsertionProb=1,
    numberOfInsertionSites=10000,
    numberOfChromosomes=6,
    baseRecombinationRate=0.01,
    baseTau=1,
    numberOfPiRNA=6,
    piPercentage=3,
    disablePiRecombination=False,
    DisablePiSelection=False,
):
    # Create piRNA dictionary storing the coordinates
    piRNAcoord = {}
    # Create piRNA array
    piRNArray = np.zeros(numberOfInsertionSites)
    # Define selection penalty for the insertion site
    if baseSelection == None:
        SelectionCoef = np.zeros(numberOfInsertionSites)
    elif baseSelection == "Random":
        SelectionCoef = np.random.normal(-0.02, 0.01, numberOfInsertionSites)
    else:
        SelectionCoef = np.full(numberOfInsertionSites, baseSelection)

    # Define insertion probability for the insertion site
    if baseInsertionProb == None:
        insertionProbability = np.zeros(numberOfInsertionSites)
    elif baseInsertionProb == "Random":
        insertionProbability = np.random.uniform(0.01, 0.99, numberOfInsertionSites)
    else:
        insertionProbability = np.full(numberOfInsertionSites, baseInsertionProb)

    RecombinationRates = np.full(
        numberOfInsertionSites, baseRecombinationRate, dtype=float
    )
    # Generate piRNA information
    totalPiRNALength = int(numberOfInsertionSites * (piPercentage / 100))
    individualPiRNALength = int(totalPiRNALength / numberOfPiRNA)
    # Make chromosomes of equal lengths
    chromosomeLocation = [
        round(x)
        for x in np.linspace(
            0, numberOfInsertionSites, numberOfChromosomes + 1
        ).tolist()
    ]
    # Make chromosomes of random length
    # chromosomeLocation = np.random.choice(
    #     np.arange(RecombinationRates.size), replace=False, size=numberOfChromosomes,
    # )
    # Only insert more chromosomes if needed and randomly insert piRNA clusters
    if numberOfChromosomes > 1:
        RecombinationRates[chromosomeLocation[1:-1]] = 0.499
        # Insert piRNA uniformly in chromosomes
        counter = 0
        while counter < numberOfPiRNA:
            for prime5, prime3 in zip(chromosomeLocation, chromosomeLocation[1:]):
                piRNALocation = np.random.choice(
                    np.arange(prime5 + 1, prime3 - individualPiRNALength - 1),
                    replace=False,
                )
                retryCounter = 0
                # Check if the piRNA is overlapping with a previous piRNA
                while (
                    piRNALocation in np.nonzero(piRNArray)[0].tolist()
                    or (piRNALocation + individualPiRNALength)
                    in np.nonzero(piRNArray)[0].tolist()
                ):
                    # get a new location
                    piRNALocation = np.random.choice(
                        np.arange(prime5 + 1, prime3 - individualPiRNALength - 1),
                        replace=False,
                    )
                    retryCounter += 1
                    if retryCounter > 5:
                        raise SystemExit(
                            "Too many chromosomes defined (relative to insertion sites)"
                        )
                        break

                piRNArray[
                    piRNALocation : piRNALocation + individualPiRNALength
                ] = baseTau
                piRNAcoord[counter] = (
                    piRNALocation,
                    piRNALocation + individualPiRNALength,
                )
                counter += 1
                if counter == numberOfPiRNA:
                    break
        counter = 0

    else:
        counter = 0
        # Insert piRNA randomly in single chromosome
        piRNAcoordinates = [
            individualPiRNALength * i + x
            for i, x in enumerate(
                sorted(
                    random.sample(
                        range(numberOfInsertionSites - individualPiRNALength),
                        numberOfPiRNA,
                    )
                )
            )
        ]
        for i in piRNAcoordinates:
            piRNArray[i : i + individualPiRNALength] = baseTau
            piRNAcoord[counter] = (
                i,
                i + individualPiRNALength,
            )
            counter += 1
        counter = 0
    piRNAindices = np.nonzero(piRNArray)[0].tolist()

    # To disable piRNA selection
    if DisablePiSelection != False:
        # Check if the value is of type float
        if isinstance(DisablePiSelection, float):
            SelectionCoef[piRNAindices] = DisablePiSelection
        else:
            # No selection
            SelectionCoef[piRNAindices] = 0.0

    # To enable disable recombination
    if disablePiRecombination != False:
        # Check if the value is of type float
        if isinstance(disablePiRecombination, float):
            RecombinationRates[piRNAindices] = disablePiRecombination
        else:
            # No recombination in piRNA cluster
            RecombinationRates[piRNAindices] = 0

    # For recombination
    rate2Map = np.insert(
        np.cumsum(-0.5 * np.log(1 - (2 * RecombinationRates))),
        0,
        0,
        axis=0,
    )

    genome = np.vstack(
        (
            SelectionCoef,
            insertionProbability,
            RecombinationRates,
            piRNArray,
        )
    ).T
    return (
        genome,
        piRNAcoord,
        piRNAindices,
        rate2Map,
    )


def generatePopulation(
    genomeMatrix,
    piRNAindices,
    NumberOfIndividual=1000,
    NumberOfTransposonTypes=2,
    NumberOfInsertionsPerType=[2, 2],
    FrequencyOfInsertions=[0.5, 0.5],
    ExcisionRates=[0, 0],
    RepairRates=[1, 1],
    InsertionRates=[1, 1],
    HardyWeinberg=False,
    numberOfPiRNA=6,
    Mig=False,
):
    # Generate arrays for managing transposon content
    transposonID = [0]
    transposonInsertionSite = [0]
    transposonSelectionPenalty = [0]
    transposonExcision = [0]
    transposonRepair = [0]
    transposonInsertion = [0]

    # Generate aux information
    NumberOfInsertionSites = len(genomeMatrix)

    # Create sets to store transposon pedigree
    TEset = {}
    for i in range(NumberOfTransposonTypes):
        TEset[i + 1] = set()

    # Empty population array
    population = np.zeros(
        (NumberOfIndividual, 3),
        dtype=np.ndarray,
    )
    if Mig:
        assert (
            NumberOfTransposonTypes == 2
        ), "Only two transposons are supported with unidirectional migration"
        populationMig = np.zeros(
            (int(NumberOfIndividual * FrequencyOfInsertions[1]), 3),
            dtype=np.ndarray,
        )
        populationMig[0:, 2] = 1

    # Set base fitness
    population[0:, 2] = 1
    assert (
        NumberOfTransposonTypes
        == len(NumberOfInsertionsPerType)
        == len(FrequencyOfInsertions)
    ), "Length mismatch for insertions per type and total number of insertions/frequency information"

    # If HardyWeinberg distribution is required
    # Do not use this option! In major regression
    if HardyWeinberg != False:
        print("Not yet implemented")
        pass

    elif Mig == False:
        # Performing this for each transposon independently
        indices = list(range(NumberOfIndividual))
        counter = 1
        shuffle(indices)
        for i in list(range(NumberOfTransposonTypes)):
            shuffle(indices)
            for k in range(int(NumberOfIndividual * FrequencyOfInsertions[i])):
                indices = indices[
                    0 : int(NumberOfIndividual * FrequencyOfInsertions[i])
                ]
                for l in range(NumberOfInsertionsPerType[i]):
                    # Create transposon
                    shuffle(indices)
                    transposonID.append(i + 1)
                    transposonInsertionSite.append(
                        random.choice(
                            [
                                o
                                for o in list(range(NumberOfInsertionSites))
                                if o not in piRNAindices
                            ]
                        )
                    )
                    transposonSelectionPenalty.append(
                        genomeMatrix[transposonInsertionSite[-1]][0]
                    )
                    transposonExcision.append(ExcisionRates[i])
                    transposonRepair.append(RepairRates[i])
                    transposonInsertion.append(InsertionRates[i])
                    transposonMatrix = np.array(
                        [
                            transposonID,
                            transposonInsertionSite,
                            transposonSelectionPenalty,
                            transposonExcision,
                            transposonRepair,
                            transposonInsertion,
                        ],
                        dtype="object",
                    ).T

                    # Add transposon to set
                    TEset[i + 1].add(counter)

                    # Insert transposon into individuals
                    allele = random.choice([0, 1])
                    if population[indices[k]][allele] == 0:
                        population[indices[k]][allele] = [counter]
                    else:
                        population[indices[k]][allele].append(counter)
                    # Calculate fitness
                    population[indices[k]][2] = calculateFitness(
                        transposonMatrix,
                        v1=population[indices[k]][0],
                        v2=population[indices[k]][1],
                    )
                    counter += 1
            indices = list(range(NumberOfIndividual))

    elif Mig:
        # Creating two sets of population
        counter = 1
        for i in list(range(NumberOfTransposonTypes)):
            if i == 0:
                populationCur = population
                indicesCur = list(range(NumberOfIndividual))
                shuffle(indicesCur)
            if i == 1:
                populationCur = populationMig
                indicesCur = list(range(len(populationMig)))
                shuffle(indicesCur)
            for k in range(int(len(populationCur) * FrequencyOfInsertions[i])):
                indices = indicesCur[
                    0 : int(len(populationCur) * FrequencyOfInsertions[i])
                ]
                for l in range(NumberOfInsertionsPerType[i]):
                    shuffle(indices)
                    transposonID.append(i + 1)
                    transposonInsertionSite.append(
                        random.choice(
                            [
                                o
                                for o in list(range(NumberOfInsertionSites))
                                if o not in piRNAindices
                            ]
                        )
                    )
                    transposonSelectionPenalty.append(
                        genomeMatrix[transposonInsertionSite[-1]][0]
                    )
                    transposonExcision.append(ExcisionRates[i])
                    transposonRepair.append(RepairRates[i])
                    transposonInsertion.append(InsertionRates[i])
                    transposonMatrix = np.array(
                        [
                            transposonID,
                            transposonInsertionSite,
                            transposonSelectionPenalty,
                            transposonExcision,
                            transposonRepair,
                            transposonInsertion,
                        ],
                        dtype="object",
                    ).T

                    # Add transposon to set
                    TEset[i + 1].add(counter)

                    # Insert transposon into individuals
                    allele = random.choice([0, 1])
                    if populationCur[indices[k]][allele] == 0:
                        populationCur[indices[k]][allele] = [counter]
                    else:
                        populationCur[indices[k]][allele].append(counter)
                    # Calculate fitness
                    populationCur[indices[k]][2] = calculateFitness(
                        transposonMatrix,
                        v1=population[indices[k]][0],
                        v2=population[indices[k]][1],
                    )
                    counter += 1

    # Some cleaning
    transposonMatrix = np.array(
        [
            transposonID,
            transposonInsertionSite,
            transposonSelectionPenalty,
            transposonExcision,
            transposonRepair,
            transposonInsertion,
        ],
        dtype="object",
    ).T
    if Mig:
        return population, populationMig, transposonMatrix, TEset

    return population, transposonMatrix, TEset


def initHGT(
    populationMatrix,
    transposonMatrix,
    genomeMatrix,
    TEset,
    piRNAindices,
    numberOfTranspositionEvents,
    FrequencyOfInsertion,
    ExcisionRate,
    RepairRate,
    InsertionRate,
):
    popSize = len(populationMatrix)
    numCarrier = int(FrequencyOfInsertion * popSize)
    # Randomly select the carriers
    populationIndices = list(range(popSize))
    shuffle(populationIndices)
    carrierIndices = populationIndices[0:numCarrier]
    for i in carrierIndices:
        numberOfTranspositionEvents += 1
        site = random.choice(
            [i for i in range(len(genomeMatrix)) if i not in piRNAindices]
        )
        transposonMatrix[numberOfTranspositionEvents, 0] = 2
        transposonMatrix[numberOfTranspositionEvents, 1] = site
        transposonMatrix[numberOfTranspositionEvents, 2] = genomeMatrix[site][0]
        transposonMatrix[numberOfTranspositionEvents, 3] = ExcisionRate
        transposonMatrix[numberOfTranspositionEvents, 4] = RepairRate
        transposonMatrix[numberOfTranspositionEvents, 5] = InsertionRate
        # Randomly select an allele
        allele = random.choice([0, 1])
        if isinstance(populationMatrix[i, allele], list):
            # Check for pre-existing transposon in a site
            TElocations = (
                transposonMatrix[populationMatrix[i, allele], 1].astype(int).tolist()
            )
            if site in TElocations:
                populationMatrix[i, allele][
                    TElocations.index(site)
                ] = numberOfTranspositionEvents
            else:
                populationMatrix[i, allele].append(numberOfTranspositionEvents)
        else:
            populationMatrix[i, allele] = [numberOfTranspositionEvents]
        # Add transposon to set
        TEset[2].add(numberOfTranspositionEvents)
    return populationMatrix, transposonMatrix, TEset, numberOfTranspositionEvents
