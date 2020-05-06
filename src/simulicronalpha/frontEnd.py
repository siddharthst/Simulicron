import click
import time
import pickle
from popSim import runBatch

import click


@click.command(no_args_is_help=True)
@click.option(
    "--ns", default=1000, help="Number of distinct simulations", show_default=True
)
@click.option(
    "--nc", default=6, help="Number of chromosomes in the genome", show_default=True
)
@click.option(
    "--ni",
    default=1000,
    help="Number of insertion sites in the genome",
    show_default=True,
)
@click.option(
    "--br",
    default=0.1,
    help="Base recombination rate in genome (recombination rate between chromosomes is always .5)",
    show_default=True,
)
@click.option(
    "--nind",
    default=1000,
    help="Number of individuals in a simulation",
    show_default=True,
)
@click.option(
    "--insrtall",
    default=False,
    help="Insert transposon in all individuals(bool)?",
    show_default=True,
)
@click.option(
    "--nt", default=2, help="Number of transposon insertions(unique)", show_default=True
)
@click.option(
    "--ng", default=100000, help="Maximum number of generations", show_default=True
)
@click.option(
    "--be", default=1, help="Base excision rate for transposons", show_default=True
)
@click.option(
    "--ifr",
    default=False,
    help="Insertion frequency for transposons (can conflict with othe options!)",
    show_default=True,
)
@click.option(
    "--hw",
    default=False,
    help="Follow HardyWeinberg distribution ?(Bool) (can conflict with othe options!)",
    show_default=True,
)
@click.option(
    "--threads", default=1, help="Number of threads used by program", show_default=True,
)
@click.option(
    "--output",
    default="sim_" + time.strftime("%Y%m%d_%H%M%S") + ".pickle",
    help="Name of the file containing result (default will use the current time)",
    show_default=True,
)
def cli(ns, nc, ni, br, nind, insrtall, nt, ng, be, ifr, hw, threads, output):
    result = runBatch(
        numberOfSimulations=ns,
        numberOfChromosomes=nc,
        numberOfInsertionSites=ni,
        baseRecombinationRate=br,
        NumberOfIndividual=nind,
        InsertIntoOne=False,
        InsertIntoAll=insrtall,
        NumberOfTransposonInsertions=nt,
        NumberOfGenerations=ng,
        baseSelection=1,
        baseExcision=be,
        baseRepair=1,
        baseInsertion=1,
        consecutiveTransposons=False,
        changeRecombination=False,
        baseTrRecombination=0.1,
        insertionFrequency=ifr,
        HardyWeinberg=hw,
        numberOfThreads=threads,
    )
    with open(output, "wb") as f:
        pickle.dump((result), f)


if __name__ == "__main__":
    cli()

"""
@click.command()
@click.option("--numberOfSimulations|ns", default=1000, help="Number of distinct simulations")
@click.option("--numberOfChromosomes|nc", default=6, help="Number of chromosomes in the genome")
@click.option("--numberOfInsertionSites|ni", default=1000, help="Number of insertion sites in the genome")
@click.option("--baseRecombinationRate|br", default=0.1, help="Base recombination rate in genome (recombination rate between chromosomes is always .5)")
@click.option("--NumberOfIndividual|nind", default=1000, help="Number of individuals in a simulation")
@click.option("--InsertIntoAll|insrtall", default=False, help="Insert transposon in all individuals(bool)?")
@click.option("--NumberOfTransposonInsertions|nt", default=2, help="Number of transposon insertions(unique)")
@click.option("--NumberOfGenerations|ng", default=100000, help="Maximum number of generations")
@click.option("--baseExcision|be", default=1, help="Base excision rate for transposons")
@click.option("--insertionFrequency|ifr", default=False, help="Insertion frequency for transposons (can conflict with othe options!)")
@click.option("--HardyWeinberg|hw", default=False, help="Follow HardyWeinberg distribution ?(Bool) (can conflict with othe options!)")"""
