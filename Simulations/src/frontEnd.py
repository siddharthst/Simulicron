# ~ Copyright or © or Copr. CNRS, contributor(s) : Siddharth Sing Tomar (2020-2023)

# ~ Contact: siddharth.egce@gmail.com , arnaud.le-rouzic@universite-paris-saclay.fr

# ~ This software is a computer program whose purpose is to [describe
# ~ functionalities and technical features of your software].

# ~ This software is governed by the CeCILL license under French law and
# ~ abiding by the rules of distribution of free software.  You can  use, 
# ~ modify and/ or redistribute the software under the terms of the CeCILL
# ~ license as circulated by CEA, CNRS and INRIA at the following URL
# ~ "http://www.cecill.info". 

# ~ As a counterpart to the access to the source code and  rights to copy,
# ~ modify and redistribute granted by the license, users are provided only
# ~ with a limited warranty  and the software's author,  the holder of the
# ~ economic rights,  and the successive licensors  have only  limited
# ~ liability. 

# ~ In this respect, the user's attention is drawn to the risks associated
# ~ with loading,  using,  modifying and/or developing or reproducing the
# ~ software by the user in light of its specific status of free software,
# ~ that may mean  that it is complicated to manipulate,  and  that  also
# ~ therefore means  that it is reserved for developers  and  experienced
# ~ professionals having in-depth computer knowledge. Users are therefore
# ~ encouraged to load and test the software's suitability as regards their
# ~ requirements in conditions enabling the security of their systems and/or 
# ~ data to be ensured and,  more generally, to use and operate it in the 
# ~ same conditions as regards security. 

# ~ The fact that you are presently reading this means that you have had
# ~ knowledge of the CeCILL license and that you accept its terms.



import click
import time
import pickle
from popSim import runBatch

import click
import ast


class PythonLiteralOption(click.Option):
    def type_cast_value(self, ctx, value):
        try:
            return ast.literal_eval(value)
        except:
            raise click.BadParameter(value)


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
@click.option("--npi", default=2, help="Number of piRNA", show_default=True)
@click.option("--tau", default=2, help="piRNA regulatory factor", show_default=True)
@click.option("--piper", default=2, help="piRNA length - as percentage of genome", show_default=True)
@click.option("--nt", default=2, help="Number of transposon familes", show_default=True)
@click.option(
    "--tif",
    default=[0.5, 0.5],
    help="Transposon insertion frequencies",
    show_default=True,
    cls=PythonLiteralOption,
)
@click.option(
    "--tin",
    default=[1, 1],
    help="Number of insertions per individual",
    show_default=True,
    cls=PythonLiteralOption,
)
@click.option(
    "--be",
    default=[0.1, 0.1],
    help="Base excision rate for transposons",
    cls=PythonLiteralOption,
    show_default=True,
)
@click.option(
    "--ng", default=10000, help="Maximum number of generations", show_default=True
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
def cli(ns, nc, ni, br, nind, npi, tau, piper, nt, tif, tin, be, ng, threads, output):
    result = runBatch(
        numberOfSimulations=ns,
        numberOfChromosomes=nc,
        numberOfInsertionSites=ni,
        baseRecombinationRate=br,
        NumberOfIndividual=nind,
        baseTau=tau,
        numberOfPiRNA=npi,
        piPercentage=piper,
        NumberOfTransposonTypes=nt,
        NumberOfInsertionsPerType=tin,
        FrequencyOfInsertions=tif,
        NumberOfGenerations=ng,
        baseSelection=1,
        ExcisionRates=be,
        RepairRates=1,
        InsertionRates=1,
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
