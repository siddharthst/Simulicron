# Research lines

## Project P1: Population genetics of piRNA-regulated transposable elements

### Background

Transposable elements are "selfish DNA" sequences which dynamics cannot be understood without specific population genetics models. Pop gen TE models generally root in Charlesworth & Charkesworth 1983, which hypothesised two (non-exclusive) mechanisms to limit the (otherwise exponential) increase in active TE numbers: transposition regulation and natural selection. Transposition regulation in the early 1980s was essentially phenomenological: transposition rate decreases during the invasion, but no specific mechanism was posited. In Charlesworth & Charlesworth 1983 (and following work, including Charlesworth 1991, Le Rouzic & Capy 2005...), the transpositon rate per copy was expressed as u(n), meaning that the transposition rate was negatively associated with the number of copies in the genome (generally with an exponential decay, as u(n) = u(0) exp(-kn)). Work on TE regulation in the late XXth century was a bit confusing, as contradictory or misleading results emerged from model species cases (e.g. P element in D. melanogaster). The leading idea was "self-regulation", based on the conjecture that TE had an evolutionary benefit in regulating its transposition rate both in soma and germline; the best known example being an alternative splicing mechanism and/or maternally-transmitted factors in the P element. Other explanations involved biochemical constraints (such as the OverProduction Inhibition in class II elements, based on the idea that when transposase becomes inactivated when expressed too much). Later, it has also been suggested that mutations in TE sequences could also participate to the decrease in TE activity, either by eliminating passively TE copies from the pool of "transposable" sequences (Kijima & Innan 2013), or by introducing non-autonomous competing copies (Le Rouzic & Capy 2007). 

The importance of small RNA in TE regulation came up in the 2010s, so that basically all the literature on TE dynamics was based on uncertain and/or wrong or dubious regulation processes. Only a very few recent papers (including Kofler 2019 + 2 earlier preprint refs) model explicitly RNA-induced regulation (double check if there was not something vaguely similar in an old CHarlesworth paper, like Charlesworth 1991 or Charlesworth & Langley 198x). RNA regulation is based on the presence of regulatory TE copies, that appear through a regular transposition process. A TE insertion becomes "regulatory" when it lands into a specific part of the genome ("cluster"). When such a copy is present in the genome, the transposition rate of all TEs from the same family decreases. This type of regulation may differ from the "traditional" copy-number based regulation models by several important aspects: (i) regulation is heritable (regulation does not stop if the number of copies decrease), (ii) the occurence of regulation is associated with the transposition rate, and not with the copy number, (iii) if TEs are deleterious, regulatory copies are under positive selection, (iv) TE regulatory sites are in linkage disequilibrium with TEs (which could matter a lot in case of migrations), (v) potential maternal inheritance of the regulatory factor. 

### Questions

In spite of the presence of prior work on this (mainly by R. Kofler, other papers tend to be of lower quality), there remains some dark spots or unanswered issues. 

#### P1.1: analytical model of RNA-regulation dynamics

So far, the only published models are individual-based simulations. Analytical approches (e.g. through difference equations) would bring a lot of undetsranding about the parameters that matter in e.g. (i) the max copy number, (ii) the number and distribution of regulatory copies in the genome, (iii) the "host" strategy (number and size of clusters), (iv) efficiency of the regulation...  

##### Pros

* This would improve the state of the art and help including recent TE dynamic models into the existing pop gen framework. 
* Definitely complementary with Kofler 2019. 
* Analytical results need to be confirmed and compared with simulations, the simulation program already exists. 

##### Cons

* Not very exciting, as most of the results were already published. 
* Could be difficult, as it can be quite technical (pop gen modelling and maths). So far TE dynamics difference equations do not have a dynamic solution, and could only be solved assuming an equilibrium. If no equilibrium (which is likely for RNA regulation models), how to approach the analysis? Collaboration necessary? (Who might be interested?) 
* No plan B if the problem cannot be simplified to ensure mathematical tractability.


#### P1.2 TE dynamics in structured populations

The RNA regulation model identifies regulatory loci that are heritable and follow their own dynamics. In practice, the invasion of a TE in a "naive" population (after a HTT event) might then be very different from the invasion from migration, in which the migrant is likely to bring regulatory sites along with active copies. This might affect the spread of the TE in subsequent populations, and the whole dynamics might depend on the population structure of the species. Very few prior work on population models in the TE context: Quesneville & Anxolabehere  Theor Pop Biol 1998, Deceliere et al. Genetics 2005. 

##### Pros

* The question is really opened, and could lead to interesting data-grounded research (e.g. in various Drosophila species which might have drastically different population stuctures). 
* The dynamics in subsequent populations will very likely be different than in the first population, little risk for a "negative result" here. 
* Investigating this question probably does not need large changes in the simulation program (just add "population" classes which individually behaves as in the current algorithm, and add the possibility to exchange migrants at a given rate). 

##### Cons

* Simulations in metapopulations might lead to important performance issues (it requires to handle a lot of individuals at the same time if each population is large). 
* Understanding the metapopulation context requires a good understanding of the mono-population dynamics, which is not established knowledge for this kind of models. 
* The risk of being scooped is substantial, given that R. Kofler is working on tracking TE during the invasion wave in a species (based on sequence data and/or experimental evolution). 

### Project P1.3 Co-invasion by several TE families

So far, what happens when two (or more) TE families invade the genome at the same time has not be seriously studied. Although this might not be a frequent event, it is likely that co-invasions happen (for instance, after an introgression across two close species, which might carry different TE families). Co-invasion will also be the focus on some experimental evolution in the lab, which makes the question relevant. 

Two different TE families may interact at several levels : (i) through classical population-genetics mechanisms (epistasis on fitness), (ii) through sharing the transposition machinery (e.g. autonomous vs non-autonomous copiesà, and (iii) through co-regulation. 

Interactions through selection are not specific to transposable elements, the same issue arise with several alleles segregating at the same time. Basically, if the effects on fitness are independent (log-additive), the dynamics of both families will be exactly the same as when they invade alone. If there is negative epistasis on fitness, deleterious effects will be synergistic, and both families will interact negatively (the opposite for positive selection, but the situation is less realistic). 

Interactions through sharing part of the transposition machinery is interesting, but it is questionable whether these are really two TE families (isn't cross-mobilization by defintion mean that the two kinds of sequences belong to the same family?)

Interactions through regulation is probably the interesting part of the question. Two completely independent families may interact through the dynamics of pi-clusters, but the effect is probably small. However, several close TE families can probably cross-repress, because pi-RNA are small sequences that may match several families. This cross-repression effect could be modeled by a parameter between 0 (independent regulation) and 1 (any TE regulate both families). 

#### Pros

* The question is original
* It might be associated with data analysis (probability for a small RNA to match different families in the genome of some species)
* Interesting side questions can be addressed (e.g. what happens if both invasions do not start at the same time?)

#### Cons

* Only the co-repression part may lead to non-trivial results
* The occurence of co-invasion in natural species is probably low, because invasion is very fast (the question may not be biologically very relevant). 

## Project P2: Modelling genome dynamics

### Background

Understanding genome evolution is in its infancy. The structure of the genome are its diversity across species and population are only known superficially since a couple of decades, and even the most recent sequencing technologies may not grant us full access to all copy number variation and duplication / translocation events. So far, it is probably fait to claim that (i) most of the structural variation in genomes has been neglected, probably because of technical issues that tend to igore them in assembled genomes, and (ii) their importance in evolution (including phenotypic evolution) has been either largely exaggerated, or in the contrary, underestimated. 

Transposable elements are one of the most prominent sequences that are involved on a regular basis in genomic rearragements and other structural polymorphisms. Strictly speaking, their repeated nature make them one of the most direct sources of copy number variation, but they are also invloved more indirectly in segmental duplications or other large-scale events. Although there are vague claims in the literature about hypothetical benefits of TE sequences due to their potential to promote "genome evolvability" (sometimes called "genome plasticity"), the theoretical and empirical bases of such claims remain largely elusive. In particular, there is little evidence that TE-related structural rearrangemets are more (or less) often under positive or negative selection than structural variants of other origins. 

Due to their large number and overwhelming presence in most genomes, TEs are also one of the major sources or DNA loss and gain, and these sequences participate directly to the evolution of genome size. It is known since a long time that species may differ substantially in the amount of DNA in their genome, in a way that is mostly uncorrelated to the number of genes (the "C-value paradox"). However, whether this "extra DNA" is really neutral from a selective point of view remains a matter of debate. 

Addressing such questions from a theoretical point of view remains challenging, because most of the claims from the literature are not backed by formal models (Koonin?). It is also difficult to build models that are not too case-specific, or on the opposite that are not too general. 

### Questions

#### P2.1 A model of co-evolution of TEs and heterochromatin

The "junk DNA" is often structured functionally, and packed into specific regions of condensed chromatin (heterochromatin) tagged by specific epigenetic marks that are never or rarely expressed. The reason why such regions exist may not be completely straightforward, as their emergence could raise a 'chicken or egg' paradox. Are TEs inserting into heterochromatin regions less deleterious (less likely to be expressed or to disrupt an existing gene), or are heterochromatin marks condensed into TE-rich regions to ensure their silencing? Do TE sequences have any interest in attracting (or avoiding) heterochromatic marks or regions? When these regions grow, isn't it just more likely that most TEs land into existing heterochromatic clusters and increase their size in a passive way? 

##### Pros

* The question looks interesting, and may require approches combining population genetics, genomics, modelling, and bioinformatics.
* It remains in the frame of TE evolution, although in this case TEs would not be "agents", but rather a forcing external factor that drives genome evolution. 
* A lot of possibilities, a lot of potential questions to address. 
* Some obvious links with data analysis

##### Cons

* A new field for us, and thus a lot of new literature to digest before being sure to make a contribution to the field.
* A lot of technical challenges to overcome with variable size genomes, in particular when handling recombination. If recombination requires to "blast" two genomes, is the computation time compatible with long-term large-population individual-based simulations? 
* The simulation algorithm would be very different from the TE simulation, and would require a new simulation program. Does the development of two independent simulation programs compatible with the frame of a PhD? 
* The potential to obtain non-tautological results is difficult to estimate (the risk is that the results may simply depend on the necessary hypotheses made when building the model). 
