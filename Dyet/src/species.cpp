/***************************************************************************
                          espece.cpp  -  description
                             -------------------
    begin                : Fri Nov 15 2002
    copyright            : (C) 2002 by Arnaud Le Rouzic
    email                : lerouzic@cnrs-gif.fr
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
 
#include "species.h"
#include "system.h"
#include "exception.h"
 
//--------------------- OBJET ESPECE ---------------------------------------
using namespace std;

void Species::Migration(void)
// Realise l'ensemble des migrations a une generation dans une espece
{
  Update_PopulationNumber();
  if (nb_pop.effectif < 2) return;
  // Code pas beau...

  MigrationByPop();
  // Migration_ind_par_ind(); 
}

Species * Species::Reproduction(void)
{
  /* But de la méthode : assurer la reproduction dans la totalité de
     l'espèce. Renvoie l'espèce issue de la reproduction */
     
  Species * descendance = new Species();
     
  // La reproduction se fait population par population.

  it_pop prems = ref_pop.begin();
  it_pop fin = ref_pop.end();

	while (prems != fin)
	{
		try {
			Population * pop (NULL);
			pop =  (*prems)->Reproduction();
			descendance->AddPopulation(pop);	
			++prems;
		} catch (exec_nomore_indiv_reproduction & e) { 
			delete *prems;
			prems = ref_pop.erase(prems); // works with a list or a vector
			fin = ref_pop.end();
		} catch (exec_toomuch & e) {
			delete descendance;
			throw exec_toomuch();
		}
	}
	if (descendance->PopulationNumber() == 0)
	{
		delete descendance;
		descendance = NULL;
		throw exec_nomore_population();
	}
  UnUpdate();
  return descendance;  
}  

void Species::AddPopulation(Population * pop)
{
  ref_pop.push_back(pop);
  UnUpdate();
}

Ptr2Element Species::PickRandomElement(void) const
{
	assert (TEnumber() > 0);
	Population * pop = ChooseRandomPopByTE();
	return pop->PickRandomElement();
}

Species * Species::VoidCopy(void) const
{
	Species * cpy = new Species();
	
	c_it_pop pop = ref_pop.begin();
	c_it_pop endpop = ref_pop.end();
	for (; pop != endpop; ++pop)
	{
		cpy->AddPopulation((*pop)->VoidCopy());
	}
	return cpy;
}

void Species::ClearElements(void)
{
	it_pop pop, endpop;
	endpop = ref_pop.end();
	for (pop = ref_pop.begin(); pop != endpop; ++pop)
		(*pop)->ClearElements();
}
	
void Species::AddElement(Ptr2Element ptr)
{
	unsigned int what_pop = generator->Uniform_int(PopulationNumber());
	#if (USE_VECTOR_POP)
		ref_pop[what_pop]->AddElement(ptr);
	#endif
	#if (USE_LIST_POP)
		it_pop pop = ref_pop.begin();
		for (unsigned int i = 0; i < what_pop; ++i)
			++pop;
		(*pop)->AddElement(ptr);
	#endif
	UnUpdate();
}

/***************************************************************************/
/*                           Méthodes private                              */
/***************************************************************************/

void Species::MigrationByPop(void)
{
  /* Migration pop par pop :
     le nombre de migrants est déterminé à partir du nombre total d'individus,
     pour chaque migration, on tire une population de départ et une population
     d'arrivée.
     Avantage : rapide
     Inconvénient : faux, car une petite pop subira autant d'événements de
     migration qu'une grande pop.
  */
  
  Update_PopulationNumber();
  Update_IndividualNumber();
  unsigned int nb_migr = MigrantsNumber(nb_ind.effectif);

  // Une seule population -> On saute la boucle.
  if (nb_pop.effectif < 2) return;

  for (unsigned int m = 0; m < nb_migr; m++)
  {
    unsigned int dep;
    unsigned int arr;

    Update_IndividualNumber();
   
    /* Pour que cette boucle ne bloque pas, il faut au moins un individu dans
       l'espèce. Il faudrait tester ça!!! */
    do {
      dep = generator->Uniform_int(nb_pop.effectif);
    } while (ref_pop[dep]->Size() == 0);
    // Si le container est une liste, ça coince ici!!!
    do {
      arr = generator->Uniform_int(nb_pop.effectif);
    } while (dep == arr);  // Il faut que la pop d'arrivée soit différente!

    Population * a (ref_pop[arr]);
    Population * d (ref_pop[dep]);

    // a->Ajoute_individu(d->Retire_individu());
		// Oui mais non : il faut connaitre le sexe de l'individu migrant!
    a->Update();
    d->Update();
  }
}

void Species::MigrationByInd(void)
{
  /* Migration ind par ind
     Chaque individu a la même probabilité de migrer. Le tirage du nombre
     de migrants par population dépend de l'effectif de chaque population.
     Avantage : Réaliste (?)
     Inconvénient : Lent.
  */

  if (nb_pop.effectif < 2) return;

  unsigned int nb_migr;
  unsigned int arr;
  it_pop e (ref_pop.end());
  unsigned int p = 0;

  for (it_pop b = ref_pop.begin(); b != e; b++)
  {
    // On maintient une double boucle (sur p et sur b) : pas beau.
    nb_migr = MigrantsNumber((*b)->Size());
    for (unsigned int m = 0; m < nb_migr; m++)
    {
      // tirage de la population d'arrivée
      do {
        arr = generator->Uniform_int(nb_pop.effectif);
      } while (arr == p);
      // ref_pop[arr]->Ajoute_individu((*b)->Retire_individu());
			// Même problème : sexe du migrant inconnu.
    }
    p++;
  }
}

unsigned int Species::MigrantsNumber(unsigned long int p) const
/* Détermine le nombre d'individus migrants parmi p individus
*/
{  
  double esper_mig = _tx * (static_cast<double>(p));

  return generator->Poisson(esper_mig);
}

Population * Species::ChooseRandomPopByTE(void) const
{
	assert (TEnumber() > 0);
	
	if (PopulationNumber() == 1)
		return *(ref_pop.begin());
	
	unsigned int tmp = 0;
	unsigned int and_the_winner_is = generator->Uniform_int(TEnumber());
	c_it_pop pop, popend;
	popend = ref_pop.end();
	for (pop = ref_pop.begin(); pop != popend; ++pop)
	{
		tmp += (*pop)->TEnumber();
		if (and_the_winner_is < tmp)
			return (*pop);
	}
	return NULL;
}
/*****************************************************************************/
/*    Construction, destruction, nettoyage, recopie                        */
/*****************************************************************************/

Species::Species(const Species & sp)
{
  Copy(sp);
  UnUpdate();
}

Species::Species(void)
{	
	assert (parameter != NULL);
	assert (generator != NULL);
	assert (genome != NULL);
  
  UnUpdate();
	my_system = NULL;
}

Species::Species(System * sys)
{
	assert (sys != NULL);
	
	assert (parameter != NULL);
	assert (generator != NULL);
	assert (genome != NULL);
	
	UnUpdate();
	
	my_system = sys;
	my_te_number = 0;
}


Species &Species::operator = (const Species &species)
{
  // Les variables _copenviro et r ont été initialisées lors de la construction
  Clean();
  Copy(species);
  return *this;
}


Species::~Species(void)
{
  Clean();
}


void Species::Clean(void)
{
  // Destruction de toutes les populations
  it_pop e = ref_pop.end();
  for (it_pop b = ref_pop.begin(); b != e; b++)
  {
    delete *b;
  }
  ref_pop.clear();
  UnUpdate();
}

void Species::Copy(const Species & modele)
{  
  Population * inter;
  
  // Recopie les variables
  nb_pop = modele.nb_pop;
  nb_ind = modele.nb_ind;
  nb_et = modele.nb_et;
  
  _tx = modele._tx;

  // Recopie les vecteurs
  c_it_pop e = modele.ref_pop.end();
  for (c_it_pop b = modele.ref_pop.begin(); b != e; b++)
  {
    // inter = (*b)->Virtual_new_Population();
		inter = Population::Copy_Constructor(**b);
    ref_pop.push_back(inter);
  }
  UnUpdate(); // Normalement : non nécesaire...
}

void Species::TE_Add(unsigned int count) // default = 1
{
	my_te_number += count;
	my_system->TE_Add(count);
}

void Species::TE_Remove(unsigned int count) // default = 1
{
	assert (count <= my_te_number);
	my_te_number -= count;
	my_system->TE_Remove(count);
}
