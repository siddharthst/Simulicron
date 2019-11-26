/***************************************************************************
                          espece.inl  -  description
                             -------------------
    begin                : Mon 11 July 2003
    copyright            : (C) 2003 by Arnaud Le Rouzic
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
 
 // Inline functions, included by species.h
 
#include "population.h"
#include "system.h"
 
 /********************* méthodes inline *********************************/
void Species::Update(void)
{
  it_pop p, pf;
  pf = ref_pop.end();
  for (p = ref_pop.begin(); p != pf; p++)
  {
    (*p)->Update();
  }
  Update_PopulationNumber();
  Update_IndividualNumber();
  Update_TEnumber();
}

unsigned int Species::PopulationNumber(void)
{
  Update_PopulationNumber();
  return nb_pop.effectif;
}

unsigned int Species::PopulationNumber(void) const
{
  if (nb_pop.mise_a_jour == false)
  {
    std::cerr << "Le nombre de populations dans l'espece n'a pas ete mis a jour." << std::endl;
    std::cerr << " Si ca plante, c'est normal..." << std::endl;
  }
  return nb_pop.effectif;
}

unsigned int Species::IndividualNumber(void)
{
  Update_IndividualNumber();
  return nb_ind.effectif;
}

unsigned int Species::TEnumber(void)
{
  Update_TEnumber();
  return nb_et.effectif;
}

unsigned int Species::TEnumber(void) const
{
	  if (nb_et.mise_a_jour == false)
  {
    std::cerr << "Le nombre l'elements dans l'espece n'a pas ete mis a jour." << std::endl;
    std::cerr << " Si ca plante, c'est normal..." << std::endl;
  }
  return nb_et.effectif;
}

double Species::AverageFitness(void) const
{
	double compte = 0.0;
	c_it_pop pop, endpop;
	endpop = ref_pop.end();
	
	for (pop = ref_pop.begin(); pop != endpop; pop++)
		compte += (*pop)->AverageFitness();

	return compte/(static_cast<double>(PopulationNumber()));
}


void Species::UnUpdate(void)
{
    nb_ind.mise_a_jour = false;
    nb_pop.mise_a_jour = false;
    nb_et.mise_a_jour = false;
}

void Species::Update_PopulationNumber(void)
{
  if (!nb_pop.mise_a_jour)
  {    
    nb_pop.effectif = ref_pop.size();
    nb_pop.mise_a_jour = true;
  }
}

void Species::Update_IndividualNumber(void)
{
  // calcule le nombre d'individus dans l'espèce. Ne pas abuser de cette
  // fonction, le calcul est long.
  if (nb_ind.mise_a_jour == false)
  {
    nb_ind.effectif = 0;
    c_it_pop fin (ref_pop.end());
    for (c_it_pop deb = ref_pop.begin(); deb != fin; deb++)
      nb_ind.effectif += (*deb)->Size();
    nb_ind.mise_a_jour = true;
  }
}

void Species::Update_TEnumber(void)
{
  // Calcule le nombre d'éléments dans l'espèce.
  if (!nb_et.mise_a_jour)
  {
    nb_et.effectif = 0;
    c_it_pop fin (ref_pop.end());
    for (c_it_pop deb = ref_pop.begin(); deb != fin; deb++)
      nb_et.effectif += (*deb)->TEnumber();
    nb_et.mise_a_jour = true;
  }
}
