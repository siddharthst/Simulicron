/***************************************************************************
                          individu.inl  -  description
                             -------------------
    begin                :  Fri July 25 2003
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
 
 /* inline functions included in individu.h */
 

/*********************** méthodes inline ******************************/
bool Individual::IsBurst(void) const
{
	return burst;
}

double Individual::Fitness(void)
 { 
		Update_Fitness(); 
	 return fitness.compte;
}

double Individual::Fitness(void) const
{ 
	assert(fitness.mise_a_jour); 
	return fitness.compte;
}

unsigned int Individual::TEnumber(void)
{ 
	Update_TEnumber(); 
	return nb_et.compte;
}

unsigned int Individual::TEnumber(void) const
{ 
	assert(nb_et.mise_a_jour); 
	return nb_et.compte;
}

unsigned int Individual::DuplicationNumber(void) const
{
	return duplication_number;
}

unsigned int Individual::PotentialDuplicationNumber(void) const
{
	return potential_duplication_number;
}

double Individual::sum_duplication_rates(void) const
{
	return gam1->sum_duplication_rates() + gam2->sum_duplication_rates();
}

double Individual::sum_deletion_rates(void) const
{
	return gam1->sum_deletion_rates() + gam2->sum_deletion_rates();
}

double Individual::sum_multiplicative_burst_factors(void) const
{
	return gam1->sum_multiplicative_burst_factors() + gam2->sum_multiplicative_burst_factors();
}

double Individual::sum_additive_burst_factors(void) const
{
	return gam1->sum_additive_burst_factors() + gam2->sum_additive_burst_factors();
}

double Individual::sum_selection_burst_factors(void) const
{
	return gam1->sum_selection_burst_factors() + gam2->sum_multiplicative_burst_factors();
}

void Individual::Update(void)
{
  Update_TEnumber();
  Update_Fitness();
  Update_Activity();
}

void Individual::UnUpdate(void)
{
  VERB;
  // Méthode appelée à chaque fois que le génome de l'individu est modifié.  
    nb_et.mise_a_jour = false;
    fitness.mise_a_jour = false;
    activite.mise_a_jour = false;
}

const char * Individual::e_indiv::what(void) const throw()
{
  return "Erreur générique sur l'objet 'Individu'.";
}

const char * Individual::e_bad_act::what(void) const throw()
{
  return "Erreur de code : activité négative, impossible.";
}
