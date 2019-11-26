/***************************************************************************
                          systeme.inl  -  description
                             -------------------
    begin                : Mon Aug 25 2003
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
 
#include "species.h" 
#include "indivgroup.h"
#include "individual.h"
#include "gamete.h"


/* inline functions for class System (system.h) */ 

bool System::IsUpToDate(void) const
{
	return (nb_esp.mise_a_jour && nb_et.mise_a_jour);
}

unsigned int System::SpeciesNumber(void)
{
  Update_SpeciesNumber();
  return nb_esp.effectif;
}

unsigned int System::SpeciesNumber(void) const
{
  if (nb_esp.mise_a_jour != true)
  {
    std::cerr << "Tentative d'acces au nombre d'especes du systeme alors qu'il n'a"
         << "pas ete mis a jour." << std::endl;
  }
  return nb_esp.effectif;
}

long unsigned int System::TEnumber(void)
{
	Update_TEnumber();
  return nb_et.effectif;
}

long unsigned int System::TEnumber(void) const
{
  if (nb_et.mise_a_jour != true)
  {
    std::cerr << "Tentative d'acces au nombre d'elements du systeme alors qu'il n'a"
         << "pas ete mis a jour." << std::endl;
  }
  return nb_et.effectif;
}

long unsigned int System::IndividualNumber(void) const
{
	unsigned int nbind = 0;
	c_it_sp spe, endspe;
	endspe = ref_esp.end();
	for (spe = ref_esp.begin(); spe != endspe; spe++)
		nbind += (*spe)->IndividualNumber();
	return nbind;
}

void System::Update(void)
{
  // Mise a jour complete du systeme.
  it_sp d, f;
  f = ref_esp.end();
  for (d = ref_esp.begin(); d != f; d++)
  {
    (*d)->Update();
  }
  Update_TEnumber();
  Update_SpeciesNumber();
  Update_PopulationNumber();
}

void System::UnUpdate(void)
{
  nb_pop.mise_a_jour = false;
  nb_esp.mise_a_jour = false;
  nb_et.mise_a_jour = false;
}


void System::Update_PopulationNumber(void)
{
  if (nb_pop.mise_a_jour == false)
  {
    nb_pop.effectif = 0;
    c_it_sp b, e;
    e = ref_esp.end();
    for (b = ref_esp.begin(); b != e; b++)
      nb_pop.effectif += (*b)->PopulationNumber();
    nb_pop.mise_a_jour = true;
  }
}

void System::Update_SpeciesNumber(void)
{
  if (nb_esp.mise_a_jour == false)
  {
    nb_esp.effectif = ref_esp.size();
    nb_esp.mise_a_jour = true;
  }
}

void System::Update_TEnumber(void)
{
  if (nb_et.mise_a_jour == false)
  {
    nb_et.effectif = 0;
    c_it_sp b, e;
    e = ref_esp.end();
    for (b = ref_esp.begin(); b != e; b++)
      nb_et.effectif += (*b)->TEnumber();
    nb_et.mise_a_jour = true;
  }
}

// HTransferts
/* 
void System::HTransfert_base::Do(double htrate)
{
	unsigned int ht_number = HTNumber(htrate);
	
	for (unsigned int ht = 0; ht < ht_number; ht++)
	{
		Species * destination_species;
		Species * origin_species = this->ChooseOriginSpecies();
		Ptr2Element jumping_te = origin_species->PickRandomElement(); // equiv to a transposition
		
		destination_species = origin_species->VoidCopy();
		destination_species->UnUpdate();
		destination_species->Update();
		destination_species->AddElement(jumping_te);
		destination_species->UnUpdate();
		destination_species->Update();
		p_system->AddSpecies(destination_species);
		p_system->Update();
	}
}

unsigned int System::HTransfert_TENumber::HTNumber(double htrate) const
{ 
	return p_system->generator->Poisson(static_cast<double>(p_system->TEnumber()) * htrate);
}

Species * System::HTransfert_TENumber::ChooseOriginSpecies(void) const
{ 
	return p_system->ChooseRandomSpeciesByTE();
}

unsigned int System::HTransfert_SpNumber::HTNumber(double htrate) const
{
	return p_system->generator->Poisson(static_cast<double>(p_system->SpeciesNumber()) * htrate);
}

Species * System::HTransfert_SpNumber::ChooseOriginSpecies(void) const
{ 
	return p_system->ChooseRandomSpecies();
}

unsigned int System::HTransfert_Const::HTNumber(double htrate) const
{
	return p_system->generator->Poisson(htrate);
}

Species * System::HTransfert_Const::ChooseOriginSpecies(void) const
{
	return p_system->ChooseRandomSpecies();
}

*/

void System::TE_Add(unsigned int count) // default = 1
{
	my_te_number += count;
}

void System::TE_Remove(unsigned int count) // default = 1
{
	assert (count <= my_te_number);
	my_te_number -= count;
}
