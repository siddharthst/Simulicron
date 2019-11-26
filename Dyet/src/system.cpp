/***************************************************************************
                          systeme.cpp  -  description
                             -------------------
    begin                : Thu Dec 26 2002
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

#include "system.h"
#include "species.h"

#include <algorithm>
#include <set>

 using namespace std;
 
//********************** OBJET SYSTEME *************************************

// public:

/************************ Constructeurs et destructeurs ********************/

System::System(void)
{
  // Initialisation du système
  // Pas grand chose à faire.
	
  UnUpdate();
  nb_et.effectif = 0;
  nb_esp.effectif = 0;
  nb_pop.effectif = 0;
	
	assert (generator != NULL);
	assert (genome != NULL);
	assert (parameter != NULL); 
}

System::System(const System & sys)
{
  Copy(sys);
  UnUpdate();
}


System & System::operator = (const System & sys)
{	
  if (this != &sys)  // On évite ainsi l'auto-affectation
  {
    // Il faut nettoyer complètement l'objet avant de le réemplir.
    Clean();
    Copy(sys);
    UnUpdate();
	}
  return *this;
}

System::~System(void)
{
  Clean();
}

/************************ Méthodes d'évolution du système *****************/

System * System::Reproduction(void)
{
  /* Assure la reproduction du système. */

  System * sys = new System();
	Update();
	
  it_sp debut = ref_esp.begin();
  it_sp fin = ref_esp.end();

	while(debut != fin)
	{
		try {
		Species * spe (NULL);
		spe = (*debut)->Reproduction();
    sys->AddSpecies(spe);
		++debut;
		}	catch (exec_nomore_population & e) {
			delete *debut;
			debut = ref_esp.erase(debut);	// equivalent to ++debut
			fin = ref_esp.end(); // end() move with erase(). 
		} catch (exec_toomuch & e) {
			delete sys;
			throw exec_toomuch();
		}
	}
  sys->Update();
	if (sys->SpeciesNumber() == 0)
	{
		delete sys;
		sys = NULL;
		throw exec_nomore_species();
	}
	return sys;
}

void System::Migration(void)
{
  /* Assure la migration entre toutes les populations de toutes les espèces
     du système. */

  c_it_sp debut, fin;

  fin = ref_esp.end();

  for (debut = ref_esp.begin(); debut != fin; debut++)
  {
    (*debut)->Migration();
  }
}

void System::HTransfert(void)
{
	
}

Species * System::ChooseRandomSpecies(void) const
{
	if (TEnumber() == 0)
		throw exec_noelements();
	
	if (SpeciesNumber() == 1)
		return *ref_esp.begin();
	
	Species * response = NULL;
	
	c_it_sp spe;
	c_it_sp endspe = ref_esp.end();
	do {		
		unsigned int tmp = 0;
		unsigned int and_the_winner_is = generator->Uniform_int(SpeciesNumber());
		for (spe = ref_esp.begin() ; spe != endspe; ++spe)
		{
			tmp ++;
			if (and_the_winner_is < tmp)
			{
				response = (*spe);
				break;
			}
		}
	} while (response->TEnumber() == 0);
	return response;
}

Species * System::ChooseRandomSpeciesByTE(void) const
{
	if (SpeciesNumber() == 1)
		return *ref_esp.begin();
	
	c_it_sp spe;
	c_it_sp endspe = ref_esp.end();
	unsigned int tmp = 0;
	unsigned int and_the_winner_is = generator->Uniform_int(TEnumber());
	for (spe = ref_esp.begin() ; spe != endspe; ++spe)
	{
		tmp += (*spe)->TEnumber();
		if (and_the_winner_is < tmp)
			return (*spe);
	}
	return NULL;
}

/*************************** Echantillonnage - mesures *********************/

void System::AutoCheckUp(void)
{
		if (SpeciesNumber() == 0)
		{
			cerr << "Nomore_species detected in System::AutoCheckUp()." << endl;
			throw exec_nomore_species(); // Extinction before TE cleaning.
			// Theoretically, this exception must be thrown during Reproduction.
		}
	
	// throws exceptions if the system is in bad conditions.
		
	if (TEnumber() > static_cast<unsigned int>(MAX_TE_LIMIT))
		throw exec_toomuch();

	
	// Destruction of species without TEs
	/* it_sp spe = ref_esp.begin();
	while (spe != ref_esp.end())
	{
		if ((*spe)->TEnumber() == 0)
		{
			delete *spe;
			spe = ref_esp.erase(spe);
		} else {
			++spe;
		}
	} */
		
	UnUpdate();
	Update();

	if (ref_esp.begin()== ref_esp.end())
		throw exec_noelements(); // Now, if the system dies, it is because of element loss.
}

#define UPDATE_SET 1

vector<Ptr2Element> System::UpdateMemory(void) const
{
	  /* 	Returns a vector containing all existing elements.
		Elements are counted only once. */
	vector<Ptr2Element> memory;
	#if UPDATE_SET
	set<Ptr2Element, IsIdBefore> set_memory;
	#endif
	
	System::c_it_sp esp, esp_end;
	Species::c_it_pop pop, pop_end;
	IndivGroup::c_it_ind ind, ind_end;
	Gamete::c_it_elem elem, elem_end;

	esp_end = ref_esp.end();
	for (esp = ref_esp.begin(); esp != esp_end; esp++)
	{
		pop_end = (*esp)->ref_pop.end();
		for (pop = (*esp)->ref_pop.begin(); pop != pop_end; pop++)
		{			
			const IndivGroup * group;
			const Gamete * gam;
			for (int sex = 0; sex < 2; sex++)
			{
				if (sex == 0) 
					group = (*pop)->males;
				else 
					group = (*pop)->females;
				
				ind_end = group->ref_ind.end();
				for (ind = group->ref_ind.begin(); ind != ind_end; ind++)
				{
					for (int g = 0; g < 2; g++)
					{
						if (g == 0) gam = (*ind)->gam1; else gam = (*ind)->gam2;
						elem_end = gam->ref_et.end();
						for (elem = gam->ref_et.begin(); elem != elem_end; elem++)
						{
							#if UPDATE_SET
							set_memory.insert(*elem); // not inserted if *elem already exists
							#else
							if (find(memory.begin(), memory.end(), *elem) == memory.end())
							{
								memory.push_back((*elem));
							}
							#endif
						}
					}
				}
			}
		}
	}
	#if UPDATE_SET
	set<Ptr2Element, IsIdBefore>::iterator set_it;
	for ( set_it = set_memory.begin(); set_it != set_memory.end(); set_it++)
	{
		memory.push_back(*set_it);
	}
	#endif
	return memory;
}

/******************************************************************************/
/*             Lecture - Ecriture (08/04/03)                                  */
/******************************************************************************/


System * System::Load(const string nomdufic)
{
  // Charge un systeme à partir du fichier nomdufic

  ImportSystem parser;
	parser.SetSilence(); 
	// No need of output here ; warnings will be displayed during the real import.
	
  // cout2 << "Lecture de la population initiale dans le fichier '";
  // cout2 << nomdufic << "'" << endl;  

	parser.parse_file(nomdufic); // throws a lot of things...

	return parser.GetSystem();
}

void System::Save(const string nomdufic) // const 
{
  /* Enregistre la totalité du système dans le fichier "nomdufic" au format XML */

  // Il faut que le système soit bien à jour.
	Update();
  assert (IsUpToDate());
	ExportSystem exportation(nomdufic, this);
  exportation.Save();
} 



// private:


void System::AddSpecies(Species * esp)
{

  ref_esp.push_back(esp);
  UnUpdate();
}

void System::RemoveSpecies(it_sp toberemoved)
{
	assert (!ref_esp.empty());
	assert(toberemoved != ref_esp.end());
	
	delete *toberemoved;
	ref_esp.erase(toberemoved);
}

/******************************************************************************/
/*   Recopies et nettoyage de mémoire  (07/01/03)                             */
/******************************************************************************/


void System::Copy(const System & modele)
{
  /* Recopie intégralement le modele dans le système courant.
     Utilisé par le constructeur de copie et l'opérateur d'affectation.
  */
  
  // recopie des variables
  nb_pop = modele.nb_pop;
  nb_esp = modele.nb_esp;
  nb_et = modele.nb_et;

  // Recopie du tableau d'especes
  c_it_sp d, f;
  f = modele.ref_esp.end();
  for (d = modele.ref_esp.begin(); d != f; d++)
  {
    ref_esp.push_back(new Species(**d));
  }
  UnUpdate(); 
	Update();
}

void System::Clean(void)
{
  /* Détruit complètement le système.
     Méthode utilisée par : - le destructeur
                            - l'opérateur d'affectation
  */
  it_sp d, f;
  f = ref_esp.end();
  
	// cout << "destruction systeme" << endl;
  for (d = ref_esp.begin(); d != f; d++)
  {
    delete *d;
  }
  ref_esp.clear(); 
  UnUpdate();
}
