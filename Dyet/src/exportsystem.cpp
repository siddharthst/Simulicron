/***************************************************************************
                          exportsystem.cpp  -  description
                             -------------------
    begin                : Wed Feb 12 2003
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

#include "exportsystem.h"

#include <sstream>
#include "genome.h"
#include "system.h"

using namespace std;

ExportSystem::ExportSystem(const System * syst)
{
	if (syst != NULL)
		sys = syst;
	else
		_e_systnull();
}

ExportSystem::ExportSystem(void)
{
}

ExportSystem::ExportSystem(const string nomfic)
{
  if (Ouv_fic(nomfic) == false)
    _e_ouvfic(nomfic);
}

ExportSystem::ExportSystem(const string nomfic, const System * syst)
{
  if (Ouv_fic(nomfic) == false)
    _e_ouvfic(nomfic);

  if (syst != NULL)
    sys = syst;
  else
    _e_systnull();
}

bool ExportSystem::Ouv_fic(const string nomdufic)
{
  assert (!ficsx.is_open()); // Impossible d'ouvrir un fichier déja ouvert

  ficsx.open(nomdufic.c_str());
  if (ficsx.is_open()) return true;
  else return false;
}

void ExportSystem::Save(void)
{
  /* The system must be updated. It is not possible to update here because sys is const
    sys->Maj_nb_esp();
    sys->Maj_nb_et();
    sys->Maj_nb_pop();
  */
  
  Print_header();
  Print_genome();
  Print_memory();
  Print_systeme(sys); // sys is not necessary, since it is a variable of class ExportSystem
  Print_ccl();

  ficsx.close();
}

string ExportSystem::DTD(void)
{
  ostringstream o;

  o << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>" << endl;
  // Il faut prendre en compte le chemin du fichier pop.dtd
  // string path_dtd (CURR_PATH + "dtd/pop.dtd");
  string path_dtd = ""; // no dtd for the moment.
	// o << "<!DOCTYPE dyet SYSTEM \"" << path_dtd << "\">" << endl;
	o << "<!DOCTYPE dyet>" << endl;

  return o.str();
}                

void ExportSystem::Print_header(void)
{  
  ficsx << ExportSystem::DTD();
  ficsx << "<dyet>" << endl;
}

void ExportSystem::Print_ccl(void)
{
  ficsx << "</dyet>" << endl; // Tout pour l'instant.
}

void ExportSystem::Print_genome(void)
{
  // Ecriture du génome commun au système
  const Genome * g = BaseBiolObject::GetGenome();

  assert (!g->EstVide());
  unsigned int nbc = g->NbChromo();

  string indent1(" ");
  string indent2("  ");
  ficsx << indent1 << "<genome>" << endl;
  for (unsigned int c = 0; c < nbc; c++)
    ficsx << indent2 << "<chromo>" << g->ChromoSz(c) << "</chromo>" << endl;
  ficsx << indent1 << "</genome>" << endl;
}
  

void ExportSystem::Print_memory(void)
{
  UpdateMemory();

  string indent(" ");

  ficsx << indent << "<memory_table>" << endl;
  for (unsigned int i = 0; i < memory.size(); i++)
  {
    Print_element(memory[i]);
  }
  ficsx << indent << "</memory_table>" << endl;
}
    
void ExportSystem::Print_systeme(const System * syst)
{
  assert (syst != NULL);
  string indent (" ");
  ficsx << indent << "<system>" << endl;
	System::c_it_sp spe, endspe;
	endspe = syst->ref_esp.end();
	for (spe = syst->ref_esp.begin(); spe != endspe; spe++)
	{
		Print_espece(*spe);
	}
  ficsx << indent << "</system>" << endl;
}

void ExportSystem::Print_espece(const Species * esp)
{
  assert (esp != NULL);
  string indent ("  ");
  ficsx << indent << "<species>" << endl;
  for (unsigned int p = 0; p < esp->PopulationNumber(); p++)
  {
    Print_population(esp->ref_pop[p]);
  }
  ficsx << indent << "</species>" << endl;
}

void ExportSystem::Print_population(const Population * pop)
{
  assert (pop != NULL);
  string indent("   ");
  ficsx << indent << "<population>" << endl;
	// Hermaphrodite pop
	if (pop->males == pop->females)
	{
		unsigned int nbind = pop->Size();
		for (unsigned int i = 0; i < nbind; i++)
			Print_individu(pop->males->ref_ind[i], "individu");
	} else {
		// Sexual pop
		unsigned int nbmal = pop->males->Size();
		unsigned int nbfem = pop->females->Size();
		for (unsigned int i = 0; i < nbmal; i++)
			Print_individu(pop->males->ref_ind[i], "male");
		for (unsigned int i = 0; i < nbfem; i++)
			Print_individu(pop->females->ref_ind[i], "female");
	}
  ficsx << indent << "</population>" << endl;
}

void ExportSystem::Print_individu(const Individual * ind, const std::string & name)
{
  assert (ind != NULL);
	assert ((name == "individu") || (name == "male") || (name == "female"));
  string indent ("    ");
  ficsx << indent << "<" << name << ">" << endl;
  Print_gamete(ind->gam1);
  Print_gamete(ind->gam2);
  ficsx << indent << "</" << name << ">" << endl;
}

void ExportSystem::Print_gamete(const Gamete * gam)
{
  Gamete::c_it_elem elem, elem_end;
  string indent ("     ");  
  
  assert (gam != NULL);
  ficsx << indent << "<gamete>";
  elem_end = gam->ref_et.end();
  for (elem = gam->ref_et.begin(); elem != elem_end; elem++)
  {
    Print_ref_element(*elem);
  }
  ficsx << "</gamete>" << endl;
}

void ExportSystem::Print_ref_element(const Ptr2Element te)
{
  ficsx << "<te>" << te->id() << "</te>";
}

void ExportSystem::Print_element(const Ptr2Element te)
{
	string indent1 ("  ");
	string indent2 ("   ");

	ficsx << indent1 << "<element id = \"" << te->id() << "\">" << endl;

	Print_ElementFeature(te->fixed_features._family, indent2);
	Print_ElementFeature(te->fixed_features._mean_fitness, indent2);
	Print_ElementFeature(te->fixed_features._sd_fitness, indent2);
	Print_ElementFeature(te->fixed_features._dup_mut, indent2);
	Print_ElementFeature(te->fixed_features._dup_sd, indent2);
	Print_ElementFeature(te->fixed_features._del_mut, indent2);
	Print_ElementFeature(te->fixed_features._del_sd, indent2);
	Print_ElementFeature(te->fixed_features._act_mut, indent2);
	Print_ElementFeature(te->fixed_features._act_sd, indent2);
	Print_ElementFeature(te->fixed_features._reg_mut, indent2);
	Print_ElementFeature(te->fixed_features._reg_sd, indent2);

	Print_ElementFeature(te->variable_features._position, indent2);
	Print_ElementFeature(te->variable_features._chromo, indent2);
	Print_ElementFeature(te->variable_features._fitness, indent2);
	Print_ElementFeature(te->variable_features._duplication_rate, indent2);
	Print_ElementFeature(te->variable_features._deletion_rate, indent2);
	Print_ElementFeature(te->variable_features._activity, indent2);
	Print_ElementFeature(te->variable_features._regulation_factor, indent2);

	ficsx << indent1 << "</element>" << endl;
}

template <typename T> void ExportSystem::Print_ElementFeature(const ElementFeature<T> & feat, const string & indent)
{
	ficsx << indent << "<" << feat.xml_tag << ">" << feat.value << "</" << feat.xml_tag << ">" << endl;
}

void ExportSystem::UpdateMemory(void)
{
  memory = sys->UpdateMemory();
}

void ExportSystem::_e_ouvfic(const std::string & filnam = "unknown") const
{
	// should be an exception!!!
  cerr << "Erreur lors de la sauvegarde du systeme." << endl;
	cerr << "Impossible de creer le fichier '" << filnam << "' sur le disque" << endl;
  exit(1);
}

void ExportSystem::_e_systnull(void) const
{
  cerr << "Erreur de programmation : le système à sauvegarder n'est pas valide" << endl;
  exit(1);
}
