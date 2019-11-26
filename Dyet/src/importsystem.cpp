/***************************************************************************
                          importsystem.cpp  -  description
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

#include <cassert>

#include "importsystem.h"
#include "system.h"
#include "genome.h"

using namespace std;

bool ImportSystem::blabla = false;


System * ImportSystem::GetSystem(void)
{
	// retourne un pointeur vers le système lu
	assert (sys != NULL);
	sys->Update();
	return sys;
}  

/******************* Surcharge des méthodes de xmlpp::SaxParser :
   on_start_element
   on_end_element
   on_characters
*/

void ImportSystem::on_start_element (const string & bal, const AttributeMap & atmap)
{
	if (blabla) cout << "Beginning tag : " << bal << flush;
  // début d'une balise
  if (_in_memory_table){
    // Traite_in_memory_start(bal, _CurrentAttributes);
    Traite_in_memory_start(bal, atmap);
    return;
  }
  if (_in_system){
    Traite_in_system_start(bal, atmap);
    return;
  }
  if (_in_genome){
    Traite_in_genome_start(bal);
    return;
  }
    
  // pas de switch sur une string
  if (bal == "memory_table"){
    _in_memory_table = true;
  } else {
  if (bal == "system"){
		// Genome has been initialized in genome_start.
    sys = new System();
    _in_system = true;
  } else {
  if (bal == "genome"){
			genome = new Genome();
    _in_genome = true;
  } else {
  if (bal == "dyet"){
    // Rien de particulier, cette balise est plutot rassurante!
    if (blabla) cout << "Lecture xml" << flush;
  } else {
    _e_xml(bal);
  }}}}
  return;  
}

void ImportSystem::on_end_element (const string & bal)
{
	if (blabla) 
		cout << " ending tag " << bal << endl;
  // Fin de la balise "bal"
  if (bal == "dyet") {
    if (blabla) cout << " fin" << endl;
    return;
  }
  if (_in_genome){
    Traite_in_genome_end(bal);
  } else {
  if (_in_memory_table){
    Traite_in_memory_end(bal);
  } else {
  if (_in_system){
    Traite_in_system_end(bal);
  } else {
    _e_xml(bal);
  }}}  
}

void ImportSystem::on_characters (const string & characters)
{
	if (blabla)
		cout << " characters : " << characters << flush;
  // Lecture des donnees
  istringstream data(characters);
	
  if (_in_genome){
    if (_in_chromo){
      double f;
      data >> f;
      genome->Nouveau_chromo(f);
    }
  }      
  	if (_in_memory_table){
    		if (_in_family){
      			unsigned int f;
      			data >> f;
      			(*ptrel)->Set_family(f);
    		}
    		if (_in_locus){
      			double l;
      			data >> l;
      			(*ptrel)->Set_position(l);
      			// Chromosome setting is done here
    		}
    		if (_in_chromosome){
      			int c;
      			data >> c;
      			if (((*ptrel)->chromo()) != c)
		  	{
				  cerr 	<< "Warning!! Element id=" << (*ptrel)->id() 
						<< " at locus " << (*ptrel)->position()
						<< " is said to be on chromosome " << c << "." << endl;
				  cerr	<< "It cannot. The real chromosome will be set to " << (*ptrel)->chromo() << "." << endl;
		  	}
      			// no affectation here; it is just a verification.
		  	// chromosome should be a simple xml attribute.
    		}
   		if (_in_fitness){
      			double f;
      			data >> f;
      			(*ptrel)->Set_fitness(f);
    		}
		if (_in_mean_fitness){
			double mf;
			data >> mf;
			(*ptrel)->Set_mean_fitness(mf);
		}
		if (_in_sd_fitness){
			double vf;
			data >> vf;
			(*ptrel)->Set_sd_fitness(vf);
		}
		if (_in_dup_rate){
			double u;
			data >> u;
			(*ptrel)->Set_dup_rate(u);
		}
		if (_in_del_rate){
			double v;
			data >> v;
			(*ptrel)->Set_del_rate(v);
		}
		if (_in_activity){
			double q;
			data >> q;
			(*ptrel)->Set_activity(q);
    		}
		if (_in_dup_mut){
			double f;
			data >> f;
			(*ptrel)->Set_mut_dup_rate(f);
		}
		if (_in_del_mut){
			double f;
			data >> f;
			(*ptrel)->Set_mut_del_rate(f);
		}
		if (_in_act_mut){
			double f;
			data >> f;
			(*ptrel)->Set_mut_activity(f);
		}
		if (_in_dup_sd){
			double r;
			data >> r;
			(*ptrel)->Set_sd_dup_rate(r);
		}
		if (_in_del_sd){
			double r;
			data >> r;
			(*ptrel)->Set_sd_del_rate(r);
		}
		if (_in_act_sd){
			double f;
			data >> f;
			(*ptrel)->Set_sd_activity(f);
		}
		if (_in_regulation){
			double f;
			data >> f;
			(*ptrel)->Set_regulation_factor(f);
		}
		if (_in_reg_mut){
			double f;
			data >> f;
			(*ptrel)->Set_mut_reg(f);
		}
		if (_in_reg_sd){
			double f;
			data >> f;
			(*ptrel)->Set_sd_reg(f);
		}
  }
  if (_in_system){
    if (_in_species){
      if (_in_population){
        if (_in_individual){
          if (_in_gamete){
            if (_in_te)
            {
		unsigned int ref;
		data >> ref;
		assert(!data.fail()); // might throw an exception
              // assert(ref <= refelem); Not needed.
		if (ref > refelem) // The warning is already displayed in ImportSystem::Traite_in_memory_start
		{ 
			/* cerr << "Warning: element " << ref << " should be numbered " << refelem << "." << endl; */ }
             		assert (memory.find(ref) != memory.end()); 
			ptrel = &(memory[ref]);
            	}
  } } } } }
}

/****************************** Traitement des <balise> ********************/

// void ImportSystem::Traite_in_memory_start(const string & bal, const XML_Char ** args)

void ImportSystem::Traite_in_memory_start(const string & bal, const AttributeMap & atmap)
{
	// Traitement du début d'une balise dans "memory_table"
	
	if (_in_element)
	{
		Traite_in_element_start(bal);
		return;
	}

	if (bal == "element")
	{
		AttributeMap::const_iterator it = atmap.find("id");
		istringstream i_reflue (it->second);
		long unsigned int reflue;
		i_reflue >> reflue;

		/*	Import system was previously designed for initial data sets :
			Elements numbered "O" to "n", indexed by 'refelem'.
			It must not be always the case. */
		// assert (reflue == refelem); No! 
		if ((reflue != refelem) && (!silence))
			cerr << "Warning: element " << reflue	<< " should be numbered "<< refelem << "." << endl;
		_in_element = true;
		try {
			Element * el = new Element();
			ptrel = new Ptr2Element(el);
		} catch (exception & e) {cerr << e.what() << endl;}

		std::pair<long int, Ptr2Element> mypair(reflue, *ptrel);
		assert(memory.insert(mypair).second);    
		assert ((*ptrel)->id() == refelem);
		refelem++;
  	} else {
    		_e_xml(bal);
  	}
  	return;
}

void ImportSystem::Traite_in_genome_start(const string & bal)
{
  // Traitement de l'intérieur de la balise <genome>
  if (bal == "chromo"){
    _in_chromo = true;
  } else {
    _e_xml(bal);
  }
}

void ImportSystem::Traite_in_element_start(const string & bal)
{
	// Traitement interne à la balise <element>
	assert (ptrel != NULL); // an element must exist.

  	if (bal == (*ptrel)->fixed_features._family.xml_tag) {
   		_in_family = true;
  	} else {
  	if (bal == (*ptrel)->variable_features._position.xml_tag) {
    		_in_locus = true;
	} else {
	if (bal == (*ptrel)->variable_features._chromo.xml_tag) {
    		_in_chromosome = true;
  	} else {
  	if (bal == (*ptrel)->variable_features._fitness.xml_tag) {
		_in_fitness = true;
	} else {
  	if (bal == (*ptrel)->variable_features._duplication_rate.xml_tag ) {
    		_in_dup_rate = true;	
	} else {
  	if (bal == (*ptrel)->variable_features._deletion_rate.xml_tag ) {
    		_in_del_rate = true;	
	} else {
  	if (bal == (*ptrel)->variable_features._activity.xml_tag ) {
    		_in_activity = true;	
	} else {
  	if (bal == (*ptrel)->variable_features._regulation_factor.xml_tag ) {
    		_in_regulation = true;	
	} else {
  	if (bal == (*ptrel)->fixed_features._reg_mut.xml_tag ) {
    		_in_reg_mut = true;	
	} else {
  	if (bal == (*ptrel)->fixed_features._reg_sd.xml_tag ) {
    		_in_reg_sd = true;	
	} else {
  	if (bal == (*ptrel)->fixed_features._mean_fitness.xml_tag) {
   		_in_mean_fitness = true;
  	} else {
  	if (bal == (*ptrel)->fixed_features._sd_fitness.xml_tag) {
   		_in_sd_fitness = true;
  	} else {
  	if (bal == (*ptrel)->fixed_features._dup_mut.xml_tag) {
   		_in_dup_mut = true;
  	} else {
  	if (bal == (*ptrel)->fixed_features._dup_sd.xml_tag) {
   		_in_dup_sd = true;
  	} else {
  	if (bal == (*ptrel)->fixed_features._del_mut.xml_tag) {
   		_in_del_mut = true;
  	} else {
  	if (bal == (*ptrel)->fixed_features._del_sd.xml_tag) {
   		_in_del_sd = true;
  	} else {
  	if (bal == (*ptrel)->fixed_features._act_mut.xml_tag){
   		_in_act_mut = true;
  	} else {
  	if (bal == (*ptrel)->fixed_features._act_sd.xml_tag) {
   		_in_act_sd = true;
	} else {
    		_e_xml(bal);
  	}}}}}}}}}}}}}}}}}}
  	return;
}

void ImportSystem::Traite_in_system_start(const string & bal, const AttributeMap & atmap)
{
  // Gère la rencontre d'une nouvelle balise dans <system>
  if (_in_species){
    Traite_in_species_start(bal, atmap);
    return;
  }

  if (bal == "species"){
    _in_species = true;
    spe = new Species();
  } else {
    _e_xml(bal);
  }
  return;
}

void ImportSystem::Traite_in_species_start(const string & bal, const AttributeMap & atmap)
{
  // Gère la rencontre d'une nouvelle balise à l'intérieur de <species>

  if (_in_population){
    Traite_in_population_start(bal, atmap);
    return;
  }

  if (bal == "population"){
    if (blabla) cout << "." << flush;
    _in_population = true;
    // pop= new Population();
		pop = NULL; 
		// Initialization cannot be made until population status (sexual or not)
		// is known.
  } else {
    _e_xml(bal);
  }
  return;
}

void ImportSystem::Traite_in_population_start(const string & bal, const AttributeMap & atmap)
{
  // Gère la rencontre d'une nouvelle balise à l'intérieur de <population>

  if (_in_individual){
    Traite_in_individual_start(bal);
    return;
  }

  if ((bal == "individual") || (bal == "male") || (bal == "female")){
		_in_individual = true;
		if (bal == "individual")
		{
    if (pop == NULL)
			pop = new Herm_Population();
		}
		if (bal == "male")
		{
			if (pop == NULL)
				pop = new Dioic_Population();
			_in_male = true;
		}
		if (bal == "female")
		{
			if (pop == NULL)
				pop = new Dioic_Population();
			_in_female = true;
		}
  gam_1ou2 = 0;
  ind = new Individual();
  }
  return;
}

void ImportSystem::Traite_in_individual_start(const string & bal)
{
  // Gère la rencontre d'une nouvelle balise à l'intérieur de <individual>

  if (_in_gamete){
    Traite_in_gamete_start(bal);
    return;
  }

  if (bal == "gamete"){
    _in_gamete = true;
    if (gam_1ou2 == 0){
      gam_1ou2 = 1;
    } else {
    if (gam_1ou2 == 1){
      gam_1ou2 = 2;
    } else {
      _e_gam();
    }}
  } else {
    _e_xml(bal);
  }
  return;
}

void ImportSystem::Traite_in_gamete_start(const string & bal)
{
  // Gère la rencontre d'une nouvelle balise à l'intérieur de <gamete>

  if (bal == "te") {
    _in_te = true;
  } else {
    _e_xml(bal);
  }
  return;
}

/****************************** Traitement des </balise> ******************/

void ImportSystem::Traite_in_memory_end(const string & bal)
{
  // gère la rencontre d'une balise de fin dans <memory_table>

  if (bal == "memory_table"){
    _in_memory_table = false;
    // Rien d'autre à faire, le tableau memory est rempli, c'est tout
    return;
  }
  if (bal == "element"){
    _in_element = false;
		ptrel = NULL;
    return;
  }
  
  assert (_in_element); // (sinon, il y a un problème...)
  Traite_in_element_end(bal);
  return;
}

void ImportSystem::Traite_in_genome_end(const string & bal)
{
  if (bal == "genome"){
    _in_genome = false;
		BaseBiolObject::SetGenome(genome);
		// no 'delete genome;' : genome is sent to BaseBiolObject.
  } else {
  if (bal == "chromo"){
    _in_chromo = false;
  } else {
    _e_xml(bal);
  }}
  return;
}  

void ImportSystem::Traite_in_element_end(const string & bal)
{
  // Gère la rencontre d'une balise de fin dans <element>

	assert (ptrel != NULL); // an element must exist.

  	if (bal == (*ptrel)->fixed_features._family.xml_tag) {
   		_in_family = false;
  	} else {
  	if (bal == (*ptrel)->variable_features._position.xml_tag) {
    		_in_locus = false;
	} else {
	if (bal == (*ptrel)->variable_features._chromo.xml_tag) {
    		_in_chromosome = false;
  	} else {
  	if (bal == (*ptrel)->variable_features._fitness.xml_tag) {
		_in_fitness = false;
	} else {
  	if (bal == (*ptrel)->variable_features._duplication_rate.xml_tag ) {
    		_in_dup_rate = false;	
	} else {
  	if (bal == (*ptrel)->variable_features._deletion_rate.xml_tag ) {
    		_in_del_rate = false;	
	} else {
  	if (bal == (*ptrel)->variable_features._activity.xml_tag ) {
    		_in_activity = false;	
	} else {
  	if (bal == (*ptrel)->variable_features._regulation_factor.xml_tag ) {
    		_in_regulation = false;	
	} else {
	if (bal == (*ptrel)->fixed_features._reg_mut.xml_tag ) {
		_in_reg_mut = false;
	} else {
	if (bal == (*ptrel)->fixed_features._reg_sd.xml_tag ) {
		_in_reg_sd = false;
	} else {
  	if (bal == (*ptrel)->fixed_features._mean_fitness.xml_tag) {
   		_in_mean_fitness = false;
  	} else {
  	if (bal == (*ptrel)->fixed_features._sd_fitness.xml_tag) {
   		_in_sd_fitness = false;
  	} else {
  	if (bal == (*ptrel)->fixed_features._dup_mut.xml_tag) {
   		_in_dup_mut = false;
  	} else {
  	if (bal == (*ptrel)->fixed_features._dup_sd.xml_tag) {
   		_in_dup_sd = false;
  	} else {
  	if (bal == (*ptrel)->fixed_features._del_mut.xml_tag) {
   		_in_del_mut = false;
  	} else {
  	if (bal == (*ptrel)->fixed_features._del_sd.xml_tag) {
   		_in_del_sd = false;
  	} else {
  	if (bal == (*ptrel)->fixed_features._act_mut.xml_tag){
   		_in_act_mut = false;
  	} else {
  	if (bal == (*ptrel)->fixed_features._act_sd.xml_tag) {
   		_in_act_sd = false;
	} else {
    		_e_xml(bal);
  	}}}}}}}}}}}}}}}}}}
  return;
}

void ImportSystem::Traite_in_system_end(const string & bal)
{
  // Gère la rencontre d'une balise de fin à l'intérieur de <systeme>

  if (bal == "system"){
    _in_system = false;
    // Qqch d'autre???
    return;
  }
  if (bal == "species"){
    _in_species = false;
    assert (spe != NULL);
    sys->AddSpecies(spe);
    return;
  }

  assert (_in_species);
  Traite_in_species_end(bal);
  return;
}

void ImportSystem::Traite_in_species_end(const string & bal)
{
  // Gère la rencontre d'une balise de fin à l'intérieur de <species>

  if (bal == "population"){
    _in_population = false;
    assert (pop != NULL);
    spe->AddPopulation(pop);
    return;
  }
  assert (_in_population);
  Traite_in_population_end(bal);
  return;
}

void ImportSystem::Traite_in_population_end(const string & bal)
{
  // Gère la rencontre d'une balise de fin à l'intérieur de <systeme>

  if ((bal == "individual") || (bal == "male") || (bal == "female")){
    assert (ind != NULL);  
		_in_individual = false;
		// Hermaphrodit
		if (!(_in_male || _in_female))
			pop->AddMale(ind); // Sex no matters

		// Male
		if (_in_male)
		{
			pop->AddMale(ind);
			_in_male = false;
		}
		// Female
		if (_in_female)
		{
			pop->AddFemale(ind);
			_in_female = false;
		}
    return;
  }
  assert (_in_individual);
  Traite_in_individual_end(bal);
  return;
}

void ImportSystem::Traite_in_individual_end(const string & bal)
{
  // Gère la rencontre d'une balise de fin à l'intérieur de <individual>

  if (bal == "gamete"){
    _in_gamete = false;
    // Rien d'autre ici : les éléments sont insérés au fur et à mesure.
    return;
  } else {
  if (bal == "te"){
    _in_te = false;
    assert ((*ptrel).pointeur() != NULL);
    ind->AddElement(*ptrel, gam_1ou2 - 1);
    return;
  } else {
    _e_xml(bal);
  }}
  return;
}

/*************************** Gestion des erreurs ***************************/

void ImportSystem::_e_xml(string bal_def) const
{
  cerr << "Une erreur de lecture du XML s'est produite" << endl;
  cerr << "La balise <" << bal_def << "> n'est pas reconnue." << endl;
  cerr << "Fin du programme." << endl;
  exit(1);
}

void ImportSystem::_e_gam(void) const
{
  cerr << "Le parser XML a détecté 3 balises <gamete> ... </gamete> adjacentes." << endl;
  cerr << "Jusqu'à preuve du contraire, un individu n'a que deux gametes." << endl;
  cerr << "Contactez le revendeur du programme, qui se fera une joie de" << endl;
  cerr << "corriger ce bug. Merci." << endl;
  exit(1);
}

/************************* Constructeur *************************************/

ImportSystem::ImportSystem(const std::string & file)
{	
	silence = false;
	
	sys = NULL;
	spe = NULL;
	pop = NULL;
	ind = NULL;
	ptrel = NULL;

  // Au début, on se situe hors de toutes les balises.
	_in_memory_table = false;
	_in_genome = false;
	_in_chromo = false;
	_in_element = false;
	_in_family = false;
	_in_locus = false;
	_in_chromosome = false;
	_in_fitness = false;
	_in_mean_fitness = false;
	_in_sd_fitness = false;
	_in_dup_rate = false;
	_in_del_rate = false;
	_in_activity = false;
	_in_dup_mut = false;
	_in_del_mut = false;
	_in_act_mut = false;
	_in_dup_sd = false;
	_in_del_sd = false;
	_in_act_sd = false;
	_in_regulation = false;
	_in_reg_mut = false;
	_in_reg_sd = false;
	_in_system = false;
	_in_species = false;
	_in_population = false;
	_in_individual = false;
	_in_male = false;
	_in_female = false;
	_in_gamete = false;
	_in_te = false;

	refelem = 0;
	 
	Element::Reset(); // reset element counter.
	// Individual::Reset(); // the parameters can change (?)
	
	if (file != "")
		parse_file(file);
}
