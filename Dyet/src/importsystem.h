/***************************************************************************
                          importsystem.h  -  description
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

#ifndef IMPORTSYSTEM_H
#define IMPORTSYSTEM_H

#include <libxml++/libxml++.h>
#include <string>
#include <map>

#include "element.h"

class System; 
class Gamete;
class Element;
class IndivGroup;
class Population;
class Dioic_Population;
class Herm_Population;
class Species;
class Individual;

class ImportSystem : public xmlpp::SaxParser
{
	public:

	ImportSystem(const std::string & file = "");
	~ImportSystem(void) {}

	System * GetSystem(void);
	
	void SetSilence(bool s = true) { silence = s; }

	private:
  
	void on_start_element (const std::string &, const AttributeMap &);
	void on_end_element (const std::string &);
 	void on_characters (const std::string &);


	void Traite_in_memory_start(const std::string &, const AttributeMap &);
  	void Traite_in_genome_start(const std::string &);
  	void Traite_in_element_start(const std::string &);
  	void Traite_in_system_start(const std::string &, const AttributeMap &);
  	void Traite_in_species_start(const std::string &, const AttributeMap &);
  	void Traite_in_population_start(const std::string &, const AttributeMap &);
  	void Traite_in_individual_start(const std::string &);
  	void Traite_in_gamete_start(const std::string &);

  	void Traite_in_memory_end(const std::string &);
  	void Traite_in_genome_end(const std::string &);
  	void Traite_in_element_end(const std::string &);
  	void Traite_in_system_end(const std::string &);
  	void Traite_in_species_end(const std::string &);
  	void Traite_in_population_end(const std::string &);
  	void Traite_in_individual_end(const std::string &);


  	void _e_xml(std::string bal_def = "--unknown--") const;
  	void _e_gam(void) const;


  	bool _in_memory_table;
  	bool _in_genome;
  	bool _in_chromo;
  	bool _in_element;
  	bool _in_family;
  	bool _in_locus;
  	bool _in_chromosome;
	bool _in_mean_fitness;
	bool _in_sd_fitness;
  	bool _in_fitness;
  	bool _in_dup_rate;
  	bool _in_del_rate;
  	bool _in_activity;
	bool _in_dup_mut;
	bool _in_dup_sd;
	bool _in_del_mut;
	bool _in_del_sd;
	bool _in_act_mut;
	bool _in_act_sd;
	bool _in_regulation;
	bool _in_reg_mut;
	bool _in_reg_sd;
  	bool _in_system;
  	bool _in_species;
  	bool _in_population;
  	bool _in_individual;
	bool _in_male;
	bool _in_female;
  	bool _in_gamete;
  	bool _in_te;

  	System * sys;
  	Species * spe;
  	Population * pop;
  	Individual * ind;
  	Ptr2Element * ptrel;
	Genome * genome;

  	std::map<long int, Ptr2Element> memory;
  	unsigned int refelem;
  	unsigned int gam_1ou2;
  	static bool blabla;  
	bool silence; //not the same as !blabla. silence == true -> No warnings at all
};

#endif // IMPORTSYSTEM_H
