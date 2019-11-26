/***************************************************************************
                          systeme.h  -  description
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
 
#ifndef SYSTEME_H
#define SYSTEME_H

#include <assert.h>
#include <vector>

#include "dyet.h"
#include "base.h"
#include "importsystem.h"
#include "exportsystem.h"
#include "element.h"

#define CONTAINER_SP std::vector
#define MAX_TE_LIMIT 1000000

class Species;
	
class System : public BaseBiolObject
{
	friend class ExportSystem;
	friend class ImportSystem;
  
	public:

  	System(void);
  	System(const System &);
	System & operator = (const System &);
 	~System(void);

	System * Reproduction(void);
 	void Migration(void);
	void HTransfert(void);
  
	void AutoCheckUp(void);

  	// Inputs / Outputs
	static System * Load(const std::string);
	void Save(const std::string); // const;

	inline void Update(void); 
	inline void UnUpdate(void);
	inline bool IsUpToDate(void) const;
	inline unsigned int SpeciesNumber(void);
	inline unsigned int SpeciesNumber(void) const;
 	inline long unsigned int TEnumber(void);
	inline long unsigned int TEnumber(void) const;
	inline long unsigned int IndividualNumber(void) const;
	std::vector<Ptr2Element> UpdateMemory(void) const;
	
	inline void TE_Add(unsigned int count = 1);
	inline void TE_Remove(unsigned int count = 1);

	private:	  

	typedef CONTAINER_SP<Species *> table_sp;
	typedef table_sp::iterator it_sp;
 	typedef table_sp::const_iterator c_it_sp;

	struct compteur_int
		{	unsigned long int effectif;
			bool mise_a_jour; };

	compteur_int nb_pop; // Nombre total de populations dans le système
	compteur_int nb_esp; // Nombre total d'espèces dans le système
	compteur_int nb_et;  // Nombre total d'éléments dans le système
	
	long unsigned int my_te_number;

	table_sp ref_esp; // Vecteur de pointeurs vers les espèces

	void AddSpecies(Species *);
	void RemoveSpecies(it_sp);
	
	Species * ChooseRandomSpeciesByTE(void) const;
	Species * ChooseRandomSpecies(void) const;
	
	void Copy(const System &);
	void Clean(void);

	inline void Update_TEnumber(void);
	inline void Update_PopulationNumber(void);
	inline void Update_SpeciesNumber(void);
};

#include "system.inl"

#endif
