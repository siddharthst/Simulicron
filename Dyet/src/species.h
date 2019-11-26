/***************************************************************************
                          species.h  -  description
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
  
#ifndef ESPECE_H
#define ESPECE_H

#include <assert.h>
#include <vector>
#include <list>

#include "dyet.h"
#include "base.h"
#include "element.h"

#define CONTAINER_POP std::vector

class Population;
class System;

class Species : public BaseBiolObject
{
	friend class ExportSystem;
	friend class System;
	
	public:
	 
	Species(const Species &);
	Species(void);
	Species(System *);
	Species &operator = (const Species&);
	~Species(void);


	Species * Reproduction(void);
	void Migration(void);

	void AddPopulation(Population *);
  
	Species * VoidCopy(void) const;
	void ClearElements(void);
	void AddElement(Ptr2Element);

	Ptr2Element PickRandomElement(void) const;
	Population * ChooseRandomPopByTE(void) const;

	inline void Update(void);
	inline void UnUpdate(void);
	inline unsigned int PopulationNumber(void);
	inline unsigned int PopulationNumber(void) const;
	inline unsigned int IndividualNumber(void);
	inline unsigned int TEnumber(void);
	inline unsigned int TEnumber(void) const;
	inline double AverageFitness(void) const;
	
	/* inline */ void TE_Add(unsigned int count = 1);
	/* inline */ void TE_Remove(unsigned int count = 1);
	
	protected :
	System * my_system;	
	
	typedef CONTAINER_POP<Population *> table_pop;
	typedef table_pop::iterator it_pop;
	typedef table_pop::const_iterator c_it_pop;

	double _tx; // taux de migration, initialisé à chaque construction.
	
	struct compteur_int
		{	unsigned long int effectif;
			bool mise_a_jour; };
			
	compteur_int nb_pop; // Nombre de populations dans l'espèce
	compteur_int nb_ind; // Nombre d'individus dans l'espèce
	compteur_int nb_et;  // Nombre d'éléments dans l'espèce
	table_pop ref_pop;
	
	long unsigned int my_te_number;

  inline void Update_PopulationNumber(void);
  inline void Update_IndividualNumber(void);
  inline void Update_TEnumber(void);

  void MigrationByPop(void);
  void MigrationByInd(void);
  unsigned int MigrantsNumber(long unsigned int) const;
  
  void Clean(void);
  void Copy(const Species &);
};

#include "species.inl"

#endif
