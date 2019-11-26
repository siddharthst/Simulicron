/***************************************************************************
                          gamete.h  -  description
                             -------------------
    begin                : Mon Dec 23 2002
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

/* ToDo :
    - Gérer correctement les exceptions (pour l'instant : que des "assert")
    - Tester si les algorithmes fonctionnent dans toutes les conditions!
    - Vérifier si le passage vector -> list est au point
    - Beaucoup de méthodes qui devraient être const ne le sont pas, à cause
    de pb sur les const_iterator.
 */
 
#ifndef GAMETE_H
#define GAMETE_H

#include <math.h>
#include <assert.h>
#include <vector>
#include <list>
#include <algorithm>

#include "dyet.h"
#include "base.h"
#include "element.h"

#define TE_CONTAINER std::list
#define TE_USE_LIST true
#define TE_USE_VECTOR false

#define DELETION_TRANSP_DEPENDANT false

class Individual;

class Gamete: public BaseBiolObject
{
  public:

	friend class ExportSystem;
	friend class System;
	friend class Recombination;
	friend class Recombination_Global;
	friend class Recombination_Local;

	typedef TE_CONTAINER<Ptr2Element> tabl_elem;
	typedef tabl_elem::iterator it_elem;
	typedef tabl_elem::const_iterator c_it_elem;

  public:
	Gamete();
	Gamete(const Gamete &);
	~Gamete();

	bool AddElement(Ptr2Element);
	void DeleteElement(it_elem);
	void ClearElements(void);
	bool Sort(void);
  	
	tabl_elem * Duplication(const double activity); //const
	void Excision(const double activity);
	
	Ptr2Element PickRandomElement(void) const;

	void Mutation_replication(void);
	// void Mutation_transposition(void);

  // Interface
  	inline unsigned int TEnumber(void) const;
  	inline double SumActivity(const long int copy_number) const; 
  	inline double SumFitness(bool burst = false) const; 
	
  private:
  	Individual * my_individual;
	
  	tabl_elem ref_et; // Tableau de pointeurs vers les éléments du gamete	
};

// is_sorted implementation (test if the container is sorted or not.
// this function is not STL standard.
template<class ForwardIterator> bool is_sorted(ForwardIterator first, ForwardIterator last);

#include "gamete.inl"

#endif
