/***************************************************************************
                          indivgroup.h  -  description
                             -------------------
    begin                : Tue July 15 2003
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

#ifndef INDIVGROUP_H
#define INDIVGROUP_H

#define CONTAINER_IND std::vector
#define USE_VECTOR_IND true
#define USE_LIST_IND false

#define EPSILON_W 0.000001

#include <assert.h>
#include <vector>
#include <list>
#include <algorithm>


#include "dyet.h"
#include "exception.h"
#include "base.h"
#include "element.h"

class Population;
class Individual;

class IndivGroup : public BaseBiolObject
{
	friend class ExportSystem;
	friend class System;
	
	public:		
	
		void AddIndividual(Individual *);
		void AddElement(Ptr2Element);
	
		Individual * RemoveRandomIndiv(void);
		Individual * ChooseRandomIndiv(void) const;
		Individual * ChooseRandomIndivByTE(void) const;
	
		Individual * GetParent(void) const;
		
		void ClearElements(void);
		
		unsigned int Size(void) const;
		unsigned int TEnumber(void) const;
		double Sum_Fitness(void) const;
		bool AreAllSterile(void) const;
	
		void Update(void);
		void Update_IndivNumber(void);
		void Update_TEnumber(void);
		void Update_Fitness(void);
		void UnUpdate(void);
		bool IsUpToDate(void) const;
		
		IndivGroup * VoidCopy(void) const;
	
		IndivGroup(void);
	  IndivGroup(const IndivGroup &);
		~IndivGroup(void);
		IndivGroup & operator = (const IndivGroup &);
	
	protected:
		Population * my_population;
		
		typedef CONTAINER_IND<Individual *> tabl_ind;
  	typedef tabl_ind::iterator it_ind;
  	typedef tabl_ind::const_iterator c_it_ind;
	
  	struct table_fitness
  	{ std::vector<double> table;	bool up_to_date; };
		struct integer_compt
		{ unsigned int value; bool up_to_date;};

  	tabl_ind ref_ind; 
  	integer_compt nb_ind;
  	integer_compt nb_et;
  	table_fitness fitness;		
		
		unsigned int ReturnIndivIndex(double) const;
		unsigned int ReturnIndivIndex_stl(double) const;
		
		void Clean(void);
		void Copy(const IndivGroup &);
};

#endif   // INDIVGROUP_H
