/***************************************************************************
                          populations.h  -  description
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

#ifndef POPULATIONS_H
#define POPULATIONS_H

#include <assert.h>

#include "base.h"
#include "indivgroup.h"

class Species;
class Individual;

class Population : public BaseBiolObject
{
	friend class ExportSystem;
	friend class System;
		
	public:
		virtual Population * Reproduction(void);
		virtual Population * VoidCopy(void) const = 0;
		Individual * RemoveRandomIndiv(void);
		Ptr2Element PickRandomElement(void) const;
	
		void AddMale(Individual *);
		void AddFemale(Individual *);
		void AddElement(Ptr2Element);
	
		inline void ClearElements(void);
	
		inline virtual unsigned int Size(void) const;
		inline virtual unsigned int TEnumber(void) const;
		inline virtual double AverageFitness(void) const;
	
		void Update(void);
	
		virtual Population & operator = (const Population & p) {Clear(); Copy(p);return *this;}
		virtual ~Population(void)  {  } 
		virtual Population * Virtual_new_Population(void) const = 0;
		static Population * Copy_Constructor(const Population &);
		
		/* inline */ void TE_Add(unsigned int count = 1);
		/* inline */ void TE_Remove(unsigned int count = 1);
		
	protected:
		Species * my_species;
	
		IndivGroup * males;
		IndivGroup * females;
		static long unsigned int te_max_number;
	
		long unsigned int my_te_number;
	
		void Virtual_Constructor(void);
		inline virtual unsigned int MaleOffspringSize(void) const;
		inline virtual unsigned int FemaleOffspringSize(void) const;
		inline virtual void Clear(void);
		inline virtual void Copy(const Population &);
		inline virtual Individual * ChooseRandomIndivByTE(void) const;
};

class Dioic_Population : public Population
{
	friend class Export_systeme;
		
	public:
		Dioic_Population(void);
		Dioic_Population(const Dioic_Population &);
		Dioic_Population(Species *);
		Dioic_Population & operator =(const Dioic_Population &);
		virtual	~Dioic_Population(void);
		virtual Dioic_Population * VoidCopy(void) const;
	
		inline Dioic_Population * Virtual_new_Population(void) const
			{ return new Dioic_Population();}
};

class Herm_Population : public Population
{
	friend class Export_systeme;
				
	public:
		Herm_Population(void);
		Herm_Population(const Herm_Population &);
		Herm_Population(Species *);
		Herm_Population & operator = (const Herm_Population &);
		virtual ~Herm_Population(void);
		void AddIndividual(Individual *);
		virtual Individual * ChooseRandomIndivByTE(void) const;
		virtual Herm_Population * VoidCopy(void) const;
	
		inline Herm_Population * Virtual_new_Population(void) const
		{return new Herm_Population();}
		unsigned int Size(void) const;
		unsigned int TEnumber(void) const;
		double AverageFitness(void) const;
		unsigned int DuplicationNumber(void) const;
		unsigned int PotentialDuplicationNumber(void) const;
		
		unsigned int MaleOffspringSize(void) const;
		unsigned int FemaleOffspringSize(void) const;
	
	protected:
		void Copy(const Population &);
		void Clear(void);
};

// inline functions
#include "population.inl"

#endif // POPULATIONS_H
