/***************************************************************************
                          individu.h  -  description
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
 
#ifndef INDIVIDU_H
#define INDIVIDU_H

#include <string>
#include <exception>

#include "dyet.h"
#include "base.h"
#include "gamete.h"
#include "element.h"

#ifndef RECOMBINATION_LOCAL 
	#define RECOMBINATION_LOCAL false
	#define RECOMBINATION_GLOBAL !RECOMBINATION_LOCAL
#endif

class IndivGroup;

class Individual : public BaseBiolObject
{
	friend class ExportSystem; 
	friend class System;
	
  public:
  
	Individual(const Individual&);
	Individual(void);
	Individual & operator=(const Individual &);
	~Individual(void);

	void Fecondation(const Individual *, const Individual *);
	void Transposition(void);

	void AddElement(Ptr2Element, short int);
	Ptr2Element PickRandomElement(void) const;
	void ClearElements(void);
  
	inline void Update(void);
	inline double Fitness(void);
	inline double Fitness(void) const;
	inline unsigned int TEnumber(void);
	inline unsigned int TEnumber(void) const;
	
	std::string Description(void);
	
	void Sort(void);

	protected:
		
	/* static const unsigned int dysg_cause_threshold = 1;
	static const unsigned int dysg_regul_threshold = 1; */
			
	IndivGroup * my_indivgroup;
	
	static	const bool verbose = false;	
	
	struct compteur_int
		{ 	long unsigned int compte;
			bool mise_a_jour; };

	struct compteur_double
		{ 	double compte;
    			bool mise_a_jour; };


	compteur_int nb_et;  
	compteur_double fitness; 
	compteur_double activite; 

	Gamete *gam1; 
	Gamete *gam2;

	void Clean(void);
	void Copy(const Individual &);

	inline void UnUpdate(void);

  // Mise à jour des variables propres au génome de l'individu
	void Update_TEnumber(void);
	void Update_Fitness(void);
	void Update_Activity(void);

  // Scan des gamètes
	double Sum_Activity(const long int copy_number) const;
	double Sum_FitnessAdd(void);
	// double Sum_FitnessMul(void) const; 
    
	  //**************** exceptions ******************
	class e_indiv : public std::exception
	{ public:
		inline virtual const char * what(void) const throw();	
		virtual ~e_indiv() throw() {}
	};

	class e_bad_act : public e_indiv
	{ public:
		inline virtual const char * what(void) const throw();
		~e_bad_act() throw() {}
	};
};

#include "individual.inl"

#endif
