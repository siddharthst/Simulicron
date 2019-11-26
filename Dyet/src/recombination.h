/***************************************************************************
                          recombination.h  -  description
                             -------------------
    begin                : Fri Jan  9 11:15:48 2004
    copyright            : (C) 2003 by Arnaud Le Rouzic
    email                : <lerouzic@pge.cnrs-gif.fr>
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/


#ifndef RECOMBINATION_H
#define RECOMBINATION_H

#include <vector>

#include "base.h"
#include "gamete.h"

class Element;

#define DIST_HALD true
#define DIST_TRUNC false

class Recombination : public BaseBiolObject
{
	public:
		virtual ~Recombination() {}
		virtual Gamete * Do (const Gamete *, const Gamete *) const {return NULL;}
	
	protected:
		// Genetic distances calculations
  	inline double _distance_tronc(const Element *, const Element *) const;
  	inline double _distance_hald(const Element*, const Element*) const;
  	inline double _distance(const Element *, const Element *) const;
};

class Recombination_Global : public Recombination
{
	public:
		Recombination_Global() {}
		virtual ~Recombination_Global() {}
		virtual Gamete * Do (const Gamete *, const Gamete *) const;
			
	protected:
		std::vector<double> RandomRecombinationPoints(void) const;
		inline void Swap(Gamete::c_it_elem *&, Gamete::c_it_elem *&) const;
		inline Gamete::c_it_elem * First(Gamete::c_it_elem *, Gamete::c_it_elem *, const Gamete::c_it_elem &, const Gamete::c_it_elem &) const;
};

class Recombination_Local : public Recombination
{
	public:
		Recombination_Local();
		virtual ~Recombination_Local() {}
		virtual Gamete * Do (const Gamete *, const Gamete *) const;
			
	protected:

  struct It_gam{
		It_gam(unsigned short int g, Gamete::c_it_elem i) : ga(g), it(i) { }
		It_gam(void) { }
    unsigned short int ga;
    Gamete::c_it_elem it;};
		
	// Méthodes privées accessoires pour la recombinaison de type "it"
  inline It_gam Nextone(Gamete::c_it_elem &, const Gamete::c_it_elem &, 
												Gamete::c_it_elem &, const Gamete::c_it_elem &) const;
	inline It_gam ReturnIt(short int, Gamete::c_it_elem &) const;
  inline void Change_gam(short unsigned int &) const;
};
		
#include "recombination.inl"

#endif	//RECOMBINATION_H
