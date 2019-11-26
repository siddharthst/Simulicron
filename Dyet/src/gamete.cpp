/***************************************************************************
                          gamete.cpp  -  description
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
#include <cassert>

#include "gamete.h"
#include "generator.h"

#define EPSILON_ACTIVITY 0.000001

using namespace std;

/******************************** OBJET GAMETE ********************************/

void Gamete::Mutation_replication(void)
{
	it_elem it_b, it_e;
	it_b = ref_et.begin();
	it_e = ref_et.end();
	
	for (; it_b != it_e; ++it_b)
	{
		bool is_mutation = false;
		Ptr2Element nouvo  = ((*it_b)->Mutation_replication(is_mutation));
		if (is_mutation)
			*it_b = nouvo;
	}
}


Gamete::tabl_elem * Gamete::Duplication(const double activity) // const
{ 
  tabl_elem * v = new tabl_elem;
  if (activity < EPSILON_ACTIVITY)
	return v;

  it_elem fin = ref_et.end();

  for (it_elem current (ref_et.begin()); current != fin ; ++current)
  {	
	if ( generator->ProbaSample(activity * (*current)->duplication_rate()) )
	{
		v->push_back((*current)->Duplication());
		// Duplication includes Mutation_duplication
	}
  }
  return v;
}

void Gamete::Excision(const double activity)
{
  
	if (TE_USE_LIST)
	{ // If the container is a list, we just have to do an "erase"
	it_elem fin(ref_et.end());

	// Strange "for" loop : erase will modify the iterator, and ++courrant
	// may put the iterator to end() + 1.
   	for (it_elem courrant = ref_et.begin(); courrant != fin ; )
    	{
		double real_deletion_rate;
		if (DELETION_TRANSP_DEPENDANT)
			real_deletion_rate = (*courrant)->deletion_rate() * activity;
		else
			real_deletion_rate = (*courrant)->deletion_rate();
			
      		if (generator->ProbaSample(real_deletion_rate))
      		{
        		courrant = ref_et.erase(courrant);
      		} else {
        		++courrant;
		}
    	}
	} // TE_USE_LIST
  /* if (TE_USE_VECTOR)
  { // It is more complex with a vector
    it_elem fin(ref_et.end());
    vector<it_elem> v_temp;

    for (it_elem courrant = ref_et.begin(); courrant != fin ; courrant++)
    {
			double deletion_rate;
			
			if (burst)
				deletion_rate = (*courrant)->deletion_rate_burst();
			else
				deletion_rate = (*courrant)->deletion_rate();
			
      if (generator->Uniform() < deletion_rate * activite)
      {
        v_temp.push_back(courrant);
      }
    }
    int sz (v_temp.size());
    for (int a = 0; a < sz; a++)
      DeleteElement(v_temp[a]);
    Sort(); // Is that necessary?
  } */
}  

bool Gamete::AddElement(Ptr2Element elem)
{
  ref_et.push_back(elem);
  return true;
}

bool Gamete::Sort(void)
{
  // Sort the elements following their position in the genome
  #if (TE_USE_VECTOR)
	if (!is_sorted(ref_et.begin(), ref_et.end())
    sort(ref_et.begin(), ref_et.end());
  #endif
  #if (TE_USE_LIST)
	if (!is_sorted(ref_et.begin(), ref_et.end()))
    ref_et.sort();
  #endif
  return true;
}


void Gamete::DeleteElement(it_elem ptr)
{
  /* Delete an element of the gamete */
  
  it_elem fin;

  #if (TE_USE_LIST)
    fin = --ref_et.end();  // does not compile with vectors (? Am I sure?)
  #endif
  #if (TE_USE_VECTOR)
    fin = ref_et.end() - 1; // does not compile with lists
  #endif


  ptr = fin;
  ref_et.pop_back();
}
   
Ptr2Element Gamete::PickRandomElement(void) const
{
	unsigned int and_the_winner_is = generator->Uniform_int(TEnumber());
	unsigned int tmp = 0;
	c_it_elem te, teend;
	teend = ref_et.end();
	for (te = ref_et.begin(); te != teend; te++)
	{
		tmp++;
		if (and_the_winner_is < tmp)
			return (*te)->Duplication();
	}
	return *teend; //bof
}

void Gamete::ClearElements(void)
{
	ref_et.clear();
}


/****************************************************************************/
/*                  Gamete : constructors and destructor                  */
/****************************************************************************/

Gamete::Gamete()
{
	assert (generator != NULL);
	assert (genome != NULL);
	assert (parameter != NULL);
}

Gamete::Gamete(const Gamete & gam)
// Copy Constructor
{
	c_it_elem fin (gam.ref_et.end());
	for (c_it_elem deb = gam.ref_et.begin(); deb != fin; ++deb)
	{
		Ptr2Element p(*deb);
		ref_et.push_back(p);
	}
}


Gamete::~Gamete(void)
{
  // Nothing to do
}
