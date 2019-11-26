/***************************************************************************
                          gamete.inl  -  description
                             -------------------
    begin                : Wed Jul 30 2003
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

// inline functions for class Gamete
// File included in "gamete.h"

template<class ForwardIterator> bool is_sorted(ForwardIterator first, ForwardIterator last)
{
   if (first != last)
   {
      ForwardIterator prev = first;
      for (++first; first != last; ++first)
      {
         if (*first < *prev)
            return false;
         prev = first;
      }
   } return true;
}

/****************** INTERFACE ****************************/

unsigned int Gamete::TEnumber(void) const
{ 
	return ref_et.size();
}

double Gamete::SumFitness(bool burst) const
{
  	double sum_fitness = 0;
  	c_it_elem beg = ref_et.begin();
	c_it_elem end = ref_et.end();
  	
	for (; beg != end; ++beg) 
		sum_fitness += (*beg)->fitness();

  return sum_fitness;
}

double Gamete::SumActivity(const long int copy_number) const
{
  /* sums the activity of all elements 
	note that, for each element, activity depends on the total genomic copy number */

	double temp = 0.;
	c_it_elem deb = ref_et.begin();
	c_it_elem fin = ref_et.end();

	if (copy_number == 0)
		assert (deb == fin);

	for (; deb != fin; ++deb) 
		temp += (*deb)->real_activity(copy_number);
 
  return temp;
}
