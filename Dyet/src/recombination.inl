/***************************************************************************
                          recombination.cpp  -  description
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

/* inline functions from recombination.h */

double Recombination::_distance(const Element * e1, const Element * e2) const
{
	assert (DIST_HALD == !DIST_TRUNC);
	assert ((e1 != NULL) && (e2 != NULL));
	double dist = -1.0;
	
	if (DIST_HALD)
		dist = _distance_hald(e1, e2);
	else
		dist = _distance_tronc(e1, e2);
	
	assert (dist >= 0.);
	return dist;
}

double Recombination::_distance_tronc(const Element * e1, const Element * e2) const
{
	  /* Calcul de la distance génétique par "troncature" :
     m <= 0.5  :  r = m
     m > 0.5   :  r = 0.5
  */
  double dist;

  if (e1->chromo() == e2->chromo())
  {
    dist = e1->position() - e2->position();
    assert (dist >= 0);
    if (dist >= 50)
      return 0.5;
    else
      return dist/100;
  }
  else
    return 0.5;
}

double Recombination::_distance_hald(const Element* e1, const Element* e2) const
{
	/* Calcul de la distance génétique par la formule de Haldane :
     r = 0/5 * (1 - exp(-2m))
  */
	
	if (e1 == e2) return 0.;

  if (e1->chromo() == e2->chromo())
    return 0.5*(1 - exp(-2*(e1->position() - e2->position())/100));
  else
    return 0.5;
}

void Recombination_Local::Change_gam(unsigned short int & curr_gam) const
{
  // Modifie le gamète courant.
  if (curr_gam == 0)
    curr_gam = 1;
  else
    curr_gam = 0;
}

Recombination_Local::It_gam Recombination_Local::ReturnIt(short int gam, Gamete::c_it_elem & iter) const
{
	Recombination_Local::It_gam itg(gam, iter++);
	return itg;
}

void Recombination_Global::Swap(Gamete::c_it_elem * &a, Gamete::c_it_elem * &b) const
{
	Gamete::c_it_elem * tmp = a;
	a = b;
	b = tmp;
}

Gamete::c_it_elem * Recombination_Global::First
		(Gamete::c_it_elem * a, Gamete::c_it_elem * b, 
			const Gamete::c_it_elem & max1, const Gamete::c_it_elem & max2) const
{
	bool end_a((*a == max1) || (*a == max2));
	bool end_b((*b == max1) || (*b == max2));
	
	if (end_a || end_b)
	{
		if (end_a)
		{
			if (end_b)
				return NULL;
			else
				return b;
		}
		else
			return a;
	}
	
	if ((**a)->position() < (**b)->position())
		return a;
	return b;
}
