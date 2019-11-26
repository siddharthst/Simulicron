/***************************************************************************
                          genome.cpp  -  description
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

#include "genome.h"

using namespace std;

/****************************** OBJET GENOME ****************************/

// non static

void Genome::Nouveau_chromo(const double longueur)
{
  chromo.push_back(longueur);
  taille += longueur;
	chromo_cum.push_back(taille);
  initOK = true;
}

unsigned int Genome::Retourne_chromo(double pos) const
{
  // retourne le chromosome correspondant à la position pos
  if (!initOK) throw e_gen_noinit();
	if (pos > taille) throw e_gen_tooshort();

	return find_if(chromo_cum.begin(), chromo_cum.end(), bind2nd(greater<double>(), pos)) - chromo_cum.begin();
}

double Genome::ChromoSz(const unsigned int i) const
{
  // Retourne la taille du chromosome i
  if (static_cast<int>(i) > NbChromo()) throw e_gen_nec();
  if (!initOK) throw e_gen_noinit();
  return chromo[i];
}

string Genome::Description(void) const
{
  if (!initOK) throw e_gen_noinit();
  ostringstream o;
  o << "Génome: taille=" << taille;
  o << ", nb_chromo=" << chromo.size() << endl;
  for (unsigned int a = 0; a < chromo.size(); a++)
    o << "K" << a << ":" << chromo[a] << "cM, ";
  return o.str();
}

// constructors // destructors:

Genome::Genome(void)
{
  taille = 0;
  initOK = false;
}

Genome::Genome(const Genome & g)
{
	Copy(g);
}

Genome & Genome::operator=(const Genome & g)
{
	Copy(g);
	return *this;
}

Genome::~Genome(void)
{
}

void Genome::Copy(const Genome & g)
{
	chromo = g.chromo;
	chromo_cum = g.chromo_cum;
	taille = g.taille;
	if (chromo.size() == 0)
		initOK = false;
	else
		initOK = true;
	assert (initOK == g.initOK);
}
