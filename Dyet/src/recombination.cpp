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

#include "recombination.h"

#include "element.h"
#include "generator.h"
#include "genome.h"

Recombination_Local::Recombination_Local(void)
{
	//nothing to do.
}

Gamete * Recombination_Local::Do(const Gamete * pater, const Gamete * mater) const
{
 /* Méthode de recombinaison. Algorithme utilisant les itérateurs
     de la STL */
	Gamete * recombinant_gam = new Gamete();	
     
  if ((mater->ref_et.empty()) && (pater->ref_et.empty()))
    return recombinant_gam;
  // Si les deux gamètes sont vides, il n'y a rien à faire.
  
  It_gam curr_it;
  It_gam next_it;
  unsigned short int curr_gam;
  const Gamete::c_it_elem end_mere = mater->ref_et.end();
  const Gamete::c_it_elem end_pere = pater->ref_et.end();

	Gamete::c_it_elem pere = pater->ref_et.begin();
  Gamete::c_it_elem mere = mater->ref_et.begin();
  
  curr_gam = generator->Uniform_int(2);
  // cout << "Currgam : " << curr_gam << endl;
  curr_it = Nextone(mere, end_mere, pere, end_pere);
  
  if (curr_gam == curr_it.ga)
  {
      recombinant_gam->AddElement(*(curr_it.it));
  }
  
  for(;;)
  {
    if ((mere == end_mere) && (pere == end_pere))
      break;
    next_it = Nextone(mere, end_mere, pere, end_pere);
    if ( generator->Uniform() < _distance((*(next_it.it)).pointeur(),
                                    (*(curr_it.it)).pointeur()))
    {
      Change_gam(curr_gam);
    }
    curr_it = next_it;
    if (curr_gam == curr_it.ga)
    {
      recombinant_gam->AddElement(*(curr_it.it));
    }

  } 
  return recombinant_gam;    
}

Recombination_Local::It_gam Recombination_Local::Nextone
  (Gamete::c_it_elem & it1, const Gamete::c_it_elem & endit1,
   Gamete::c_it_elem & it2, const Gamete::c_it_elem & endit2) const
{
	
	if (it1 == endit1)
	{
		if (it2 == endit2)
		{ 
			// if (end1 && end2)
			assert(false);
		} else {
			// if (end1 && !end2)
			return ReturnIt(1, it2);
		}
	} else {
		if (it2 == endit2)
		{
			// if (!end1 && end2)
			return ReturnIt(0, it1);
		} else {
			// if (!end1 && !end2)
			if ((*it1)->position() < (*it2)->position())
			{
				return ReturnIt(0, it1);
			} else {
				return ReturnIt(1, it2);				
			}
		}
	}
}


/******************* RECOMBINATION_GLOBAL ***************************/

Gamete * Recombination_Global::Do(const Gamete * pater, const Gamete * mater) const
{
	const bool verbose = false;
	
	Gamete * recombinant_gam = new Gamete();
	
	std::vector<double> recomb_points = RandomRecombinationPoints();
	
	if (verbose){ for (unsigned int i = 0; i < recomb_points.size(); ++i)	std::cout << recomb_points[i] << " ";	std::cout << std::endl;	}
	
	Gamete::c_it_elem end_male = pater->ref_et.end();
	Gamete::c_it_elem end_female = mater->ref_et.end();
	std::vector<double>::const_iterator it_rec = recomb_points.begin();
	
	// Initialization (pater in first, no matters)
	Gamete::c_it_elem * matrix(&(pater->ref_et.begin()));
	Gamete::c_it_elem * nomatrix(&(mater->ref_et.begin()));
	
	Gamete::c_it_elem * first (First(matrix, nomatrix, end_male, end_female));
	
	if (generator->Uniform_int(2) == 0)
		Swap(matrix, nomatrix);

	if (first != NULL)
		do {
		
			while(*it_rec < (**first)->position())
				++it_rec;
			unsigned int chromo = (**first)->chromo();			
			
			if (*first == *matrix)
			{
				recombinant_gam->AddElement(**matrix);
				// if (verbose) std::cout << "Element Added!!!!!!!!!!!" << std::endl;
			}
		
			++(*first);
		
			first = First(matrix, nomatrix, end_male, end_female);
		
			if (first == NULL)
				break;
		
			if (chromo != (**first)->chromo())
			{
				if (generator->Uniform_int(2) == 0)
					Swap(matrix, nomatrix);
			} else {
				if (*it_rec < (**first)->position())
				{
					bool recomb = false;
					do 
					{
						++it_rec;
						recomb = !recomb;
					} while (*it_rec < (**first)->position());
					// if (verbose) std::cout << "Recombination" << std::endl;
					if (recomb)
						Swap(matrix, nomatrix);
				}
			}
		} while (true);
	return recombinant_gam;	
}

std::vector<double> Recombination_Global::RandomRecombinationPoints(void) const
{
	unsigned int nb_rec = generator->Poisson(genome->Taille() / 100.);
	std::vector<double> rec_points(nb_rec);
	for (unsigned int r = 0; r < nb_rec; ++r)
	{
		// memory is already reserved
		rec_points[r] = genome->Taille() * generator->Uniform();
	}
	// std::cout << "nbrec : " << nb_rec << " size : " << rec_points.size() << std::endl;
	rec_points.push_back(genome->Taille());
	std::sort(rec_points.begin(), rec_points.end());
	return rec_points;
}
