/***************************************************************************
                          indivgroup.cpp  -  description
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

#include "indivgroup.h"

#include "individual.h"
#include "population.h"

using namespace std;

void IndivGroup::Update(void)
{
	it_ind end = ref_ind.end();
	for (it_ind beg = ref_ind.begin(); beg != end; beg++)
		(*beg)->Update();
	
	UnUpdate();
	Update_IndivNumber();
	Update_TEnumber();
	Update_Fitness();
}

void IndivGroup::Update_IndivNumber(void)
{
	if (nb_ind.up_to_date)
		return;
	else
	{
		nb_ind.value = ref_ind.size();
		nb_ind.up_to_date = true;
	}
}

void IndivGroup::Update_TEnumber(void)
{
		if (nb_et.up_to_date)
		return;
	else
	{
		c_it_ind end = ref_ind.end();
		nb_et.value = 0;
		for (c_it_ind cur = ref_ind.begin(); cur != end; ++cur)
			nb_et.value += (*cur)->TEnumber();
		nb_et.up_to_date = true;
	}
}

void IndivGroup::UnUpdate(void)
{
	nb_et.up_to_date = false;
	nb_ind.up_to_date = false;
	fitness.up_to_date = false;
}

bool IndivGroup::IsUpToDate(void) const
{
	return (nb_ind.up_to_date && nb_et.up_to_date);
}

unsigned int IndivGroup::Size(void) const
{
	assert(nb_ind.up_to_date);
	return nb_ind.value;
}

bool IndivGroup::AreAllSterile(void) const
{
	assert(fitness.up_to_date);
	if (Sum_Fitness() <= EPSILON_W)
		return true;
	return false;
}

unsigned int IndivGroup::TEnumber(void) const
{
	assert(nb_et.up_to_date);
	return nb_et.value;
}


void IndivGroup::AddIndividual(Individual * ind)
{
  ref_ind.push_back(ind);
  UnUpdate();
}

Individual * IndivGroup::RemoveRandomIndiv(void)
{
 	Update_IndivNumber();
	assert (generator != NULL);
	unsigned int numero_indiv (generator->Uniform_int(nb_ind.value));
	Individual * inter = NULL;	
	
  /* Il va falloir écrire un algo pour les vector et un pour les list
     * vector :
     - on stocke le pointeur
     - on échange le pointeur avec le dernier du vecteur
     - on supprime le dernier pointeur du vecteur   */
  #if (USE_VECTOR_IND)
    inter  = ref_ind[numero_indiv];
    if (numero_indiv != nb_ind.value - 1)
      ref_ind[numero_indiv] = ref_ind[nb_ind.value - 1];
    ref_ind.pop_back();
  #endif
    /* * list :
    - on stocke le pointeur
    - on recherche l'individu (long!)
    - on supprime l'individu de la liste (+ rapide?).*/

  #if (USE_LIST_IND)
    it_ind iterat (ref_ind.begin());
    for (int a = 0; a < numero_indiv; a++)
      iterat++;
    inter = *iterat;
    ref_ind.erase(iterat);
  #endif
  
	UnUpdate();
	assert (inter != NULL);
	return inter;
}

Individual * IndivGroup::GetParent(void) const
{
	// Returns a pointer towards an individual randomly chosen following
	// his fitness value.
	 double fittot = Sum_Fitness();
	
		if (fittot <= EPSILON_W)
			throw exec_nomore_indiv_reproduction_fitness();
		// theoretically, this exception must have been thrown before.
	
	unsigned int indice_parent = ReturnIndivIndex_stl(generator->Uniform()	* fittot);
	
	#if (USE_VECTOR_IND)
	  return ref_ind[indice_parent];
  #endif  // USE_VECTOR_IND
	#if (USE_LIST_IND)
		c_it_parent p (ref_ind.begin();
    for (int a = 0; a < indice_parent; a++)
      p++;
		return *p;
	#endif // USE_LIST_IND
}

Individual * IndivGroup::ChooseRandomIndivByTE(void) const
{
	unsigned int tmp = 0;
	c_it_ind ind, endind;
	endind = ref_ind.end();
	unsigned int and_the_winner_is = generator->Uniform_int(TEnumber());
	for (ind = ref_ind.begin(); ind != endind; ++ind)
	{
		tmp += (*ind)->TEnumber();
		if (and_the_winner_is < tmp)
			return (*ind);
	}
	return NULL;
}

Individual * IndivGroup::ChooseRandomIndiv(void) const
{
	unsigned int what_ind = generator->Uniform_int(Size());
	#if (USE_VECTOR_IND)
		return ref_ind[what_ind];
	#endif
	#if (USE_LIST_IND)
		c_it_ind ind = ref_ind.begin();
		for (unsigned int i; i < what_ind; ++i)
			++ind;
		return (*ind);
	#endif
}

void IndivGroup::ClearElements(void)
{
	it_ind ind, endind;
	endind = ref_ind.end();
	for (ind = ref_ind.begin(); ind != endind; ++ind)
		(*ind)->ClearElements();
}

void IndivGroup::AddElement (Ptr2Element ptr)
{
	Individual * destination = ChooseRandomIndiv();
	destination->AddElement(ptr, -1);
	UnUpdate();
}

void IndivGroup::Update_Fitness(void)
{
	double fitcum (0.);
	c_it_ind fin (ref_ind.end());
    
	fitness.table.clear();   
	for (c_it_ind cur = ref_ind.begin(); cur != fin; cur++)
	{
		fitcum += (*cur)->Fitness();
		fitness.table.push_back(fitcum);
	}
	fitness.up_to_date = true;
}

double IndivGroup::Sum_Fitness(void) const
{
	assert (nb_ind.up_to_date);
	assert (fitness.up_to_date);
	return fitness.table[nb_ind.value - 1];
}
	
void IndivGroup::Copy(const IndivGroup & ig)
{
	c_it_ind fin (ig.ref_ind.end());
  for (c_it_ind cur = ig.ref_ind.begin(); fin != cur; cur++)
  {
    Individual * ind = new Individual(**cur);
    ref_ind.push_back(ind);
  }

  nb_ind = ig.nb_ind;
	nb_et = ig.nb_et;
  fitness = ig.fitness;
  UnUpdate(); 
}

void IndivGroup::Clean(void)
{
	Update_IndivNumber();
  it_ind fin (ref_ind.end());
  for (it_ind cur = ref_ind.begin(); cur != fin; cur++)
    delete *cur;

  ref_ind.clear();
  UnUpdate();
}

IndivGroup * IndivGroup::VoidCopy(void) const
{
	IndivGroup * cpy = new IndivGroup();
	int nind = Size();
	for (int i = 0; i < nind; i++)
	{
		cpy->AddIndividual(new Individual());
	}
	return cpy;
}


IndivGroup::IndivGroup(void)
{
	assert (generator != NULL);
	assert (genome != NULL);
	assert (parameter != NULL);
	UnUpdate();
}

IndivGroup::IndivGroup(const IndivGroup & ig)
{
	Copy(ig);
}

IndivGroup & IndivGroup::operator = (const IndivGroup & ig)
{
	Clean();
	Copy(ig);
	return *this;
}

IndivGroup::~IndivGroup(void)
{
	Clean();
}

unsigned int IndivGroup::ReturnIndivIndex(double result) const
{
	// Dichotomic algorithm.
	assert (nb_ind.up_to_date);
	assert (fitness.up_to_date);
  assert (nb_ind.value == fitness.table.size());
  assert (nb_ind.value == ref_ind.size());

  int min = 0;
  int max = nb_ind.value - 1;
  int mil;
  for (unsigned int a = 0; a < fitness.table.size(); a++)
  do {
    mil = (min + max) / 2;
    if (fitness.table[mil] <= result)
      min = mil + 1;
    else
      max = mil;
  } while (min != max);
  return max;
}

unsigned int IndivGroup::ReturnIndivIndex_stl(double result) const
{
	// optimized search algorithm.
	assert (nb_ind.up_to_date);
	assert (fitness.up_to_date);
  assert (nb_ind.value == fitness.table.size());
  assert (nb_ind.value == ref_ind.size());

	return find_if(fitness.table.begin(), fitness.table.end(), bind2nd(greater<double>(), result)) - fitness.table.begin();
}
