/***************************************************************************
                          populations.cpp  -  description
                             -------------------
    begin                : Wed July 16 2003
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

#include "population.h"
#include "exception.h"

using namespace std;

unsigned long int Population::te_max_number = 0; 
// static variable initialized in Virtual_Constructor() 

/************************ Population ***********************/

Population * Population::Reproduction(void) // not const because Update_fitness
{
	males->Update_Fitness();
	females->Update_Fitness();
	
		if ((males->Size() == 0) || (females->Size() == 0))
			throw exec_nomore_indiv_reproduction_sexratio();
		if ((males->AreAllSterile()) || (females->AreAllSterile()))
			throw exec_nomore_indiv_reproduction_fitness();
	
	Population * offspring = Virtual_new_Population();
		
	unsigned int moffsize = MaleOffspringSize();
	unsigned int foffsize = FemaleOffspringSize();
	unsigned int offsize = moffsize + foffsize;
	long unsigned te_nb_in_offspring_pop = 0;

	for (unsigned int a = 0; a < offsize; a++)
	{
		Individual * i = new Individual();
		i->Fecondation(males->GetParent(), females->GetParent());
		i->Transposition();
		
		if (a < moffsize)
			offspring->AddMale(i);			
		else
			offspring->AddFemale(i);
		
		te_nb_in_offspring_pop += i->TEnumber();
		assert (te_max_number > 0);
		if (te_nb_in_offspring_pop > te_max_number)
		/* If _LIMITE_NB_ET is passed over in one population, then 
		the simulation run have to be stopped. It is not a garantee
		that the limit is not reached when all the populations are considered. */
		{
			delete offspring;
			throw exec_toomuch();
		}
	}
	return offspring;
}

Individual * Population::RemoveRandomIndiv(void)
{
	if (generator->Uniform_int(2) == 0)
		return males->RemoveRandomIndiv();
	else
		return females->RemoveRandomIndiv();
	
	// It works with Herm_Population, even if it is not optimal.
	// There is no way to know sex of the removed individual... !
}

void Population::AddMale(Individual * m)
{
	assert (m != NULL);
	males->AddIndividual(m);
}

void Population::AddFemale(Individual * f)
{
	assert (f != NULL);
	females->AddIndividual(f);
}

void Population::AddElement(Ptr2Element ptr)
{
	unsigned what_ind = generator->Uniform_int(Size());
	if (what_ind < males->Size())
		males->AddElement(ptr);
	else
		females->AddElement(ptr);
}

Ptr2Element Population::PickRandomElement(void) const
{
	assert(TEnumber() > 0);
	Individual * ind = ChooseRandomIndivByTE();
	return ind->PickRandomElement();
}

void Population::Virtual_Constructor(void)
{
	assert (generator != NULL);
	assert (parameter != NULL);
	assert (genome != NULL);
	
	if (te_max_number == 0)
	{ 
		// Problem in Parameter lib : only integers can be parsed.
		// I fear an overflow -> assert.
		te_max_number = MAX_TE_LIMIT;
		assert ((te_max_number >= 1000) && (te_max_number <= 1000*1000));
	}
}

Population * Population::Copy_Constructor(const Population & pop)
{
	Population * ppp = pop.Virtual_new_Population();
	ppp->Copy(pop);
	return ppp;
}

void Population::Update(void)
{
	males->Update();
	females->Update();
}

void Population::TE_Add(unsigned int count) // default = 1
{
	my_te_number += count;
	my_species->TE_Add(count);
}

void Population::TE_Remove(unsigned int count) // default = 1
{
	assert (count <= my_te_number);
	my_te_number -= count;
	my_species->TE_Remove(count);
}


/***************************** Dioic_Population ******************/

Dioic_Population * Dioic_Population::VoidCopy(void) const
{
	assert (males != females);
	Dioic_Population * cpy = new Dioic_Population();
	delete cpy->males;
	delete cpy->females;
	cpy->males = this->males->VoidCopy();
	cpy->females = this->females->VoidCopy();
	return cpy;
}

Dioic_Population::Dioic_Population(void)
{
	Virtual_Constructor();
	males = new IndivGroup();
	females = new IndivGroup();
}

Dioic_Population::Dioic_Population(const Dioic_Population & p)
{
	Copy(p);
}

Dioic_Population & Dioic_Population::operator = (const Dioic_Population & p)
{
	Clear();
	Copy(p);
	return *this;
}

Dioic_Population::~Dioic_Population(void)
{
	Clear();
}

/**************************** Herm_Population *************************/

Herm_Population::Herm_Population(void)
{
	Virtual_Constructor();
	// males and females are the same IndivGroup!
	males = new IndivGroup();
	females = males;
}

Herm_Population::Herm_Population(const Herm_Population & p)
{
	Copy(p);
}

Herm_Population & Herm_Population::operator = (const Herm_Population & p)
{
	Clear();
	Copy(p);
	return *this;
}

Herm_Population::~Herm_Population(void)
{
	Clear();
}

void Herm_Population::AddIndividual(Individual * i)
{
	assert (i != NULL);
	males->AddIndividual(i);
}

// virtual methods redefinition

Individual * Herm_Population::ChooseRandomIndivByTE(void) const
{
	assert (males == females);
	return males->ChooseRandomIndivByTE();
}

Herm_Population * Herm_Population::VoidCopy(void) const
{
	assert (males == females);
	Herm_Population * cpy = new Herm_Population();
	cpy->males = this->males->VoidCopy();
	cpy->females = cpy->males;
	return cpy;
}

void Herm_Population::Clear(void)
{
	assert (males == females);
	delete males;
	males = NULL;
	females = NULL;
}

void Herm_Population::Copy(const Population & pop)
{
	const Herm_Population & p (dynamic_cast<const Herm_Population &>(pop));
	assert (p.males == p.females);
	males = new IndivGroup(*p.males);
	females = males;
}

unsigned int Herm_Population::Size(void) const
{
	assert (males->IsUpToDate());
	return males->Size();
}

unsigned int Herm_Population::TEnumber(void) const
{
	assert (males->IsUpToDate());
	return males->TEnumber();
}

double Herm_Population::AverageFitness(void) const
{
	assert(males->IsUpToDate());
	return (males->Sum_Fitness() / static_cast<double>(males->Size()));
}

unsigned int Herm_Population::MaleOffspringSize(void) const
{
	return males->Size();
}
	
unsigned int Herm_Population::FemaleOffspringSize(void) const
{ 
	return 0; // females are also males in an Herm_Population
}
