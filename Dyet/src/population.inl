/***************************************************************************
                          populations.inl  -  description
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

/* This is a .inl file, it is included in population.h and contains
					inline methods definitions. */

#include "species.h"
#include "generator.h"

/* inline virtual*/ 
void Population::Clear(void)
{
	assert (males != females);
	delete males;
	delete females;
}

void Population::Copy(const Population & p)
{
	// Warning : 'this' must be empty.
	assert (p.males != p.females);
	males = new IndivGroup(*p.males);
	females = new IndivGroup(*p.females);
}

unsigned int Population::Size(void) const
{
	assert (males != females);
	assert (males->IsUpToDate() && females->IsUpToDate());
	return males->Size() + females->Size();
}

unsigned int Population::TEnumber(void) const
{
	assert (males != females);
	assert (males->IsUpToDate() && females->IsUpToDate());
	return males->TEnumber() + females->TEnumber();
}

double Population::AverageFitness(void) const
{
	assert (males!= females);
	assert (males->IsUpToDate() && females->IsUpToDate());
	return (males->Sum_Fitness() + females->Sum_Fitness())/(static_cast<double>(males->Size() + females->Size()));
}	

Individual * Population::ChooseRandomIndivByTE(void) const
{
	assert (males != females);
	unsigned int and_the_winner_is = generator->Uniform_int(TEnumber());
	if (and_the_winner_is < males->TEnumber())
		return males->ChooseRandomIndivByTE();
	else
		return females->ChooseRandomIndivByTE();
	return NULL;
}

void Population::ClearElements(void)
{
	males->ClearElements();
	females->ClearElements();
}

unsigned int Population::MaleOffspringSize(void) const
{
	// The size of the populations are considered as fixed.
	return males->Size();
}

unsigned int Population::FemaleOffspringSize(void) const
{
	return females->Size();
}
