/***************************************************************************
                          base.cpp  -  description
                             -------------------
    begin                : Tue Aug 12 17:42:59 2003
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

#include <cassert>

#include "base.h"

const Parameters * BaseBiolObject::parameter = NULL;
const Generator * BaseBiolObject::generator = NULL;
const Genome * BaseBiolObject::genome = NULL;

void BaseBiolObject::Reset_static(void)
{
	generator = NULL;
	genome = NULL;
	parameter = NULL;
}

void BaseBiolObject::SetGenerator(const Generator * g)
{
	assert(BaseBiolObject::generator == NULL);
	BaseBiolObject::generator = g;
}

const Generator * BaseBiolObject::GetGenerator()
{
	return BaseBiolObject::generator;
}

void BaseBiolObject::SetGenome(const Genome * g)
{
	assert (BaseBiolObject::genome == NULL);
	BaseBiolObject::genome = g;
}

const Genome * BaseBiolObject::GetGenome(void)
{
	return BaseBiolObject::genome;
}

void BaseBiolObject::SetParameter(const Parameters * p)
{
	assert (BaseBiolObject::parameter == NULL);
	BaseBiolObject::parameter = p;
}

const Parameters * BaseBiolObject::GetParameter() 
{
	return BaseBiolObject::parameter;
}

bool BaseBiolObject::IsGeneratorInit(void)
{
	if (generator == NULL)
		return false;
	else
		return true;
}

bool BaseBiolObject::IsGenomeInit(void)
{
	if (genome == NULL)
		return false;
	else
		return true;
}

bool BaseBiolObject::IsParameterInit(void)
{
	if (parameter == NULL)
		return false;
	else
		return true;
}

bool BaseBiolObject::IsInit(void)
{
	return (BaseBiolObject::IsGeneratorInit() && BaseBiolObject::IsGenomeInit() && BaseBiolObject::IsParameterInit());
}
