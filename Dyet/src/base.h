/***************************************************************************
                          base.h  -  description
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

#ifndef _BASEBIOLOBJECT_H_
#define _BASEBIOLOBJECT_H_

#include "dyet.h"

class Generator;
class Genome;
class Parameters;

class BaseBiolObject
{
	public:
		virtual ~BaseBiolObject() { }
	
		static void Reset_static(void);
	
		static void SetGenerator(const Generator *);
		static const Generator * GetGenerator();
		static void SetGenome(const Genome *);
		static const Genome * GetGenome();
		static void SetParameter(const Parameters *);
		static const Parameters * GetParameter();
	
		static bool IsGeneratorInit(void);
		static bool IsGenomeInit(void);
		static bool IsParameterInit(void);
		static bool IsInit(void);
	
	protected:
		static const Generator * generator;
		static const Genome * genome;
		static const Parameters * parameter;	
};


#endif	//_BASEBIOLOBJECT_H_
