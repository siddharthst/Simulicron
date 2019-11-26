/***************************************************************************
                          simul.h  -  description
                             -------------------
    begin                : Jul 13 2004
    copyright            : (C) 2005 by Arnaud Le Rouzic
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
 
#ifndef SIMUL_H
#define SIMUL_H

#include <string>

class System;

class Simulation
{
	public:
	Simulation(System *, const std::string & out = "");
	void Go();  
	
	private:
	System * system;
	std::string prefix;
	
	bool verbose;
}; 
  
#endif // SIMUL_H
