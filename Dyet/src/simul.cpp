/***************************************************************************
                          simul.cpp  -  description
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

#include <ostream>
#include <cassert> 
   
#include "simul.h"

#include "system.h"
#include "parameters.h"

using namespace std;

Simulation::Simulation(System * s, const string & out)
{
	assert (s != NULL);
	system = s;
	if (out == "")	
		{ prefix = "pop"; }
	else
		{ prefix = out;}
	verbose = false;
}

void Simulation::Go()
{
	int gen_nb = BaseBiolObject::GetParameter()->generations;
	int gen_st = BaseBiolObject::GetParameter()->steps;

 	for (int gen = 1; gen <= gen_nb; gen++)
	{
		// transposition and reproduction are simultaneous
		try {
			if ((gen == 1) || (gen == gen_nb) || (gen % gen_st == 0))
			{ 
				ostringstream o;
				o << prefix;
				o << "G" << (gen<10?"0":"") << (gen<100?"0":"") << (gen<1000?"0":"") << gen << ".xml";
				system->Save(o.str()); 
			}
		
			System * newsys = system->Reproduction();
			delete system;
			system = newsys;
    			// system->HTransfert();
    			// system->Migration();
    		
			system->AutoCheckUp();
		} 
		catch (exec_nomore_species & e) 
			{ if (verbose)
				// cerr << err_msg_dead(gen); 
				// OutCollection::SendFlag("EX");
			break; }
		catch (exec_toomuch & e) 
			{ 
				delete system; 
				system = NULL; 
				throw e;			
			break; }		
		catch (exec_noelements & e) 
			{ if (verbose)
				// cerr << err_msg_noel(gen); 
			break; }
  	}
	delete system;
	system = NULL;
}
