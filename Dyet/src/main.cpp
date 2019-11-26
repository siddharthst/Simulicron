/***************************************************************************
                          main.cpp  -  description
                             -------------------
    begin                : lun sep 23 15:42:23 CEST 2002
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
 
#include <iostream>
#include <string>
#include <sstream>
#include <exception>

#ifndef DYET
  #define DYET
#endif

#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include "dyet.h"

#include "argstream.h"
#include "parameters.h"
#include "generator.h"
#include "importsystem.h"
#include "simul.h"

using namespace std;

int main(int argc, char *argv[])
{
	//Intro();
	try {  
		string parfile = "";
		string initfile = "";
		string outfile = "";
		string seedfile = "seed.txt";
	
		argstream as(argc, argv);
		as 	>> parameter('p', "param", parfile, "General parameters file name", false)
			>> parameter('i', "init", initfile, "Initial population XML file", false)
			>> parameter('o', "output", outfile, "Output files prefix", false)
			>> parameter('s', "seed", seedfile, "Random number generator seed file", false)
			// >> option('v', "version", vers, "Display current version")
			>> help();
		as.defaultErrorHandling();
		
		BaseBiolObject::SetParameter(new Parameters(parfile));
		BaseBiolObject::SetGenerator(new Generator(seedfile));
		
		ImportSystem import; // BaseBiolObject::genome will be initialized here
		import.parse_file(initfile);
		
		Simulation simul(import.GetSystem(), outfile);
	
  		simul.Go();
  		
	} catch (exception & e)
		{ cout << "Execution Error. Exception thrown." << endl;
			cout << e.what() << endl;
			exit(-1);
	} catch (...)
		{ cout << "Error: unexpected exception." << endl;
			exit(-1);
	}
	// Conclusion();
	return EXIT_SUCCESS;
}
