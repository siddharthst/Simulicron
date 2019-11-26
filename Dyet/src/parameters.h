/***************************************************************************
                          parameters.h  -  description
                             -------------------
    begin                : ven jul 02 2004
    copyright            : (C) 2004 by Arnaud Le Rouzic
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

#ifndef PARAMETERS_H
#define PARAMETERS_H

/* struct Parameters contains all the values of the parameters needed for a simulation.
	It can be considered as a "global" struct.
*/

#include <string>
#include <exception>
#include <sstream>

struct Parameters
{
	public:
	long int generations; 		// Number of generations
	long int steps;			// Number of steps between two consecutive outputs
	bool bursts;			// Whether transposition bursts are enable
	int burst_type;			// Burst type : 
						// 0: no burst
						// 1: random bursts
						// 2: hybrid dysgenesis
						// 3: continuous self-regulation
						// 4: threshold self-regulation
	private: static const int max_type = 4;
	
	public:
	long double burst_value;	// Differs following the burst type:
						// 0: No significance
						// 1: burst frequency
						// 2: No significance
						// 3: Half-life (or period) of the curve
						// 4: Threshold value
						
	unsigned int activity_matrix;

	Parameters(const std::string & filename);
	void SetValue(const std::string & label, const std::string & value);
	~Parameters() { }
	
	// Exceptions
	class e_parameters: public std::exception
	{ public:
		virtual const char * what() const throw() = 0;
		virtual ~e_parameters() throw() { };
	};
	
	class e_order: public e_parameters
	{ public:
		e_order(const std::string & l = "unknown") : label(l) { }
		virtual const char * what() const throw() { std::ostringstream o; o << "Order error on label " << label << "." << std::endl; return o.str().c_str();}
		virtual ~e_order() throw() { };
	  private:
		std::string label;
	};
	
	class e_unused: public e_parameters
	{ public:
		e_unused(const std::string & l = "unknown") : label(l) { }
		virtual const char * what() const throw() { std::ostringstream o; o << "Label " << label << " will be ignored." << std::endl; return o.str().c_str();}
		virtual ~e_unused() throw(){ };
	  private:
		std::string label;
	};
};

#endif // PARAMETERS_H
