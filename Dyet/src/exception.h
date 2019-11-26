/***************************************************************************
                          exception.h  -  description
                             -------------------
    begin                : Mon June 23 2003
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
 
#ifndef DYET_EXCEPTION_H
#define DYET_EXCEPTION_H

#include <exception>
#include <string>
#include <sstream>

#include "dyet.h"

class dyet_exception : public std::exception
{ public:
	virtual const char * what(void) const throw()
	{ return "Dyet general exception";}
	virtual ~dyet_exception(void) throw() { }
};

class execution_error : public dyet_exception
{ public:
	virtual const char * what(void) const throw() 
	{ return "General execution exception!";}
	virtual ~execution_error(void) throw() { }
};

class exec_toomuch : public execution_error
{ public:
	virtual const char * what(void) const throw ()
	{ return "Too many elements in the system.";}
	virtual ~exec_toomuch(void) throw() { }
};

class exec_noelements : public execution_error
{ public:
	virtual const char * what(void) const throw ()
	{ return "No more elements in the system.";}
	virtual ~exec_noelements(void) throw() { }
};

class exec_nomore : public execution_error
{ public:
	virtual const char * what(void) const throw()
	{ return "Biological object extinction";}
	virtual ~exec_nomore(void) throw() { }
};

class exec_nomore_indiv : public exec_nomore
{ public:
	virtual const char * what(void) const throw()
	{ return "No more individuals in the population";}
	virtual ~exec_nomore_indiv(void) throw() { }
};

class exec_nomore_indiv_reproduction : public exec_nomore_indiv
{ public:
	virtual const char * what (void) const throw()
	{ return "No reproduction possible.";}
	virtual ~exec_nomore_indiv_reproduction(void) throw() { }
}; 

class exec_nomore_indiv_reproduction_fitness : public exec_nomore_indiv_reproduction
{ public :
	virtual const char * what (void) const throw()
	{ return "No reproduction possible : population fitness drop to 0.";}
	virtual ~exec_nomore_indiv_reproduction_fitness(void) throw() { }
};

class exec_nomore_indiv_reproduction_sexratio : public exec_nomore_indiv_reproduction
{ public:
	virtual const char * what (void) const throw()
	{ return "No reproduction possible : the two sexes are not present";}
	virtual ~exec_nomore_indiv_reproduction_sexratio(void) throw() { }
};

class exec_nomore_population : public exec_nomore
{ public:
	virtual const char * what(void) const throw()
	{ return "No more populations in the species";}
	virtual ~exec_nomore_population(void) throw() { }
};

class exec_nomore_species : public exec_nomore
{ public:
	virtual const char * what(void) const throw()
	{ return "No more species in the system";}
	virtual ~exec_nomore_species(void) throw() { }
};

class e_switch_par : public dyet_exception
{ public:
	e_switch_par(const std::string & cod = "unknown") : code(cod) { }
	virtual const char * what(void) const throw()
	{ std::ostringstream o;
		o << "Error reading parameter " << code << " : bad value.";
		std::string s = o.str();
		return s.c_str(); }
	virtual ~e_switch_par(void) throw() { }
	protected:
		std::string code;
};

#endif // DYET_EXCEPTION_H
