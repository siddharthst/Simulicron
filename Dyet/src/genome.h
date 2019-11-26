/***************************************************************************
                          genome.h  -  description
                             -------------------
    begin                : Mon Dec 23 2002
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
 
#ifndef GENOME_H
#define GENOME_H

#include "dyet.h"

#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <exception>
#include <algorithm>

#include <assert.h>

class Genome
{
  public:
		  
  Genome(void);
	Genome(const Genome &);
	Genome & operator=(const Genome &);
  ~Genome(void);  

  void Nouveau_chromo(double);
  unsigned int Retourne_chromo(double) const;
  std::string Description(void) const;
  double ChromoSz(const unsigned int) const;

  inline bool EstVide(void) const
    {return !initOK;}
  inline double Taille(void) const
    { if (!initOK) throw e_gen_noinit();
      return taille;}
  inline int NbChromo(void) const
    { if (!initOK) throw e_gen_noinit();
      return chromo.size();}
  

  protected:
  
  std::vector<double> chromo;
	std::vector<double> chromo_cum;
	double taille;
  bool initOK;

	void Copy(const Genome &);
  // exceptions
  public:
  class e_gen : public std::exception
  { public:
		e_gen(void) { std::cerr << what() << std::endl;}
    virtual inline const char * what(void) const throw();
    virtual ~e_gen(void) throw() {}
  };

  class e_gen_tooshort : public e_gen
  { public:
		e_gen_tooshort(void) { std::cerr << what() << std::endl;}
    virtual inline const char * what(void) const throw();
    virtual ~e_gen_tooshort(void) throw() {}
  };

  class e_gen_noinit : public e_gen
  { public:
		e_gen_noinit(void) { std::cerr << what() << std::endl;}
    virtual inline const char * what(void) const throw();
    virtual ~e_gen_noinit(void) throw() {}
  };

  class e_gen_nec : public e_gen
  { public:
		e_gen_nec(void) { std::cerr << what() << std::endl;}
    virtual inline const char * what(void) const throw();
    virtual ~e_gen_nec(void) throw() {}
  };
};


const char * Genome::e_gen::what(void) const throw()
{
  return "Generic error in object \"Genome\"";
}

const char * Genome::e_gen_tooshort::what(void) const throw()
{
  return "Access attempt to a locus exceeding the genome size";
}

const char * Genome::e_gen_noinit::what(void) const throw()
{
  return "Access attempt to a non-initialized genome (0 chromosomes)";
}

const char * Genome::e_gen_nec::what(void) const throw()
{
  return "Access attempt to a non-existing chromosome.";
}

#endif
