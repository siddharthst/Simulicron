/***************************************************************************
                          generator.h  -  description
                             -------------------
    begin                : Tue Mar 25 2003
    copyright            : (C) 2003 by Arnaud Le Rouzic
    email                : lerouzic@pge.cnrs-gif.fr
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

/***************************************************************************
 Utilisation de la classe Generator.

 string Display(void) : renvoie la description du générateur (nom - graine).
 double Uniform(void) : renvoie une valeur tirée aléatoirement dans ]0;1].
 unsigned long int Uniform_int(unsigned long int) :
                      : renvoie un entier entre 0 (inclu) et l'argument (exclu).
 unsigned long int Poisson(double) : tirage dans une loi de Poisson.

 */
 
#ifndef GENERATOR_H
#define GENERATOR_H

#define DEBUG_GENERATOR false

#define DEFAULT_SEED_FILE "seed.txt"

// C++ headers
#include <string>
#include <sstream>
#include <iostream>
#include <exception>
#include <vector>

// C headers
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <cassert>

class Generator
{

  public:
  	Generator(const std::string & fil = "");
	Generator(const Generator &);
	Generator & operator=(const Generator &);
  	~Generator();	

 	std::string Display(void) const;

 	// Random samples
	inline bool ProbaSample(const long double & proba) const;
	inline double Uniform(void) const;
	inline unsigned long int Uniform_int(unsigned long int) const;
	inline unsigned long int Poisson(double) const;
	inline double Gaussian(double mu, double sigma) const;
	inline void Shuffle(std::vector<int> &) const;

  private:
	
	gsl_rng * r;
	  
  	std::string fic;
	unsigned long int graine;

  	void Initialize();
  
  	void GetSeed(void);
  	void ReadSeed(void);
  	void WriteSeed(const unsigned long int) const;
  	void CreateFile(void) const;
  	unsigned long int MakeSeed(void) const;
  	unsigned long int NewSeed(void) const;

  	bool ReadOK(void) const;
  	bool SeedOK(unsigned long int) const;
  	bool FileNameOK(const std::string &) const;
	
	void Clean(void);
	void Copy(const Generator &);

  public:
  // Exceptions
  class e_generator : public std::exception
  { public:
    inline virtual const char * what(void) const throw();
    virtual ~e_generator() throw() {}
  };
  class e_generator_fatal : public e_generator
  { public:
    inline virtual const char * what(void) const throw();
    virtual ~e_generator_fatal() throw() {}
  };
  class e_generator_nofatal : public e_generator
  { public:
    inline virtual const char * what(void) const throw();
    virtual ~e_generator_nofatal() throw() {}
  };
  class e_no_name : public e_generator_nofatal
  { public:
    inline const char * what(void) const throw();
  };
  class e_const_void : public e_generator_fatal
  { public:
    inline const char * what(void) const throw();
  };
  class e_bad_name : public e_generator_nofatal
  { public:
    inline const char * what(void) const throw();
  };
  class e_cant_open_file : public e_generator_nofatal
  {
    public:
    inline const char * what(void) const throw();
  };
};

// Méthodes d'échantillonnage inline

bool Generator::ProbaSample(const long double & proba) const
{
  /*
	Returns "true" with a frequency proba 
  */
  if (gsl_rng_uniform(r) < proba)
	return true;
  else
	return false;
}

double Generator::Uniform(void) const
{
  /*
  "This function returns a double precision floating point number
  uniformly distributed in the range [0,1). The range includes 0.0
  but excludes 1.0. The value is typically obtained by dividing the
  result of gsl_rng_get(r) by gsl_rng_max(r) + 1.0 in double precision.
  Some generators compute this ratio internally so that they can provide
  floating point numbers with more than 32 bits of randomness
  (the maximum number of bits that can be portably represented in a
  single unsigned long int)."
  (http://sources.redhat.com/gsl/ref/gsl-ref_17.html)
  */
  return gsl_rng_uniform(r);
}

unsigned long int Generator::Uniform_int(unsigned long int max) const
{
  /*
  This function returns a random integer from 0 to n-1 inclusive.
  All integers in the range [0,n-1] are equally likely, regardless
  of the generator used. An offset correction is applied so that zero
  is always returned with the correct probability, for any minimum
  value of the underlying generator.
  (http://sources.redhat.com/gsl/ref/gsl-ref_17.html)
  */
  return gsl_rng_uniform_int(r, max);
}

unsigned long int Generator::Poisson(double mu) const
{
  /*
  Random: unsigned int gsl_ran_poisson (const gsl_rng * r, double mu)
  This function returns a random integer from the Poisson distribution with
  mean mu. The probability distribution for Poisson variates is,
  p(k) = {\mu^k \over k!} \exp(-\mu)      for k >= 0.
  (http://sources.redhat.com/gsl/ref/gsl-ref_19.html#SEC310)
  */

  return gsl_ran_poisson(r, mu);
}



double Generator::Gaussian(double mu, double sigma) const
{
	/*
  This function returns a Gaussian random variate, with mean mu 
  and standard deviation sigma. The probability distribution for 
  Gaussian random variates is,

  p(x) dx = {1 \over \sqrt{2 \pi \sigma^2}} \exp (-x^2 / 2\sigma^2) dx

  for x in the range -\infty to +\infty. Use the transformation 
  z = \mu + x on the numbers returned by gsl_ran_gaussian to obtain 
  a Gaussian distribution with mean \mu. This function uses the Box-Mueller 
  algorithm which requires two calls the random number generator r. 
  (http://sources.redhat.com/gsl/ref/gsl-ref_19.html#SEC283)
*/
	return (gsl_ran_gaussian(r, sigma) + mu);
}

void Generator::Shuffle(std::vector<int> & vectoshuffle) const
{
  gsl_ran_shuffle(r, &vectoshuffle[0], vectoshuffle.size(), sizeof(int));
  // Communication avec une bibliothèque C délicate. Normalement, ça marche,
  // même si ça n'est pas garanti par la norme.
}                                     

// exceptions.

const char * Generator::e_generator::what(void) const throw()
{
  return "Generic error in the random number generator.";
}

const char * Generator::e_generator_fatal::what(void) const throw()
{
  return "Fatal error using the random number generator. The program will stop.";
}

const char * Generator::e_generator_nofatal::what(void) const throw()
{
  return "Non-fatal error using the random number generator.";
}

const char * Generator::e_no_name::what(void) const throw()
{
  std::ostringstream os;
  os << "Random generator seed file name not specified. Default name: "
     << DEFAULT_SEED_FILE << ".";
  std::string rep = os.str();
  return rep.c_str();
}  

const char * Generator::e_const_void::what(void) const throw()
{
  return "Programming error: the default constructor must not be called.";
}

const char * Generator::e_bad_name::what(void) const throw()
{
  std::ostringstream os;
  os << "Bad random generator seed file. Default name: "
     << DEFAULT_SEED_FILE << ".";
  std::string rep = os.str();
  return rep.c_str();
}

const char * Generator::e_cant_open_file::what(void) const throw()
{
  std::ostringstream os;
  os << "Error opening the specified random generator seed file." << std::endl;
  os << "A new one will be created.";
  std::string rep = os.str();
  return rep.c_str();
}

#endif
