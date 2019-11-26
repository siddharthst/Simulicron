/***************************************************************************
                          generator.cpp  -  description
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

#include <fstream>

#include "generator.h"
 
 using namespace std;
 
/******************************************************************************/
/*                                CLASS GENERATOR                            */
/******************************************************************************/

/***************************** PUBLIC FUNCTIONS ******************************/

// Constructors, Destructor

Generator::Generator(const string & fil )
{
	r = NULL;
  try {
    if (FileNameOK(fil))
      fic = fil;
    else
      throw e_bad_name();
  }
  catch (e_bad_name & e)
  {
    if (DEBUG_GENERATOR)
      cerr << e.what();
    fic = DEFAULT_SEED_FILE;
  }
  GetSeed();
  Initialize();
}

Generator::Generator(const Generator & gen)
{
	Copy(gen);
}

Generator & Generator::operator=(const Generator & gen)
{
	Clean();
	Copy(gen);
	return *this;
}

Generator::~Generator(void)
{
	Clean();
}


// Other Public Functions

string Generator::Display(void) const
{
  // Display the features of the generator (algorithm and seed)

  ostringstream oflu;

  oflu << "Generateur : " <<  gsl_rng_name (r) << ", Graine : " << graine;
  return oflu.str();
}


// Sample functions are inline and included in "generator.h"

/********************************* FUNCTIONS PRIVATE **********************/


// Initialisation

void Generator::Initialize(void)
{
  /* Init of the generator
     Modify the value of the seed in the seed file.
  */
  if (DEBUG_GENERATOR)
  {
	  cerr << "Random Number Generator Initialization.";
	  cerr << " Seed in file: " << fic << endl;
  }
  assert (r == NULL);
  r = gsl_rng_alloc (gsl_rng_mt19937);
  gsl_rng_set (r, graine);
  WriteSeed(NewSeed());
}

void Generator::GetSeed(void)
{
  // Get a valid seed
  try {
    if (!ReadOK())
      throw e_cant_open_file();
  }
  catch (e_cant_open_file & e)
  {
    if (DEBUG_GENERATOR)
      cerr << e.what() << endl;
    CreateFile();
  }
  ReadSeed();  
}

void Generator::ReadSeed(void)
{
  // Read the seed in the seed file
  ifstream rf;
  unsigned long int readseed;

  rf.open(fic.c_str(), ios::in);

  rf >> readseed;
  rf.close();

  if (SeedOK(readseed))
    graine = readseed;
  else
    graine = MakeSeed();
}

void Generator::WriteSeed(const unsigned long int gr) const
{
  /* Write the seed gr in the seed file */
  if (DEBUG_GENERATOR)
    cerr << "Seed " << gr << " has been written in the file " << fic << "." << endl;

  ofstream ficsor;

  ficsor.open(fic.c_str(), ios::trunc | ios::out);
  ficsor << gr;
  ficsor.close();
}

void Generator::CreateFile(void) const
{
  // Create a new seed file, with a seed from the timer
  if (DEBUG_GENERATOR)
   cerr << "New seed file. The seed is computed from the computer timer." << endl;  
  
  WriteSeed(MakeSeed());
}

// Get

unsigned long int Generator::MakeSeed(void) const
{
  time_t top;

  time(&top);
  return top;
}

unsigned long int Generator::NewSeed(void) const
{
  // Returns a random seed.
  if (DEBUG_GENERATOR)
    cerr << "New random seed." << endl;

  return gsl_rng_get(r);
}

// Test functions

bool Generator::ReadOK(void) const
{
  // Tests if the file can be read
  ifstream rf;

  rf.open(fic.c_str(), ios::in);
  if (! rf.is_open())
  {
    if (DEBUG_GENERATOR)
      cerr << "Impossible de lire le fichier de graines " << fic << "." << endl;
    return false;
  }
  rf.close();
  if (DEBUG_GENERATOR)
    cerr << "Lecture du fichier de graines possible." << endl;

  return true;
}

bool Generator::SeedOK(unsigned long int seedtotest) const
{
  // Returns true if the seed is OK

  if (seedtotest <= 0)
  {
    if (DEBUG_GENERATOR) {
      cerr << "The seed read (" << seedtotest << ") in the file ";
      cerr << fic << " is not valid!" << endl;}
    return false;
  } else {
    if (DEBUG_GENERATOR)
      cerr << "Seed read " << seedtotest << " OK." << endl;
    return true;
  }
}

bool Generator::FileNameOK(const string & filetotest) const
{
  // Returns true if the file name looks like a valid file name
	
  if (filetotest.size() < 3)
  {
    if (DEBUG_GENERATOR){
      cerr << "The file name " << filetotest << " is not " << endl;
      cerr << "a valid file name. It will be replaced by ";
      cerr << DEFAULT_SEED_FILE << "." << endl;}
    return false;
  } else {
    if (DEBUG_GENERATOR)
      cerr << "File name '" << filetotest << "' OK." << endl;
    return true;
  }
}

void Generator::Copy(const Generator & templ)
{
	graine = templ.graine;
	fic = templ.fic;
}

void Generator::Clean(void)
{
		gsl_rng_free (r);
}
