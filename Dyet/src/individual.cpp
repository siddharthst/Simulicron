/***************************************************************************
                          individu.cpp  -  description
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

using namespace std; 

#include "individual.h"
#include "recombination.h"
#include "generator.h"
#include "parameters.h"
#include "exception.h"

// ----------------------- OBJECT INDIVIDUAL ------------------------------

void Individual::Fecondation(const Individual * p, const Individual *m)
{
  UnUpdate();
	
	assert (RECOMBINATION_LOCAL == !RECOMBINATION_GLOBAL);
	Recombination * recombinate;
	
	if (RECOMBINATION_LOCAL)
		recombinate = new Recombination_Local();
	else
		recombinate = new Recombination_Global();
	
	// current gametes have to be deleted!
	delete gam1;
	delete gam2;
	gam1 = NULL;
	gam2 = NULL;
	
	gam1 = recombinate->Do(p->gam1, p->gam2);
	gam2 = recombinate->Do(m->gam1, m->gam2);
 
	delete recombinate;
	recombinate = NULL;
	
	gam1->Mutation_replication();
	gam2->Mutation_replication();
	Update();
}  

void Individual::Transposition(void)
{
	Gamete::tabl_elem * vect_elem1;
	Gamete::tabl_elem * vect_elem2;

	Update_Activity();
	Update_TEnumber();
	double activity = activite.compte;
	long int copy_nb = nb_et.compte;

	if (activity < 0.) activity = 0.0;

	// Duplication
   	vect_elem1 = gam1->Duplication(activity);
   	vect_elem2 = gam2->Duplication(activity);

	// Excision
	gam1->Excision(activity);
	gam2->Excision(activity);
	  
	// Duplicated elements are reinserted
	Gamete::it_elem end = vect_elem1->end();
    	for (Gamete::it_elem beg = vect_elem1->begin(); beg != end; beg++)
    	{
     		AddElement((*beg), -1); // -1 : random gamete
    	}
    
    	end = vect_elem2->end();
    	for (Gamete::it_elem beg = vect_elem2->begin(); beg != end; beg++)
    	{
      		AddElement((*beg), -1);
    	}

	// not very clean. news have been done in Duplication()
    	delete vect_elem1;
    	delete vect_elem2;         

  	UnUpdate();
	// is it necessary here?
	gam1->Sort();
	gam2->Sort();
  	Update();
}

void Individual::AddElement(Ptr2Element el, short int ga)
{
  /* Insere l'élément el dans le gamete ga. 0 : mere, 1 : pere, -1 : hasard */

  assert ((ga >= -1) && (ga <= 1));
  switch (ga)
  {
    case -1 : // hasard
      if (generator->Uniform_int(1) == 0)
        gam1->AddElement(el);
      else
        gam2->AddElement(el);
    break;
    case 0 : // mère (gamete 1)
      gam1->AddElement(el);
    break;
    case 1: // père
      gam2->AddElement(el);
    break;
  } // Pas de default : c'est testé au début.
	UnUpdate();
}


Ptr2Element Individual::PickRandomElement(void) const
{
	assert (TEnumber() > 0);
	unsigned int and_the_winner_is = generator->Uniform_int(TEnumber());
	if (and_the_winner_is < gam1->TEnumber())
		return gam1->PickRandomElement();
	else
		return gam2->PickRandomElement();
	return NULL;
}

void Individual::ClearElements(void)
{
	gam1->ClearElements();
	gam2->ClearElements();
}

	
// ********************** Mises à jour et calculs ************************

void Individual::Sort(void)
{
  gam1->Sort();
  gam2->Sort();
}

void Individual::Update_TEnumber(void)
{
  // Met à jour la quantité d'éléments portés par l'individu.
  if (!nb_et.mise_a_jour)
  {
    nb_et.compte = gam1->TEnumber() + gam2->TEnumber();
    nb_et.mise_a_jour = true;
  }
}

void Individual::Update_Fitness(void)
{
  // Met à jour la fitness de l'individu.
  if (!fitness.mise_a_jour)
  {
	fitness.compte = Sum_FitnessAdd();
    fitness.mise_a_jour = true;
  }
  // cout3 << "Fitness de l'individu : " << fitness.compte << endl;
}

void Individual::Update_Activity(void)
{
  // Met à jour l'activité.
  if (!activite.mise_a_jour)
  {
    Update_TEnumber();
		if (nb_et.compte == 0)
			activite.compte = 0;
		else
    	activite.compte = (Sum_Activity(nb_et.compte)/(static_cast<double>(nb_et.compte)));
    activite.mise_a_jour = true;
  }
}
    
double Individual::Sum_FitnessAdd(void)
{
        /* WARNING : FITNESS MUST NOT BE < 0! */
 	double w;
        
        // only one "return" : clean, but a little slow...
        
        Update_TEnumber();
        Update_Activity();
        if (nb_et.compte == 0)
                w = 1.0;
        else
        {
                w = gam1->SumFitness() + gam2->SumFitness() + 1.0;
        }
        if (w < 0.) 
                w = 0;
        return w;
}

/* double Individual::Sum_FitnessMul(void) const
{
	// not implemented
  return 1;
} */

double Individual::Sum_Activity(const long int copynumber) const
{
   return gam1->SumActivity(copynumber) + gam2->SumActivity(copynumber);
}


/****************************************************************************/
/*   Constructeurs, destructeur, nettoyage, recopie (08/01/03)              */
/****************************************************************************/

Individual::Individual()
{	
	assert (parameter != NULL);
	assert (genome != NULL);
	assert (generator != NULL);
		
  gam1 = NULL;
  gam2 = NULL;  
  // Constructeur d'individu vide.
  gam1 = new Gamete();  // Construit deux gamètes vides
  gam2 = new Gamete();
	
  UnUpdate(); // Mieux vaut tout initialiser
}

Individual::Individual(const Individual &indiv)
// Constructeur de copie
{
  Copy(indiv);
}

Individual & Individual::operator = (const Individual & indiv)
{
  if (&indiv != this)
  {
    Clean();
    Copy(indiv);
  }
  return *this;
}

Individual::~Individual(void)
{
  Clean();
}


void Individual::Clean(void)
{
  // !!! Il ne faut surtout pas deleter les éléments! Ils sont communs
  // a tous les individus!!

  delete gam1;
  delete gam2;
  UnUpdate(); // est-ce vraiment utile?
}

void Individual::Copy(const Individual & modele)
{
  nb_et = modele.nb_et;
  fitness = modele.fitness;
  activite = modele.activite;

  UnUpdate(); // Il faut mieux dans un premier temps...
  
  gam1 = new Gamete(*(modele.gam1));
  gam2 = new Gamete(*(modele.gam2));
}

std::string Individual::Description(void)
{
	std::ostringstream outp;
	Update();
	outp << "Copy nb:" << TEnumber();
	outp	<< " fitness:" << Fitness();
	outp	<< " activity:" << activite.compte;
	return outp.str();
}
