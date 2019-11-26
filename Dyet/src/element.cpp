/***************************************************************************
                          element.cpp  -  description
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

#include <cmath>

#include "element.h" 
#include "generator.h"
#include "parameters.h"
#include "genome.h"

using namespace std;

//---------------------- OBJET ELEMENT --------------------------------

// reference count initialization
unsigned long int Element::_idref = 0;


Ptr2Element Element::Duplication(void) const
{  
	/* 	Creates a new element after a transposition.
		Identical to the "father" element, except for the variable features */
	Ptr2Element nouvo (new Element(*this));

	nouvo->SetRandomPosition();
	nouvo->Mutation_transposition();
	return nouvo;
}

Ptr2Element Element::Mutation_replication(bool & mutation) const
{
	/* 	Mutations during the DNA replication process are a rare phenomenon.
		Except for the fitness change (which has not really any biological meaning), 
		all mutations occur with the probability defined by the mutation rate for
		each feature. */
	
	Ptr2Element nouvo (new Element(*this));
	mutation = false;
	
	if (MUTATION_FITNESS_DURING_REPLICATION)
	{
		mutation = true;
		nouvo->SetRandomFitness();
	}
	if (MUTATION_DUPRATE_DURING_REPLICATION)
	{
		if (generator->ProbaSample(nouvo->duplication_mutation_rate()))
		{
			mutation = true;
			nouvo->SetRandomDupRate();
		}
	}
	if (MUTATION_DELRATE_DURING_REPLICATION)
	{
		if (generator->ProbaSample(nouvo->deletion_mutation_rate()))
		{
			mutation = true;
			nouvo->SetRandomDelRate();
		}
	}
	if (MUTATION_ACTIVITY_DURING_REPLICATION)
	{
		if (generator->ProbaSample(nouvo->activity_mutation_rate()))
		{
			mutation = true;
			nouvo->SetRandomActivity();
		}
	}		
	if (MUTATION_REGULATION_DURING_REPLICATION)
	{ 
		if (generator->ProbaSample(nouvo->regulation_mutation_rate()))
		{
			mutation = true;
			nouvo->SetRandomRegulation();
		}
	}

	return nouvo;
}

void Element::Mutation_transposition(void)
{
	/* 	Mutations during the transposition process always occur.
		The corresponding mutation rate is not used. */
	
	if (MUTATION_FITNESS_DURING_TRANSPOSITION)
	{
		SetRandomFitness();
	}
	if (MUTATION_DUPRATE_DURING_TRANSPOSITION)
	{
		SetRandomDupRate();
	}	
	if (MUTATION_DELRATE_DURING_TRANSPOSITION)
	{
		SetRandomDelRate();
	}
	if (MUTATION_ACTIVITY_DURING_TRANSPOSITION)
	{
		SetRandomActivity();
	}
	if (MUTATION_REGULATION_DURING_TRANSPOSITION) 
	{ 
		SetRandomRegulation();
	}
}

/***************************************************************************/
/*                   PRIVATE VARIABLES SETTINGS                         */
/***************************************************************************/

/*************************************/
/*      Random settings              */
/*************************************/

void Element::SetRandomFitness(void)
{
	if (fixed_features._sd_fitness.value != 0.)
		variable_features._fitness.value = 
			generator->Gaussian(fixed_features._mean_fitness.value, fixed_features._sd_fitness.value);
	else
		variable_features._fitness.value = fixed_features._mean_fitness.value;
}

void Element::SetRandomPosition(void)
{
	Set_position(genome->Taille() * (generator->Uniform()));
}

void Element::SetRandomDupRate(void)
{
	int sign = 0;
	generator->ProbaSample(0.5)? sign = 1 : sign = -1;

	double log_deviation = std::abs(generator->Gaussian(0., duplication_mutation_sd()));
	
	double new_log_dup_rate = sign*log_deviation + log10(duplication_rate());
	Set_dup_rate(pow(10.0, new_log_dup_rate));
	
	if (duplication_rate() <= low_limit_rates)
		Set_dup_rate(low_limit_rates);
}

void Element::SetRandomDelRate(void)
{	
	int sign = 0;
	generator->ProbaSample(0.5)? sign = 1 : sign = -1;

	double log_deviation = std::abs(generator->Gaussian(0., deletion_mutation_sd()));
	
	double new_log_del_rate = sign*log_deviation + log10(deletion_rate());
	
	if (new_log_del_rate > 0.0)
		new_log_del_rate = 0.0;
	
	Set_del_rate(pow(10.0, new_log_del_rate));
	
	if (deletion_rate() <= low_limit_rates)
		Set_del_rate (low_limit_rates);
}

void Element::SetRandomActivity(void)
{
	#if (ACTIVITT_CONT)
		double delta_act = generator->Gaussian(0., activity_mutation_sd());

		if (!ACTIVITY_MUTATION_UP)
			delta_act = - std::abs(delta_act);

		if ((ACTIVITY_ABSORBING_STATE_0) && (activity() == 0.0))
			delta_act = 0.0;

		Set_activity(activity() + delta_act);
	#else 
		ActivityMatrix * p_matrix = ActivityMatrix::GetMatrix(BaseBiolObject::GetParameter()->activity_matrix);
		// not very clear. activity_matrix is the index of the matrix used, set in the parameter file.

		Set_activity(p_matrix->RandomActivity(activity(), generator->Uniform()));
	#endif
}

void Element::SetRandomRegulation(void)
{
	Set_regulation_factor(regulation_factor() + generator->Gaussian(0, regulation_mutation_sd()));
}

/**********************************/
/*      Fixed settings            */
/**********************************/

void Element::Set_position(const double p)
{
  /* 	1) change the value of the position
		2) determine the new chromosome */
     
  assert (p >= 0);
  assert (p <= genome->Taille() );

  variable_features._position.value = p;
  variable_features._chromo.value = genome->Retourne_chromo(variable_features._position.value);
}

void Element::Set_fitness(const double fit)
{
	// no conditions...
	variable_features._fitness.value = fit;
}

void Element::Set_mean_fitness(const double mean_w)
{
	// no special conditions...
	fixed_features._mean_fitness.value = mean_w;
}

void Element::Set_sd_fitness(const double var_w)
{
	assert(var_w >= 0.);
	fixed_features._sd_fitness.value = var_w;
}

void Element::Set_family(const short int f)
{
  // assert (f < 2); 
	fixed_features._family.value = f;
}

// ##### Duplication rate

void Element::Set_dup_rate(const double dr)
{
	assert (dr >= 0.0);
	if (dr < low_limit_rates)
		variable_features._duplication_rate.value = 0.0;
	else
		variable_features._duplication_rate.value = dr;
}

void Element::Set_mut_dup_rate(const double mdr)
{
	assert ((mdr >= 0.0) && (mdr <= 1.0));
	if ((mdr < low_limit_rates) && (mdr != 0.0))
	{
		cerr << "Warning: dup_mut_rate < LOW_LIMIT_RATES!" << endl;
		fixed_features._dup_mut.value = 0.0;
	} else {
		fixed_features._dup_mut.value = mdr;
	}
}

void Element::Set_sd_dup_rate(const double sddup)
{
	assert (sddup >= 0.0);
	fixed_features._dup_sd.value = sddup;
}

// ##### Deletion rate

void Element::Set_del_rate(const double dr)
{
	assert (dr >= 0.0);
	assert (dr <= 1.0);
	if (dr < low_limit_rates)
		variable_features._deletion_rate.value = 0.0;
	else
		variable_features._deletion_rate.value = dr;
}

void Element::Set_mut_del_rate(const double mdr)
{
	assert ((mdr >= 0.0) && (mdr <= 1.0));
	if ((mdr < low_limit_rates) && (mdr != 0.0))
	{
		cerr << "Warning: del_mut_rate < LOW_LIMIT_RATES!" << endl;
		fixed_features._del_mut.value = 0.0;
	} else {
		fixed_features._del_mut.value = mdr;
	}
}

void Element::Set_sd_del_rate(const double sddel)
{
	assert (sddel >= 0.0);
	fixed_features._del_sd.value = sddel;
}

// ##### Activity

void Element::Set_activity(double act)
{
	#if ACIVITY_CONT
		if ((ACTIVITY_THRESHOLD_0) && (act < 0.0))
			{ act = 0; }
		if ((ACTIVITY_THRESHOLD_1) && (act > 1.0))
			{ act = 1; }
		variable_features._activity.value = act;
	#else
		variable_features._activity.value = act;
	#endif
}

void Element::Set_mut_activity(const double muta)
{
	assert (muta >= 0.0);
	assert (muta <= 1.0);
	if ((muta < low_limit_rates) && (muta != 0.0))
	{
		cerr << "Warning: act_mut_rate < LOW_LIMIT_RATES!" << endl;
		fixed_features._act_mut.value = 0.0;
	} else {
		fixed_features._act_mut.value = muta;
	}	
}

void Element::Set_sd_activity(const double actmdec)
{
	assert (actmdec >= 0.0);
	fixed_features._act_sd.value = actmdec;
}

// ##### Regulation factor

void Element::Set_regulation_factor(const double regfac)
{
	if (regfac < 0.0)
		variable_features._regulation_factor.value = 0.0;
	else 
		variable_features._regulation_factor.value = regfac;	
}

void Element::Set_mut_reg(const double regmut)
{
	assert (regmut >= 0.0);
	assert (regmut <= 1.0);
	fixed_features._reg_mut.value = regmut;
}

void Element::Set_sd_reg(const double regsd)
{
	assert (regsd >= 0.0);
	fixed_features._reg_sd.value = regsd;
}

/***************************************************************************/
/*                CONSTRUCTORS AND DESTRUCTORS                           */
/***************************************************************************/

Element::Element()
{
  _id = _idref++;	
  fixed_features._family.value = 0;
}

Element::Element(const Element& et)
{  
  Copy(et);
}

Element & Element::operator = (const Element & et)
{
	Clean();
	Copy(et);
	return *this;
}



Element::~Element(void) {
  Clean();
}


void Element::Clean(void)
{
  // Nothing to do : only scalar variables
}

void Element::Copy(const Element & modele)
{

  /* copies the element "modele" in this
     Warning : this is completely erased */
	
	_id = _idref++; // The reference is not the same!
	fixed_features = modele.fixed_features;
	variable_features = modele.variable_features;
}

/*################################# Internal classes constructors ##########################*/

Element::FixedFeatures::FixedFeatures(void) : 
				_family		("family", 	0),
				_mean_fitness	("mean_fitness",0.0),
				_sd_fitness	("sd_fitness", 	0.0),
				_dup_mut	("dup_mut_rate",0.0),
				_dup_sd		("dup_mut_sd", 	0.0),
				_del_mut	("del_mut_rate",0.0), 
				_del_sd		("del_mut_sd", 	0.0),
				_act_mut	("act_mut_rate",0.0),
				_act_sd		("act_mut_sd", 	0.0), 
				_reg_mut	("reg_mut_rate",0.0),
				_reg_sd		("reg_mut_sd", 	0.0)
{

}

Element::VariableFeatures::VariableFeatures(void) :
				_position	("locus", 	0.0), 
				_chromo		("chromosome", 	0),
				_fitness	("fitness", 	0.0),
				_duplication_rate("dup_rate", 	0.0),
				_deletion_rate	("del_rate", 	0.0),
				_activity	("activity", 	1.0),
				_regulation_factor("regulation",0.0)
{

}

/*############################### ActivityMatrix ################################*/

ActivityMatrix * ActivityMatrix::instance = NULL;

ActivityMatrix::ActivityMatrix(int matrix)
{
	if ((matrix < 0) || (matrix > 3))
		throw e_MatrixIndex(matrix);

	double revert_rate = 0.0;

	if ((matrix == 0)||(matrix == 2))
	{
		if (matrix == 2) revert_rate = 0.1;

		// No regulation, mutations only Act -> No
			transitionMatrix[Regulator][Regulator] 	= 1.0; 	
			transitionMatrix[Regulator][NoActivity]	= 0.0;
			transitionMatrix[Regulator][Active] 	= 0.0;	
			transitionMatrix[NoActivity][Regulator]	= 0.0; 	
			transitionMatrix[NoActivity][NoActivity]= 1.0 - revert_rate;	
			transitionMatrix[NoActivity][Active]	= revert_rate;
			transitionMatrix[Active][Regulator] 	= 0.0; 
			transitionMatrix[Active][NoActivity]	= 1.0;
			transitionMatrix[Active][Active] 	= 0.0;	
	}
	if ((matrix == 1)||(matrix == 3))
	{
		if (matrix == 3) revert_rate = 0.1;
		// Regulation allowed
			transitionMatrix[Regulator][Regulator] 	= 0.0; 	
			transitionMatrix[Regulator][NoActivity] = 0.5;
			transitionMatrix[Regulator][Active] 	= 0.5;	
			transitionMatrix[NoActivity][Regulator] = revert_rate; 	
			transitionMatrix[NoActivity][NoActivity]= 1.0 - 2.0*revert_rate;	
			transitionMatrix[NoActivity][Active] 	= revert_rate;
			transitionMatrix[Active][Regulator]	= 0.5; 
			transitionMatrix[Active][NoActivity] 	= 0.5;
			transitionMatrix[Active][Active] 	= 0.0;	
	}

	StatesValues[Regulator] = -1.0;
	StatesValues[NoActivity] = 0.0;
	StatesValues[Active] = 1.0;

	epsilon = 0.0001;
}

double ActivityMatrix::GetProba(double cur_activity, double new_activity)
{
	return transitionMatrix[act2enum(cur_activity)][act2enum(new_activity)];
}

double ActivityMatrix::RandomActivity(double activity, double random_double)
{
	double sum = 0.0;
	int enum_act = act2enum(activity);
	for (unsigned int state = 0; state < ActivityMatrix::NbrStates; state++)
	{
		sum += transitionMatrix[enum_act][state];
		if (random_double < sum) // < or <=?
			return enum2act(state);
	}
	throw e_SumProba(sum, random_double);
}

int ActivityMatrix::act2enum(double activity)
{
	for (unsigned int state = 0; state < NbrStates; state++)
	{
		if (TestValuePlusorMinusEpsilon(static_cast<int>(StatesValues[state]), activity))
			return state;
	}
	throw e_Activity(activity);
}

double ActivityMatrix::enum2act(int state)
{
	assert (state >= 0);
	assert (state < NbrStates);

	return StatesValues[state];
}

bool ActivityMatrix::TestValuePlusorMinusEpsilon(int value, double activity)
{
	if ((activity > static_cast<double>(value) - epsilon) && (activity < static_cast<double>(value) + epsilon))
		return true;
	else
		return false;
}
