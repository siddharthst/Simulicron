/***************************************************************************
                          element.h  -  description
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

   
#ifndef ELEMENT_H
#define ELEMENT_H

#include <assert.h>
#include <iostream>
#include <sstream>
#include <string>

#include "dyet.h" 
#include "base.h"

#define low_limit_rates 0.0000001

#define MUTATION_FITNESS_DURING_TRANSPOSITION true
#define MUTATION_FITNESS_DURING_REPLICATION false

#define MUTATION_DUPRATE_DURING_TRANSPOSITION false
#define MUTATION_DUPRATE_DURING_REPLICATION true

#define MUTATION_DELRATE_DURING_TRANSPOSITION false
#define MUTATION_DELRATE_DURING_REPLICATION true

#define MUTATION_ACTIVITY_DURING_TRANSPOSITION false
#define MUTATION_ACTIVITY_DURING_REPLICATION true

#define MUTATION_REGULATION_DURING_TRANSPOSITION false
#define MUTATION_REGULATION_DURING_REPLICATION true

#define ACTIVITY_CONT true
#define ACTIVITY_DISC !ACTIVITY_CONT

// Activity_discrete options


// Activity continuous options
#define ACTIVITY_MUTATION_UP true
#define ACTIVITY_THRESHOLD_1 true
#define ACTIVITY_THRESHOLD_0 true
#define ACTIVITY_ABSORBING_STATE_0 true



class Ptr2Element; 

template <typename T> struct ElementFeature
{
	ElementFeature(const std::string & s, const T & v) : xml_tag(s), value(v) { }
	std::string xml_tag;
	T value;
};

class Element : public BaseBiolObject
{
	// Import and Export become easier...
	friend class ImportSystem;
	friend class ExportSystem;
		
	public:
  	// Constructors and destructor
	Element();
  	Element(const Element&);
  	~Element();
  	Element & operator = (const Element &);
	
	// This is a conception error : Erase function might be in System instead of here
	static void Reset(void) 
		{Element::_idref = 0;}
  
  	Ptr2Element Duplication() const;
	Ptr2Element Mutation_replication(bool &) const;
	void Mutation_transposition();
		
	// Frequently needed features
	inline unsigned long int id(void) const {return _id;}
  	inline short int family(void) const {return fixed_features._family.value;}
	
  	inline double position(void) const {return variable_features._position.value;}
  	inline unsigned int chromo(void) const {return variable_features._chromo.value;}
		
  	inline double fitness(void) const {return variable_features._fitness.value;}	
	inline double mean_fitness(void) const {return fixed_features._mean_fitness.value;}
	inline double sd_fitness(void) const {return fixed_features._sd_fitness.value;}

  	inline double duplication_rate(void) const {return variable_features._duplication_rate.value;}
	inline double duplication_mutation_rate(void) const {return fixed_features._dup_mut.value;}
	inline double duplication_mutation_sd(void) const {return fixed_features._dup_sd.value;}
  
	inline double deletion_rate(void) const {return variable_features._deletion_rate.value;}
	inline double deletion_mutation_rate(void) const {return fixed_features._del_mut.value;}
	inline double deletion_mutation_sd(void) const {return fixed_features._del_sd.value;}

  	inline double real_activity(const long int copy_number) const 
		{ return (activity())/((static_cast<double>(copy_number - 1)) * regulation_factor() + 1.0); }
	
	inline double activity(void) const {return variable_features._activity.value;}
	inline double activity_mutation_rate(void) const {return fixed_features._act_mut.value;}
	inline double activity_mutation_sd(void) {return fixed_features._act_sd.value;}

	inline double regulation_factor(void) const {return variable_features._regulation_factor.value;}
	inline double regulation_mutation_rate(void) const {return fixed_features._reg_mut.value;}
	inline double regulation_mutation_sd(void) const {return fixed_features._reg_sd.value;}

  	// Comparison of the position of 2 elements
  	// No Ptr2Element : too long
  	inline static bool IsBefore(const Element* et1, const Element* et2)
  		{ return (et1->position() < et2->position());}
		
	protected:	
	// ******************** static variables *******************
  	// Absolute numbering of the element
  	static unsigned long int _idref;
	
	// **************** Element's data *************************
	unsigned long int _id;  // Element's Id 
	
	protected :
		class FixedFeatures {
			friend class Element;
			friend class ExportSystem; // does not modify the data
			friend class ImportSystem; // idem
			FixedFeatures();
			protected :
  				ElementFeature<short int> _family; 

				ElementFeature<double> _mean_fitness;
				ElementFeature<double> _sd_fitness;

				ElementFeature<double> _dup_mut;
				ElementFeature<double> _dup_sd;

				ElementFeature<double> _del_mut;
				ElementFeature<double> _del_sd;

				ElementFeature<double> _act_mut;
				ElementFeature<double> _act_sd;

				ElementFeature<double> _reg_mut;
				ElementFeature<double> _reg_sd;
		};
		class VariableFeatures {
			friend class Element;
			friend class ExportSystem; // only in read mode
			friend class ImportSystem; // only in read mode
			public:
			VariableFeatures();
			protected :
				ElementFeature<double> _position;
				ElementFeature<unsigned int> _chromo;
				ElementFeature<double> _fitness;
				ElementFeature<double> _duplication_rate;
				ElementFeature<double> _deletion_rate;
				ElementFeature<double> _activity;
				ElementFeature<double> _regulation_factor;
		};

	 
	FixedFeatures fixed_features;
	VariableFeatures variable_features;
		
  	// Random settings
	void SetRandomFitness(void);
	void SetRandomPosition(void);
	void SetRandomDupRate(void);
	void SetRandomDelRate(void);
	void SetRandomActivity(void);
	void SetRandomRegulation(void);
		
	// Not random settings
	void Set_family(const short int);
	void Set_position(const double);
	void Set_fitness(const double);
	void Set_mean_fitness(const double);
	void Set_dup_rate(const double);
	void Set_del_rate(const double);
	void Set_activity(double);  
	
	void Set_sd_fitness(const double);

	void Set_sd_dup_rate(const double);
	void Set_mut_dup_rate(const double);
	void Set_sd_del_rate(const double);
	void Set_mut_del_rate(const double);

	void Set_mut_activity(const double);
	void Set_sd_activity(const double);
	void Set_regulation_factor(const double);
	void Set_mut_reg(const double);
	void Set_sd_reg(const double);

	// *********************** Clean and copy ************************
	void Clean(void);
	void Copy(const Element &);
};


/* The following class is a "SmartPointer". It manages automatically
   the memory needed by the elements. Each element 'knows' how much times
   it is pointed, and is able to self-destroy when there is no more pointer 
   on it.
   Element* must just be replaced in the code by Ptr2Element.


   This class is stringly similar to the example given in the following book 
   (in french):
   "Pour mieux développer avec C++"
   Aurélien Géron, Fatmé Tawbi
   Dunod, Collection "Informatiques".
*/

class Ptr2Element
{
	private:
	static const bool details = false;
	struct Compteur {
		Element * ptr;
		long int compte;
		Compteur(Element * lePtr = NULL) : ptr(lePtr), compte(1) { }
    		~Compteur() { if (ptr != NULL)  { delete ptr;}}
		int operator++() { return ++compte;}
		int operator--() {
			if (--compte!=0) return compte;
				delete this;
      			return 0;}
  	};
	Compteur * leCompteur;
  
	public:
	Ptr2Element(void) :  leCompteur(new Compteur(NULL))
		{ std::cout << "Ne doit pas être appelé!!" << std::endl;}
	Ptr2Element(Element * lePtr) : leCompteur(new Compteur(lePtr)) 
		{ if (details) std::cout << "Ptr2Element(Element *) compt:" << leCompteur->compte << std::endl; }
	Ptr2Element (const Ptr2Element & smt) : leCompteur(smt.leCompteur) 
	{ 
		++(*leCompteur);
		 if (details) std::cout << "Ptr2Element(const Ptr2Element &) compt:" << leCompteur->compte << std::endl; 
  	}
	~Ptr2Element() 
	{
		--(*leCompteur);
		if (details) std::cout << "~Ptr2Element() compt:" <<  leCompteur->compte << std::endl; 	
	}
	Ptr2Element & operator =  (const Ptr2Element & smt) 
	{
		if (details) std::cout << "operator=(Ptr2Element&)"  << std::endl; 
		if (leCompteur != smt.leCompteur) {
      			-- (*leCompteur);
      			leCompteur = smt.leCompteur;
      			++(*leCompteur);
    		}
    		return *this;
  	}
  	Element * operator -> () const {return leCompteur->ptr;}
	Element & operator *  () const {return *(leCompteur->ptr);}
	Element * pointeur() const {return leCompteur->ptr;}

	inline static bool Est_avant(const Ptr2Element & et1, const Ptr2Element & et2)
		{ return et1->position() < et2->position(); }
		 //Element::Est_avant(et1.pointeur(), et2.pointeur());}
	friend inline bool operator<(const Ptr2Element &, const Ptr2Element &);
	friend inline bool operator ==(const Ptr2Element &, const Ptr2Element &);
}; 

bool operator < (const Ptr2Element & p1, const Ptr2Element & p2)
{
  return (p1->position() < p2->position());
}

bool operator == (const Ptr2Element & p1, const Ptr2Element & p2)
{
	return (p1.pointeur() == p2.pointeur());
}

struct IsIdBefore
{
	bool operator () (const Ptr2Element & p1, const Ptr2Element & p2) const
		{ return p1->id() < p2->id(); }
};

/*############ Discrete Activity ############*/

class ActivityMatrix
{
	private:
		static ActivityMatrix * instance;
		enum ActivityStates {Regulator, NoActivity, Active, NbrStates};
		// 3 states : -1, 0, 1
		double transitionMatrix[NbrStates][NbrStates];
		double StatesValues[NbrStates];
		double epsilon;

		int act2enum(double activity);
		double enum2act(int state);
		bool TestValuePlusorMinusEpsilon(int value, double activity);
		ActivityMatrix(int matrix);
		~ActivityMatrix() { }
	public:
		static ActivityMatrix * GetMatrix(unsigned int matrix)
			{ if (instance == NULL)	instance = new ActivityMatrix(matrix); return instance; }
		double GetProba(double cur_activity, double new_activity);
		double RandomActivity(double activity, double random_double);

		class ActivityMatrixException : public std::exception
		{
			public: 
				ActivityMatrixException() { }
				virtual const char * what() const throw() { return "Undefined ActivityMatrix Exception\n";}
				virtual ~ActivityMatrixException() throw () { }
		};
				

		class e_SumProba : public ActivityMatrixException
		{
			public: 
				e_SumProba(double s, double u) : sum(s), uni(u) { }
				const char * what(void) const throw() 
					{ std::ostringstream o; 
					o << "Error in ActivityMatrix : sum " << sum << "> 1 (unif=" << uni << ").\n";
					return (o.str()).c_str(); }
				~e_SumProba() throw() { }
			private: double sum, uni;
		};
		class e_Activity : public ActivityMatrixException
		{
			public:
				e_Activity(double a) : activity(a) { }
				const char * what(void) const throw()
					{ std::ostringstream o;
					o << "Error in ActivityMatrix : activity " << activity << " does not fit in accepted values.\n";
					return (o.str()).c_str(); }
				~e_Activity() throw() { }
			private:
				double activity;
		};
		class e_MatrixIndex : public ActivityMatrixException
		{
			public:
				e_MatrixIndex(int i) : index(i) { }
				const char * what() const throw()
					{ std::ostringstream o;
					o << "Error in Matrix Index : index = " << index << " is not a valid matrix index.\n";
					return (o.str()).c_str(); }
				~e_MatrixIndex() throw() { }
			private:
				int index;
		};
};


#endif
