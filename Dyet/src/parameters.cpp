/***************************************************************************
                          parameters.cpp  -  description
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

#include <fstream>
 
#include "parameters.h"
#include "tools.h"

using namespace std;

Parameters::Parameters(const string & filename)
{
	
	/* The file must be perfect here */
	
	ifstream my_file;
	my_file.open(filename.c_str(), ios::in);
	
	char line[256];
	while (my_file.getline(line, 256))
	{
		istringstream i(line);
		string label;
		i >> label;
		try
		{
			if (label[0] != '#')
			{
				string value;
				i >> value;
				SetValue(label, value);
			}
		} catch (e_unused & e) { /* cerr << e.what() << endl; */ }
	}
}

void Parameters::SetValue(const string & label, const string & value)
{
	if (label == "GEN_NB")
	{
		generations = Tools::extract_int(value, 0);
	} else
	if (label == "GEN_ST")
	{
		if (generations < 0) throw e_order(label);
		steps = Tools::extract_int(value, 0, generations);
	} else
	if (label == "BURSTS")
	{
		bursts = Tools::extract_bool(value);
	} else
	if (label == "BU_TYP")
	{
		burst_type = Tools::extract_int(value, 0, max_type);
	} else
	if (label == "BU_VAL")
	{
		if (burst_type < 0) throw e_order(label);
		switch (burst_type)
		{
			case 0:
			case 2:
				throw e_unused(label);
			break;
			case 1:
				burst_value = Tools::extract_double(value, 0.0, 1.0);
			break;
			case 3:
			case 4:
				burst_value = static_cast<double>(Tools::extract_int(value, 0));
			break;
			default:
				throw e_order(label);
			break;
		}
	}
	if (label == "ACT_MA")
	{
		activity_matrix = Tools::extract_int(value, 0);
	}
}
