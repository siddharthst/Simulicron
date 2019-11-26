/***************************************************************************
 *            tools.cpp
 *
 *  Thu Jul 15 15:56:53 2004
 *  Copyright  2004  lerouzic
 *  lerouzic@pge.cnrs-gif.fr
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "tools.h"

#include <sstream>

using namespace std;

namespace Tools
{
	
int extract_int(const std::string & value, int min, int max)
{
	istringstream i(value);
	int result;
	i >> result;
	if ((result < min) || (result > max))
		throw NotInRange<int>(result, min, max);
	return result;
}

double extract_double(const std::string & value, double min, double max)
{
	istringstream i(value);
	double result;
	i >> result;
	if ((result < min) || (result > max))
		throw NotInRange<double>(result, min, max);
	return result;
}
	
bool extract_bool(const std::string & value)
{
	string upvalue;
	for (unsigned int c = 0; c < value.size(); c++)
		upvalue += toupper(value[c]);
	if (   (upvalue == "Y")
		|| (upvalue == "YES")
		|| (upvalue == "O")
		|| (upvalue == "OUI")
		|| (upvalue == "TRUE")
		|| (upvalue == "OK")
		|| (upvalue == "1"))
		return true;
	else
		return false;
	// no control
}

} // namespace Tools
