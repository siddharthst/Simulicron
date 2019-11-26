/***************************************************************************
 *            tools.h
 *
 *  Thu Jul 15 15:47:41 2004
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
 
#ifndef TOOLS_H
#define TOOLS_H

#include <limits>
#include <exception>
#include <string>

namespace Tools
{
	int extract_int(const std::string & value, 
					int min = std::numeric_limits<int>::min(),
					int max = std::numeric_limits<int>::max());
	double extract_double(const std::string & value, 
					double min = std::numeric_limits<double>::min(),
					double max = std::numeric_limits<double>::max());
	bool extract_bool(const std::string & value);

	template <class T>
	class NotInRange : public std::exception
	{
		public:
		NotInRange(T v, T i, T m)
			: val(v), min(i), max(m) { }
		~NotInRange() throw() { }
		const char * what() const throw()
		{ std::ostringstream o; o << "Value " << val << " is not in the range [" << min << ":" << max << "]." << std::endl; return o.str().c_str(); }
		private:
		T val, min, max;
	};
} //namespace Tools

#endif // TOOLS_H
