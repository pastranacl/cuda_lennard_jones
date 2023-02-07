/****************************************************************************

    This is part of cpputils
    Copyright (C) 2019  Cesar L. Pastrana

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

****************************************************************************/


#ifndef _CPPUTILS_
#define _CPPUTILS_

#include <string>
#include <sstream>


/*******************************************************/
/* Numeric functions                                   */
/*******************************************************/

int imin(int *v, int n);
int imax(int *v, int n);
double dmin(double *v, int n);
double dmax(double *v, int n);

int iargmin(int *v, int n);
int iargmax(int *v, int n);
int dargmin(double *v, int n);
int dargmax(double *v, int n);




/*******************************************************/
/* String Functions                                    */
/*******************************************************/
std::string strcat(std::string str1, std::string str2);
std::string int2str(int n);














#endif
