/****************************************************************************

    This is part of bdvesenz
    Copyright (C) 2021  Cesar L. Pastrana

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


#include <iostream>
#include <fstream>
#include <sstream> 
#include <algorithm>
#include <string>
#include <ctime>
#include <vector>
#include <math.h>
#include <unistd.h>

#ifndef PI
#define PI 3.141592653589793238462643383279502884197169
#endif

#define   AVOGADRO_PREFACTOR              6.0221415              // Avogadro number without the 10^23 factor

#define INIT_COORDS_FNAME       "enzyme_coords.dat"              //  Path of the initial file to load enzyme coardinates
#define INIT_COORDS_SAVE       "enzyme_coords_used.xyz"          // Generated or loaded coordiantes
#define PARAMETERS_FILE         "params.conf"                    //  Parameters file




/*********************************   Namespaces and type definitions   *******************************/



// Standard namespace
using std::cout;
using std::cin;
using std::endl;
using std::vector;
using std::ofstream; 
using std::ifstream; 
using std::string;
using std::getline;
using std::stringstream;
using std::nothrow; // Do not send exception but a null pointer in new


typedef vector<vector<double>> dvector_vl;




void enzyme_positions(struct enzymes *enz, struct tube *tb);
double det_dmin(struct enzymes *enz, int id);

