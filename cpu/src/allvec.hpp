/****************************************************************************

    This is part of allvec
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


#ifndef _ALLVEC_
#define _ALLVEC_

#include <iostream>

using namespace std;


int *ivector(int M);
double *dvector(int M);

void val_ivector(double *ivec, double val, int N);
void val_dvector(double *dvec, double val, int N);

void zeros_ivector(int *ivec, int M);
void zeros_dvector(double *dvec, int M);

void free_ivector(int *v);
void free_dvector(double *v);

/*****************************************************/

int **imatrix(int M, int N);
double **dmatrix(int M, int N);

void val_imatrix(int **imtx, int val, int M, int N);
void val_dmatrix(double **dmtx, double val, int M, int N);

void zeros_imatrix(int **imtx, int M, int N);
void zeros_dmatrix(double **dmtx, int M, int N);

void free_imatrix(int **dmtx, int M);
void free_dmatrix(double **dmtx, int M);

/*****************************************************/

int ***i3matrix (int M, int N, int K);
double ***d3matrix (int M, int N, int K);

void val_i3matrix(int ***imtx3, int val, int M, int N, int K);
void val_d3matrix(double ***dmtx3, double val, int M, int N, int K);
void zeros_i3matrix(int ***imtx3, int M, int N, int K);
void zeros_d3matrix(double ***dmtx3, int M, int N, int K);

void free_i3matrix(int ***mtx, int M, int N);
void free_d3matrix(double ***mtx, int M, int N);


#endif
