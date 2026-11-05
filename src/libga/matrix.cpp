/*
 * matrix.cpp
 *
 *  Created on: November 18, 2011
 *      Author: Pieter Chys
 */

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <vector>
#include <cmath>

#include "matrix.h"


Matrix::Matrix(int i,int j) {
	// TODO Auto-generated constructor stub

	size_x = i;
	size_y = j;


	//	allocate coordinate pointer **co
	el = (double**) malloc(sizeof(double *)* size_x);

	for(int k = 0; k < size_x; k++) {
		el[k] = (double *) malloc(sizeof(double) * size_y);
	}
	int temp;	/*advice parallellisation routine*/
	temp = size_y;

	for(int k = 0; k < size_x; k++) {
		//#pragma ivdep	/*preprocessor directive parallellisation*/
		for(int l = 0; l < temp; l++) {
			el[k][l] =  0.0;
		}
	}

	/*for(int k = 0; k < size_x; k++) {
		for(int l = 0; l < size_y; l++) {
			el[k][l] =  0.0;
		}
	}*/
}

Matrix::~Matrix() {
	// TODO Auto-generated destructor stub

	//	Free pointer **el
	for(int i = 0; i < size_x; i++)
	{
		free(el[i]);
	}
	free(el);
}

void Matrix::copy(int x,int y,double **input,int diff_x,int diff_y)
{
	for(int i = 0; i < diff_x; i++)
	{
		//#pragma ivdep
		for(int j = 0; j < diff_y; j++)
		{
			el[x + i][y + j] 		=  input[i][j];
		}
	}
}

void Matrix::dihedral(int i,int j){	//	some problems possible
	if (el[i][j] > M_PI)
	{										//	correct dihedral value according to conventions
		el[i][j] = el[i][j] - 2*M_PI;
	} else if (el[i][j] < -M_PI)
	{
		el[i][j] = el[i][j] + 2*M_PI;
	}
}
