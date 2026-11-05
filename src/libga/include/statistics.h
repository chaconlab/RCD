#ifndef GUARD_statistics_h_
#define GUARD__statistics_h_

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "matrix.h"


using namespace std;

double total_rmsd(double** co_i,double** co,double **pco_i,double **pco,int number_atoms_loop);

// Compute the RMSD of the loop backbone (N, CA, C)
double total_rmsd_co(double** co_i,double** co,int number_atoms_loop);


double rmsd_total(double **,
				  double **,
				  int);

void rmsd_atom(double **,
		       double **,
		       int,
		       double *);

void rmsd_squared(double **,
		          double **,
		          int,
		          double **);

double max(double *,
		   int,
		   int *);

double min(double *,
		   int,
		   int *);

double avg(double *,
		   int);

double stdev(double *,
		     int);

double skewness(double *,
		        int);

double kurtosis(double *,
		        int);

void makehistogram(double *,
		           int ,
		           double ,
		           double ,
		           double *,
		           double *,
		           int);

double randomangle(void);

#endif /* GUARD_statistics_h_ */

/*
 * STATISTICS.H
 *
 *  Created on: Jan 21, 2011
 *      Author: Pieter Chys
 */

