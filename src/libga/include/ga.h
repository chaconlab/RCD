#ifndef GUARD_ga_h
#define GUARD_ga_h

//	ga.h header file
#include <cstdio>
#include <cstdlib>
#include <cmath>


/*standard geometric algebra*/
void geom11(double *,
			double *,
		    double *);

double geom11in(double *,
				double *,
				int);

void geom11out(double *,
			   double *,
			   double *);

void geom22(double *,
			double *,
			double *);

void  geom23(double *,
			 double *,
			 double *);

void  geom32(double *,
			 double *,
			 double *);

void geom23opt(double *,
			   double *,
			   double *);

void geom32opt(double *,
			   double *,
			   double *);

void geomi2(double *,
			double *);

/*spinor functions*/
void spinor(double *,
			double,
			double *);

void revspinor(double *,
			   double *);

void spinorproduct(double *,
				   double *,
				   double *);

void spinorproduct(double *,
				   double *);

void spinor3prod1(double *,
				  double *,
				  double *);

void spinor3prod2(double *,
			      double *,
				  double *);

void spinormat(double *,
			   double,
			   double *);

void mapspinor(double* bivector1,
			   double* bivector2);

void spinormatrix(double *,
				  double **);

void matrixspinor(double **,
				  double *);

//void 	spinorproduct(double*,double*,double *);
//void 	spinorproduct(double*,double*,double,double*);

/*quaternion functions*/
void quaternion(double *,
				double,
				double *);

/* void conjquat(double *,
			  double *);
			  */

void quatmatrix(double *,
		        double **);

void quatmatrix(double *,
		        double [][3]);



#endif /* GUARD_ga_h*/

/*==============================================================================
 * GA.H
 *
 *  Created on: Dec 7, 2010
 *      Author: Pieter Chys
 *===============================================================================
 */
