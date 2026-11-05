#ifndef GUARD_own_math_h
#define GUARD_own_math_h

//own_math.h
#include <cmath>

/*vector functions*/
void quatrotate(double **co,double **co1,double **co2,int nratoms, double *quaternion);
void quatrotate_cg(double **co,double **co1,double **co2,int nratoms, double *quaternion, int mode);

void quatproduct(double *qt2, double *qt1, double *quaternion);

double distance(double *,double *);
double distance2(double *,double *);
bool isnear_fast(double *p1, double *p2, double clash_factor);

double dotprod(double *,
			   double *);

void crossprod(double *,
			   double *,
			   double *);

double norm(double *);

void normalize(double *, double *);

void normalise(double *);

void diffvector(double **, int, int,
				double *);

void diffvectorseq(double **,
		 	 	   int,
		 	 	   double *);

void flipud(double **,
			int);

/*matrix functions*/
void transpose(double **,
			   double **);

double det33(double **);

void matrixXvector(double **,
				   double *,
				   double *);

void matrixXvector(double [][3],
		 	 	   double *);

void matrixXmatrix(double **,
				   double **);

void mapmatrix(double **,
			   double **);

void matrixXmatrix(double **,
				   double **,
				   double **);	//	other functions

//void 	matrixXmatrix(double **,double **, double **,double **);
void matrixprod1(double **,
				 double **);

void matrixprod2(double **,
				 double **);

void min(double *);

void randnvector(double *);

void subtract(double *,
			  double *,
			  double *);

void add(double *,
		 double *,
		 double *);



#endif /* OWN_MATH_H_ */

/*
 * own_math.h
 *
 *  Created on: Dec 7, 2010
 *      Author: Pieter Chys
 */
