/*
 * STATISTICS.CPP
 *
 *  Created on: Jan 21, 2011
 *      Author: Pieter Chys
 */

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>


#include "statistics.h"
#include "matrix.h"

using namespace std;

/*===============================================================================
 *
 *  RMSD_TOTAL.CPP
 *
 *	Created on: Jan 21, 2011
 *      Author: Pieter Chys
 *
 *===============================================================================
 * RMSD_TOTAL.CPP calculates the root-mean-square-deviation (RMSD)
 * between the coordinates of two superimposed chain backbones.
 * This is a single double number that is returned. The chain
 * backbones coordinates are a N x 3 matrix with 3 columns (X,Y,Z).
 * The backbone coordinates are referenced by double pointers.
 *
 *
 *--------------------------------------------------------------------------------
 *
 *  INPUT: - double** pointer to initial coordinates co_i (N x 3)
 *  	   - double** pointer to active coordinates co (N x 3)
 *  	   - int number which gives the second dimension size of
 *  		 the coordinates matrices = number of atoms in backbone
 *
 *-------------------------------------------------------------------------------
 *
 *  OUTPUT: - rmsd value is returned as a double number
 *
 *
 *
 *-------------------------------------------------------------------------------
 *
 *  DEPENDENCIES: <cmath> needed
 *
 *-------------------------------------------------------------------------------
 *
 *	REMARKS:	- should you divide by N or (N-1) in the calculation?
 *				  at present -> N-1 because largest value....
 *
 *===============================================================================
 */

double rmsd_total(double** co_i,double** co,int number_atoms_loop) 
{
	double sum_total = 0; 
    for(int i = 0; i < number_atoms_loop; i++) 
	{
		double sum = 0.0;
		for(int j = 0; j < 3; j++) 
		{
			sum += pow(co[i][j]-co_i[i][j],2);
		}

		sum_total += sum;
	}

	return sqrt((sum_total/number_atoms_loop));				// divide by N or (N-1)?
}


double total_rmsd(double** co_i,double** co,double **pco_i,double **pco,int number_atoms_loop)	//:TODO: Well enough debugged?
{
	int number; 					/*total number of atoms including the oxygen atoms	of the peptide backbone (interating and increasing)	*/
	double sum_total; 				/*total sum of RMSD contributions (iterating and increasing)											*/
	double sum;						/*sub sum floating point number																			*/

	sum_total = 0.0;				/*initialization for each computation									*/
	number    = number_atoms_loop;	/*initialization: number of atoms without oxygen atoms included			*/

	//DEBUGGING PROTOCOL
	//printf("\n");
	//printf("RMSD computation on total loop\n" );
	//printf("\n");
	//printf("Total number of atoms without oxygen atoms: %i\n" ,number);
	//printf("\n");
	//printf("Final number should be: %i\n",number + (number_atoms_loop/3) - 1 );

	for(int i = 0; i < number_atoms_loop; i++) 
	{
		//printf("Current atom index: %i\n",i);
		sum = 0.0;					/*initialization to zero*/

		for(int j = 0; j < 3; j++) 
		{
			sum += pow( co[i][j] - co_i[i][j], 2 );
		}

		if ( (i + 1 ) % 3 == 0 && i != (number_atoms_loop - 1) )   /*select oxygen atoms and not the last one -> can rotate freely because of next bond to get into correct position*/
		{
			for(int j = 0; j < 3; j++) 
			{
				sum += pow(pco[i][j] - pco_i[i][j],2);
			}

			//printf("Current oxygen coordinate: %8.3f %8.3f %8.3f\n",pco[i][0],pco[i][1],pco[i][2]);
			//printf("Native oxygen coordinate %8.3f %8.3f %8.3f\n",pco_i[i][0],pco_i[i][1],pco_i[i][2]);

			number++;
		}

		sum_total += sum;
	}

	//DEBUGGING PROTOCOL
	//printf("\n");
	//printf("Total number of atoms counted: %i\n",number);


	return sqrt((sum_total/number));
}

// Compute the RMSD of the loop backbone (N, CA, C)
double total_rmsd_co(double** co_i,double** co,int number_atoms_loop)
{
	int number; 					/*total number of atoms including the oxygen atoms	of the peptide backbone (interating and increasing)	*/
	double sum_total; 				/*total sum of RMSD contributions (iterating and increasing)											*/
	double sum;						/*sub sum floating point number																			*/

	sum_total = 0.0;				/*initialization for each computation									*/
	number    = number_atoms_loop;	/*initialization: number of atoms without oxygen atoms included			*/

	for(int i = 0; i < number_atoms_loop; i++)
	{
		//printf("Current atom index: %i\n",i);
		sum = 0.0;					/*initialization to zero*/

		for(int j = 0; j < 3; j++)
		{
			sum += pow( co[i][j] - co_i[i][j], 2 );
		}

//		if ( (i + 1 ) % 3 == 0 && i != (number_atoms_loop - 1) )   /*select oxygen atoms and not the last one -> can rotate freely because of next bond to get into correct position*/
//		{
//			for(int j = 0; j < 3; j++)
//			{
//				sum += pow(pco[i][j] - pco_i[i][j],2);
//			}
//
//			//printf("Current oxygen coordinate: %8.3f %8.3f %8.3f\n",pco[i][0],pco[i][1],pco[i][2]);
//			//printf("Native oxygen coordinate %8.3f %8.3f %8.3f\n",pco_i[i][0],pco_i[i][1],pco_i[i][2]);
//
//			number++;
//		}

		sum_total += sum;
	}

	//DEBUGGING PROTOCOL
	//printf("\n");
	//printf("Total number of atoms counted: %i\n",number);


	return sqrt((sum_total/number));
}



/*===============================================================================
 *
 *  RMSD_ATOM.CPP
 *
 *	Created on: Jan 21, 2011
 *      Author: Pieter Chys
 *
 *===============================================================================
 * RMSD_ATOM.CPP calculates the root-mean-square-deviation (RMSD)
 * between the coordinates of two superimposed chain backbones atom
 * per atom and is returned in a preallocated single array vector
 * with the same size as the number of atoms in the backbone.
 * Backbones coordinates are in a N x 3 matrix with 3 columns (X,Y,Z).
 * The backbone coordinates are referenced by double pointers.
 *
 *
 *--------------------------------------------------------------------------------
 *
 *  INPUT: - double** pointer to initial coordinates co_i (N x 3)
 *  	   - double** pointer to active coordinates co (N x 3)
 *  	   - int number which gives the second dimension size of
 *  		 the coordinates matrices = number of atoms in backbone
 *  	   - double* pointer preallocated to store (number of atoms
 *  		 in backbone)  coordinates
 *
 *-------------------------------------------------------------------------------
 *
 *  OUTPUT: - rmsd values are returned as a double number in array[number_atoms_loop]
 *
 *-------------------------------------------------------------------------------
 *
 *  DEPENDENCIES: <cmath> needed
 *
 *-------------------------------------------------------------------------------
 *
 *	REMARKS:	- should you divide by N or (N-1) in the calculation?
 *				  at present -> N-1 because largest value....
 *
 *===============================================================================
 */


void rmsd_atom(double** co_i,double** co,int number_atoms_loop,double* rmsd) {
	for(int i = 0; i < number_atoms_loop; i++) {
		rmsd[i] = 0.0;
		for(int j = 0; j < 3; j++) {
			rmsd[i] += pow(co[i][j]-co_i[i][j],2);
		}

		rmsd[i] = sqrt(rmsd[i]);
	}
}

/*===============================================================================
 *
 *  RMSD_SQUARED.CPP
 *
 *	Created on: Jan 21, 2011
 *      Author: Pieter Chys
 *
 *===============================================================================
 * RMSD_SQUARED.CPP calculates the root-mean-square-deviation (RMSD)
 * between the coordinates of two superimposed chain backbones between
 * any possible pair of atoms that exists.  This resembles a 2D
 * correlation like protocol and results in a 2-dimensional rmsd
 * profile.  Only half of this profile is unique and thus two
 * triangles exist which are equivalent. The result contains
 *(N x N) rmsd values in which N represents the number of atoms
 * in the considered loop backbone. The chain
 * backbones coordinates are a N x 3 matrix with 3 columns (X,Y,Z).
 * The backbone coordinates are referenced by double pointers.
 *
 *
 *--------------------------------------------------------------------------------
 *
 *  INPUT: - double** pointer to initial coordinates co_i (N x 3)
 *  	   - double** pointer to active coordinates co (N x 3)
 *  	   - int number which gives the second dimension size of
 *  		 the coordinates matrices = number of atoms in backbone
 *  	   - double** pointer to preallocated rmsd matrix with size
 *  		 N x N with N = number_atoms_loop
 *
 *-------------------------------------------------------------------------------
 *
 *  OUTPUT: - rmsd values are returned as a double numbers in the matrix
 *
 *-------------------------------------------------------------------------------
 *
 *  DEPENDENCIES: <cmath> needed
 *
 *-------------------------------------------------------------------------------
 *
 *	REMARKS:	- should you divide by N or (N-1) in the calculation?
 *				  at present -> N-1 because largest value....
 *
 *===============================================================================
 */


void rmsd_squared(double** co_i,double** co,int number_atoms_loop,double** rmsd) {
	for(int i = 0; i < number_atoms_loop; i++) {
		for(int j = 0; j > i && j < number_atoms_loop; j++) {
			double sum = 0.0;
			for(int k = 0; k < 3; k++) {
				sum += pow(co[i][k]-co_i[i][k],2);
			}

			rmsd[i][j] = sqrt(sum);
		}
	}
}

/*===============================================================================
 *
 *  MAX.CPP
 *
 *	Created on: February 11, 2011
 *      Author: Pieter Chys
 *
 *===============================================================================
 *
 * MAX.CPP finds the largest value in a set of data points. It operates
 * on data of type double.
 *
 *
 *--------------------------------------------------------------------------------
 *
 *  INPUT: - double *pointer to the data points (DATA)
 *  	   - int size is the number of elements in DATA
 *  	   - int *pointer to store the position of the minimum element (in vector)
 *
 *-------------------------------------------------------------------------------
 *
 *  OUTPUT: - double which is the maximum value found
 *  		- dereferenced marker contains position maximum value (C++ style)
 *
 *-------------------------------------------------------------------------------
 *
 *  DEPENDENCIES:
 *
 *-------------------------------------------------------------------------------
 *
 *	REMARKS:
 *
 *===============================================================================
 */

double max(double *data,int size, int *marker) {
	double 	max = data[0];
	*marker	=	0;
	for (int i = 0; i < size; i++) {
		if (data[i] > max) {
			max = data[i];
			*marker =	i;
		} else {	}
	}

	return max;

}

/*===============================================================================
 *
 *  MIN.CPP
 *
 *	Created on: February 11, 2011
 *      Author: Pieter Chys
 *
 *===============================================================================
 *
 * MIN.CPP finds the smallest value in a set of data points. It operates
 * on data of type double.
 *
 *
 *--------------------------------------------------------------------------------
 *
 *  INPUT: - double *pointer to the data points (DATA)
 *  	   - int size is the number of elements in DATA
 *  	   - int *pointer to store the position of the minimum element (in vector)
 *
 *-------------------------------------------------------------------------------
 *
 *  OUTPUT: - double which is the minimum value found
 *  		- dereferenced marker contains position minimum value (C++ style)
 *
 *-------------------------------------------------------------------------------
 *
 *  DEPENDENCIES:
 *
 *-------------------------------------------------------------------------------
 *
 *	REMARKS:
 *
 *===============================================================================
 */

double min(double *data,int size,int *marker) {
	double 	min = data[0];
	*marker = 0;
	for (int i = 0; i < size; i++) {
		if (data[i] < min) {
			min = data[i];
			*marker = i;
		} else { 	}
	}

	return min;
}

/*===============================================================================
 *
 *  AVG.CPP
 *
 *	Created on: March 2, 2011
 *      Author: Pieter Chys
 *
 *===============================================================================
 *
 * AVG.CPP calculates the average value of a vector array
 *
 *
 *--------------------------------------------------------------------------------
 *
 *  INPUT: - double *pointer to the data points (DATA)
 *  	   - int size is the number of elements in DATA
 *
 *-------------------------------------------------------------------------------
 *
 *  OUTPUT: - double which is the calculated average
 *-------------------------------------------------------------------------------
 *
 *  DEPENDENCIES:
 *
 *-------------------------------------------------------------------------------
 *
 *	REMARKS:
 *
 *===============================================================================
 */

double avg(double *data,int size) {
	double 	sum = 0.0;
	for (int i = 0; i < size; i++) {
		sum += data[i];
	}

	return sum/size;
}

/*===============================================================================
 *
 *  STDEV.CPP
 *
 *	Created on: March 2, 2011
 *      Author: Pieter Chys
 *
 *===============================================================================
 *
 * STDEV.CPP calculates the standard deviation on a given data set
 *
 *--------------------------------------------------------------------------------
 *
 *  INPUT: - double *pointer to the data points (DATA)
 *  	   - int size is the number of elements in DATA
 *
 *-------------------------------------------------------------------------------
 *
 *  OUTPUT: - double which is the calculated standard deviation
 *-------------------------------------------------------------------------------
 *
 *  DEPENDENCIES:
 *
 *-------------------------------------------------------------------------------
 *
 *	REMARKS:
 *
 *===============================================================================
 */

double stdev(double *data, int size) {
	double sum = 0.0;
	double mean;
	mean     = avg(data,size);
	for (int i = 0; i < size; i++) {
		sum  	+=  pow((double) ( data[i] - mean),2);
	}
	return  sqrt((double) (sum/(size - 1)));
}

/*===============================================================================
 *
 *  SKEWNESS.CPP
 *
 *	Created on: March 2, 2011
 *      Author: Pieter Chys
 *
 *===============================================================================
 *
 * SKEWNESS.CPP calculates the skewness of an experimental data set
 * as defined by standard statistical analysis.  It uses the third
 * moments of the statistics
 *
 *--------------------------------------------------------------------------------
 *
 *  INPUT: - double *pointer to the data points (DATA)
 *  	   - int size is the number of elements in DATA
 *
 *-------------------------------------------------------------------------------
 *
 *  OUTPUT: - double which is the calculated skewness of the experimental
 *  		  data
 *-------------------------------------------------------------------------------
 *
 *  DEPENDENCIES:
 *
 *-------------------------------------------------------------------------------
 *
 *	REMARKS:
 *
 *===============================================================================
 */

double skewness(double *data, int size) {
	double m3 = 0.0;
	double m2 = 0.0;
	double mean;
	mean     = avg(data,size);
	for (int i = 0; i < size; i++) {
		m3 	+=  pow((double) ( data[i] - mean),3);
		m2  	+=  pow((double) ( data[i] - mean),2);
	}
	m3 /= size;
	m2 /= size;
	m2 = pow(sqrt(m2),3);
	return  m3/m2;
}

/*===============================================================================
 *
 *  KURTOSIS.CPP
 *
 *	Created on: March 2, 2011
 *      Author: Pieter Chys
 *
 *===============================================================================
 *
 * KURTOSIS.CPP calculates the kurtosis of an experimental data set
 * as defined by standard statistical analysis.  It uses the fourth
 * moments of the statistics
 *
 *--------------------------------------------------------------------------------
 *
 *  INPUT: - double *pointer to the data points (DATA)
 *  	   - int size is the number of elements in DATA
 *
 *-------------------------------------------------------------------------------
 *
 *  OUTPUT: - double which is the calculated skewness of the experimental
 *  		  data
 *-------------------------------------------------------------------------------
 *
 *  DEPENDENCIES:
 *
 *-------------------------------------------------------------------------------
 *
 *	REMARKS:
 *
 *===============================================================================
 */


double kurtosis(double *data, int size) {
	double m4 = 0.0;
	double m2 = 0.0;
	double mean;
	mean     = avg(data,size);
	for (int i = 0; i < size; i++) {
		m4 	+=  pow((double) ( data[i] - mean),4);
		m2  	+=  pow((double) ( data[i] - mean),2);
	}
	m4 /= size;
	m2 /= size;
	m2 = pow(m2,2);
	return  m4/m2;
}

/*===============================================================================
 *
 *  KURTOSIS.CPP
 *
 *	Created on: March 2, 2011
 *      Author: Pieter Chys
 *
 *===============================================================================
 *
 * KURTOSIS.CPP calculates the kurtosis of an experimental data set
 * as defined by standard statistical analysis.  It uses the fourth
 * moments of the statistics
 *
 *--------------------------------------------------------------------------------
 *
 *  INPUT: - double *pointer to the data points (DATA)
 *  	   - int size is the number of elements in DATA
 *
 *-------------------------------------------------------------------------------
 *
 *  OUTPUT: - double which is the calculated skewness of the experimental
 *  		  data
 *-------------------------------------------------------------------------------
 *
 *  DEPENDENCIES:
 *
 *-------------------------------------------------------------------------------
 *
 *	REMARKS:
 *
 *===============================================================================
 */


void makehistogram(double *data,int size,double min,double max,double *xarray,double *yarray, int size_hist){
	double step;
	double xlim[size_hist+1];
	step	=	(max - min)/ size_hist;
	//printf("%f\n",max);
	//printf("%f\n",min);
	//printf("%f\n",step);

	for(int i = 0; i < size_hist; i++) {
		xarray[i] 	=	min  + (step/2) + i*step;
		//printf("%f\n",xarray[i]);
	}

	for(int i = 0; i < size_hist+1; i++) {
		xlim[i] 	=	min  - (step/2) + (i+1)*step;
	}

	for(int i = 0; i < size_hist; i++) {
		double counter = 0.0;
		for (int j = 0; j < size; j++) {
			if (data[j] > xlim[i] && data[j] < xlim[i+1]) {
				counter += 1.0;
			}
			/*
			} else if (i == 0  && data[j] < xarray[i+1]) {
				counter +=1.0;
			} else if (i == (size_hist)  && data[j] > xarray[i-1]) {
				counter += 1.0;
			}
			*/
		}

		yarray[i] = counter/size;

	}
}

double randomangle(void)
{
	return ( 1.0  - 2.0*(double (rand())/double (RAND_MAX)) )* M_PI;
	//return ( 2*(double (rand())/double (RAND_MAX)) -1)* M_PI;	// OLD CODE BUT NOT CORRECT IN PRINCIPLE
}

