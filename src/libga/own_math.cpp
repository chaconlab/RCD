/*==============================================================================
 * OWN_MATH.CPP
 *
 *  Created on: Dec 10, 2010
 *      Author: Pieter Chys
 *
 *==============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "own_math.h"





void quatrotate(double **co,double **co1,double **co2,int nratoms, double *quaternion)
{
	//perror("1st sub-bypass\n");

	//ALLOCATE ARRAYS
	double **difference;												/*difference vectors with respect to the first atom*/
	difference = (double**) malloc( sizeof(double*) * nratoms );
	for(int i = 0; i < nratoms; i++ )
	{
		difference[i]	= (double*) malloc( sizeof(double) * 4);
	}

	for(int i = 0; i < nratoms; i++ )
	{	for(int j = 0; j < 4; j++)
		{
			difference[i][j]	= 0.0;
		}
	}

	double *between;													/*auxiliary quaternion*/
	between 	=	(double *) malloc ( sizeof(double) * 4);
	for(int i = 0; i < 4; i++ )
	{
		between[i] = 0.0;
	}

	//perror("2nd sub-bypass\n");

	double *cquaternion;												/*conjugate rotation quaternion*/
	cquaternion 	=	(double *) malloc ( sizeof(double) * 4);
	for(int i = 0; i < 4; i++ )
	{
		cquaternion[i] = 0.0;
	}




	//MAKE CONJUGATE QUATERNION
	cquaternion [0] 	=	 quaternion[0];
	cquaternion [1] 	=	-quaternion[1];
	cquaternion [2] 	=	-quaternion[2];
	cquaternion [3] 	=	-quaternion[3];




	//MAIN CHAIN COORDINATES
	//CONSTRUCT DIFFERENCE VECTORS
	for(int i = 0; i < nratoms ; i++ )
	{	for(int j = 0; j < 3; j++)
		{
			difference[i][j+1]	= co[i][j] - co[0][j];					/*difference with first atom*/
		}																/*shift assignment since multiplied as a quaternion*/
	}

	//EXECUTE THE ROTATION
	for (int i = 0; i < nratoms ; i++)
	{
		quatproduct(difference[i],cquaternion,between);
		quatproduct(quaternion,between,difference[i]);

		for(int i = 0; i < 4; i++ )
		{
			between[i] = 0.0;
		}
	}

	//CONSTRUCT COORDINATES AGAIN
	for(int i = 0; i < nratoms ; i++ )
	{	for(int j = 0; j < 3; j++)
		{
			co[i][j]	= difference[i][j+1] + co[0][j];					/*difference with first atom*/
		}																	/*shift assignment since constructed as a quaternion with 4 elements*/
	}


	//SIDE CHAIN 1 COORDINATES
	//CONSTRUCT DIFFERENCE VECTORS
	for(int i = 0; i < nratoms ; i++ )
	{	for(int j = 0; j < 3; j++)
		{
			difference[i][j+1]	= co1[i][j] - co[0][j];					/*difference with first atom*/
		}																/*shift assignment since multiplied as a quaternion*/
	}

	//EXECUTE THE ROTATION
	for (int i = 0; i < nratoms ; i++)
	{
		quatproduct(difference[i],cquaternion,between);
		quatproduct(quaternion,between,difference[i]);

		for(int i = 0; i < 4; i++ )
		{
			between[i] = 0.0;
		}
	}

	//CONSTRUCT COORDINATES AGAIN
	for(int i = 0; i < nratoms ; i++ )
	{	for(int j = 0; j < 3; j++)
		{
			co1[i][j]	= difference[i][j+1] + co[0][j];					/*difference with first atom*/
		}																	/*shift assignment since constructed as a quaternion with 4 elements*/
	}



	//SIDE CHAIN 1 COORDINATES
	//CONSTRUCT DIFFERENCE VECTORS
	for(int i = 0; i < nratoms ; i++ )
	{	for(int j = 0; j < 3; j++)
		{
			difference[i][j+1]	= co2[i][j] - co[0][j];					/*difference with first atom*/
		}																/*shift assignment since multiplied as a quaternion*/
	}

	//EXECUTE THE ROTATION
	for (int i = 0; i < nratoms ; i++)
	{
		quatproduct(difference[i],cquaternion,between);
		quatproduct(quaternion,between,difference[i]);

		for(int i = 0; i < 4; i++ )
		{
			between[i] = 0.0;
		}
	}

	//CONSTRUCT COORDINATES AGAIN
	for(int i = 0; i < nratoms ; i++ )
	{	for(int j = 0; j < 3; j++)
		{
			co2[i][j]	= difference[i][j+1] + co[0][j];					/*difference with first atom*/
		}																	/*shift assignment since constructed as a quaternion with 4 elements*/
	}





	//DEALLOCATE ARRAYS
	for (int i = 0; i < nratoms ; i++)
	{
		free(difference[i]);
	}
	free(difference);

	free(cquaternion);
	free(between);

}


void quatrotate_cg(double **co,double **co1,double **co2,int nratoms, double *quaternion, int mode)
{
	//perror("1st sub-bypass\n");

	//ALLOCATE ARRAYS
	double **difference;												/*difference vectors with respect to the first atom*/
	difference = (double**) malloc( sizeof(double*) * nratoms );
	for(int i = 0; i < nratoms; i++ )
	{
		difference[i]	= (double*) malloc( sizeof(double) * 4);
	}

	for(int i = 0; i < nratoms; i++ )
	{	for(int j = 0; j < 4; j++)
		{
			difference[i][j]	= 0.0;
		}
	}

	double *between;													/*auxiliary quaternion*/
	between 	=	(double *) malloc ( sizeof(double) * 4);
	for(int i = 0; i < 4; i++ )
	{
		between[i] = 0.0;
	}

	//perror("2nd sub-bypass\n");

	double *cquaternion;												/*conjugate rotation quaternion*/
	cquaternion 	=	(double *) malloc ( sizeof(double) * 4);
	for(int i = 0; i < 4; i++ )
	{
		cquaternion[i] = 0.0;
	}




	//MAKE CONJUGATE QUATERNION
	cquaternion [0] 	=	 quaternion[0];
	cquaternion [1] 	=	-quaternion[1];
	cquaternion [2] 	=	-quaternion[2];
	cquaternion [3] 	=	-quaternion[3];




	//MAIN CHAIN COORDINATES
	//CONSTRUCT DIFFERENCE VECTORS
	for(int i = 0; i < nratoms ; i++ )
	{	for(int j = 0; j < 3; j++)
		{
			difference[i][j+1]	= co[i][j] - co[0][j];					/*difference with first atom*/
		}																/*shift assignment since multiplied as a quaternion*/
	}

	//EXECUTE THE ROTATION
	for (int i = 0; i < nratoms ; i++)
	{
		quatproduct(difference[i],cquaternion,between);
		quatproduct(quaternion,between,difference[i]);

		for(int i = 0; i < 4; i++ )
		{
			between[i] = 0.0;
		}
	}

	//CONSTRUCT COORDINATES AGAIN
	for(int i = 0; i < nratoms ; i++ )
	{	for(int j = 0; j < 3; j++)
		{
			co[i][j]	= difference[i][j+1] + co[0][j];					/*difference with first atom*/
		}																	/*shift assignment since constructed as a quaternion with 4 elements*/
	}

	if(mode>0)
	{
		//SIDE CHAIN 1 COORDINATES
		//CONSTRUCT DIFFERENCE VECTORS
		for(int i = 0; i < nratoms ; i++ )
			for(int j = 0; j < 3; j++)
				difference[i][j+1]	= co1[i][j] - co[0][j];					/*difference with first atom*/
																		/*shift assignment since multiplied as a quaternion*/
		//EXECUTE THE ROTATION
		for (int i = 0; i < nratoms ; i++)
		{
			quatproduct(difference[i],cquaternion,between);
			quatproduct(quaternion,between,difference[i]);

			for(int i = 0; i < 4; i++ )
			{
				between[i] = 0.0;
			}
		}

		//CONSTRUCT COORDINATES AGAIN
		for(int i = 0; i < nratoms ; i++ )
			for(int j = 0; j < 3; j++)
				co1[i][j]	= difference[i][j+1] + co[0][j];					/*difference with first atom*/
	}

	//DEALLOCATE ARRAYS
	for (int i = 0; i < nratoms ; i++)
	{
		free(difference[i]);
	}
	free(difference);

	free(cquaternion);
	free(between);

}


void quatproduct(double *qt2, double *qt1, double *quaternion) /* calculates product quaternion = qt2*qt1*/
{
	double temp[4];
	temp[0] =	qt2[0]*qt1[0] - qt2[1]*qt1[1] - qt2[2]*qt1[2] - qt2[3]*qt1[3];
	temp[1] =	qt2[0]*qt1[1] + qt2[1]*qt1[0] + qt2[2]*qt1[3] - qt2[3]*qt1[2];
	temp[2] =	qt2[0]*qt1[2] - qt2[1]*qt1[3] + qt2[2]*qt1[0] + qt2[3]*qt1[1];
	temp[3] =	qt2[0]*qt1[3] + qt2[1]*qt1[2] - qt2[2]*qt1[1] + qt2[3]*qt1[0];
	quaternion[0] = temp[0];
	quaternion[1] = temp[1];
	quaternion[2] = temp[2];
	quaternion[3] = temp[3];
}



/*==============================================================================
 * DISTANCE.CPP
 *
 *  Created on: Nov 25, 2010
 *      Author: Pieter Chys
 *
 *==============================================================================
 *
 *
 *DISTANCE.C calculates the distance between two points in space.
 *This is done by using two pointers to arrays which contain the
 *point positions. The distance is returned as a double.
 *
 *----------------------------------------------------------------------------
 *
 *INPUT:-double* pointer to vector1[3]
 *      -double* pointer to vector2[2]
 *
 *-------------------------------------------------------------------------
 *
 *OUTPUT:-double number which is the distance
 *
 *===============================================================================
 */

double distance(double* point1,double* point2)
{
	double v[3];
	for(int i = 0; i < 3; i++)
	{
		v[i] = point2[i] - point1[i];
	}

	return sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
}

double distance2(double* point1,double* point2)
{
	double v[3];
	v[0] = point2[0] - point1[0];
	v[1] = point2[1] - point1[1];
	v[2] = point2[2] - point1[2];
	return (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

bool isnear_fast(double *p1, double *p2, double clash_factor)
{
	double d;
	d = pow(p1[0]-p2[0],2);
	if( d < clash_factor )
	{
		d += pow(p1[1]-p2[1],2);
		if( d < clash_factor )
		{
			d += pow(p1[2]-p2[2],2);
			if( d < clash_factor )
				return true;
		}
	}
	return false;
}

/*==============================================================================
 * DOTPROD.CPP
 *
 *  Created on: Nov 25, 2010
 *      Author: Pieter Chys
 *
 *==============================================================================
 *
 *
 *DOTPROD.C calculates the dot product between two vectors which
 *have three components (Cartesian components).  The data are
 *transferred as pointers to the corresponding arrays.
 *
 *----------------------------------------------------------------------------
 *
 *INPUT:-double* pointer to vector1[3]
 *      -double* pointer to vector2[3]
 *
 *-------------------------------------------------------------------------
 *
 *OUTPUT:-double number which equals the dot product
 *
 *--------------------------------------------------------------------------
 */

double dotprod(double* vector1,double* vector2) {
	double dot = 0;
	for(int i = 0; i < 3; i++) {
		dot += vector1[i]*vector2[i];
	}

	return dot;
}

/*==============================================================================
 * CROSSPROD.CPP
 *
 *  Created on: March 8, 2010
 *      Author: Pieter Chys
 *
 *==============================================================================
 *
 *
 *CROSSPROD.C calculates the cross or vector product between two Cartesian
 *vectors with components (X,Y,Z).  The data is transferred as pointer to
 *the resulting cross vector. The vector or cross product is formally:
 *
 *   C = A x B
 *
 *in which the order is relevant
 *----------------------------------------------------------------------------
 *
 *INPUT: -double *pointer to vector_a[3] = vector A
 *       -double *pointer to vector_b[3] = vector B
 *
 *-------------------------------------------------------------------------
 *
 *OUTPUT: - explicitly void
 *        - implicitly double *pointer to vector C
 *
 *
 *--------------------------------------------------------------------------
 */

//inline
void crossprod(double *vector_a,double *vector_b,double *vector_c) {
	vector_c[0] = vector_a[1]*vector_b[2] - vector_a[2]*vector_b[1];
	vector_c[1] = vector_a[2]*vector_b[0] - vector_a[0]*vector_b[2];
	vector_c[2] = vector_a[0]*vector_b[1] - vector_a[1]*vector_b[0];
}
	//	Other implementation --- slower

	//	double result[3][3];											//CODE IS SLOWER
	//	for(int i = 0; i < 3; i++){
	//		for(int j = 0; j < 3; j++){
	//			result[i][j] = vector_a[i]*vector_b[j];
	//		}
	//	}
	//	vector_c[0] = result[1][2] - result[2][1];
	//	vector_c[1] = result[2][0] - result[0][2];
	//	vector_c[2] = result[0][1] - result[1][0];

/*==============================================================================
 * norm.cpp
 *
 *  Created on: Nov 24, 2010
 *      Author: Pieter Chys
 *
 *==============================================================================
 *
 *NORM.C calculates the norm of a vector. This is its length
 *starting from the origin as base.
 *
 *----------------------------------------------------------------------------
 *
 *INPUT:-double* pointer to bond[3]
 *
 *-------------------------------------------------------------------------
 *
 *OUTPUT:-double number which is the norm of the bond
 *
 *===============================================================================
 */

double norm(double* bond) {
	double norm = sqrt(dotprod(bond,bond));
	return norm;
}

/*===============================================================================
 * normalize.cpp
 *
 *  Created on: Nov 24, 2010
 *      Author: Pieter Chys
 *
 *===============================================================================
 *
 *
 *NORMALIZE.C takes a vector and computes the corresponding
 *unit vector and its separate components. The transferral
 *of the components is done by a double* pointer to the
 *component array of the original vector.
 *
 *----------------------------------------------------------------------------
 *
 *INPUT:	- double* pointer to bond[3]
 *			- double* pointer to vector[3] to be filled
 *
 *-------------------------------------------------------------------------
 *
 *OUTPUT:	- explicitly void return
 *			- double* pointer to normalize[3]
 *
 *--------------------------------------------------------------------------
 */

void normalize(double* bond,double* normalize) {
	double res = norm(bond);
	for(int i = 0; i < 3; i++) {
		normalize[i] = bond[i]/res;
	}
}

void normalise(double* bond) {
	double res = norm(bond);
	for(int i = 0; i < 3; i++) {
		bond[i] = bond[i]/res;
	}
}


/*===============================================================================
 * DIFFVECTOR.CPP
 *
 *  Created on: Nov 24, 2010
 *      Author: Pieter Chys
 *
 *===============================================================================
 *
 *DIFFVECTOR.C computes the difference vector between two vectors
 *in 3D space.The input data is at present a double pointer
 *to the all coordinates of the involved molecule. The sequence
 *number of atom 1 and atom2 are as well given.
 *
 *--------------------------------------------------------------------------------
 *
 *INPUT: -double double** pointer to the Cartesian coordinates
 *       -int number referring to first atom position in molecule
 *       -int number referring to second atom position in molecule
 *--------------------------------------------------------------------------------
 *
 *OUTPUT:-double* pointer to  vector[3]
 *
 *================================================================================
 */

void diffvector(double** co, int number_atom_1, int number_atom_2,double* vector) {
	vector[0] = co[number_atom_2][0] - co[number_atom_1][0];
	vector[1] = co[number_atom_2][1] - co[number_atom_1][1];
	vector[2] = co[number_atom_2][2] - co[number_atom_1][2];
}

/*
 * diffvectorseq.cpp
 *
 *  Created on: Nov 24, 2010
 *      Author: Pieter Chys
 *
 *
 *===============================================================================
 *DIFFVECTORSEQ.CPP calculates the bond vector based on Cartesian coordinates for
 *a given molecule. The coordinates are passed by a double pointer and
 *used to construct bond b by making a basic subtraction of atom positions.
 *
 *-------------------------------------------------------------------------------
 *
 *INPUT: -double double** pointer to a matrix-like array of data
 *		  which refer to the Cartesian coordinates
 *	     -int number which is the bond that needs to be constructed
 *	     -bond which will contain the resulting bond vector
 *
 *--------------------------------------------------------------------------------
 *
 *OUTPUT:-pointer double* that refers to bond[3] which is the
 *        resulting bond vector
 *
 *================================================================================
 */

void diffvectorseq(double** co, int atom_number,double* bond) {						//	atom_number ends at a certain bond
	bond[0] = co[atom_number-1][0] - co[atom_number-2][0];							//	Modified at 27/01/2011: indices minus 1
	bond[1] = co[atom_number-1][1] - co[atom_number-2][1];
	bond[2] = co[atom_number-1][2] - co[atom_number-2][2];
}

/*================================================================================
 * flipud.cpp
 *
 *  Created on: Dec 9, 2010
 *      Author: Pieter Chys
 *
 *================================================================================
 *
 *FLIPUD.CPP changes the order of the coordinates in the
 *pointer (or bidimensional array). The first atom coordinate
 *becomes the last one and vice versa. The second last one
 *the second and vice versa. And so on for all atoms. The
 *function is an analog to the MATLAB function flipud.
 *
 *--------------------------------------------------------------------------------
 *
 *INPUT:-	double** pointer to coordinates
 *      -	int which denotes the number of coordinates
 *
 *--------------------------------------------------------------------------------
 *
 *OUTPUT: 	- void
 *			- input double** pointer is changed for the new order
 *
 *================================================================================
 */


void flipud(double** co, int number_atoms) {
	double temp[3];
	for(int i = 0; i < ceil(number_atoms/2); i++) { 								//no problems with  division by two
		for(int j = 0; j < 3; j++) {
			temp[j]  = co[number_atoms-i-1][j];
			co[number_atoms-i-1][j] = co[i][j];
			co[i][j] = temp[j];
		}
	}

}

/*================================================================================
 * MATRIXXVECTOR.CPP
 *
 *  Created on: Dec 9, 2010
 *      Author: Pieter Chys
 *
 *================================================================================
 *
 *MATRIXXVECTOR.CPP executes matrix multiplication  of
 *a 3x3 square matrix with a 3-dimensional vector.
 *
 *--------------------------------------------------------------------------------
 *
 *INPUT:- double** pointer to an allocated 3x3 matrix
 *      - double* pointer to the allocated vector
 *		- double* pointer to an allocated array[3] to contain result
 *
 *--------------------------------------------------------------------------------
 *
 *OUTPUT: 	- explicitly void
 *			- input double* pointer which is filled array[3]
 *
 *================================================================================
 */

void matrixXvector(double **matrix, double *vector, double *result)
{
	double temp[3];
	for(int i = 0; i < 3; i++)
	{
		temp[i] = 0.0;
		for(int j = 0; j < 3; j++)
			temp[i] += matrix[i][j]*vector[j];
	}
	for(int i = 0; i < 3; i++)
		result[i] = temp[i];
}

// Mon: faster
void matrixXvector(double matrix[][3], double *vector)
{
	double result[3];
	result[0] = matrix[0][0]*vector[0] + matrix[0][1]*vector[1] + matrix[0][2]*vector[2];
	result[1] = matrix[1][0]*vector[0] + matrix[1][1]*vector[1] + matrix[1][2]*vector[2];
	result[2] = matrix[2][0]*vector[0] + matrix[2][1]*vector[1] + matrix[2][2]*vector[2];
	vector[0] = result[0];
	vector[1] = result[1];
	vector[2] = result[2];
}


/*================================================================================
 * MATRIXXMATRIX.CPP
 *
 *  Created on: February 9, 2011
 *      Author: Pieter Chys
 *
 *================================================================================
 *
 *MATRIXXMATRIX.CPP executes matrix multiplication of 2
 *3x3 square matrices.
 *
 *--------------------------------------------------------------------------------
 *
 *INPUT:- double **pointer to first allocated 3x3 matrix
 *      - double **pointer to second allocated 3x3 matrix
 *		- double **pointer to allocated 3x3 result matrix
 *
 *--------------------------------------------------------------------------------
 *
 *OUTPUT: 	- explicitly void
 *			- input double **pointer of resulting 3x3 matrix
 *			-
 *REMARK:	- order of input matrices matters!!!
 *
 *================================================================================
 */


void mapmatrix(double **temp,double **matrix1){
	matrix1[0][0] = temp[0][0];
	matrix1[0][1] = temp[0][1];
	matrix1[0][2] = temp[0][2];
	matrix1[1][0] = temp[1][0];
	matrix1[1][1] = temp[1][1];
	matrix1[1][2] = temp[1][2];
	matrix1[2][0] = temp[2][0];
	matrix1[2][1] = temp[2][1];
	matrix1[2][2] = temp[2][2];
}


//	Implementation with for loops
void transpose(double **matrix1,double **matrix2)
{
	double temp[3][3];						//	allow that matrix1 = result_matrix --- this code is a bit slower: 1/7
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			temp[i][j] = matrix1[j][i];
		}
	}

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			matrix2[i][j] = temp[i][j];
		}
	}
}


void matrixXmatrix(double **matrix1,double **matrix2,double **output)
{
	double temp[3][3];						//	allow that matrix1 = result_matrix --- this code is a bit slower: 1/7
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			temp[i][j] = 0.0;
			for (int k = 0; k < 3; k++) {
				temp[i][j] += matrix1[i][k]*matrix2[k][j];
			}
		}
	}
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			output[i][j] = temp[i][j];
		}
	}
}


void matrixprod1(double **tot,double **input){
	double temp[3][2];
	temp[0][0]	=	tot[0][1]*input[0][0] + tot[0][2]*input[1][0];
	temp[1][0]  =	tot[1][1]*input[0][0] + tot[1][2]*input[1][0];
	temp[2][0]	=	tot[2][1]*input[0][0] + tot[2][2]*input[1][0];
	temp[0][1]	=	tot[0][2]*input[1][1] - tot[0][1]*input[0][1];
	temp[1][1]  =	tot[1][2]*input[1][1] - tot[1][1]*input[0][1];
	temp[2][1]	=	tot[2][2]*input[1][1] - tot[2][1]*input[0][1];
	tot[0][1]   = 	temp[0][0];
	tot[1][1]  	= 	temp[1][0];
	tot[2][1]  	= 	temp[2][0];
	tot[0][2] 	= 	temp[0][1];
	tot[1][2] 	= 	temp[1][1];
	tot[2][2] 	= 	temp[2][1];
}

/*=====================================================================
 * MATRIXPROD2.CPP
 *
 *  Created on: May 9, 2011
 *      Author: Pieter Chys
 *
 *=====================================================================
 *
 *
 *MATRIXPROD2.CPP performs the second matrix product to construct a
 *protein backbone chain without side chains. In principle the method
 *is applicable to any kind of polymer backbone.
 *
 *----------------------------------------------------------------------------
 *
 *INPUT: - double **pointer to TOT which is the accumulated matrix (3 x 3)
 *       - double **pointer to input which is contracted (2x2) and only
 *         contains four elements elements as in the reference article below.
 *         Both the COSINE AND SINE of the valence angle are used here.
           The COSINE needs to have negative sign in the allocated input
           matrix, the SINE hower must be positive in sign
 *
 *----------------------------------------------------------------------------
 *
 *OUTPUT:	- explicitly VOID return
 *			- double **pointer to TOT which is as well the output
 *
 *----------------------------------------------------------------------------
 *
 *REMARK:	- before MATRIXPROD2 the function MATRIXPROD1 is used
 *            after which the bond orientation can be calculated
 *          - the implementation is after Seok C. and Coutsias E.A. :
 *            "Efficiency of rotational operators for geometric manipulation
 *             of chain molecules"
 *          - the framework is pure 3x3 matrix algebra
 *          - the algorithm is written to be able to iterate by
 *            writing to the input variable by a temporary one.
 *
 *
 *============================================================================
 */

void matrixprod2(double **tot,double **input){
	double temp[3][2];
	temp[0][0]	=	tot[0][0]*input[0][0] + tot[0][1]*input[1][0];
	temp[1][0]  =	tot[1][0]*input[0][0] + tot[1][1]*input[1][0];
	temp[2][0]	=	tot[2][0]*input[0][0] + tot[2][1]*input[1][0];
	temp[0][1]	=	tot[0][1]*input[1][1] - tot[0][0]*input[0][1];
	temp[1][1]  =	tot[1][1]*input[1][1] - tot[1][0]*input[0][1];
	temp[2][1]	=	tot[2][1]*input[1][1] - tot[2][0]*input[0][1];
	tot[0][0]   = 	temp[0][0];
	tot[1][0]  	= 	temp[1][0];
	tot[2][0]  	= 	temp[2][0];
	tot[0][1] 	= 	temp[0][1];
	tot[1][1] 	= 	temp[1][1];
	tot[2][1] 	= 	temp[2][1];
}


/*================================================================================
 * MINUS.CPP
 *
 *  Created on: Dec 9, 2010
 *      Author: Pieter Chys
 *
 *================================================================================
 *
 *MINUS.CPP takes a vector and inverts all of its components;
 *this is changing the value from + to - values with same magnitude.
 *Thus vector[j] becomes -vector[j].
 *
 *--------------------------------------------------------------------------------
 *
 *INPUT:- DOUBLE* pointer to a vector (VECTOR)
 *
 *
 *
 *--------------------------------------------------------------------------------
 *
 *OUTPUT: 	- explicitly void
 *			- DOUBLE* pointer to vector which is replaced
 *
 *
 *================================================================================
 */

void min(double* vector) {
	for(int j = 0; j < 3; j++) {
		vector[j] *= -1;
	}

}

/*================================================================================
 * RANDNVECTOR.CPP
 *
 *  Created on: Dec 13, 2010
 *      Author: Pieter Chys
 *
 *================================================================================
 *
 *RANDNVECTOR.CPP create a normalized random vector. The vector
 *is normalized to magnitude one.
 *
 *--------------------------------------------------------------------------------
 *
 *INPUT:- DOUBLE* pointer to a vector (VECTOR)
 *
 *
 *--------------------------------------------------------------------------------
 *
 *OUTPUT: 	- explicitly VOID
 *			- DOUBLE* pointer to created vector
 *
 *================================================================================
 */

void randnvector(double* vector) {
	for(int j=0;j<3;j++) {
		vector[j] = ((double)rand()/(double)RAND_MAX);
	}

	normalize(vector,vector); //  WILL THIS STATEMENT WORK WELL?...CURIOUS THING...
							  //  SHOULD WORK FINE..NACHO
}


/*================================================================================
 * SUBTRACT.CPP
 *
 *  Created on: Dec 13, 2010
 *      Author: Pieter Chys
 *
 *================================================================================
 *
 *SUBTRACT.CPP subtracts two 3-dimensional vectors and returns the
 *result as a new pointer.
 *
 *
 *--------------------------------------------------------------------------------
 *
 *INPUT: -DOUBLE* pointer to vector 1
 *       -DOUBLE* pointer to vector 2
 *       -DOUBLE* pointer to subtraction vector
 *
 *
 *--------------------------------------------------------------------------------
 *
 *OUTPUT: 	- explicitly VOID
 *			- DOUBLE* pointer to vector resulting from subtraction
 *			  from the two vectors
 *
 *================================================================================
 */

//void subtract(double* vector1,double* vector2,double* subtraction) {
//	double temp[3];
//	for(int j = 0; j < 3; j++) {
//		temp[j] = vector2[j] - vector1[j];
//	}
//
//	for(int j = 0; j < 3; j++) {
//		subtraction[j] = temp[j];
//	}
//}

// vector2 - vector1 = subtraction
void subtract(double* vector1,double* vector2,double* subtraction)
{
	for(int j = 0; j < 3; j++)
		subtraction[j] = vector2[j] - vector1[j];
}

void add(double* vector1,double* vector2,double* added) {
	double temp[3];
	for(int j = 0; j < 3; j++) {
		temp[j] = vector2[j] + vector1[j];
	}

	for(int j = 0; j < 3; j++) {
		added[j] = temp[j];
	}
}



double 	det33(double** matrix) {									//	seems Ok!
	return matrix[0][0]*matrix[1][1]*matrix[2][2] +
			matrix[0][1]*matrix[1][2]*matrix[2][0] +
			 matrix[0][2]*matrix[1][0]*matrix[2][1] -
			  matrix[0][0]*matrix[1][2]*matrix[2][1] -
			  matrix[0][1]*matrix[1][0]*matrix[2][2] -
			  matrix[0][2]*matrix[1][1]*matrix[2][0];


}

