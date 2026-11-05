/*============================================================================
 * GA.CPP
 *
 *  Created on: Dec 10, 2010
 *      Author: Pieter Chys
 *
 *===============================================================================
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cassert>
#include "ga.h"
#include "own_math.h"
#include "own_io.h"


#define QUAT_TRACE_QZERO_TOLERANCE		0.1





/*====================================================================
 * GEOM11.CPP
 *
 *  Created on: November 25, 2010
 *      Author: Pieter Chys
 *
 *
 *=====================================================================
 *
 *GEOM11.CPP computes the geometric product ab of two 1-vectors.  The dimension
 *of a and b must be the same and is of Cartesian format (dimension 3).
 *The result is a bivector C of the form C = c0 + [c1s1s2 + c2s2s3 + c3s1s3]
 *The with the last 3 terms the pure bivector terms in unit blades I1,I2 and
 *I3 (s1s2,s2s3 ands1s3).
 *
 *----------------------------------------------------------------------------
 *
 *INPUT: - double *pointer to vector1[3]
 *       - double *pointer to vector2[3]
 *
 *-------------------------------------------------------------------------
 *
 *OUTPUT: - explicitly void
 *        - double *pointer to bivector[4]
 *
 *=======================================================================
 */

void geom11(double *vector1,double *vector2,double *bivector)
{
	bivector[0]	=	dotprod(vector1,vector2);
	bivector[1]	=	vector1[0]*vector2[1]-vector1[1]*vector2[0];
	bivector[2]	=	vector1[1]*vector2[2]-vector1[2]*vector2[1];
	bivector[3]	=	vector1[0]*vector2[2]-vector1[2]*vector2[0];
}


/*=====================================================================
 * GEOM11IN.CPP
 *
 *  Created on: October 15, 2010
 *      Author: Pieter Chys
 *
 *
 *=====================================================================
 *
 *GEOM11IN.CPP takes 3-dimensional vectors a and b
 *and calculates their inner geometric product. This basically
 *identical to the dot product a.b = a1s1*b1s1 + a2s2*b2s2 + a3s3*b3s3
 *and is returned as a double explicitly.
 *
 *
 *-----------------------------------------------------------------
 *	INPUT:  - double *pointer to vector a
 *	        - double *pointer to vector b
 *          - integer length which refers to the length of the vector;
 *            this is mostly dimension 3 (Cartesian coordinates).
 *
 *------------------------------------------------------------------
 *  OUTPUT: - double number which is the inner geometric product
 *          - of two vectors
 *
 *====================================================================
 */
double geom11in(double *a,double *b,int length)
{
	double runningSum = 0.0;
	for(int i = 0; i < length; i++) {
		runningSum += a[i] * b[i];
	}

	return runningSum;
}


/*=====================================================================
 * GEOM11OUT.CPP
 *
 *  Created on: November 25, 2010
 *      Author: Pieter Chys
 *
 *=====================================================================
 *
 *
 *GEOM11OUT.CPP returns the outer geometric product of
 *1-vectors a en b.  It uses therefore GEOM11.CPP by using
 *it 2 times in a row. As such this function is actually
 *a real function in GEOM11.CPP.
 *
 *----------------------------------------------------------------------------
 *
 *INPUT: -double *pointer to vector1[3]
 *       -double *pointer to vector2[3]
 *
 *-------------------------------------------------------------------------
 *
 *OUTPUT: - double *pointer to wedge[4]
 *        - explicitly returns void
 *
 *=====================================================================
 */
void geom11out(double *vector1,double *vector2,double *wedge)
{
	double product1[4],product2[4];
	geom11(vector2,vector1,product1);
	geom11(vector1,vector2,product2);
	for(int i = 0; i < 4; i++)
	{
		wedge[i] = 0.50*(product2[i] - product1[i]);
	}
}

/*====================================================================
 * GEOM22.CPP
 *
 *  Created on: March 11, 2010
 *      Author: Pieter Chys
 *
 *
 *=====================================================================
 *
 *GEOM22.CPP computes the geometric product AB of 2 bivectors.  The dimension
 *of A and B must be the same.  The result is a bivector C of the form
 *C = c0 + [c1s1s2 + c2s2s3 + c3s1s3] with the latter three terms the pure
 *bivector part in unit blades I1,I2 and I3.  This is the bivector product
 *for a Cartesian basis.
 *
 *----------------------------------------------------------------------------
 *
 *INPUT: - double *pointer to bivector1[4]
 *       - double *pointer to bivector2[4]
 *
 *-------------------------------------------------------------------------
 *
 *OUTPUT: - explicitly void
 *		  - double *pointer to bivector3[4]
 *
 *----------------------------------------------------------------------------
 *
 *REMARKS: - this function is to be used with function SPINOR.CPP; after
 *           that 2 spinors are created a resulting one can be constructed
 *           by a spinorproduct. In fact function SPINORPRODUCT.CPP is simply
 *           GEOM22.CPP but modified for nomenclature (basically identical).
 *
 *=======================================================================
 */


void geom22(double  *bivector1,double *bivector2,double *bivector3)
{
	bivector3[0]	=	bivector1[0]*bivector2[0] - dotprod(&bivector1[1],&bivector2[1]);
	bivector3[1]	=	bivector1[0]*bivector2[1] + bivector1[1]*bivector2[0] - bivector1[3]*bivector2[2] + bivector1[2]*bivector2[3];
	bivector3[2]	=	bivector1[0]*bivector2[2] + bivector1[2]*bivector2[0] - bivector1[1]*bivector2[3] + bivector1[3]*bivector2[1];
	bivector3[3]	=	bivector1[0]*bivector2[3] + bivector1[3]*bivector2[0] - bivector1[2]*bivector2[1] + bivector1[1]*bivector2[2];
}


/*====================================================================
 * GEOM23.CPP
 *
 *  Created on: November 24, 2010
 *      Author: Pieter Chys
 *
 *====================================================================
 *
 *
 *GEOM23.CPP returns the geometric product BT (different from TB)
 *in which B is an even bivector and T an odd trivector and this yields a
 *trivector Q of the form Q = [q1s1 + q2s2 + q3s3] +q4s1s2s3 in which
 *the first triple term is of Cartesian vector format and the last
 *one a pseudocalar i term.
 *
 *
 *--------------------------------------------------------------------
 *
 *INPUT:	-double *pointer to spinor[4]
 *      	-double *pointer to trivector[4]
 *
 *---------------------------------------------------------------------
 *
 *OUTPUT:	- explicitly void return
 *			- double *pointer to result[4]
 *
 *=====================================================================
 */

void geom23(double *spinor,double *trivector,double *result)
{
	result[0] 	= spinor[0]*trivector[0] + spinor[1]*trivector[1]
	          	     -spinor[2]*trivector[3] +spinor[3]*trivector[2];
	result[1]	= spinor[0]*trivector[1] - spinor[1]*trivector[0]
	         	     + spinor[2]*trivector[2] +spinor[3]*trivector[3];
	result[2]	= spinor[0]*trivector[2] - spinor[1]*trivector[3]
					- spinor[2]*trivector[1] -spinor[3]*trivector[0];
	result[3]	= spinor[0]*trivector[3] + spinor[1]*trivector[2]
						+ spinor[2]*trivector[0] -spinor[3]*trivector[1];
}


/*======================================================================
 * GEOM32.CPP
 *
 *  Created on: November 24, 2010
 *      Author: Pieter Chys
 *
 *======================================================================
 *
 *
 *GEOM32.CPP returns the geometric product TB in which
 *T is an odd trivector and B an even bivector and yields
 *an odd trivector in form P = [p1s1 + p2s2s + p3s3] +p4s1s2s3
 *as in the case for geomi2(B). The first three terms of P
 *form a vector coordinate in Cartesian space; the last term
 *is the pseudoscalar term.
 *
 *----------------------------------------------------------------------------
 *
 *INPUT:	- double* pointer to trivector[4]
 *      	- double* pointer to spinor[4]
 *
 *-------------------------------------------------------------------------
 *
 *OUTPUT:	- explicitly void
 *			- double* pointer to result[4]
 *
 *======================================================================
 */

void geom32(double *trivector,double *spinor,double *result)
{
	result[0]	= trivector[0]*spinor[0] - trivector[1]*spinor[1]
	         	    -trivector[2]*spinor[3] -trivector[3]*spinor[2];
	result[1]	= trivector[0]*spinor[1] + trivector[1]*spinor[0]
				    -trivector[2]*spinor[2] -trivector[3]*spinor[3];
	result[2]	= trivector[0]*spinor[3] + trivector[1]*spinor[2]
					+ trivector[2]*spinor[0] - trivector[3]*spinor[1];
	result[3]	= trivector[0]*spinor[2] - trivector[1]*spinor[3]
					+ trivector[2]*spinor[1] +trivector[3]*spinor[0];
}


/*====================================================================
 * GEOM23OPT.CPP
 *
 *  Created on: June 14, 2010
 *      Author: Pieter Chys
 *
 *====================================================================
 *
 *
 *GEOM23.CPP returns the geometric product BT (different from TB)
 *in which B is an even bivector and T an odd trivector and this yields a
 *trivector Q of the form Q = [q1s1 + q2s2 + q3s3] +q4s1s2s3 in which
 *the first triple term is of Cartesian vector format and the last
 *one a pseudocalar i term.
 *
 *
 *--------------------------------------------------------------------
 *
 *INPUT:	-double *pointer to spinor[4]
 *      	-double *pointer to trivector[4]
 *
 *---------------------------------------------------------------------
 *
 *OUTPUT:	- explicitly void return
 *			- double *pointer to result[4]
 *
 *=====================================================================
 */

void geom23opt(double *spinor,double *trivector,double *result)
{		//	4-th component must be zero in coordinate applications
	result[0] 	= spinor[0]*trivector[0] + spinor[1]*trivector[1] - spinor[2]*trivector[3] +spinor[3]*trivector[2];
	result[1]	= spinor[0]*trivector[1] - spinor[1]*trivector[0] + spinor[2]*trivector[2] +spinor[3]*trivector[3];
	result[2]	= spinor[0]*trivector[2] - spinor[1]*trivector[3] - spinor[2]*trivector[1] -spinor[3]*trivector[0];
	//result[3]  = 0;
	//result[3]	= spinor[0]*trivector[3] + spinor[1]*trivector[2] + spinor[2]*trivector[0] -spinor[3]*trivector[1];
}


/*======================================================================
 * GEOM32OPT.CPP
 *
 *  Created on: June 15, 2011
 *      Author: Pieter Chys
 *
 *======================================================================
 *
 *
 *GEOM32.CPP returns the geometric product TB in which
 *T is an odd trivector and B an even bivector and yields
 *an odd trivector in form P = [p1s1 + p2s2s + p3s3] +p4s1s2s3
 *as in the case for geomi2(B). The first three terms of P
 *form a vector coordinate in Cartesian space; the last term
 *is the pseudoscalar term.
 *
 *----------------------------------------------------------------------------
 *
 *INPUT:	- double* pointer to trivector[4]
 *      	- double* pointer to spinor[4]
 *
 *-------------------------------------------------------------------------
 *
 *OUTPUT:	- explicitly void
 *			- double* pointer to result[4]
 *
 *======================================================================
 */

void geom32opt(double *trivector,double *spinor,double *result)
{		//	4-th component is zero for trivector in coordinate applications
	result[0]	= trivector[0]*spinor[0] - trivector[1]*spinor[1] - trivector[2]*spinor[3];
	result[1]	= trivector[0]*spinor[1] + trivector[1]*spinor[0] - trivector[2]*spinor[2];
	result[2]	= trivector[0]*spinor[3] + trivector[1]*spinor[2] + trivector[2]*spinor[0];
	result[3]	= trivector[0]*spinor[2] - trivector[1]*spinor[3] + trivector[2]*spinor[1];
}



 /*======================================================================
 * GEOMI2.CPP
 *
 *  Created on: November 25, 2010
 *      Author: Pieter Chys
 *
 *======================================================================
 *
 *
 *GEOMI2.CPP returns the geometric product of unit trivector
 *iota (also called pseudoscalar s1s2s3 or i) with an arbitrary
 *bivector B = b0 + [b1s1s2 + b2s2s3 +b3s1s3] (unit blades in
 *bivector the end) in format of trivector
 *[i1s1 + i2s2 + i3s3] +  i4s1s2s3.
 *
 *----------------------------------------------------------------------------
 *
 *INPUT:  - double *pointer to bivector[4]
 *
 *-------------------------------------------------------------------------
 *
 *OUTPUT: - double *pointer to trivector[4]
 *
 *------------------------------------------------------------------------
 *
 *REMARKS: - at present this functions calculates the minus geomi2
 *           to bypass the effect of making a function for - for
 *           pointers to arrays
 *
 *=====================================================================
 */

void geomi2(double* bivector,double* trivector)
{
	trivector[0] =     bivector[2]; // formulas to yield -geomi2
	trivector[1] =   - bivector[3];
	trivector[2] =     bivector[1];
	trivector[3] =   - bivector[0];
}


/*=====================================================================
 * SPINOR.CPP
 *
 *  Created on: November 24, 2010
 *      Author: Pieter Chys
 *
 *=====================================================================
 *
 *
 *SPINOR.CPP constructs a spinor as a pointer (double*) to a 4-array which rotates
 *a 1-vector through an angle phi radians (double) around the unit axis vector
 *as specified by an input pointer (double*) to an array[3]. The resulting spinor
 *has form B = b0 + [b1 I1 + b2 I2 + b3 I3] with I1, I2 and I3 being unit blades.
 *The unit blades are I1 = s1s2, I2 = s2s3 and I3 = s1s3.
 *
 *----------------------------------------------------------------------------
 *
 *INPUT: - double *pointer to unit_axis[3]
 *       - double angle of rotation
 *
 *----------------------------------------------------------------------------
 *
 *OUTPUT:	- explicitly VOID return
 *			- double *pointer to spinor[4]
 *
 *----------------------------------------------------------------------------
 *
 *REMARK:	- following the construction of 2 SPINORS one can use SPINORPRODUCT
 *            to create a  third new spinor which represent the result of multiplying
 *            both spinors with each other.
 *          - an implementation with approximations for small angles of the
 *            trigonometric functions is slower than the real trigonometric functions;
 *            also without the if statement... (3 times slower)
 *
 *
 *============================================================================
 */

void spinor(double *unit_axis,double angle,double *spinor)
{
	double sin_value = sin(angle/2.0);
	spinor[0] = cos(angle/2.0);
	spinor[1] =  unit_axis[2]*sin_value;
	spinor[2] =  unit_axis[0]*sin_value;
	spinor[3] = -unit_axis[1]*sin_value;
}


/*=====================================================================
 * REVSPINOR.CPP
 *
 *  Created on: November 24, 2010
 *      Author: Pieter Chys
 *
 *=====================================================================
 *
 *
 *REVSPINOR.CPP returns the reverse spinor as a pointer (double*) to a 4-array
 *which is the complement to the SPINOR formed by SPINOR.CPP which rotates a 1-vector
 *which through an angle phi radians (double) around the unit axis vector
 *a specified by a pointer (double*) to an array[3]. The revspinor has final form
 *B = b0 + [b1 I1 + b2 I2 + b3 I3] with I1, I2 and I3 being the unit blades
 *and has only modified signs in its terms as compared to the original spinor.
 *
 *----------------------------------------------------------------------------
 *
 *INPUT: -double *pointer to spinor[4]
 *
 *-------------------------------------------------------------------------
 *
 *OUTPUT: 	- explicitly VOID
 *		    - double *pointer to revspinor[4]
 *
 *======================================================================
 */

void revspinor(double *spinor,double *revspinor)
{
	revspinor[0] = 	spinor[0];
	revspinor[1] = -spinor[1];
	revspinor[2] = -spinor[2];
	revspinor[3] = -spinor[3];
}

/*=====================================================================
 * SPINORPRODUCT.CPP
 *
 *  Created on: March 9, 2010
 *      Author: Pieter Chys
 *
 *=====================================================================
 *
 *
 *SPINORPRODUCT.CPP calculates the spinorproduct R3 = R1*R2 in that
 *specific order in which R1 is the spinor represented by pointer *a
 *and R2 the one by *b. The result is passed as a double pointer *c
 *which is the total spinor R3
 *
 *----------------------------------------------------------------------------
 *
 *INPUT: - pointer to spinor R1 (double *a, a[4])
 *		 - pointer to spinor R2 (double *b, b[4])
 *-------------------------------------------------------------------------
 *
 *OUTPUT: 	- explicitly VOID
 *		    - pointer to total spinor R3 (double *c, c[4])
 *
 *======================================================================
 */

//	Raw implementation avoiding for loops and mapping to its own allowed

void spinorproduct(double *bivector1,double *bivector2,double *result)
{
	double temp[4];
//	for(int i = 0; i < 4; i++){
//	  bivector2[i] = temp[i];
//	}
	temp[0]	=	bivector1[0]*bivector2[0] - bivector1[1]*bivector2[1] - bivector1[2]*bivector2[2] - bivector1[3]*bivector2[3];
	temp[1]	=	bivector1[0]*bivector2[1] + bivector1[1]*bivector2[0] - bivector1[3]*bivector2[2] + bivector1[2]*bivector2[3];
	temp[2]	=	bivector1[0]*bivector2[2] + bivector1[2]*bivector2[0] - bivector1[1]*bivector2[3] + bivector1[3]*bivector2[1];
	temp[3]	=	bivector1[0]*bivector2[3] + bivector1[3]*bivector2[0] - bivector1[2]*bivector2[1] + bivector1[1]*bivector2[2];
	result[0] = temp[0];
	result[1] = temp[1];
	result[2] = temp[2];
	result[3] = temp[3];
}


void spinor3prod1(double* tot,double* input,double *output)
{
	output[0]	=	tot[0]*input[0] - tot[1]*input[1];
	output[1]	=	tot[1]*input[0] + tot[0]*input[1];
	output[2]	=	tot[2]*input[0] + tot[3]*input[1];
	output[3]	=	tot[3]*input[0] - tot[2]*input[1];
}

   /*TRIAL FOR FASTER TRIGONOMETRIC FUNCTIONS AT SMALL ANGLES (does not work: code is slower)
	*double sin_value;
	*if ( angle > -0.244 && angle < 0.244){ 			//	error of = 1% implemented at present
	*	sin_value = angle/2.0- angle*angle*angle/48.0;
	*	spinor[0] = 1 - angle*angle/8.0;
	*} else {
	*	sin_value = sin(angle);
	*	spinor[0] = cos(angle);
	*}
	*/




void spinor3prod2(double* tot,double* input,double *output)
{
	output[0]	=	tot[0]*input[0] - tot[3]*input[1];
	output[1]	=	tot[1]*input[0] + tot[2]*input[1];
	output[2]	=	tot[2]*input[0] - tot[1]*input[1];
	output[3]	=	tot[3]*input[0] + tot[0]*input[1];
}

/*=====================================================================
 * SPINORMATRIX.CPP
 *
 *  Created on: May 10, 2010
 *  Author: Pieter Chys
 *
 *=====================================================================
 *
 *
 *SPINORMATRIX.CPP calculates the corresponding 3x3 matrix to a quaternion
 *(q0,q1,q2,q3) or alternatively a spinor (s0,s1,s2,s3) which is characterized
 *by a normalized rotation axis and a dihedral rotation. The relevance of this
 *algorithm is that in successive application of the same quaternion/rotor/spinor
 *to coordinates the constituting matrix elements are recalculated. Therefore,
 *a soon as we have a few coordinates to manipulate with the same
 *quaternion/spinor/rotor the conversion to a matrix and subsequent matrix
 *multiplication become more efficient.
 *
 *----------------------------------------------------------------------------
 *
 *INPUT: - pointer to spinor/quaternion (double *spinor, spinor[4])
 *
 *-------------------------------------------------------------------------
 *
 *OUTPUT: 	- explicitly VOID
 *		    - pointer to corresponding 3x3 matrix (double **matrix, matrix[3][3])
 *
 *-------------------------------------------------------------------------
 *
 *REMARKS;	- code is implemented as defined by quaternion conversion; this needs
 *            to be checked carefully
 *
 *======================================================================
 */
 void spinormatrix(double *spinor, double **matrix)
 {
	 //	Straight implementation
	 matrix[0][0] =  	2*(spinor[0]*spinor[0] + spinor[1]*spinor[1]) - 1.0;
	 matrix[0][1] =  	2*(spinor[1]*spinor[2] - spinor[0]*spinor[3]);
	 matrix[0][2] = 	2*(spinor[1]*spinor[3] - spinor[0]*spinor[2]);
	 matrix[1][0] = 	2*(spinor[1]*spinor[2] + spinor[0]*spinor[3]);
	 matrix[1][1] = 	2*(spinor[0]*spinor[0] + spinor[2]*spinor[2]) - 1.0;
	 matrix[1][2] =     2*(spinor[2]*spinor[3] - spinor[0]*spinor[1]);
	 matrix[2][0] =     2*(spinor[1]*spinor[3] - spinor[0]*spinor[2]);
	 matrix[2][1] = 	2*(spinor[2]*spinor[3] + spinor[0]*spinor[1]);
	 matrix[2][2] = 	2*(spinor[0]*spinor[0] + spinor[3]*spinor[3])- 1.0;

	 //	Loop implementation
	 /*double temp[10];
	 int k = 0;
	 for (int i = 0; i < 3; i++){
		 for (int j = 0; j < i + 1; j++){
			 temp[k] = 2*spinor[i]*spinor[j];
			 k++;
		 }
	 }
	 matrix[0][0] = temp[0] + temp[2] - 1;
	 matrix[0][1] = temp[4] - temp[6];
	 matrix[0][2] = temp[7] - temp[3];
	 matrix[1][0] = temp[4] + temp[6];
	 matrix[1][1] = temp[0] + temp[5] - 1;
	 matrix[1][2] = temp[8] - temp[1];
	 matrix[2][0] = temp[7] - temp[3];
	 matrix[2][1] = temp[8] + temp[1];
	 matrix[2][2] = temp[0] + temp[9] - 1;*/

 }




void spinormat(double *unit_axis,double angle,double *spinor)
{
	double sin_value = sin(angle/2.0);
	spinor[0] = cos(angle/2.0);
	spinor[1] =  unit_axis[0]*sin_value;
	spinor[2] =  unit_axis[1]*sin_value;
	spinor[3] =  unit_axis[2]*sin_value;
}


/*======================================================================
  *MAPSPINOR.CPP
  *
  *Created on: July 16, 2012
  *Author:     Pieter Chys
  *
  *========================================================================
  *
  *
  *
  *MAPSPINOR.CPP simply copies the input spinor coordinates into
  *the output spinor coordinates. This simply corresponds to a
  *copy+paste action of an existing spinor to a second array.
  *------------------------------------------------------------------------
  *INPUT:  - double *pointer to BIVECTOR1 which is the spinor to be copied
  *------------------------------------------------------------------------
  *OUTPUT: - explicitly VOID return
  *		   - double *pointer BIVECTOR2 which is the output spinor
  *------------------------------------------------------------------------
  *REMARK:-  this function can in fact be used to copy-paste any
  *          4-D vector; it can thus be used to copy-paste a quaternion
  *          which is most logical in the context here
  *======================================================================*/
void mapspinor(double* bivector1,double* bivector2)
{
	bivector2[0] = bivector1[0];
	bivector2[1] = bivector1[1];
	bivector2[2] = bivector1[2];
	bivector2[3] = bivector1[3];
}

//	Raw/loop implementation allowing mapping on itself

//void spinorproduct(double* bivector1,double* bivector2){
//	double temp[4];
//	temp[0]	=	bivector1[0]*bivector2[0] - bivector1[1]*bivector2[1] - bivector1[2]*bivector2[2] - bivector1[3]*bivector2[3];
//	temp[1]	=	bivector1[0]*bivector2[1] + bivector1[1]*bivector2[0] - bivector1[3]*bivector2[2] + bivector1[2]*bivector2[3];
//	temp[2]	=	bivector1[0]*bivector2[2] + bivector1[2]*bivector2[0] - bivector1[1]*bivector2[3] + bivector1[3]*bivector2[1];
//	temp[3]	=	bivector1[0]*bivector2[3] + bivector1[3]*bivector2[0] - bivector1[2]*bivector2[1] + bivector1[1]*bivector2[2];
//	for(int i = 0; i < 4; i++){
//		bivector2[i] = temp[i];
//	}
//}

// Loop	implementation trying to take advantage of for loops

//void spinorproduct(double* spinor1,double* spinor2,double *temp) {
//	double result[3][3];
//	for(int i = 0; i < 3; i++){
//		for(int j = 0; j < 3; j++){
//			result[j][i] = spinor2[j]*spinor1[i];
//		}
//	}
////	double temp2[4];
////	for(int i = 1; i < 3; i++){
////		for(int j = 0; j < 3; j++){
////			temp[i] = spinor2[j]*spinor1[i];
////		}
////	}
//
//	spinor2[0]	= result[0][0] -result[1][1] -  result[2][2]-  spinor1[3]*spinor2[3];
//	spinor2[1]	= result[0][1] + result[1][0] - spinor1[3]*spinor2[2] + spinor1[2]*spinor2[3];
//	spinor2[2]	= result[0][2] + result[2][0] - spinor1[1]*spinor2[3] + spinor1[3]*spinor2[1];
//	spinor2[3]	= spinor1[0]*spinor2[3] + spinor1[3]*spinor2[0] - result[2][1] + result[1][2];
//
//}
//	original above
//void spinorproduct(double* spinor1,double* spinor2,double *temp) {
//	double result[4][4];
//	for(int i = 0; i < 4; i++){
//		for(int j = 0; j < 4; j++){
//			result[j][i] = spinor2[j]*spinor1[i];
//		}
//	}
//	spinor2[0]	= result[0][0] - result[1][1] - result[2][2] - result[3][3];
//	spinor2[1]	= result[0][1] + result[1][0] - result[3][2] + result[2][3];
//	spinor2[2]	= result[0][2] + result[2][0] - result[1][3] + result[3][1];
//	spinor2[3]	= result[0][3] + result[3][0] - result[2][1] + result[1][2];
//}


	//	Work code spinorproduct

	//  Does not allow self referencing: spinor3 = spinor1

	//	spinor3[0]	=	spinor1[0]*spinor2[0] - spinor1[1]*spinor2[1]  -spinor1[2]*spinor2[2] -spinor1[3]*spinor2[3];
	//	spinor3[1]	=	spinor1[0]*spinor2[1] + spinor1[1]*spinor2[0] - spinor1[3]*spinor2[2] + spinor1[2]*spinor2[3];
	//	spinor3[2]	=	spinor1[0]*spinor2[2] + spinor1[2]*spinor2[0] - spinor1[1]*spinor2[3] + spinor1[3]*spinor2[1];
	//	spinor3[3]	=	spinor1[0]*spinor2[3] + spinor1[3]*spinor2[0] - spinor1[2]*spinor2[1] + spinor1[1]*spinor2[2];

	//	geom22(spinor1,spinor2,spinor3);  //	secondary function calling slows down?

	//	spinor3[0]	=	spinor1[0]*spinor2[0] - dotprod(&spinor1[1],&spinor2[1]);
	//	spinor3[0]	=	spinor1[0]*spinor2[0] - spinor1[1]*spinor2[1]  -spinor1[2]*spinor2[2] -spinor1[3]*spinor2[3];		// much slower with GNU compiler
	//	geom22(spinor1,spinor2,spinor3);									//	secondary function calling slows down?

	//	Does allow self referencing but is not fastest---above code is faster

	//	double result[4][4];
	//	for(int i = 0; i < 4; i++){
	//		for(int j = 0; j < 4; j++){
	//			result[j][i] = spinor1[i]*spinor2[j];
	//		}
	//	}
	//	spinor3[0]	= result[0][0] -result[1][1] -  result[2][2]-  result[3][3];
	//	spinor3[1]	= result[0][1] + result[1][0] - result[3][2] + result[2][3];
	//	spinor3[2]	= result[0][2] + result[2][0] - result[1][3] + result[3][1];
	//	spinor3[3]	= result[0][3] + result[3][0] - result[2][1] + result[1][2];





//void spinorproduct(double *spinor1,double *axis2, double angle2,double *spinor3) {
//	double sin_value = sin(angle2/2.0);
//	double cos_value = cos(angle2/2.0);
//	double cross[3];
//	normalize(axis2,axis2);
//	for (int i = 0; i < 3; i++) {
//		axis2[i] *= sin_value;
//	}
//
//	double axis1[3] = {spinor1[2],-spinor1[3],spinor1[1]};
//	crossprod(axis1,axis2,cross);
//	spinor3[0] = cos_value*spinor1[0] - dotprod(axis1,axis2);
//	spinor3[2] = 	spinor1[0]*axis2[0] 		+ cos_value*axis1[0] 		-cross[0];
//	spinor3[3] =  -(spinor1[0]*axis2[1] 		+ cos_value*axis1[1] 		-cross[1]);
//	spinor3[1] = 	spinor1[0]*axis2[2] 		+ cos_value*axis1[2] 		-cross[2];
//}


/*======================================================================
  *QUATERNION.CPP
  *
  *Created on: July 16, 2012
  *Author:     Pieter Chys
  *
  *========================================================================
  *
  *
  *
  *QUATERNION.CPP constructs a standard quaternion based on the input
  *rotation axis and the rotation angle that are given.  The construction
  *of a rotation quaternion is isomorphic to that of a spinor and it
  *is usefull to take a look at spinor code as well.
  *------------------------------------------------------------------------
  *INPUT:  - double *pointer to UNIT_AXIS which is the rotation axis
  *        - double ANGLE2 which is the rotation angle in radians
  *------------------------------------------------------------------------
  *OUTPUT: - explicitly VOID return
  *		   - double *pointer to QUAT which is the constructed quaternion
  *------------------------------------------------------------------------
  *REMARK:- the unit_axis is a 3-D vector which must be normalized
  *         (normalize rotation axis vector)
  *======================================================================*/
 void quaternion(double *unit_axis,double angle2,double *quat)
 {
 	double sin_value = sin(angle2);									//	code modified for half angle input
 	quat[0] =  cos(angle2);
 	quat[1] =  sin_value*unit_axis[0];
 	quat[2] =  sin_value*unit_axis[1];
 	quat[3] =  sin_value*unit_axis[2];
 }


 /*======================================================================
  *CONJQUAT.CPP
  *
  * Created on: July 16, 2012
  * Author: Pieter Chys
  *
  *========================================================================
  *
  *
  *
  *CONJQUAT.CPP constructs a conjugate quaternion based on an input
  *quaternion which constitutes a couple of rotation quaternions.
  *------------------------------------------------------------------------
  *INPUT:  - double *pointer to QUAT1 which is the input quaternion
  *------------------------------------------------------------------------
  *OUTPUT: - explicitly VOID return
  *		   - double *pointer to QUAT2 which is the output conjugated spinor
  *------------------------------------------------------------------------
  *REMARK:-
  *======================================================================*/
/* void conjquat(double *quat1,double *quat2)
 {
 	quat2[0] = 	 quat2[0];
 	quat2[1] = - quat2[1];
 	quat2[2] = - quat2[2];
 	quat2[3] = - quat2[3];
 }
*/

 /*=======================================================================
  *QUATMATRIX.CPP
  *
  * Created on: July 16, 2012
  * Author: Pieter Chys
  *
  *========================================================================
  *
  *QUATMATRIX.CPP constructs a 3-by-3 matrix on the basis of the input
  *quaternion. The resulting matrix is a rotation matrix which can
  *be used to rotate a given 3-dimensional vector.
  *-----------------------------------------------------------------------
  *INPUT:   - double *pointer to QUAT which is the input quaternion
  *-----------------------------------------------------------------------
  *OUTPUT:	- explicitly VOID return
  *	        - double **pointer to MATRIX which is the rotation matrix (3-by-3)
  *----------------------------------------------------------------------
  *REMARK:	-
  *======================================================================*/
 void quatmatrix(double *quat, double **matrix)
  {
  	 matrix[0][0] =  	2*(quat[0]*quat[0] + quat[1]*quat[1]) - 1;
  	 matrix[0][1] =  	2*(quat[1]*quat[2] - quat[0]*quat[3]);
  	 matrix[0][2] = 	2*(quat[1]*quat[3] + quat[0]*quat[2]);
  	 matrix[1][0] = 	2*(quat[1]*quat[2] + quat[0]*quat[3]);
  	 matrix[1][1] = 	2*(quat[0]*quat[0] + quat[2]*quat[2]) - 1;
  	 matrix[1][2] =     2*(quat[2]*quat[3] - quat[0]*quat[1]);
  	 matrix[2][0] =     2*(quat[1]*quat[3] - quat[0]*quat[2]);
  	 matrix[2][1] = 	2*(quat[2]*quat[3] + quat[0]*quat[1]);
  	 matrix[2][2] = 	2*(quat[0]*quat[0] + quat[3]*quat[3])- 1;
   }

 void quatmatrix(double *quat, double matrix[][3])
  {
  	 matrix[0][0] =  	2*(quat[0]*quat[0] + quat[1]*quat[1]) - 1;
  	 matrix[0][1] =  	2*(quat[1]*quat[2] - quat[0]*quat[3]);
  	 matrix[0][2] = 	2*(quat[1]*quat[3] + quat[0]*quat[2]);
  	 matrix[1][0] = 	2*(quat[1]*quat[2] + quat[0]*quat[3]);
  	 matrix[1][1] = 	2*(quat[0]*quat[0] + quat[2]*quat[2]) - 1;
  	 matrix[1][2] =     2*(quat[2]*quat[3] - quat[0]*quat[1]);
  	 matrix[2][0] =     2*(quat[1]*quat[3] - quat[0]*quat[2]);
  	 matrix[2][1] = 	2*(quat[2]*quat[3] + quat[0]*quat[1]);
  	 matrix[2][2] = 	2*(quat[0]*quat[0] + quat[3]*quat[3])- 1;
   }


 /*=====================================================================
  * SPINORPRODUCT.CPP
  *
  *  Created on:  March 8, 2010
  *      Author: Pieter Chys
  *
  *=====================================================================
  *
  *
  *SPINORPRODUCT.CPP is a different implementation that calculates
  *SPINORPRODUCT.CPP the spinorproduct R3 = R1*R2 in specific order
  *SPINORPRODUCT.CPP in which R1 is the spinor represented by pointer *a
  *and R2 the one by *b. The result is passed as a double pointer *c
  *which is the total spinor R3. Instead of using however the second
  *spinor, its unit axis and dihedral angle are directly used in the
  *implementation
  *
  *----------------------------------------------------------------------------
  *
  *INPUT: - pointer to spinor R1 (double *a, a[4])
  *		 - pointer to axis R2 (double *b, b[3])
  *		 - dihedral angle value (double)
  *-------------------------------------------------------------------------
  *
  *OUTPUT: 	- explicitly VOID
  *		    - pointer to total spinor R3 (double *c, c[4])
  *
  *----------------------------------------------------------------------------
  *
  *REMARKS;	- code is not finished at present/ not tested/ not debugged
  *			- code could be faster maybe than approach with two spinors
  *======================================================================
  */
