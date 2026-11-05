/*
 *CHAIN_MANIP.CPP
 *
 *  Created on: Dec 10, 2010
 *      Author: Pieter Chys
 *
 *
 */

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "libpdb/include/Macromolecule.h"
#include "libpdb/include/pdbIter.h"
#include "libpdb/include/Fragment.h"


#include "own_math.h"
#include "ga.h"
#include "statistics.h"
#include "own_io.h"
#include "chain_manip.h"
#include  "matrix.h"
#include  "geometry.h"

using namespace std;

// Values taken from the most populated values (manually) of Dumbrack's maps.
// http://dunbrack.fccc.edu/nmhrcm/KernRegr_CDL_v1.2.1_Oct14-2009.txt




//#define PHIMIN1  -89.2*180/M_PI
//#define PHIMIN2   60.0*180/M_PI
//#define PHIGLY1  -90.9*180/M_PI
//#define PHIGLY2   92.2*180/M_PI
//#define PHIPRO   -65.8*180/M_PI
//#define PSIMIN1  -27.0*180/M_PI
//#define PSIMIN2  138.1*180/M_PI
//#define PSIGLY1   -7.0*180/M_PI
//#define PSIGLY2  179.2*180/M_PI


/*CONSTANTS*/

/********************************************************************************
 * 																				*
 *			1st set of internal coordinates based on reference:					*
 *																				*
 *																				*
 *		 Fersht Alan, Enzyme Structure and Mechanism, 2nd edition,				*
 *		 1985, W.H. Freeman and Company, New York								*
 *		 																		*
 *		 																		*
 *		 Remark: used in original implementation (also predecessor ISAF) 		*
 *		 																		*
 ********************************************************************************/


/*bond lengths (Angstrom)*/
//#define bondNCA2 	 		1.46			/*N-CA bond*/
//#define bondCACO2    		1.51			/*CA-CO bond*/
//#define bondCON2 	 		1.33	  		/*CO-N bond*/
//
//#define bondCC				1.54
// #define bondCO				1.24
//#define bondNH				1.02
//#define bondCH				1.10
//#define bondCN				1.488

/*valence angles (rad)*/
//#define valCONCA2   			122*M_PI/180		/*valence angle CO-N-CA*/
//#define valNCACO2   			111*M_PI/180		/*valence angle N-CA-CO*/
//#define valCACON2   			116*M_PI/180		/*valence angle CA-CO-N*/

//#define valNCH	 			116*M_PI/180		/*correct values?---more or less standard average value chosen: 116 degrees*/
//#define valNCC	 			116*M_PI/180
//#define valCCO	 			116*M_PI/180
//#define valCNH	 			116*M_PI/180
//#define valCCC	 			0.667*M_PI
//#define valCCH	 			0.667*M_PI
//#define valCCN	 			0.667*M_PI			/*based on CCH and CCC; must be around 0.65 rad*/

/*dihedral angles (rad)*/
//#define dihCON				M_PI				/*dihedral for peptide plane = 180 degrees ------correct?*/








/********************************************************************************
 * 																				*
 *			2nd set of internal coordinates based on references:				*
 *																				*
 *			Engh Richard A. and Huber Robert, Accurate bond and angle parameters*
 * 			for X-ray protein structure refinement, Acta Cryst, 1991, A47,		*
 *          392-400																*
 *		                                                                        *
 *		 	Berkholz D. S., Shapovalov, M.V., Dunbrack, R.L. and Karplus P.A.,  *
 *		    Conformation dependence of backbone geometry in proteins, Structure,*
 *		    2009, 17, 1316-1325                                                 *
 *		 																		*
 *		 	Remark:	- Used as result of need to compare with Soto paper			*
 *		 			- Sigma's also used as perturbations later on				*
 *		 			- Check also Dunbrack values and variations					*
 *																				*
 ********************************************************************************/





////REGULAR RESIDUES (EXCEPT GLYCINE AND PROLINE)
//
//
//
//
////BOND LENGTHS (Angstrom)
//#define bondNCA 	 		1.458			/*N-CA bond      sigma	+- 0.019	*/
//#define bondCACO    		1.525			/*CA-CO bond	 sigma  +- 0.021	*/
//#define bondCON 	 		1.329	  		/*CO-N bond		 sigma	+- 0.014	*/
//
//#define sigma_bondNCA 	 	0.019
//#define sigma_bondCACO    	0.021
//#define sigma_bondCON 	 	0.014
//
//
////VALENCE ANGLES (rad)
//#define valCONCA   			121.7*M_PI/180		/*valence angle CO-N-CA	       */
//#define valNCACO   			111.2*M_PI/180		/*valence angle N-CA-CO	       */
//#define valCACON   			116.2*M_PI/180		/*valence angle CA-CO-N	       */
//
//#define sigma_valCONCA_EH   	1.8*M_PI/180	// VALUES FROM ENGH AND HUBER FOR SIGMA
//#define sigma_valNCACO_EH   	2.8*M_PI/180
//#define sigma_valCACON_EH   	2.0*M_PI/180
//
//#define sigma_valCONCA_B   		1.7*M_PI/180	//	VALUES FROM BERKHOLZ RELATED TO ENGH AND HUBER
//#define sigma_valNCACO_B   		1.5*M_PI/180    //  (SMALLER SIGMA AND SMALLER RANGE, OPTIMIZED VALUES)
//#define sigma_valCACON_B   		1.3*M_PI/180
//
//
//
//
//
////GLYCINE
//
////BOND LENGTHS (Angstrom)
//#define GLYbondNCA 	 			1.451			/*N-CA bond      	*/
//#define GLYbondCACO    			1.516			/*CA-CO bond	 	*/
//#define GLYbondCON 	 			1.329	  		/*CO-N bond		 	*/	// REGULAR VALUE
//
//#define sigma_GLYbondNCA 	 	0.016
//#define sigma_GLYbondCACO    	0.018
//#define sigma_GLYbondCON 	 	0.014									//  REGULAR VALUE
//
//
////VALENCE ANGLES (rad) REGULAR PEPTIDES
//#define GLYvalCONCA   			120.6*M_PI/180		/*valence angle CO-N-CA	     */
//#define GLYvalNCACO   			112.5*M_PI/180		/*valence angle N-CA-CO	     */
//#define GLYvalCACON   			116.4*M_PI/180		/*valence angle CA-CO-N	     */
//
//#define sigma_GLYvalCONCA   	1.7*M_PI/180
//#define sigma_GLYvalNCACO   	2.9*M_PI/180
//#define sigma_GLYvalCACON   	2.1*M_PI/180
//
//
//
//
////PROLINE
//
////BOND LENGTHS (Angstrom)
//#define PRObondNCA 	 			1.466			/*N-CA bond  	*/
//#define PRObondCACO    			1.525			/*CA-CO bond	*/	// REGULAR VALUE
//#define PRObondCON 	 			1.341	  		/*CO-N bond		*/
//
//#define sigma_PRObondNCA 	 	0.015
//#define sigma_PRObondCACO    	0.021								// REGULAR VALUE
//#define sigma_PRObondCON 	 	0.016
//
//
////VALENCE ANGLES (rad) REGULAR PEPTIDES
//#define PROvalCONCA   			122.6*M_PI/180		/*valence angle CO-N-CA	    */
//#define PROvalNCACO   			111.8*M_PI/180		/*valence angle N-CA-CO	    */
//#define PROvalCACON   			116.9*M_PI/180		/*valence angle CA-CO-N	    */
//
//#define sigma_PROvalCONCA   	5.0*M_PI/180
//#define sigma_PROvalNCACO   	2.5*M_PI/180
//#define sigma_PROvalCACON   	1.5*M_PI/180


// Import standard values for internal coordinate values but not for anchor extremities*/
void standic(double *bond_length,double *valence_angle,int nratoms,char *seq, char seq_ct1)
{
	int i;
	for (i = 0; i < (nratoms/3)-1; i++) // Screen all loop residues
	{
		if( !(seq[i+1] == 'P' || seq[i+1] == 'p') ) // NOT pre-PRO
		{
			if (seq[i] == 'I' || seq[i] == 'V') // Iso-Leucine or Valine
			{
				// Bond lengths
				bond_length[3*i] 		= 	IVbondCON;  // CN
				bond_length[3*i+1] 	   	= 	IVbondNCA;  // NA
				bond_length[3*i+2]     	= 	IVbondCACO; // AC
				// Bond angles
				valence_angle[3*i+3]     = 	IVvalCACON; // ACN (of the next residue)
				valence_angle[3*i+1]    = 	valCONCA;   // CNA
				valence_angle[3*i+2]    = 	valNCACO;   // NAC
			}
			else if (seq[i] == 'G') // Glycine
			{
				// Bond lengths
				bond_length[3*i] 		= 	GLYbondCON;  // CN
				bond_length[3*i+1] 	   	= 	GLYbondNCA;  // NA
				bond_length[3*i+2]     	= 	GLYbondCACO; // AC
				// Bond angles
				valence_angle[3*i+3]     = 	GLYvalCACON; // ACN (of the next residue)
				valence_angle[3*i+1]    = 	GLYvalCONCA;   // CNA
				valence_angle[3*i+2]    = 	GLYvalNCACO;   // NAC
			}
			else if (seq[i] == 'P' || seq[i] == 'p') // Proline
			{
				// Bond lengths
				bond_length[3*i] 		= 	PRObondCON;  // CN
				bond_length[3*i+1] 	   	= 	PRObondNCA;  // NA
				bond_length[3*i+2]     	= 	PRObondCACO; // AC
				// Bond angles
				valence_angle[3*i+3]     = 	PROvalCACON; // ACN (of the next residue)
				valence_angle[3*i+1]    = 	PROvalCONCA;   // CNA
				valence_angle[3*i+2]    = 	PROvalNCACO;   // NAC
			}
			else // NOT I, V, G or P (by default)
			{
				// Bond lengths
				bond_length[3*i] 		= 	bondCON;  // CN
				bond_length[3*i+1] 	   	= 	bondNCA;  // NA
				bond_length[3*i+2]     	= 	bondCACO; // AC
				// Bond angles
				valence_angle[3*i+3]     = 	valCACON; // ACN (of the next residue)
				valence_angle[3*i+1]    = 	valCONCA;   // CNA
				valence_angle[3*i+2]    = 	valNCACO;   // NAC
			}
		}
		else // NOT pre-PRO
		{
			if (seq[i] == 'I' || seq[i] == 'V') // Iso-Leucine or Valine
			{
				// Bond lengths
				bond_length[3*i] 		= 	pPRO_IVbondCON;  // CN
				bond_length[3*i+1] 	   	= 	pPRO_IVbondNCA;  // NA
				bond_length[3*i+2]     	= 	pPRO_IVbondCACO; // AC
				// Bond angles
				valence_angle[3*i+3]     = 	pPRO_IVvalCACON; // ACN (of the next residue)
				valence_angle[3*i+1]    = 	pPRO_valCONCA;   // CNA
				valence_angle[3*i+2]    = 	pPRO_valNCACO;   // NAC

			}
			else if (seq[i] == 'G') // Glycine
			{
				// Bond lengths
				bond_length[3*i] 		= 	pPRO_GLYbondCON;  // CN
				bond_length[3*i+1] 	   	= 	pPRO_GLYbondNCA;  // NA
				bond_length[3*i+2]     	= 	pPRO_GLYbondCACO; // AC
				// Bond angles
				valence_angle[3*i+3]     = 	pPRO_GLYvalCACON; // ACN (of the next residue)
				valence_angle[3*i+1]    = 	pPRO_GLYvalCONCA;   // CNA
				valence_angle[3*i+2]    = 	pPRO_GLYvalNCACO;   // NAC
			}
			else if (seq[i] == 'P' || seq[i] == 'p') // Proline
			{
				// Bond lengths
				bond_length[3*i] 		= 	pPRO_PRObondCON;  // CN
				bond_length[3*i+1] 	   	= 	pPRO_PRObondNCA;  // NA
				bond_length[3*i+2]     	= 	pPRO_PRObondCACO; // AC
				// Bond angles
				valence_angle[3*i+3]     = 	pPRO_PROvalCACON; // ACN (of the next residue)
				valence_angle[3*i+1]    = 	pPRO_PROvalCONCA;   // CNA
				valence_angle[3*i+2]    = 	pPRO_PROvalNCACO;   // NAC
			}
			else // NOT I, V, G or P (by default)
			{
				// Bond lengths
				bond_length[3*i] 		= 	pPRO_bondCON;  // CN
				bond_length[3*i+1] 	   	= 	pPRO_bondNCA;  // NA
				bond_length[3*i+2]     	= 	pPRO_bondCACO; // AC
				// Bond angles
				valence_angle[3*i+3]     = 	pPRO_valCACON; // ACN (of the next residue)
				valence_angle[3*i+1]    = 	pPRO_valCONCA;   // CNA
				valence_angle[3*i+2]    = 	pPRO_valNCACO;   // NAC
			}
		}
	}

	// Ct anchor residue
	i=(nratoms/3)-1;
	if( !(seq_ct1 == 'P' || seq_ct1 == 'p') ) // NOT pre-PRO
	{
		if (seq[i] == 'I' || seq[i] == 'V') // Iso-Leucine or Valine
		{
			// Bond lengths
			bond_length[3*i] 		= 	IVbondCON;  // CN
			bond_length[3*i+1] 	   	= 	IVbondNCA;  // NA
			bond_length[3*i+2]     	= 	IVbondCACO; // AC
			// Bond angles
			valence_angle[3*i+1]    = 	IVvalCONCA;   // CNA
			valence_angle[3*i+2]    = 	IVvalNCACO;   // NAC
		}
		else if (seq[i] == 'G') // Glycine
		{
			// Bond lengths
			bond_length[3*i] 		= 	GLYbondCON;  // CN
			bond_length[3*i+1] 	   	= 	GLYbondNCA;  // NA
			bond_length[3*i+2]     	= 	GLYbondCACO; // AC
			// Bond angles
			valence_angle[3*i+1]    = 	GLYvalCONCA;   // CNA
			valence_angle[3*i+2]    = 	GLYvalNCACO;   // NAC
		}
		else if (seq[i] == 'P' || seq[i] == 'p') // Proline
		{
			// Bond lengths
			bond_length[3*i] 		= 	PRObondCON;  // CN
			bond_length[3*i+1] 	   	= 	PRObondNCA;  // NA
			bond_length[3*i+2]     	= 	PRObondCACO; // AC
			// Bond angles
			valence_angle[3*i+1]    = 	PROvalCONCA;   // CNA
			valence_angle[3*i+2]    = 	PROvalNCACO;   // NAC
		}
		else // NOT I, V, G or P (by default)
		{
			// Bond lengths
			bond_length[3*i] 		= 	bondCON;  // CN
			bond_length[3*i+1] 	   	= 	bondNCA;  // NA
			bond_length[3*i+2]     	= 	bondCACO; // AC
			// Bond angles
			valence_angle[3*i+1]    = 	valCONCA;   // CNA
			valence_angle[3*i+2]    = 	valNCACO;   // NAC
		}
	}
	else // pre-PRO
	{
		if (seq[i] == 'I' || seq[i] == 'V') // Iso-Leucine or Valine
		{
			// Bond lengths
			bond_length[3*i] 		= 	pPRO_IVbondCON;  // CN
			bond_length[3*i+1] 	   	= 	pPRO_IVbondNCA;  // NA
			bond_length[3*i+2]     	= 	pPRO_IVbondCACO; // AC
			// Bond angles
			valence_angle[3*i+1]    = 	pPRO_IVvalCONCA;   // CNA
			valence_angle[3*i+2]    = 	pPRO_IVvalNCACO;   // NAC

		}
		else if (seq[i] == 'G') // Glycine
		{
			// Bond lengths
			bond_length[3*i] 		= 	pPRO_GLYbondCON;  // CN
			bond_length[3*i+1] 	   	= 	pPRO_GLYbondNCA;  // NA
			bond_length[3*i+2]     	= 	pPRO_GLYbondCACO; // AC
			// Bond angles
			valence_angle[3*i+1]    = 	pPRO_GLYvalCONCA;   // CNA
			valence_angle[3*i+2]    = 	pPRO_GLYvalNCACO;   // NAC
		}
		else if (seq[i] == 'P' || seq[i] == 'p') // Proline
		{
			// Bond lengths
			bond_length[3*i] 		= 	pPRO_PRObondCON;  // CN
			bond_length[3*i+1] 	   	= 	pPRO_PRObondNCA;  // NA
			bond_length[3*i+2]     	= 	pPRO_PRObondCACO; // AC
			// Bond angles
			valence_angle[3*i+1]    = 	pPRO_PROvalCONCA;   // CNA
			valence_angle[3*i+2]    = 	pPRO_PROvalNCACO;   // NAC
		}
		else // NOT I, V, G or P (by default)
		{
			// Bond lengths
			bond_length[3*i] 		= 	pPRO_bondCON;  // CN
			bond_length[3*i+1] 	   	= 	pPRO_bondNCA;  // NA
			bond_length[3*i+2]     	= 	pPRO_bondCACO; // AC
			// Bond angles
			valence_angle[3*i+1]    = 	pPRO_valCONCA;   // CNA
			valence_angle[3*i+2]    = 	pPRO_valNCACO;   // NAC
		}
	}

	//INTERNAL COORDINATES OF ANCHORS
	bond_length[0]				= 	0.0;
	valence_angle[0]			= 	0.0;
	valence_angle[1]			= 	0.0;

	// Positions 1 and 2 of bond_lengh[] are overwritten outside this function (with anchor values)
	// Positions 2 of valence_angle[] is overwritten outside this function (with anchor values)
	// These position should be filled with standard values here...???? Or inserted from input PDB...(like in randomizeIC function)
	// Idem for C-terminal end...

}



void standic_old(double *bond_length,double *valence_angle,int nratoms,int* residuemarker)			/*import standard values for internal coordinate values but not for anchor extremities*/
{
	//BULK OF THE INTERNAL COORDINATES
	for (int i = 1; i < (nratoms/3); i++) // screen non-Nt-anchor aminoacids
	{
		bond_length[3*i] 		= 	bondCON;  // CN
		bond_length[3*i+1] 	   	= 	bondNCA;  // NA
		bond_length[3*i+2]     	= 	bondCACO; // AC

		valence_angle[3*i]      = 	valCACON; // ACN
		valence_angle[3*i+1]    = 	valCONCA; // CNA
		valence_angle[3*i+2]    = 	valNCACO; // NAC
	}

	//INTERNAL COORDINATES OF ANCHORS
	bond_length[0]				= 	0.0;
	valence_angle[0]			= 	0.0;
	valence_angle[1]			= 	0.0;

	// Positions 1 and 2 of bond_lengh[] are filled outside this function (with anchor values)
	// Positions 2 of valence_angle[] is filled outside this function (with anchor values)
	// These possition should be filled with standard values here...???? Or inserted from input PDB...(like in randomizeIC function)
	// Idem for C-terminal end...

	// MODIFICATION FOR GLYCINES AND PROLINES
	for (int i = 1; i < (nratoms/3) - 1; i++ )
	{
		if (residuemarker[i] == 1)							//FINDING GLYCINES
		{
			// Bug: 8/4/2016 "Some small inaccuracy with Gly and Pro..."
			bond_length[3*i] 		= 	GLYbondCON;
			bond_length[3*i+1] 	   	= 	GLYbondNCA;
			bond_length[3*i+2]     	= 	GLYbondCACO;
			//			bond_length[3*i+3] 		= 	GLYbondCON;

			valence_angle[3*i]    = 	GLYvalCACON;
			valence_angle[3*i+1]    = 	GLYvalCONCA;
			valence_angle[3*i+2]    = 	GLYvalNCACO;
			//			valence_angle[3*i+3]    = 	GLYvalCACON;
		}
		else if (residuemarker[i] == 2)					//FINDING PROLINES
		{
			bond_length[3*i] 		= 	PRObondCON;
			bond_length[3*i+1] 	   	= 	PRObondNCA;
			bond_length[3*i+2]     	= 	PRObondCACO;
			//			bond_length[3*i+3] 		= 	PRObondCON;

			valence_angle[3*i]     = 	PROvalCACON;
			valence_angle[3*i+1]    = 	PROvalCONCA;
			valence_angle[3*i+2]    = 	PROvalNCACO;
			//			valence_angle[3*i+3]     = 	PROvalCACON;
		}
	}
}



//// Obsolete?
//void randic(double *bond_length,double *valence_angle,int nratoms,double lengthptb,double valenceptb,int *residuemarker,int randbond,int randval)	/*randomize internal coordinate values but not for anchor extremities*/
//{
//
//
//	//BULK OF THE INTERNAL COORDINATES
//	for (int i = 1; i < (nratoms/3) - 1; i++)
//	{
//		if (randbond == 1)												//:TODO:	potential performance issue due to branching? ---- also everywhere below...
//		{
//			bond_length[3*i] 		= 	bondCON   +   lengthptb*2.0 * sigma_bondCON  * ( ((double) rand()/(double) RAND_MAX) - 1.0);		//:TRICKY: BRACKET FORMATTING IS CRITICAL
//			bond_length[3*i+1] 	   	= 	bondNCA   +   lengthptb*2.0 * sigma_bondNCA  * ( ((double) rand()/(double) RAND_MAX) - 1.0);		//:TRICKY: using 2 x SIGMA (95% confidence interval)
//			bond_length[3*i+2]     	= 	bondCACO  +   lengthptb*2.0 * sigma_bondCACO * ( ((double) rand()/(double) RAND_MAX) - 1.0);
//		}
//
//		if (randval == 1)
//		{
//			valence_angle[3*i]      = 	valCACON  +   valenceptb*2.0 * sigma_valCACON_B *( ((double) rand()/(double) RAND_MAX) - 1.0);		//:TRICKY: values from Berkholz which are improved values upon Engh and Huber
//			valence_angle[3*i+1]    = 	valCONCA  +   valenceptb*2.0 * sigma_valCONCA_B *( ((double) rand()/(double) RAND_MAX) - 1.0);
//			valence_angle[3*i+2]    = 	valNCACO  +   valenceptb*2.0 * sigma_valNCACO_B *( ((double) rand()/(double) RAND_MAX) - 1.0);
//		}
//	}
//
//
//
//
//
//
//	//NEW CODE
//	//INTERNAL COORDINATES OF ANCHORS
//	if (randbond == 1)
//	{
//		bond_length[0]				= 	0.0;																									//ROW 1
//		bond_length[nratoms-3]		= 	bondCON   +  lengthptb*2.0 * sigma_bondCON  *  ( ((double) rand()/(double) RAND_MAX) - 1.0);			//3rd LAST LINE
//	}
//
//
//	if (randval == 1)
//	{
//		valence_angle[0]			= 	0.0;																								//ROW 1
//		valence_angle[1]			= 	0.0;																								//ROW 2
//		valence_angle[nratoms-2]	= 	valCONCA  +  valenceptb*2.0 * sigma_valCONCA_B * ( ((double) rand()/(double) RAND_MAX) - 1.0);		//2nd LAST LINE
//		valence_angle[nratoms-3]	= 	valCACON  +  valenceptb*2.0 * sigma_valCACON_B * ( ((double) rand()/(double) RAND_MAX) - 1.0);		//3rd LAST LINE
//	}
//
//
//
//
//
//	//ORIGINAL CODE
//	//INTERNAL COORDINATES OF ANCHORS
////	bond_length[0]				= 	0.0;	//ROW 1
////	valence_angle[0]			= 	0.0;	//ROW 1
////
////	valence_angle[1]			= 	0.0;	//ROW 2
////
////	bond_length[nratoms-3]		= 	bondCON   +  lengthptb*2.0 * sigma_bondCON  *  ( ((double) rand()/(double) RAND_MAX) - 1.0);		//3rd LAST LINE
////	valence_angle[nratoms-3]	= 	valCACON  +  valenceptb*2.0 * sigma_valCACON * ( ((double) rand()/(double) RAND_MAX) - 1.0);			//3rd LAST LINE
////
////	valence_angle[nratoms-2]	= 	valCONCA  +  valenceptb*2.0 * sigma_valCONCA * ( ((double) rand()/(double) RAND_MAX) - 1.0);		//2nd LAST LINE
//
//
//
//
//
//	//MODIFICATION FOR GLYCINES AND PROLINES
//	for (int i = 1; i < (nratoms/3) - 1; i++ )
//	{
//		if (residuemarker[i] == 1)						//FINDING GLYCINES
//		{
//			if (randbond == 1)
//			{
//				bond_length[3*i+1] 	   	= 	GLYbondNCA		+	lengthptb*2.0 * sigma_GLYbondNCA   * ( ((double) rand()/(double) RAND_MAX) - 1.0 );			//:TRICKY:	different order and indexing than above
//				bond_length[3*i+2]     	= 	GLYbondCACO		+	lengthptb*2.0 * sigma_GLYbondCACO  * ( ((double) rand()/(double) RAND_MAX) - 1.0 );
//				bond_length[3*i+3] 		= 	GLYbondCON		+	lengthptb*2.0 * sigma_GLYbondCON   * ( ((double) rand()/(double) RAND_MAX) - 1.0 );
//			}
//
//			if (randval == 1)
//			{
//				valence_angle[3*i+1]    = 	GLYvalCONCA		+	valenceptb*2.0 * sigma_GLYvalCONCA  * ( ((double) rand()/(double) RAND_MAX) - 1.0 );		//:TRICKY:	different order and indexing than above
//				valence_angle[3*i+2]    = 	GLYvalNCACO		+	valenceptb*2.0 * sigma_GLYvalNCACO  * ( ((double) rand()/(double) RAND_MAX) - 1.0 );		//:TODO:    check and discuss with Pablo if correct?
//				valence_angle[3*i+3]    = 	GLYvalCACON		+	valenceptb*2.0 * sigma_GLYvalCACON  * ( ((double) rand()/(double) RAND_MAX) - 1.0 );
//			}
//
//
//		} else if (residuemarker[i] == 2)				//FINDING PROLINES
//		{
//			if (randbond == 1)
//			{
//				bond_length[3*i+1] 	   	= 	PRObondNCA		+ 	lengthptb*2.0 * sigma_PRObondNCA  * ( ((double) rand()/(double) RAND_MAX) - 1.0 );			//:TRICKY:	different order and indexing than above
//				bond_length[3*i+2]     	= 	PRObondCACO		+ 	lengthptb*2.0 * sigma_PRObondCACO * ( ((double) rand()/(double) RAND_MAX) - 1.0 ); 			//:TODO:    check and discuss with Pablo if correct?
//				bond_length[3*i+3] 		= 	PRObondCON		+	lengthptb*2.0 * sigma_PRObondCON  * ( ((double) rand()/(double) RAND_MAX) - 1.0 );
//			}
//
//			if (randval == 1)
//			{
//				valence_angle[3*i+1]    = 	PROvalCONCA		+	valenceptb*2.0 * sigma_PROvalCONCA * ( ((double) rand()/(double) RAND_MAX) - 1.0 );
//				valence_angle[3*i+2]    = 	PROvalNCACO		+	valenceptb*2.0 * sigma_PROvalNCACO * ( ((double) rand()/(double) RAND_MAX) - 1.0 );
//				valence_angle[3*i+3]     = 	PROvalCACON		+	valenceptb*2.0 * sigma_PROvalCACON  *( ((double) rand()/(double) RAND_MAX) - 1.0 );
//			}
//
//		}
//	}
//
//
//}
//
//
//
//
//
//// Obsolete?
//void randic(double *bond_length,double *valence_angle,double *dihedral_angle_i,int nratoms,double lengthptb,double valenceptb,double *prearray1,double *prearray2,double *prearray3,double *prearray4)
//{
//
//	for (int i = 1; i < (nratoms/3) - 1; i++)	//:REMARK: must i start from 1 or 0?
//	{
//		bond_length[3*i] 		= 	bondCON   +  lengthptb*2.0*( ((double) rand()/(double) RAND_MAX) - 1.0);
//		bond_length[3*i+1] 	   	= 	bondNCA   +  lengthptb*2.0*( ((double) rand()/(double) RAND_MAX) - 1.0);
//		bond_length[3*i+2]     	= 	bondCACO  +  lengthptb*2.0*( ((double) rand()/(double) RAND_MAX) - 1.0);
//
//		valence_angle[3*i]      = 	valCACON  +  valenceptb*2.0*( ((double) rand()/(double) RAND_MAX) - 1.0);
//		valence_angle[3*i+1]    = 	valCONCA  +  valenceptb*2.0*( ((double) rand()/(double) RAND_MAX) - 1.0);
//		valence_angle[3*i+2]    = 	valNCACO  +  valenceptb*2.0*( ((double) rand()/(double) RAND_MAX) - 1.0);
//	}
//
//
//
//	/*Resetting array elements to zero*/
//	bond_length[nratoms-3]		= 	bondCON   + lengthptb*2.0* ( ((double) rand()/(double) RAND_MAX) - 1.0);
//	valence_angle[nratoms-3]	= 	valCACON  + valenceptb*2.0* ( ((double) rand()/(double) RAND_MAX) - 1.0);
//	//dihedral_angle[nr_atoms_loop-3]		= 	((double) rand()/(double) RAND_MAX)*M_PI*2 - M_PI;
//
//
//
//	/*Resetting array elements to zero*/
//	bond_length[0]		= 0.0;
//	valence_angle[0]	= 0.0;
//	//dihedral_angle[0]	= 0.0;
//
//	//bond_length[1]	= 0.0;
//	valence_angle[1]	= 0.0;
//	//dihedral_angle[1]	= 0.0;
//
//	//bond_length[2]		= 0.0;
//	//valence_angle[2]	= 0.0;
//	//dihedral_angle[2]	= 0.0;
//
//	for (int i = 0; i < nratoms; i++)
//	{
//		prearray1[i] 	    = cos( (M_PI - valence_angle[i])/2.0 );
//		prearray2[i] 	    = sin( (M_PI - valence_angle[i])/2.0 );
//		prearray3[i]		= cos ( dihedral_angle_i[i]/2.0);
//		prearray4[i]		= sin ( dihedral_angle_i[i]/2.0);
//	}
//
//}


//	dihedral_angle[3*i]     = ((double) rand()/(double) RAND_MAX)*M_PI*2 - M_PI;
//	dihedral_angle[3*i+1]   = dihCON  + dihptb*(2.0*((double) rand()/(double) RAND_MAX) - 1.0);
//	dihedral_angle[3*i+2]   = ((double) rand()/(double) RAND_MAX)*M_PI*2 - M_PI;
//
//	correctdihedral(&dihedral_angle[3*i]);
//	correctdihedral(&dihedral_angle[3*i+1]);
//	correctdihedral(&dihedral_angle[3*i+2]);



int checkself(double** co,double **p1co,int check, int index,int nr_atoms_loop,int atomshift,double clash_factor,int reverse)
{


	index 		= 0;									/*pure backbone collisions*/
	while (check == 0 && index < nr_atoms_loop)
	{
		for (int j = index + atomshift; check != 1 && j < nr_atoms_loop; j++)	/*(i + 4) to avoid atoms close to each other*/
		{
			if ( distance(co[index],co[j]) < clash_factor*(co[index][4] + co[j][4]) )	//||			/*Line only for backbone*/
				//									distance(p1co->el[kindex2],p1co->el[j]) < clash_factor*(p1co->el[kindex2][4] + p1co->el[j][4]) ||		/*Parallel chain 1 clashes*/
				//										(distance(p2co->el[kindex2],p2co->el[j]) < clash_factor*(p2co->el[kindex2][4] + p2co->el[j][4]) && p2co->el[kindex2][4] != 0 && p2co->el[j][4] != 0) ||		/*Parallel chain 2 clashes*/
				//											distance(co->el[k*nr_atoms_loop+kindex2],p1co->el[k*nr_atoms_loop+j]) < clash_factor*(co->el[k*nr_atoms_loop+kindex2][4] + p1co->el[k*nr_atoms_loop+j][4]) ||	/*Line for 1st parallel chain*/
				//												(distance(co->el[k*nr_atoms_loop+kindex2],p2co->el[k*nr_atoms_loop+j]) < clash_factor*(co->el[k*nr_atoms_loop+kindex2][4] + p2co->el[k*nr_atoms_loop+j][4]) &&  p2co->el[k*nr_atoms_loop+j][4] != 0)  || /*Line for 2nd parallel chain*/
				//													(distance(p1co->el[k*nr_atoms_loop+kindex2],p2co->el[k*nr_atoms_loop+j]) < clash_factor*(p1co->el[k*nr_atoms_loop+kindex2][4] + p2co->el[k*nr_atoms_loop+j][4]) &&  p1co->el[k*nr_atoms_loop+kindex2][4] != 0 &&  p2co->el[k*nr_atoms_loop+j][4] != 0) ) /*Line for 2nd parallel chain*/
			{
				check		   	= 1;
			}
		}

		index++;
	}


	if (check == 0)
	{
		index = 5;							/*oxygen-backbone collisions*/
		if (reverse == 1)					//:TODO: Is correct index value?
		{
			index = 3;						/*position starting from carboxy-terminus for oxygen*/
		}

		while (check == 0 && index < nr_atoms_loop)
		{

			for (int j = index + atomshift; check != 1 && j < nr_atoms_loop; j++)	/*(i + 4) to avoid atoms close to each other*/
			{
				if ( distance(p1co[index],co[j]) < clash_factor*(p1co[index][4] + co[j][4]) )

				{
					check		   	= 1;
				}
			}

			if (check == 0)
			{
				index += 3;
			}
		}
	}

	if (check == 0)
	{
		index = 5;									/*oxygen-oxygen collisions*/
		if (reverse == 1) //:TODO: Is correct index value?
		{
			index = 3;		/*position starting from carboxy-terminus for oxygen*/
		}

		while (check == 0 && index < nr_atoms_loop)
		{

			for (int j = index + 3; check != 1 && j < nr_atoms_loop; )	/*(i + 4) to avoid atoms close to each other*/
			{
				if ( distance(p1co[index],p1co[j]) < clash_factor*(p1co[index][4] + p1co[j][4]) )

				{
					check		   	= 1;
				} else
				{
					j  += 3;
				}
			}

			if (check == 0)
			{
				index += 3;
			}
		}
	}

	return check;
}





int checkoxyg(double** co, short int ***grid, double **limits,double stepsize,int sizegrid,int check, int index,int nr_atoms_loop,int cutoff,int reverse)
{
	int position[3];

	index 		= 5;
	if (reverse == 1)	//:TODO: FOR REVERSE PROTOCOL
	{
		index	= 3;
	}

	while (check == 0 && index < nr_atoms_loop - cutoff)	//:TRICKY: - 3 is ending at loop itself
	{
		for (int i = 0; i < 3; i++)
		{
			//posi2[i] = ceil( (co->el[kindex2][i] - limits[i][0])/stepsize );
			position[i] = int( (co[index][i] - limits[i][0])/stepsize ); /*for sufficiently large grid size*/
		}

		if( 0 <= position[0] && position[0]  < sizegrid &&
				0 <= position[1] && position[1]  < sizegrid &&
				0 <= position[2] && position[2]  < sizegrid )
		{
			if (grid[position[0]][position[1]][position[2]] == 1 )
			{
				check   		= 1;
			}
			else
			{
				index 	+=	3;	/*No collision and within range*/
			}

		}	else	/*should in principle not be used*/
		{
			index 	+=	3;	/*Out of range*/
		}
	}

	return check;
}


int checkbb(double** co, short int ***grid, double **limits,double stepsize,int sizegrid,int check, int index,int nr_atoms_loop,int cutoff)
{
	int position[3];

	check      	= 0;	/*reset that no collision occurs*/
	index 		= 3;	/*start on end atom anchor*/

	while (check == 0 && index < nr_atoms_loop - cutoff)	//:TRICKY: - 3 is ending at loop itself
	{
		for (int i = 0; i < 3; i++)
		{
			//posi2[i] = ceil( (co->el[kindex2][i] - limits[i][0])/stepsize );
			position[i] = int( (co[index][i] - limits[i][0])/stepsize ); /*for sufficiently large grid size*/
		}

		if( 0 <= position[0] && position[0]  < sizegrid &&
				0 <= position[1] && position[1]  < sizegrid &&
				0 <= position[2] && position[2]  < sizegrid )
		{
			if (grid[position[0]][position[1]][position[2]] == 1 )
			{
				check   		= 1;
			}
			else
			{
				index++;	/*No collision and within range*/
			}

		}	else	/*should in principle not be used*/
		{
			index++;	/*Out of range*/
		}
	}

	return check;
}




//double rachprob(double dihedral, int phipsi, int marker)
//{
//	double prob;	/*probability to be calculated*/
//	prob = 1.0;
//
//	if (phipsi == 0)	//phi bonds
//	{
//		if (marker == 0)	//:INFO: All residues except glycine and proline
//		{
//			prob = 1.0 - 0.67 - (4 - 2*(cos(2*(dihedral - PHIMIN1)) + cos(2*(dihedral - PHIMIN2))) )/
//							(8*M_PI +  ( sin(2*(M_PI-PHIMIN1)) - sin(2*(-M_PI-PHIMIN1)) + sin(2*(M_PI-PHIMIN2)) - sin(2*(-M_PI-PHIMIN2)) ) );
//			prob /= 0.33;
//
//		} else if (marker == 1)  //:INFO: Glycine residues
//		{
//			prob = 1.0 - 0.76 - (4 - 2*(cos(2*(dihedral - PHIGLY1)) + cos(2*(dihedral - PHIGLY2))) )/
//							(8*M_PI +  ( sin(2*(M_PI-PHIGLY1)) - sin(2*(-M_PI-PHIGLY1)) + sin(2*(M_PI-PHIGLY2)) - sin(2*(-M_PI-PHIGLY2)) ) );
//			prob /= 0.16;
//
//
//		} else if (marker == 2)  //:INFO: Proline residues
//		{
//			prob = 1.0 - 0.68 - (1 - (cos((dihedral - PHIPRO)) ) )/
//									(2*M_PI +  ( sin(2*(M_PI-PHIPRO)) - sin(2*(-M_PI-PHIPRO))  ) );
//			prob /= 0.42;
//		}
//
//	} else if (phipsi == 1 ) //psi bonds
//	{
//		if (marker == 0)	//:INFO: All residues except glycine and proline
//		{
//			prob = 1.0 - 0.67 - (4 - 2*(cos(2*(dihedral - PSIMIN1)) + cos(2*(dihedral - PSIMIN2))) )/
//								(8*M_PI +  ( sin(2*(M_PI-PSIMIN1)) - sin(2*(-M_PI-PSIMIN1)) + sin(2*(M_PI-PSIMIN2)) - sin(2*(-M_PI-PSIMIN2)) ) );
//			prob /= 0.33;
//
//
//		} else if (marker == 1)  //:INFO: Glycine residues
//		{
//			prob = 1.0 - 0.68 - (4 - 2*(cos(2*(dihedral - PSIGLY1)) + cos(2*(dihedral - PSIGLY2))) )/
//											(8*M_PI +  ( sin(2*(M_PI-PSIGLY1)) - sin(2*(-M_PI-PSIGLY1)) + sin(2*(M_PI-PSIGLY2)) - sin(2*(-M_PI-PSIGLY2)) ) );
//			prob /= 0.32;
//
//		} else if (marker == 2)  //:INFO: Proline residues
//		{
//
//		}
//	}
//
//	return prob;
//}


void copyco(double **co1,double **co2,int bond, int size)
{
	for (int i = bond - 1; i < size; i++)	//:REMARK: include last atom chain for rotation side chains
	{	for (int j = 0; j < 3; j++)
	{
		co2[i][j] = co1[i][j];
	}
	}
}



void copyco(double **co1,double **co2, int size)
{
	for (int i = 0; i < size; i++)
	{	for (int j = 0; j < 3; j++)
	{
		co2[i][j] = co1[i][j];
	}
	}
}



void pose(double **co1,double **co2,double *quaternion)
{
	/*variables*/
	double *a1;			/*first bond vector  --- coordinates 1*/
	double *b1;			/*second bond vector --- coordinates 1*/
	double *a2;			/*first bond vector  --- coordinates 2*/
	double *b2;			/*second bond vector --- coordinates 2*/
	double *normal1;	/*first normal --- coordinates 1*/
	double *normal2; 	/*second normal --- coordinates 2*/
	double *binormal;	/*normal to plane normal1-normal2*/
	double *tangent;	/*extra vector to check sign angle*/
	double *tangent2;	/*extra vector to check sign angle*/
	double *between;	/*auxiliary quaternion*/
	double angle1;		/*rotation angle quaternion 1*/
	double angle2;		/*rotation angle quaternion 2*/
	double *qt1;		/*1st rotation quaternion*/
	double *cqt1;		/*conjugate 1st rotation quaternion*/
	double *qt2;		/*2nd rotation quaternion*/
	double sine;		/*sine in quaternion*/
	double dotnorm; // Dummy variable to avoid bug...

	//ALLOCATION VECTORS
	a1 = (double*) malloc(sizeof(double)*4);
	for(int i = 0; i < 4; i++)
	{
		a1[i] = 0.0;
	}

	b1 = (double*) malloc(sizeof(double)*3);
	for(int i = 0; i < 3; i++)
	{
		b1[i] = 0.0;
	}

	a2 = (double*) malloc(sizeof(double)*3);
	for(int i = 0; i < 3; i++)
	{
		a2[i] = 0.0;
	}

	b2 = (double*) malloc(sizeof(double)*3);
	for(int i = 0; i < 3; i++)
	{
		b2[i] = 0.0;
	}

	normal1 = (double*) malloc(sizeof(double)*3);
	for(int i = 0; i < 3; i++)
	{
		normal1[i] = 0.0;
	}

	normal2 = (double*) malloc(sizeof(double)*3);
	for(int i = 0; i < 3; i++)
	{
		normal2[i] = 0.0;
	}

	binormal = (double*) malloc(sizeof(double)*3);
	for(int i = 0; i < 3; i++)
	{
		binormal[i] = 0.0;
	}

	tangent = (double*) malloc(sizeof(double)*3);
	for(int i = 0; i < 3; i++)
	{
		tangent[i] = 0.0;
	}

	tangent2 = (double*) malloc(sizeof(double)*3);
	for(int i = 0; i < 3; i++)
	{
		tangent2[i] = 0.0;
	}


	//ALLOCATION QUATERNIONS
	between = (double*) malloc(sizeof(double)*4);
	for(int i = 0; i < 4; i++)
	{
		between[i] = 0.0;
	}

	qt1 = (double*) malloc(sizeof(double)*4);
	for(int i = 0; i < 4; i++)
	{
		qt1[i] = 0.0;
	}

	cqt1 = (double*) malloc(sizeof(double)*4);
	for(int i = 0; i < 4; i++)
	{
		cqt1[i] = 0.0;
	}

	qt2 = (double*) malloc(sizeof(double)*4);
	for(int i = 0; i < 4; i++)
	{
		qt2[i] = 0.0;
	}


	//CREATE BONDVECTORS 2 PLANES
	for(int i = 0; i < 3; i++)
	{
		a1[i+1] = co1[1][i]  - co1[0][i];		/*bond vector 1 --- plane 1*/
		b1[i]   = co1[2][i]  - co1[1][i];		/*bond vector 2 --- plane 1*/

		a2[i]   = co2[1][i]  - co2[0][i]; 		/*bond vector 1 --- plane 2*/
		b2[i]   = co2[2][i]  - co2[1][i]; 		/*bond vector 2 --- plane 2*/
	}


	//CREATE UNIT NORMAL VECTORS
	crossprod(a1+1,b1,normal1);					/*normal plane 1*/
	crossprod(a2,b2,normal2);					/*normal plane 2*/
	normalise(normal1);
	normalise(normal2);
	crossprod(normal1,normal2,binormal); 		/*binormal to normal planes*/
	normalise(binormal);

	// BUG (10/6/2016): Instability de cojones
	dotnorm = dotprod(normal1,normal2) / ( norm(normal1) * norm(normal2) );
	if(dotnorm > 1.0) // if slightly over 1.0000001... or something like this
		angle1 = 0.0;
	else if (dotnorm < -1.0) // if slightly below -1.0000001... or something like this
		angle1 = M_PI; // 180 deg
	else
		//CALCULATE ROTATION ANGLE BETWEEN TWO PLANES
		//	angle1 = acos(dotprod(normal1,normal2)/(norm(normal1)*norm(normal2)));  /*angle between 2 normals*/
		angle1 = acos(dotnorm);  /*angle between 2 normals*/

	//	if( isnan(angle1) )
	//	{
	//		fprintf(stderr,"angle1= %f is NAN,  norm1= %f  norm2= %f\n", angle1, norm(normal1), norm(normal2) );
	//	}

	crossprod(normal1,binormal,tangent);
	normalise(tangent);


	//CORRECT OBTAINED ROTATION ANGLE
	if (dotprod(normal2,tangent) >= 0)
	{
		angle1 = -fabs(angle1);
	} else
	{
		angle1 = fabs(angle1);
	}




	//CALCULATE REPRESENTATIVE QUATERNION
	sine   = sin(angle1/2.0);
	qt1[0] = cos(angle1/2.0);
	qt1[1] =  binormal[0]*sine;
	qt1[2] =  binormal[1]*sine;
	qt1[3] =  binormal[2]*sine;


	//CONJUGATE QUATERNION
	cqt1[0] =  qt1[0];
	cqt1[1] = -qt1[1];
	cqt1[2] = -qt1[2];
	cqt1[3] = -qt1[3];


	//UPDATE BONDVECTOR A1
	quatproduct(a1,cqt1,between);
	quatproduct(qt1,between,a1);

	dotnorm = dotprod(a1+1,a2) / (norm(a1+1)*norm(a2));
	if(dotnorm > 1.0) // if slightly over 1.0000001... or something like this
		angle2 = 0.0;
	else if (dotnorm < -1.0) // if slightly below -1.0000001... or something like this
		angle2 = M_PI; // 180 deg
	else
		//CALCULATE SECOND ROTATION ANGLE (IN-PLANE)
		//:TODO: VERIFY HYPOTHESIS
		//	angle2 = acos(dotprod(a1+1,a2)/(norm(a1+1)*norm(a2)));						/*BUT plane 1 must be rotated first*/
		angle2 = acos(dotnorm);						/*BUT plane 1 must be rotated first*/

	//	if( isnan(angle2) )
	//	{
	//		fprintf(stderr,"angle2= %f is NAN,  norm(a1+1)= %f  norm(a2)= %f\n", angle2, norm(a1+1), norm(a2) );
	//	}

	crossprod(normal2,a1+1,tangent2);;
	normalise(tangent2);

	//CORRECT OBTAINED ROTATION ANGLE
	if(dotprod(tangent2,a2) >= 0)
	{
		angle2 = +fabs(angle2);
	} else
	{
		angle2 = -fabs(angle2);
	}


	sine   = sin(angle2/2.0);
	qt2[0] = cos(angle2/2.0);
	qt2[1] = normal2[0]*sine;
	qt2[2] = normal2[1]*sine;
	qt2[3] = normal2[2]*sine;

	//CALCULATE FINAL QUATERNION
	quatproduct(qt2,qt1,quaternion);


	free(a1);
	free(a2);
	free(b1);
	free(b2);
	free(normal1);
	free(normal2);
	free(tangent);
	free(tangent2);
	free(binormal);
	free(qt1);
	free(qt2);
	free(cqt1);
	free(between);


}


//	Old code
/*===============================================================================
 *
 *	VALENCEFIND.CPP
 *
 *	Created on: April 26, 2011
 *      Author: Pieter Chys
 *
 *===============================================================================
 *
 * FINDVALENCE.CPP calculates the valence angles in an atom system and this array
 * has a size 2 less than the original atom coordinate system.
 *
 *--------------------------------------------------------------------------------
 *
 *  INPUT: - double **coordinates are the absolute atom coordinates in (number_atoms x 3/4)
 *  		 format
 *	 	   - int number_atoms = number atoms of the atom system
 *
 *
 *-------------------------------------------------------------------------------
 *
 *  OUTPUT: - double *valence contains the calculated valence angles based on the
 *           locally created **bonds matrix and has a size 2 less as compared
 *           to **coordinates
 *
 *-------------------------------------------------------------------------------
 *
 *  DEPENDENCIES:	- mind you that double **coordinates is often (size x 4) in
 *  				  format
 *
 *
 *===============================================================================



void valencefind(double **coordinates,double *valence,int number_atoms){
	double bonds[number_atoms-1][3];
	for(int i = 0; i < number_atoms -1; i++) {
		bonds[i][0] = coordinates[i+1][0] - coordinates[i][0];
		bonds[i][1] = coordinates[i+1][1] - coordinates[i][1];
		bonds[i][2] = coordinates[i+1][2] - coordinates[i][2];
	}
	for(int i = 0; i < number_atoms -2; i++){
		valence[i] = 180*(PI- acos(dotprod(bonds[i],bonds[i+1])/(norm(bonds[i])*norm(bonds[i+1]))))/PI;
	}
}*/


/*=====================================================================
 * VALENCEFIND.CPP
 *
 *  Created on: May 9, 2011
 *      Author: Pieter Chys
 *
 *=====================================================================
 *
 *
 *VALENCEFIND.CPP finds the valence angles between consecutive protein
 *bonds of a protein chain and returns them in array format.
 *
 *----------------------------------------------------------------------------
 *
 *INPUT: - double **pointer COORDINATES which are format [N x 4]
 *       - int NUMBER_ATOMS which is equal to the number of
 *         atoms in the protein chain...it is the size in the
 *         second first dimension COORDINATES
 *
 *----------------------------------------------------------------------------
 *
 *OUTPUT:	- explicitly VOID return
 *			- double *pointer to VALENCE which contains the values
 *			  of the valence angles
 *
 *----------------------------------------------------------------------------
 *
 *REMARK:	- uses function DOTPROD and NORM
 *
 *============================================================================
 */

void valencefind(double **coordinates,double *valence,int number_atoms){
	double bonds[number_atoms][3];

	//printf("Entering valencefind... \n");
	//fflush(stdout);

	//	First atom has no bond coordinates
	bonds[0][0] = 0.0;
	bonds[0][1] = 0.0;
	bonds[0][2] = 0.0;

	//printf("First coordinates assigned ...\n");
	//fflush(stdout);

	//	Calculate all bond coordinates
	for(int i = 0; i < number_atoms - 1; i++) {
		bonds[i+1][0] = coordinates[i+1][0] - coordinates[i][0];
		bonds[i+1][1] = coordinates[i+1][1] - coordinates[i][1];
		bonds[i+1][2] = coordinates[i+1][2] - coordinates[i][2];
		//printf("Index i %i\n",i);
		//fflush(stdout);
	}

	//printf("Bond calculation passed ...\n");
	//fflush(stdout);

	//	First two atoms have no real valence angle
	valence[0] = 0.0;
	valence[1] = 0.0;

	//printf("First valence angles assigned...\n");
	//fflush(stdout);



	//	Calculate actual valence angles
	for(int i = 2; i < number_atoms; i++){
		valence[i] = M_PI- acos(dotprod(bonds[i-1],bonds[i])/(norm(bonds[i-1])*norm(bonds[i])));				//	radians
		//valence[i] = 180*(M_PI- acos(dotprod(bonds[i-1],bonds[i])/(norm(bonds[i-1])*norm(bonds[i]))))/M_PI;	//	degrees
		//printf("Index i %i\n",i);
		//fflush(stdout);
	}
}

/*=====================================================================
 * FINDLENGTH.CPP
 *
 *  Created on: May 9, 2011
 *      Author: Pieter Chys
 *
 *=====================================================================
 *
 *
 *FINDLENGTH.CPP finds the lengths of protein bonds in a
 *protein backbone and returns them in a 1-dimensional array
 *as a pointer
 *
 *----------------------------------------------------------------------------
 *
 *INPUT: - double **pointer COORDINATES which are format [N x 4]
 *       - int NUMBER_ATOMS which is equal to the number of
 *         atoms in the protein chain...it is the size in the
 *         second first dimension COORDINATES
 *
 *----------------------------------------------------------------------------
 *
 *OUTPUT:	- explicitly VOID return
 *			- double *pointer to BONDLENGTH which contains the values
 *			  of the bond lengths
 *			- double **pointer to BONDS which contain the bond vector
 *			  coordinates and is of size [(N-1) x 4] (can be implemented
 *			  simply as [N x 4])
 *
 *----------------------------------------------------------------------------
 *
 *REMARK:	- uses function NORM normally and further secondary functions
 *            but compilation gave some problems and hence code was changed
 *          - current function only needs standar C++ libraries
 *
 *============================================================================
 */

void findlength(double **coordinates,double *bondlength,int number_atoms){
	double bonds[number_atoms][3];			//	internal variable introduce instead of **bonds

	//	First atom has no bond coordinates
	bonds[0][0] = 0.0;
	bonds[0][1] = 0.0;
	bonds[0][2] = 0.0;

	bondlength[0] = 0.0;

	// Calculate all bond coordinates
	for(int i = 0; i < number_atoms - 1; i++) {
		bonds[i+1][0] = coordinates[i+1][0] - coordinates[i][0];
		bonds[i+1][1] = coordinates[i+1][1] - coordinates[i][1];
		bonds[i+1][2] = coordinates[i+1][2] - coordinates[i][2];

		//	Part adapted for problems iterative use of +=
		bondlength[i+1] = 0.0;

		for(int j = 0; j < 3; j++){
			bondlength[i+1] += bonds[i+1][j]*bonds[i+1][j] ;			//norm(bonds[i-1]);
		}
		bondlength[i+1] = sqrt(bondlength[i+1]);
	}
}

/*===============================================================================
 *
 *	FINDLENGTH.CPP
 *
 *	Created on: April 26, 2011
 *      Author: Pieter Chys
 *
 *===============================================================================
 *
 * FINDLENGTH.CPP calculates the bond lengths in an atom system and this array
 * has a size 1 less than the original atom coordinate system.
 *
 *--------------------------------------------------------------------------------
 *
 *  INPUT: - double **coordinates are the absolute atom coordinates in (number_atoms x 3/4)
 *  		 format
 *	 	   - int number_atoms = number atoms of the atom system
 *
 *
 *-------------------------------------------------------------------------------
 *
 *  OUTPUT: - double **bonds are the bond coordinates (with respect
 *            to the preceding end atom) in (number_atoms x 3/4) format and with
 *            size 1 less than the original **coordinates
 *          - double *length contains the calculated lengths based on the
 *            **bonds matrix and has the same size as that array
 *
 *-------------------------------------------------------------------------------
 *
 *  DEPENDENCIES:	- mind you that double **coordinates is often (size x 4) in
 *  				  format
 *
 *
 *===============================================================================
 */

void    findlength(double **coordinates, double **bonds,double *length,int number_atoms){
	for(int i = 1; i < number_atoms; i++) {
		bonds[i-1][0] = coordinates[i][0] - coordinates[i-1][0];
		bonds[i-1][1] = coordinates[i][1] - coordinates[i-1][1];
		bonds[i-1][2] = coordinates[i][2] - coordinates[i-1][2];
		length[i-1] = norm(bonds[i-1]);
	}
}


/*=====================================================================
 * DIHEDRALFIND.CPP
 *
 *  Created on: May 9, 2011
 *      Author: Pieter Chys
 *
 *=====================================================================
 *
 *
 *DIHEDRALFIND.CPP finds the successive dihedral angles in a
 *protein backbone and returns them in a 1-dimensional array
 *as a pointer. The function derives from the Nacho's library
 *and a selfwritten function to determine one dihedral angle
 *with this code.
 *
 *----------------------------------------------------------------------------
 *
 *INPUT: - double **pointer COORDINATES which are format [N x 4]
 *       - int NUMBER_ATOMS which is equal to the number of
 *         atoms in the protein chain...it is the size in the
 *         second first dimension COORDINATES
 *
 *----------------------------------------------------------------------------
 *
 *OUTPUT:	- explicitly VOID return
 *			- double *pointer to DIHEDRAL which contains the values
 *			  of the dihedrals
 *
 *----------------------------------------------------------------------------
 *
 *REMARK:	- uses functions DOTPROD and NORM
 *
 *============================================================================
 */

void dihedralfind(double **coordinates,double *dihedral,int number_atoms) {
	double angle_dihedral;
	double s;
	double ct;
	double bond1[3];
	double bond2[3];
	double bond3[3];
	double normal1[3];
	double normal2[3];

	double **co;												//	new pointer to keep code in this format
	co  = coordinates;


	//	First three atoms have formally no dihedral angle
	dihedral[0] = 0;
	dihedral[1] = 0;
	dihedral[2] = 0;

	for(int j = 0; j < number_atoms - 3;j++ ){
		/* Calculate vector C,B,C */
		for (int i = 0; i < 3; i++) {
			bond1[i] = co[0][i] - co[1][i];
			bond2[i] = co[2][i] - co[1][i];
			bond3[i] = co[2][i] - co[3][i];
		}

		/* Normal to plane 1 */
		normal1[0] = bond1[1] * bond2[2] - bond1[2] * bond2[1];
		normal1[1] = bond1[2] * bond2[0] - bond1[0] * bond2[2];
		normal1[2] = bond1[0] * bond2[1] - bond1[1] * bond2[0];
		/* Normal to plane 2 */
		normal2[0] = bond2[2] * bond3[1] - bond2[1] * bond3[2];
		normal2[1] = bond2[0] * bond3[2] - bond2[2] * bond3[0];
		normal2[2] = bond2[1] * bond3[0] - bond2[0] * bond3[1];

		ct = dotprod(normal1,normal2)/(norm(normal1)*norm(normal2));

		//	In case normalisation fails
		if ( ct > 1.0 ) {
			ct = 1.0;
		} else if ( ct < ( -1.0 ) ) {
			ct = -1.0;
		}

		angle_dihedral = acos( ct );

		s = bond2[0] * ( normal1[2] * normal2[1] - normal1[1] * normal2[2] ) +
				bond2[1] * ( normal1[0] * normal2[2] - normal1[2] * normal2[0]  ) +
				bond2[2] * ( normal1[1] * normal2[0] - normal1[0] * normal2[1] );

		if ( s < 0.0 ) {
			angle_dihedral = -angle_dihedral;
		}

		angle_dihedral = ( angle_dihedral > 0.0 ) ? M_PI - angle_dihedral : -( M_PI + angle_dihedral );
		dihedral[j+3] = angle_dihedral;

		co++;
	}

	//delete co;			//	WRONG DON'T DELETE ORIGINAL OBJECT
	//return angle_dihedral ;
}

/*===============================================================================
 *
 *	DIHEDRAL.CPP
 *
 *	Created on: March 4, 2010
 *      Author: Pieter Chys
 *
 *===============================================================================
 *
 * DIHEDRAL.CPP calculates the dihedral angle in a 4 atom sequence
 * as the angle between the 2 planes formed by each 3-point collections.
 *
 *--------------------------------------------------------------------------------
 *
 *  INPUT: - double **pointer to the 4 atom coordinates which contain the x,y,z
 *  		 coordinates. The 2 planes formed by three consecutive points are used
 *           to construct normals which then are used to calculate the angles.
 *
 *-------------------------------------------------------------------------------
 *
 *  OUTPUT: - double explicitly returned which is the dihedral angle in radians
 *
 *-------------------------------------------------------------------------------
 *
 *  DEPENDENCIES:	- code take from C++ lab library and modified
 *
 *===============================================================================
 */


double dihedral(double **co) {
	double angle_dihedral;
	double s;
	double ct;
	double bond1[3];
	double bond2[3];
	double bond3[3];
	double normal1[3];
	double normal2[3];

	/* Calculate vector C,B,C */
	for (int i = 0; i < 3; i++) {
		bond1[i] = co[0][i] - co[1][i];
		bond2[i] = co[2][i] - co[1][i];
		bond3[i] = co[2][i] - co[3][i];
	}

	/* Normal to plane 1 */
	normal1[0] = bond1[1] * bond2[2] - bond1[2] * bond2[1];
	normal1[1] = bond1[2] * bond2[0] - bond1[0] * bond2[2];
	normal1[2] = bond1[0] * bond2[1] - bond1[1] * bond2[0];
	/* Normal to plane 2 */
	normal2[0] = bond2[2] * bond3[1] - bond2[1] * bond3[2];
	normal2[1] = bond2[0] * bond3[2] - bond2[2] * bond3[0];
	normal2[2] = bond2[1] * bond3[0] - bond2[0] * bond3[1];

	ct = dotprod(normal1,normal2)/(norm(normal1)*norm(normal2));

	//	In case normalisation fails
	if ( ct > 1.0 ) {
		ct = 1.0;
	} else if ( ct < ( -1.0 ) ) {
		ct = -1.0;
	}

	angle_dihedral = acos( ct );

	s = bond2[0] * ( normal1[2] * normal2[1] - normal1[1] * normal2[2] ) +
			bond2[1] * ( normal1[0] * normal2[2] - normal1[2] * normal2[0]  ) +
			bond2[2] * ( normal1[1] * normal2[0] - normal1[0] * normal2[1] );

	if ( s < 0.0 ) {
		angle_dihedral = -angle_dihedral;
	}

	angle_dihedral = ( angle_dihedral > 0.0 ) ? M_PI - angle_dihedral : -( M_PI + angle_dihedral );

	return angle_dihedral ;
}


/*=====================================================================
 * DIHEDRALFIND.CPP
 *
 *  Created on: May 9, 2011
 *      Author: Pieter Chys
 *
 *=====================================================================
 *
 *
 *DIHEDRALFIND.CPP finds the successive dihedral angles in a
 *protein backbone and returns them in a 1-dimensional array
 *as a pointer. The function derives from the Nacho's library
 *and a self-written function to determine one dihedral angle
 *with this code.
 *
 *----------------------------------------------------------------------------
 *
 *INPUT: - double **pointer COORDINATES which are format [N x 4]
 *       - int NUMBER_ATOMS which is equal to the number of
 *         atoms in the protein chain...it is the size in the
 *         second first dimension COORDINATES
 *
 *----------------------------------------------------------------------------
 *
 *OUTPUT:	- explicitly VOID return
 *			- double *pointer to DIHEDRAL which contains the values
 *			  of the dihedrals
 *
 *----------------------------------------------------------------------------
 *
 *REMARK:	- uses functions DOTPROD and NORM
 *
 *============================================================================
 */

void finddihedral(double **coordinates,double *dihedral,int number_atoms) {
	double angle_dihedral;
	double s;
	double ct;
	double bond1[3];
	double bond2[3];
	double bond3[3];
	double normal1[3];
	double normal2[3];

	double **co;												//	new pointer to keep code in this format
	co  = coordinates;


	//	First three atoms have formally no dihedral angle
	dihedral[0] = 0;
	dihedral[1] = 0;
	dihedral[2] = 0;

	for(int j = 0; j < number_atoms - 3;j++ )
	{
		/* Calculate vector C,B,C */
		for (int i = 0; i < 3; i++)
		{
			bond1[i] = co[0][i] - co[1][i];
			bond2[i] = co[2][i] - co[1][i]; // Bond that contains current dihedral angle
			bond3[i] = co[2][i] - co[3][i];

			//	  if (i==0) printf("->b1 %8.3f %8.3f b2 %8.3f %8.3f b3 %8.3f %8.3f\n",co[0][i],co[1][i],co[2][i],co[1][i],co[2][i],co[3][i]);
		}

		/* Normal to plane 1 */
		normal1[0] = bond1[1] * bond2[2] - bond1[2] * bond2[1];
		normal1[1] = bond1[2] * bond2[0] - bond1[0] * bond2[2];
		normal1[2] = bond1[0] * bond2[1] - bond1[1] * bond2[0];
		/* Normal to plane 2 */
		normal2[0] = bond2[2] * bond3[1] - bond2[1] * bond3[2];
		normal2[1] = bond2[0] * bond3[2] - bond2[2] * bond3[0];
		normal2[2] = bond2[1] * bond3[0] - bond2[0] * bond3[1];

		ct = dotprod(normal1,normal2)/(norm(normal1)*norm(normal2));

		//	In case normalisation fails
		if ( ct > 1.0 )
		{
			ct = 1.0;
		}
		else if ( ct < ( -1.0 ) )
		{
			ct = -1.0;
		}

		angle_dihedral = acos( ct );

		s = bond2[0] * ( normal1[2] * normal2[1] - normal1[1] * normal2[2] ) +
				bond2[1] * ( normal1[0] * normal2[2] - normal1[2] * normal2[0]  ) +
				bond2[2] * ( normal1[1] * normal2[0] - normal1[0] * normal2[1] );

		if ( s < 0.0 )
		{
			angle_dihedral = -angle_dihedral;
		}

		angle_dihedral = ( angle_dihedral > 0.0 ) ? M_PI - angle_dihedral : -( M_PI + angle_dihedral );

		//	Correction algorithm to get value in [-Pi,+Pi] range --- BUT IS THIS REALLY NEEDED?
		if (angle_dihedral > M_PI)
		{
			angle_dihedral = angle_dihedral - 2*M_PI;
		}
		else if (angle_dihedral < -M_PI)
		{
			angle_dihedral = angle_dihedral + 2*M_PI;
		}

		dihedral[j + 3] = angle_dihedral;

		co++;
	}

	//delete co;			//	WRONG DON'T DELETE ORIGINAL OBJECT
	//return angle_dihedral ;
}

void correctdihedral(double *angle)
{
	if (*angle  >  M_PI)
	{
		*angle -= 2*M_PI;
	}
	else if (*angle < -M_PI)
	{
		*angle += 2*M_PI;
	}
}

// ESto es una mierda
void correctdihedral(double angle)
{
	if (angle  >  M_PI)
	{
		angle -= 2*M_PI;
	}
	else if (angle < -M_PI)
	{
		angle += 2*M_PI;
	}
}


/*===============================================================================
 *
 *  ROTATEMATRIX.CPP
 *
 *	Created on: Dec 10, 2010
 *      Author: Pieter Chys
 *
 *===============================================================================
 * ROTATEMATRIX.CPP constructs a 3x3 rotation matrix to rotate
 * a chosen bond vector into the coordinates [1 0 0] with the
 * start point of the vector located at the origin. Thus the
 * bond vector is put into -as if it is- an X axis. The x and
 * y axes are now constructed in a random fashion but always
 * forming an orthogonal plane to the X axis.
 *
 *--------------------------------------------------------------------------------
 *
 *  INPUT: - double* pointer to the involved bond vector
 *  	   - double** pointer to a preallocated 3x3 array
 *  		 which will be the rotation matrix
 *
 *-------------------------------------------------------------------------------
 *
 *  OUTPUT: - function return is explicitly VOID
 *  		- double** pointer which refers to the constructed rotation
 *  		  matrix
 *
 *-------------------------------------------------------------------------------
 *
 *  DEPENDENCIES:
 *
 *===============================================================================
 */

void rotatematrix(double* bond,double** matrix)
{
	// static
	double rvecn[3],bondn[3]; 									//	static needed? -> later also?
	randnvector(rvecn);
	normalize(bond,bondn);

	// static
	double wedge1[4],wedge2[4],trivec1[4],trivec2[4],trivec1n[3],trivec2n[3];
	geom11out(bondn,rvecn,wedge1);								//	this version works faster...order natural and excellent
	// geom11out(rvecn,bondn,wedge1);							//	old version...works fine
	geomi2(wedge1,trivec1); 									//	minus needs to be programmed....
	normalize(trivec1,trivec1n);

	geom11out(bondn,trivec1n,wedge2);
	geomi2(wedge2,trivec2); 									//	minus needs to be programmed....
	normalize(trivec2,trivec2n);

	for(int i = 0; i < 3; i++) {
		matrix[0][i] = bondn[i];
		matrix[1][i] = -trivec1n[i]; 							//	minus from above are introduced here
		matrix[2][i] = -trivec2n[i]; 							//	minus from above are introduced here
	}

}


void rotatematrix2(double** co,double** matrix){

	double bond1[3];							//	First bond in the plane to be rotated to XY-plane in universal XYZ frame
	double bond2[3];							//	Second bond in the plane to be rotated to XY-plane in universal XYZ frame
	double xbond[3];							//	X vector in local frame aligned with the first bond vector
	double zbond[3];							//	Constructed orthogonal vector on top of plane to be rotated (making RH local frame)
	double ybond[3];							//	Constructed y vector in the plane of the two bonds to be rotated

	//double	check[3];

	//	Constructing the two bond vectors
	for(int i = 0; i < 3; i++){
		bond1[i]	=	co[1][i]	-	co[0][i];
		bond2[i]	=	co[2][i]	-	co[1][i];
		xbond[i]	=	bond1[i];
	}

	//	Constructing frame vectors
	crossprod(bond1,bond2,zbond);
	crossprod(zbond,bond1,ybond);

	//	Normalize frame vectors: self-mapped---original vector components lost
	normalize(xbond,xbond);
	normalize(ybond,ybond);
	normalize(zbond,zbond);


	//	Creating the transformation matrix
	//	The columms of this matrix are the X,Y,Z vectors in the local frame but written out in terms of
	//	the universal frame: See book of Craig on robotics for very simple and good explanation
	//	Introduction to Robotics- Mechanics and control (3 Edition), Craig J.
	for(int i = 0; i < 3; i++) {
		matrix[0][i] = xbond[i];
		matrix[1][i] = ybond[i]; 							//	minus from above are introduced here
		matrix[2][i] = zbond[i]; 							//	minus from above are introduced here
	}

	//	Check whether the code works
	//	MATRIX IS PERFECTLY CORRECT NOW
	/*
	matrixXvector(matrix,xbond,check);
	show1("Original second bond","X Y Z components",xbond,3);
	show1("Resulting bond","X Y Z components",check,3);
	getchar();*/

}


/*===============================================================================
 *
 *  ROTATEFRAME.CPP
 *
 *	Created on: Dec 13, 2010
 *      Author: Pieter Chys
 *
 *===============================================================================
 * ROTATEFRAME.CPP takes the coordinates of a polymer chain, takes
 * a given bond and makes this bond align with an X axis with absolute
 * unit vector [1 0 0].  A rotation matrix is thereby constructed
 * and applied to the terminal chain segment involved.
 *
 *THIS FUNCTION DOES NOT SEEM TO BE ABSOLUTELY NECESSARY
 *
 *--------------------------------------------------------------------------------
 *
 *  INPUT: - int referring to the bond which needs to be rotated
 *  	   - double** pointer to the coordinates of the whole chain
 *
 *-------------------------------------------------------------------------------
 *
 *  OUTPUT: - function return is explicitly VOID
 *  		- double** pointer which refers to the new polymer
 *  		  coordinates (RESULT)
 *
 *-------------------------------------------------------------------------------
 *
 *  DEPENDENCIES: ROTATEMATRIX, MATRIXXVECTOR
 *
 *===============================================================================
 */

void rotateframe(double** co, int bond_selected, double** result,int number_atoms)
{
	// 	Double** Pointer Created + Initialized
	double** matrix;
	pointeralloc(&matrix,3,3);

	//	Actual algorithm
	double bond[3],bondn[3],temp[3];
	diffvectorseq(co,bond_selected+1,bond);					//	Modified 27/01/2011
	normalize(bond,bondn);
	rotatematrix(bondn,matrix); 							//	trick used for double* matrix[3]
	for(int i = 0; i < number_atoms; i++)	{
		matrixXvector(matrix,co[i],temp);
		result[i]=temp;
	}

	freepointer(matrix,3);
}



/*===============================================================================
 *
 *  ROTATE.CPP
 *
 *	Created on: Dec 13, 2010
 *      Author: Pieter Chys
 *
 *===============================================================================
 * ROTATE.CPP takes the coordinates of a specified
 * segment of the polymer chain and rotates them through
 * an angle of radians as specified.
 *
 *--------------------------------------------------------------------------------
 *
 *  INPUT: - int referring to the number of atoms in the chain
 *  	   - int referring to the bond number
 *  	   - double** pointer to the coordinates of the whole chain
 *  	   - double value for angle to rotate
 *
 *-------------------------------------------------------------------------------
 *
 *  OUTPUT: - function return is explicitly VOID
 *  		- double** pointer to the coordinates of the polymer
 *  		  chain
 *
 *-------------------------------------------------------------------------------
 *
 *  DEPENDENCIES:
 *
 *===============================================================================
 */

/*void rotate(double** co, int bond_number, int number_atoms, double angle)
{

	bool debug=false; //to switch off commands which take time...e.g.
	if(debug){
		//some code...
	}

	double result[4];

	// 	Initialization of part matrix
	double** part; 											//	NULL for pointers..is pointing to nothing.
	pointeralloc(&part,number_atoms-bond_number-1,4);

	//	Actual algorithm
	double rotor[4],revrotor[4];
	double bond[3],bondn[3];

	//printf("\n%s\n\n","Running the rotate function...1");
	selectpart(co,bond_number+2,number_atoms,part); 										//	CORRECT ATOM_NUMBER FOR STARTING
	//printf("%-30s %8.3f %8.3f %8.3f\n","Part end initially",part[number_atoms-bond_number-2][0], part[number_atoms-bond_number-2][1],part[number_atoms-bond_number-2][2]);
	movechain(part,co,bond_number+1,number_atoms-bond_number-1,2);
	//printf("%-30s %8.3f %8.3f %8.3f\n","Translation coordinates",co[bond_number][0], co[bond_number][1],co[bond_number][2]);
	//printf("%-30s %8.3f %8.3f %8.3f\n","Translated coordinates",part[number_atoms-bond_number-2][0], part[number_atoms-bond_number-2][1],part[number_atoms-bond_number-2][2]);
	diffvectorseq(co,bond_number+1,bond);													//	Modified 27/01/2011
	normalize(bond,bondn);
	spinor(bondn,angle,rotor);
	revspinor(rotor,revrotor);
	for(int i=0;i<number_atoms-bond_number-1;i++) {
		geom32(part[i],rotor,result);
		geom23(revrotor,result,part[i]);
	}
	//printf("%s\n","Geometric algebra loop finished");
	//printf("%-30s %8.3f %8.3f %8.3f\n","Part end after rotation",part[number_atoms-bond_number-2][0], part[number_atoms-bond_number-2][1],part[number_atoms-bond_number-2][2]);
	movechain(part,co,bond_number+1,number_atoms-bond_number-1,1);
	//printf("%-30s %8.3f %8.3f %8.3f\n","Part end after translation",part[number_atoms-bond_number-2][0], part[number_atoms-bond_number-2][1],part[number_atoms-bond_number-2][2]);
	//printf("%-30s %8.3f %8.3f %8.3f\n","Loop end before substitution",co[number_atoms-1][0], co[number_atoms-1][1],co[number_atoms-1][2]);
	for(int i=0;i<number_atoms-bond_number-1;i++) {
		for(int j = 0; j < 3; j++) {
			co[i+bond_number+1][j] = part[i][j];
		}
	}
	//printf("%-30s %8.3f %8.3f %8.3f\n","Loop end after substitution",co[number_atoms-1][0], co[number_atoms-1][1],co[number_atoms-1][2]);
	//printf("%s\n","Rotate function finished");

	//	Free allocation
	freepointer(part,number_atoms-bond_number-1);
}*/

void rotatetree_co(double **co, int bond_number, int number_atoms, double angle)
{
	double rotor[4];
	double bond[3],part[3];
	double matrix[3][3];
	double sum;

	/*Calculate normalized bond to rotate*/
	bond[0] = co[bond_number][0] - co[bond_number-1][0];
	bond[1] = co[bond_number][1] - co[bond_number-1][1];
	bond[2] = co[bond_number][2] - co[bond_number-1][2];
	sum = sqrt(bond[0]*bond[0] + bond[1]*bond[1] + bond[2]*bond[2]);
	bond[0] /= sum;
	bond[1] /= sum;
	bond[2] /= sum;
	//	faster quaternion/matrix=spinor/matrix code ----- seems to work also---slightly faster
	quaternion(bond,angle/2,rotor);
	quatmatrix(rotor,matrix);

	for(int i = bond_number+1; i < number_atoms; i++)
	{
		part[0] = co[i][0] - co[bond_number-1][0];
		part[1] = co[i][1] - co[bond_number-1][1];
		part[2] = co[i][2] - co[bond_number-1][2];
		matrixXvector(matrix,part);
		co[i][0] = part[0] + co[bond_number-1][0];
		co[i][1] = part[1] + co[bond_number-1][1];
		co[i][2] = part[2] + co[bond_number-1][2];
	}
}

void rotatetree_cg(double **co, double **p1co,int bond_number, int number_atoms, double angle, int reverse, int mode)
{
	double rotor[4];
	double bond[3],part[3];
	double matrix[3][3];
	double sum;
	double x,y,z;

	/*Calculate normalized bond to rotate*/
	bond[0] = co[bond_number][0] - co[bond_number-1][0];
	bond[1] = co[bond_number][1] - co[bond_number-1][1];
	bond[2] = co[bond_number][2] - co[bond_number-1][2];
	sum = sqrt(bond[0]*bond[0] + bond[1]*bond[1] + bond[2]*bond[2]);
	bond[0] /= sum;
	bond[1] /= sum;
	bond[2] /= sum;
	//	faster quaternion/matrix=spinor/matrix code ----- seems to work also---slightly faster
	quaternion(bond,angle/2,rotor);
	quatmatrix(rotor,matrix);


	x=co[bond_number-1][0];
	y=co[bond_number-1][1];
	z=co[bond_number-1][2];

	if(reverse == 0)
	{
		// VALID FOR FORWARD...
		for(int i = bond_number; i < number_atoms; i++)
		{
			part[0] = co[i][0] - x;
			part[1] = co[i][1] - y;
			part[2] = co[i][2] - z;
			matrixXvector(matrix,part);
			co[i][0] = part[0] + x;
			co[i][1] = part[1] + y;
			co[i][2] = part[2] + z;

			if(mode > 0)
			{
				// (mode 3 and O or CB) or (mode 1 and O) or (mode 2 and CB)
				if( (i%3 != 0 && mode == 3) || (mode == 1 && i%3 == 2)  || (mode == 2 && i%3 == 1) )
				{
					part[0] = p1co[i][0] - x;
					part[1] = p1co[i][1] - y;
					part[2] = p1co[i][2] - z;
					matrixXvector(matrix,part);
					p1co[i][0] = part[0] + x;
					p1co[i][1] = part[1] + y;
					p1co[i][2] = part[2] + z;
				}
			}
		}
	}
	else
	{
		// VALID FOR REVERSE
		for(int i = bond_number; i < number_atoms; i++)
		{
			part[0] = co[i][0] - x;
			part[1] = co[i][1] - y;
			part[2] = co[i][2] - z;
			matrixXvector(matrix,part);
			co[i][0] = part[0] + x;
			co[i][1] = part[1] + y;
			co[i][2] = part[2] + z;

			if(mode > 0)
			{
				// (mode 3 and O or CB) or (mode 1 and O) or (mode 2 and CB)
				if( (i%3 != 2 && mode == 3) || (mode == 1 && i%3 == 0)  || (mode == 2 && i%3 == 1) )
				{
					part[0] = p1co[i][0] - x;
					part[1] = p1co[i][1] - y;
					part[2] = p1co[i][2] - z;
					matrixXvector(matrix,part);
					p1co[i][0] = part[0] + x;
					p1co[i][1] = part[1] + y;
					p1co[i][2] = part[2] + z;
				}
			}
		}
	}
}

// slow version PAblo 2020
int rotatetree_cg_check(double **co, double **p1co,int bond_number, int number_atoms,
		double angle, int reverse, int mode, vlVolume *vlm,
		float posV[3], int stepV[3], float unitV[3], int limitV[3], int *marker2)
{
	double rotor[4];
	double bond[3],part[3];
	double matrix[3][3];
	double sum;
	double x0,y0,z0;

	//vlPoint3f pos;
	float *m_pData = (float*)vlm->getVoxelVoidPtr(vlPoint3ui(0,0,0));
	//vlm->getPosition(&pos);
	int x,y,z;
	int clashes = 0; // counting clashes



	/*Calculate normalized bond to rotate*/
	bond[0] = co[bond_number][0] - co[bond_number-1][0];
	bond[1] = co[bond_number][1] - co[bond_number-1][1];
	bond[2] = co[bond_number][2] - co[bond_number-1][2];
	sum = sqrt(bond[0]*bond[0] + bond[1]*bond[1] + bond[2]*bond[2]);
	bond[0] /= sum;
	bond[1] /= sum;
	bond[2] /= sum;
	//	faster quaternion/matrix=spinor/matrix code ----- seems to work also---slightly faster
	quaternion(bond,angle/2,rotor);
	quatmatrix(rotor,matrix);


	x0=co[bond_number-1][0];
	y0=co[bond_number-1][1];
	z0=co[bond_number-1][2];

	if(reverse == 0)
	{
		// VALID FOR FORWARD...
		for(int i = bond_number; i < number_atoms; i++)
		{
			part[0] = co[i][0] - x0;
			part[1] = co[i][1] - y0;
			part[2] = co[i][2] - z0;
			matrixXvector(matrix,part);
			co[i][0] = part[0] + x0;
			co[i][1] = part[1] + y0;
			co[i][2] = part[2] + z0;

			x=(int)roundf( ((co[i][0]- posV[0]) / unitV[0]) );
			y=(int)roundf( ((co[i][1]- posV[1]) / unitV[1]) );
			z=(int)roundf( ((co[i][2]- posV[2]) / unitV[2]) );

			if( (x < limitV[0] && y < limitV[1] && z < limitV[2] && x >= 0 && y >= 0 && z >= 0) && ( *(m_pData + (z*stepV[2]+y*stepV[1]+x)) >0.00001) )
				return 1;

			if(mode > 0)
			{
				// (mode 3 and O or CB) or (mode 1 and O) or (mode 2 and CB)
				if(
				    ( i%3 == 1 && marker2[i] != 2 && mode == 3) || // CB
					( i%3 == 2  && mode == 3 )  // 0
					|| (mode == 1 && i%3 == 2)  || ((mode == 2 && i%3 == 1 && (marker2[i] != 2))) )
				{
					part[0] = p1co[i][0] - x0;
					part[1] = p1co[i][1] - y0;
					part[2] = p1co[i][2] - z0;
					matrixXvector(matrix,part);
					p1co[i][0] = part[0] + x0;
					p1co[i][1] = part[1] + y0;
					p1co[i][2] = part[2] + z0;
					x=(int)roundf( ((p1co[i][0]- posV[0]) / unitV[0]) );
					y=(int)roundf( ((p1co[i][1]- posV[1]) / unitV[1]) );
					z=(int)roundf( ((p1co[i][2]- posV[2]) / unitV[2]) );

					if( (x < limitV[0] && y < limitV[1] && z < limitV[2] && x >= 0 && y >= 0
							&& z >= 0) && (*(m_pData + (z*stepV[2]+y*stepV[1]+x))>0.00001) )
						return 1;

				}
			}
		}
	}
	else
	{
		// VALID FOR REVERSE
		for(int i = bond_number; i < number_atoms; i++)
		{
			part[0] = co[i][0] - x0;
			part[1] = co[i][1] - y0;
			part[2] = co[i][2] - z0;
			matrixXvector(matrix,part);
			co[i][0] = part[0] + x0;
			co[i][1] = part[1] + y0;
			co[i][2] = part[2] + z0;
			x=(int)roundf( ((co[i][0]- posV[0]) / unitV[0]) );
			y=(int)roundf( ((co[i][1]- posV[1]) / unitV[1]) );
			z=(int)roundf( ((co[i][2]- posV[2]) / unitV[2]) );

			if( (x < limitV[0] && y < limitV[1] && z < limitV[2] && x >= 0 && y >= 0
					&& z >= 0) && (*(m_pData + (z*stepV[2]+y*stepV[1]+x))>0.00001) )
				return 1;


			if(mode > 0)
			{
				// (mode 3 and O or CB) or (mode 1 and O) or (mode 2 and CB)
				if( (i%3 == 2 && mode == 3 && marker2[i] != 2) // CB
					 || (i%3 == 0   && mode == 3) // 0
				     || (mode == 1  && i%3 == 0)
					 || ((mode == 2 && i%3 == 1) && (marker2[i] != 2)) )
				{
					part[0] = p1co[i][0] - x0;
					part[1] = p1co[i][1] - y0;
					part[2] = p1co[i][2] - z0;
					matrixXvector(matrix,part);
					p1co[i][0] = part[0] + x0;
					p1co[i][1] = part[1] + y0;
					p1co[i][2] = part[2] + z0;
					x=(int)roundf( ((p1co[i][0]- posV[0]) / unitV[0]) );
					y=(int)roundf( ((p1co[i][1]- posV[1]) / unitV[1]) );
					z=(int)roundf( ((p1co[i][2]- posV[2]) / unitV[2]) );

					if( (x < limitV[0] && y < limitV[1] && z < limitV[2] && x >= 0 && y >= 0
							&& z >= 0) && (*(m_pData + (z*stepV[2]+y*stepV[1]+x))>0.00001) )
						return 1;

				}
			}
		}
	}
	return 0;
}





// Mon: New rotate & check function (faster)
// Coordinates may clash before "bond_number" rotable bond.
int rotatetree_check(vlVolume *vlm, double **co, int bond_number, double angle, int natoms)
{
	double rotor[4];
	double bond[3],part[3];
	double matrix[3][3];
	double sum;

	vlPoint3f pos;
	float *m_pData = (float*)vlm->getVoxelVoidPtr(vlPoint3ui(0,0,0));
	vlm->getPosition(&pos);
	vlStep m_step = vlm->stepping();
	float unitsx = vlm->units().x();
	float unitsy = vlm->units().y();
	float unitsz = vlm->units().z();
	vlDim m_dimLimit = vlm->dim();
	int limitx = m_dimLimit.x();
	int limity = m_dimLimit.y();
	int limitz = m_dimLimit.z();
	int x,y,z;
	int stepy = m_step.y();
	int stepz = m_step.z();

	// Calculate normalized bond to rotate
	bond[0] = co[bond_number][0] - co[bond_number-1][0];
	bond[1] = co[bond_number][1] - co[bond_number-1][1];
	bond[2] = co[bond_number][2] - co[bond_number-1][2];
	sum = sqrt(bond[0]*bond[0] + bond[1]*bond[1] + bond[2]*bond[2]);
	bond[0] /= sum;
	bond[1] /= sum;
	bond[2] /= sum;
	// Spinor/Quaternion stuff...
	quaternion(bond,angle/2,rotor);
	quatmatrix(rotor,matrix);

	for(int i = 3; i < natoms; i++) // Exclude first anchor
	{
		if(i>bond_number) // Moving only those atoms after the rotable bond ("bond_number")
		{
			part[0] = co[i][0] - co[bond_number-1][0];
			part[1] = co[i][1] - co[bond_number-1][1];
			part[2] = co[i][2] - co[bond_number-1][2];
			matrixXvector(matrix,part);
			co[i][0] = part[0] + co[bond_number-1][0];
			co[i][1] = part[1] + co[bond_number-1][1];
			co[i][2] = part[2] + co[bond_number-1][2];
		}

		if(i<natoms-2) // Exclude last anchor
		{
			x = (int)roundf( ((co[i][0]- pos.x()) / unitsx) );
			y = (int)roundf( ((co[i][1]- pos.y()) / unitsy) );
			z = (int)roundf( ((co[i][2]- pos.z()) / unitsz) );

			//			if( (x < limitx && y < limity && z < limitz && x >= 0 && y >= 0 && z >= 0) && (*(m_pData + (z*stepz+y*stepy+x))>0.00001))
			if( (x < limitx && y < limity && z < limitz && x >= 0 && y >= 0 && z >= 0) )
				if (*(m_pData + (z*stepz+y*stepy+x))>0.00001)
					return 1; // early exit
		}
	}

	return 0;
}

// Mon: New rotate & check function (faster)
// WARNING: Coordinates CAN'T clash before "bond_number" rotable bond.
int rotatetree_check2(vlVolume *vlm, double **co, int bond_number, double angle, int natoms)
{
	double rotor[4];
	double bond[3],part[3];
	double matrix[3][3];
	double sum;

	vlPoint3f pos;
	float *m_pData = (float*)vlm->getVoxelVoidPtr(vlPoint3ui(0,0,0));
	vlm->getPosition(&pos);
	vlStep m_step = vlm->stepping();
	float unitsx = vlm->units().x();
	float unitsy = vlm->units().y();
	float unitsz = vlm->units().z();
	vlDim m_dimLimit = vlm->dim();
	int limitx = m_dimLimit.x();
	int limity = m_dimLimit.y();
	int limitz = m_dimLimit.z();
	int x,y,z;
	int stepy = m_step.y();
	int stepz = m_step.z();

	// Calculate normalized bond to rotate
	bond[0] = co[bond_number][0] - co[bond_number-1][0];
	bond[1] = co[bond_number][1] - co[bond_number-1][1];
	bond[2] = co[bond_number][2] - co[bond_number-1][2];
	sum = sqrt(bond[0]*bond[0] + bond[1]*bond[1] + bond[2]*bond[2]);
	bond[0] /= sum;
	bond[1] /= sum;
	bond[2] /= sum;
	// Spinor/Quaternion stuff...
	quaternion(bond,angle/2,rotor);
	quatmatrix(rotor,matrix);

	for(int i = bond_number+1; i < natoms; i++) // Exclude not-moved atoms
	{
		// Moving only those atoms after the rotable bond ("bond_number")
		part[0] = co[i][0] - co[bond_number-1][0];
		part[1] = co[i][1] - co[bond_number-1][1];
		part[2] = co[i][2] - co[bond_number-1][2];
		matrixXvector(matrix,part);
		co[i][0] = part[0] + co[bond_number-1][0];
		co[i][1] = part[1] + co[bond_number-1][1];
		co[i][2] = part[2] + co[bond_number-1][2];

		if(i<natoms-2) // Exclude last anchor
		{
			x = (int)roundf( ((co[i][0]- pos.x()) / unitsx) );
			y = (int)roundf( ((co[i][1]- pos.y()) / unitsy) );
			z = (int)roundf( ((co[i][2]- pos.z()) / unitsz) );

			if( (x < limitx && y < limity && z < limitz && x >= 0 && y >= 0 && z >= 0) && (*(m_pData + (z*stepz+y*stepy+x))>0.00001))
				return 1; // early exit
		}
	}

	return 0;
}

int rotatetree_check3(vlVolume *vlm, double **co, int bond_number, double angle, int natoms)
{
	if(bond_number == natoms-2)
		return 0; // Last anchor atoms never clash... (early exit)

	double rotor[4];
	double bond[3],part[3];
	double matrix[3][3];
	double sum;
	double back[3];

	vlPoint3f pos;
	float *m_pData = (float*)vlm->getVoxelVoidPtr(vlPoint3ui(0,0,0));
	vlm->getPosition(&pos);
	vlStep m_step = vlm->stepping();
	float unitsx = vlm->units().x();
	float unitsy = vlm->units().y();
	float unitsz = vlm->units().z();
	vlDim m_dimLimit = vlm->dim();
	int limitx = m_dimLimit.x();
	int limity = m_dimLimit.y();
	int limitz = m_dimLimit.z();
	int x,y,z;
	int stepy = m_step.y();
	int stepz = m_step.z();
	int index0 = bond_number-1;

	// Calculate normalized bond to rotate
	bond[0] = co[bond_number][0] - co[index0][0];
	bond[1] = co[bond_number][1] - co[index0][1];
	bond[2] = co[bond_number][2] - co[index0][2];
	sum = sqrt(bond[0]*bond[0] + bond[1]*bond[1] + bond[2]*bond[2]);
	bond[0] /= sum;
	bond[1] /= sum;
	bond[2] /= sum;
	// Spinor/Quaternion stuff...
	quaternion(bond,angle/2,rotor);
	quatmatrix(rotor,matrix);

	// Moved atom backup (First atom after "bond_number")
	back[0] = co[bond_number+1][0];
	back[1] = co[bond_number+1][1];
	back[2] = co[bond_number+1][2];

	// Moving First atom after "bond_number"
	part[0] = co[bond_number+1][0] - co[index0][0];
	part[1] = co[bond_number+1][1] - co[index0][1];
	part[2] = co[bond_number+1][2] - co[index0][2];
	matrixXvector(matrix,part);
	co[bond_number+1][0] = part[0] + co[index0][0];
	co[bond_number+1][1] = part[1] + co[index0][1];
	co[bond_number+1][2] = part[2] + co[index0][2];

	x = (int)roundf( ((co[bond_number+1][0]- pos.x()) / unitsx) );
	y = (int)roundf( ((co[bond_number+1][1]- pos.y()) / unitsy) );
	z = (int)roundf( ((co[bond_number+1][2]- pos.z()) / unitsz) );

	if( (x < limitx && y < limity && z < limitz && x >= 0 && y >= 0 && z >= 0) && (*(m_pData + (z*stepz+y*stepy+x))>0.00001))
	{
		// Restore moved atom position
		co[bond_number+1][0] = back[0];
		co[bond_number+1][1] = back[1];
		co[bond_number+1][2] = back[2];
		return 1; // early exit
	}

	for(int i = bond_number+2; i < natoms; i++) // Moving remaining atoms
	{
		part[0] = co[i][0] - co[index0][0];
		part[1] = co[i][1] - co[index0][1];
		part[2] = co[i][2] - co[index0][2];
		matrixXvector(matrix,part);
		co[i][0] = part[0] + co[index0][0];
		co[i][1] = part[1] + co[index0][1];
		co[i][2] = part[2] + co[index0][2];
	}

	return 0;
}


void rotatetree(double** co,double **pco1,double **pco2, int bond_number, int number_atoms, double angle, double **part,
		double **part1,double **part2,double **matrix)
{
	//double result[4];
	double rotor[4];
	//double revrotor[4];
	double bond[3];
	double sum;

	/*Calculate normalized bond to rotate*/
	sum = 0.0;
	for(int i = 0; i < 3; i++)
	{
		bond[i] = co[bond_number][i] - co[bond_number-1][i];
		sum		+= bond[i]*bond[i];
	}
	sum = sqrt(sum);
	for(int i = 0; i < 3; i++)
	{
		bond[i] /= sum;
	}

	for(int i = bond_number; i < number_atoms; i++)
	{	for (int j = 0; j < 3; j++)
	{
		part[ i-bond_number][j]	=   co[i][j]  -  co[bond_number-1][j];
		part1[i-bond_number][j] = pco1[i][j]  -  co[bond_number-1][j];
		if(pco2[i][0] != 0)	/*consider that some atoms do not exist*/
		{
			part2[i-bond_number][j] = pco2[i][j] - co[bond_number-1][j];
		}
	}
	}

	//	derived spinor code ---- optimized functions geom32 and geom23opt --- seems to work also
	/*spinor(bond,angle,rotor);
	revspinor(rotor,revrotor);
	for(int i = 0; i < number_atoms-bond_number - 1;i++) {
		geom32opt(part[i],rotor,result);
		geom23opt(revrotor,result,part[i]);
	}*/

	//	faster spinor/matrix code
	/*spinormat(bond,angle,rotor);	//	this code does not work properly
	spinormatrix(rotor,matrix);
	for(int i = 0; i < number_atoms-bond_number - 1;i++) {
		matrixXvector(matrix,part[i]);
	}*/



	//	faster quaternion/matrix=spinor/matrix code ----- seems to work also---slightly faster
	quaternion(bond,angle/2,rotor);
	quatmatrix(rotor,matrix);
	for(int i = 0; i < number_atoms - bond_number + 1 ;i++) /*SHOULD -1/+1/0 BE ADDED OR NOT?*/ /*NOW SEEMS FIXED... NOT LARGE ENOUGH RANGE*/
	{
		matrixXvector(matrix,part[i],part[i]);
		matrixXvector(matrix,part1[i],part1[i]);
		matrixXvector(matrix,part2[i],part2[i]);
	}

	/*Map back the coordinates*/
	for(int i = bond_number; i < number_atoms; i++)
	{	for (int j = 0; j < 3; j++)
	{
		co[i][j] 	=  part[ i-bond_number][j]   + co[bond_number-1][j];
		pco1[i][j] 	=  part1[i-bond_number][j]   + co[bond_number-1][j];
		if(part2[i-bond_number][0] != 0)
		{
			pco2[i][j] = part2[i-bond_number][j]  + co[bond_number-1][j];
		}
	}
	}

	//printf("%50s %8.3f %8.3f %8.3f\n","Renewed end coordinate",co[number_atoms - 1][0],co[number_atoms - 1][1],co[number_atoms - 1][2]);
}


/*EXPERIMENT TO PUT ABOVE IN PARALLEL --- NOT MUCH FASTER*/
//	quaternion(bond,angle/2,rotor);
//	quatmatrix(rotor,matrix);
//
//	/*Retrieve atoms and calculate relative positions*/
//	for(int i = bond_number; i < number_atoms; i++)
//	{	for (int j = 0; j < 3; j++)
//		{
//			part[i -bond_number][j]	=   co[i][j]  -  co[bond_number-1][j];
//			part1[i-bond_number][j] = pco1[i][j]  -  co[bond_number-1][j];
//			if(pco2[i][0] != 0)	/*consider that some atoms do not exist*/
//			{
//				part2[i-bond_number][j] = pco2[i][j] - co[bond_number-1][j];
//			}
//		}
//
//		matrixXvector(matrix,part[i-bond_number],part[i-bond_number]);
//		matrixXvector(matrix,part1[i-bond_number],part1[i-bond_number]);
//		matrixXvector(matrix,part2[i-bond_number],part2[i-bond_number]);
//
//		for (int j = 0; j < 3; j++)
//		{
//			co[i][j] 	=  part[i-bond_number][j]   + co[bond_number-1][j];
//			pco1[i][j] 	= part1[i-bond_number][j]   + co[bond_number-1][j];
//			if(part2[i-bond_number][0] != 0)
//			{
//				pco2[i][j] = part2[i-bond_number][j]  + co[bond_number-1][j];
//			}
//		}
//	}


//	pure spinor code ---- code that works
/*	spinor(bond,angle,rotor);
	revspinor(rotor,revrotor);
	for(int i = 0; i < number_atoms-bond_number - 1;i++) {
		geom32(part[i],rotor,result);
		geom23(revrotor,result,part[i]);
	}*/



void rotate(double** co, int bond_number, int number_atoms, double angle, double **part,double **matrix){
	//double result[4];
	double rotor[4];
	//double revrotor[4];
	double bond[3];

	selectpart(co,bond_number + 2,number_atoms,part); 										//	CORRECT ATOM_NUMBER FOR STARTING
	movechain(part,co,bond_number + 1,number_atoms-bond_number - 1,2);
	diffvectorseq(co,bond_number + 1,bond);													//	Modified 27/01/2011
	normalise(bond);


	//	pure spinor code ---- code that works
	/*	spinor(bond,angle,rotor);
	revspinor(rotor,revrotor);
	for(int i = 0; i < number_atoms-bond_number - 1;i++) {
		geom32(part[i],rotor,result);
		geom23(revrotor,result,part[i]);
	}*/



	//	derived spinor code ---- optimized functions geom32 and geom23opt --- seems to work also
	/*spinor(bond,angle,rotor);
	revspinor(rotor,revrotor);
	for(int i = 0; i < number_atoms-bond_number - 1;i++) {
		geom32opt(part[i],rotor,result);
		geom23opt(revrotor,result,part[i]);
	}*/

	//	faster spinor/matrix code
	/*spinormat(bond,angle,rotor);	//	this code does not work properly
	spinormatrix(rotor,matrix);
	for(int i = 0; i < number_atoms-bond_number - 1;i++) {
		matrixXvector(matrix,part[i]);
	}*/


	//	faster quaternion/matrix=spinor/matrix code ----- seems to work also---slightly faster
	quaternion(bond,angle/2,rotor);
	quatmatrix(rotor,matrix);
	for(int i = 0; i < number_atoms - bond_number - 1;i++) {
		matrixXvector(matrix,part[i],part[i]);
	}


	movechain(part,co,bond_number+1,number_atoms-bond_number-1,1);

	for(int i = 0; i < number_atoms-bond_number-1;i++) {
		for(int j = 0; j < 3; j++) {
			co[i+bond_number+1][j] = part[i][j];
		}
	}

	//printf("%50s %8.3f %8.3f %8.3f\n","Renewed end coordinate",co[number_atoms - 1][0],co[number_atoms - 1][1],co[number_atoms - 1][2]);
}

void rotate_angles(double **co, double *angles,int *bondnumbers, int bondrange,int number_atoms, double *valence,double *bond,double *dihedral){
	double term1;											//	intermediate term in calculations

	double *spinor1;										//	first spinor for 2-step scheme: dihedral spinor
	spinor1 = (double*) malloc(sizeof(double)*2);
	for(int i = 0; i < 2; i++){
		spinor1[i] = 0.0;
	}

	double *spinor2;										//	second spinor for 2-step scheme: valence angle spinor
	spinor2 = (double*) malloc(sizeof(double)*2);
	for(int i = 0; i < 2; i++){
		spinor2[i] = 0.0;
	}

	double *spinortot;										//	accumulated spinor for 2-step scheme
	spinortot = (double*) malloc(sizeof(double)*4);
	for(int i = 0; i < 4; i++){
		spinortot[i] = 0.0;
	}

	double *between;
	between = (double*) malloc(sizeof(double)*4);
	for(int i = 0; i < 4; i++){
		between[i] = 0.0;
	}

	//	Determination anchor dihedrals
	double *localdihedral;
	localdihedral = (double*) malloc(sizeof(double)*10);
	for(int i = 0; i < 10; i++){
		localdihedral[i] = 0.0;
	}

	//	Determination valence angles
	double *localval;
	localval = (double*) malloc(sizeof(double)*10);
	for(int i = 0; i < 10; i++){
		localval[i] = 0.0;
	}

	printf("Allocation variables\n");
	fflush(stdout);

	//	Copy anchor coordinates
	double **copy;
	pointeralloc(&copy,10,3);


	//	Copying anchor coordinates
	//	In last three slots
	for(int i = 0; i < 8; i++){
		for(int j = 0; j < 3; j++){
			copy[i+2][j] = co[i][j];
		}
	}

	//	Creating virtual slot for first two atoms to generate correct internal parameters
	for(int i = 0; i < 3; i++){
		copy[1][i] = co[0][i];
		copy[0][i] = co[0][i];
	}
	copy[1][0] -= 1.0;
	copy[1][1] -= 1.0;
	copy[0][0] -= 2.0;
	copy[0][1] -=-1.0;

	printf("Retrieving anchor coordinates\n");
	fflush(stdout);

	//	Obtain dihedrals and valence angles
	finddihedral(copy,localdihedral,10);
	printf("Retrieving anchor dihedrals\n");
	fflush(stdout);
	valencefind(copy,localval,10);

	printf("Retrieving anchor valence angles \n");
	fflush(stdout);

	// initialisation start spinor
	spinortot[0] = 1.0;
	spinortot[1] = 0.0;
	spinortot[2] = 0.0;
	spinortot[3] = 0.0;

	//	Generation correct startspinor
	for (int i = 2; i < 10; i++) {
		//	THIS CODE SEEMS TO BE CLOSER TO REAL STRUCTURE EXPERIMENTALLY --- in this mode algorith seems to start working except but first
		spinor1[0]    = 	cos((localdihedral[i])/2);		//	THE INPUT FOR THE DIHEDRAL MUST BE IN RANGE [-PI,PI]
		spinor1[1]    = 	sin((localdihedral[i])/2);
		spinor2[0]    = 	cos((M_PI - localval[i])/2);	//	THIS IS CORRECT VALENCE ANGLE FORM: PI - real valence angle (rad)
		spinor2[1]    = 	sin((M_PI - localval[i])/2);	//		*/

		spinor3prod1(spinortot,spinor1,between);
		spinor3prod2(between,spinor2,spinortot);

		//  add to previous atom coordinate
		/*term1 = 2*bond[i];
				co[i][0] = co[i-1][0] + term1*(spinortot[1]*spinortot[1] + spinortot[0]*spinortot[0] - 0.5);
				co[i][1] = co[i-1][1] + term1*(spinortot[1]*spinortot[2] + spinortot[0]*spinortot[3]);
				co[i][2] = co[i-1][2] + term1*(spinortot[1]*spinortot[3] - spinortot[0]*spinortot[2]);*/
	}

	printf("Obtaining start spinor");
	fflush(stdout);


	//	Variables to be used
	/*double valangle;									//  valence angle of bond 2 with respect to X axis
	double *xunit;										//	unit vector along X axis
	xunit = new double[3];
	double *bond2;										//	unit vector bond 2 of protein loop
	bond2 = new double[3];
	double *cross;										//	unit vector orthogonal to xunit and bond2
	cross = new double[3];
	double *spinorpose;									//	unit vector orthogonal to xunit and bond2
	spinorpose = new double[4];*/

	//	INDEXING IS VERY CRUCIAL
	//	CORRECTION DIHEDRALS IN CORRECT [-PI,PI] RANGE
	for (int i = 0; i < bondrange; i++){
		int j;
		j = bondnumbers[i];
		//printf("Bond numbers %i\n",j);
		//fflush(stdout);	-------> Debugging is Ok ... correct numbers passed on
		//	BOND STANDS IN CORRECT INDEXING RELATION: BOND [i] versus DIHEDRAL[i]
		dihedral[j+1] += angles[i];						//	update dihedral angle
		//	Correct dihedral angle in -PI,+PI range
		if (dihedral[j+1] > M_PI) {
			dihedral[j+1] = dihedral[j+1] - 2*M_PI;
		} else if (dihedral[j+1] < -M_PI) {
			dihedral[j+1] = dihedral[j+1] + 2*M_PI;
		}
	}



	// THIRD APPROACH
	/*for(int i = 0; i < 3; i++){
		xunit[i] = 0.0;
		bond2[i] = co[2][i] - co[1][i];
	}
	xunit[0] += 1.0;


	crossprod(xunit,bond2,cross);
	normalize(bond2,bond2);
	normalize(cross,cross);
	valangle = acos(dotprod(xunit,bond2));
	spinor(cross,valangle,spinorpose);
	spinor(cross,valangle,spinortot);*/


	// initialisation start spinor
	/*spinortot[0] = 1.0;
	spinortot[1] = 0.0;
	spinortot[2] = 0.0;
	spinortot[3] = 0.0;*/

	//spinorproduct(spinorpose,spinortot,spinortot);	// WORKS BEST OF TWO
	//spinorproduct(spinortot,spinorpose,spinortot);	// AND SWITCH ORDER OF OPERATION


	//	EVERYTHING WORKS EXCEPT BUT THE INITIAL ORIENTATION -> THAT IS THE MAIN PROBLEM


	//spinorproduct(spinorpose,spinortot,spinortot);


	//	ATOM [i] IS RELATED TO INTERNAL COORDINATE INDEX  [i]

	//	STARTING ATOM [i-1] for BOND [i] is RELATED TO INTERNAL COORDINATE INDEX [i]
	for (int i = 3; i < number_atoms; i++) {
		//	THIS CODE SEEMS TO BE CLOSER TO REAL STRUCTURE EXPERIMENTALLY --- in this mode algorith seems to start working except but first
		spinor1[0]    = 	cos((dihedral[i])/2);		//	THE INPUT FOR THE DIHEDRAL MUST BE IN RANGE [-PI,PI]
		spinor1[1]    = 	sin((dihedral[i])/2);
		spinor2[0]    = 	cos((M_PI - valence[i])/2);	//	THIS IS CORRECT VALENCE ANGLE FORM: PI - real valence angle (rad)
		spinor2[1]    = 	sin((M_PI - valence[i])/2);	//		*/

		spinor3prod1(spinortot,spinor1,between);
		spinor3prod2(between,spinor2,spinortot);

		//  add to previous atom coordinate
		term1 = 2*bond[i];
		co[i][0] = co[i-1][0] + term1*(spinortot[1]*spinortot[1] + spinortot[0]*spinortot[0] - 0.5);
		co[i][1] = co[i-1][1] + term1*(spinortot[1]*spinortot[2] + spinortot[0]*spinortot[3]);
		co[i][2] = co[i-1][2] + term1*(spinortot[1]*spinortot[3] - spinortot[0]*spinortot[2]);
	}

	free(spinor1);
	free(spinor2);
	free(spinortot);
	free(between);
	free(localdihedral);
	free(localval);

	freepointer(copy,10);

	//delete[] xunit;
	//delete[] cross;
	//delete[] bond2;

}


void rotate_angles2(double **co, double *angles,int *bondnumbers, int bondrange,int number_atoms, double *valence,double *bond,double *dihedral){
	double term1;											//	intermediate term in calculations
	int j;													//	running index
	int switch_debug = 0;									//	debugging switch

	double *spinor1;										//	first spinor for 2-step scheme: dihedral spinor
	spinor1 = (double*) malloc(sizeof(double)*2);
	for(int i = 0; i < 2; i++){
		spinor1[i] = 0.0;
	}

	double *spinor2;										//	second spinor for 2-step scheme: valence angle spinor
	spinor2 = (double*) malloc(sizeof(double)*2);
	for(int i = 0; i < 2; i++){
		spinor2[i] = 0.0;
	}

	double *spinortot;										//	accumulated spinor for 2-step scheme
	spinortot = (double*) malloc(sizeof(double)*4);
	for(int i = 0; i < 4; i++){
		spinortot[i] = 0.0;
	}

	double *between;
	between = (double*) malloc(sizeof(double)*4);
	for(int i = 0; i < 4; i++){
		between[i] = 0.0;
	}

	double *localdihedral;									//	Determination anchor dihedrals
	localdihedral = (double*) malloc(sizeof(double)*10);
	for(int i = 0; i < 10; i++){
		localdihedral[i] = 0.0;
	}

	double *localval;										//	Determination valence angles
	localval = (double*) malloc(sizeof(double)*10);
	for(int i = 0; i < 10; i++){
		localval[i] = 0.0;
	}

	double **copy;											//	Copy anchor coordinates
	pointeralloc(&copy,10,3);
	message("Variables allocated",switch_debug);


	//	--------------------------------------------------------------------------------------------------------------------------------------


	if (switch_debug == 1){
		show1("Dihedrals before addition","Angle (rad)",dihedral,number_atoms);
		show1("Valences before addition","Angle (rad)",valence,number_atoms);
		getchar();
	}

	//	CORRECTION DIHEDRALS IN CORRECT [-PI,PI] RANGE
	for (int i = 0; i < bondrange; i++){
		j = bondnumbers[i];
		dihedral[j + 1] += angles[i];
		correctdihedral(dihedral+j+1);
	}

	if (switch_debug == 1){
		show1("Obtained dihedral variations","Angle (rad)",dihedral,number_atoms);
		show1("Valences after dihedral addition","Angle (rad)",valence,number_atoms);
		getchar();
	}
	message("Dihedrals corrected",switch_debug);



	//	Copying anchor coordinates and keep first two positions free
	for(int i = 0; i < 8; i++){
		for(int j = 0; j < 3; j++){
			copy[i + 2][j] = co[i][j];
		}
	}


	//	Creating virtual slot for first two atoms to generate correct internal parameters
	for(int i = 0; i < 3; i++){
		copy[1][i] = co[0][i];									//	copy first atom coordinate in slot 1
		copy[0][i] = co[0][i];									//	copy first atom coordinate in slot 2
	}
	//	BEST WORKING SET
	copy[1][0] -= 	1.0;
	copy[1][2] -= 	1.0;
	copy[0][0] -= 	2.0;
	copy[0][2] -=	1.0;

	/*copy[1][0] -= 	0.7;
	copy[1][1] -= 	1.2;
	copy[1][2] -= 	2.0;
	copy[0][0] -= 	1.1;
	copy[1][1] -= 	0.4;
	copy[0][2] -=	2.0;*/


	message("Anchor coordinates retrieved",switch_debug);

	finddihedral(copy,localdihedral,10);
	valencefind(copy,localval,10);

	//	Initialize j-value to correct for fact that dihedral variations in localdihedral/localval are not taken into account
	j	= 0;

	for (int i = 0; j != 5; i++){	//5 chosen to allow that the all loops can be treated--- all loop size allow this number to be valid
		j = bondnumbers[i];
		localdihedral[j + 3] += angles[i];
		correctdihedral(localdihedral +j + 3);
	}


	if (switch_debug == 1){
		/*show1("Dihedrals locally","Angle (rad)",localdihedral,10);
				show1("Valences locally","Angle (rad)",localval,10);
				getchar();*/
	}
	message("Anchor dihedrals and valences retrieved",switch_debug);

	if(switch_debug == 1){
		printf("Results inside rotate function\n");
		for (int i = 0; i < 10; i++){
			printf("Obtained dihedrals and valence angles %8.3f %8.3f\n",localdihedral[i],localval[i]);
		}

		printf("\n");

		for (int i = 0; i < 12; i++){
			printf("Read in dihedrals and valence angles %8.3f %8.3f\n",dihedral[i],valence[i]);
		}

		printf("\n");
	}

	// initialisation start spinor
	spinortot[0] = 1.0;
	spinortot[1] = 0.0;
	spinortot[2] = 0.0;
	spinortot[3] = 0.0;

	//	Generation correct startspinor
	for (int i = 2; i < 4; i++) {

		spinor1[0]    = 	cos((localdihedral[i])/2);		//	THE INPUT FOR THE DIHEDRAL MUST BE IN RANGE [-PI,PI]
		spinor1[1]    = 	sin((localdihedral[i])/2);
		spinor2[0]    = 	cos((M_PI - localval[i])/2);	//	THIS IS CORRECT VALENCE ANGLE FORM: PI - real valence angle (rad)
		spinor2[1]    = 	sin((M_PI - localval[i])/2);

		spinor3prod1(spinortot,spinor1,between);
		spinor3prod2(between,spinor2,spinortot);
	}

	message("Startspinor obtained",switch_debug);

	//	STARTING ATOM [i-1] for BOND [i] is RELATED TO INTERNAL COORDINATE INDEX [i]
	for (int i = 3; i < number_atoms; i++) {
		//	THIS CODE SEEMS TO BE CLOSER TO REAL STRUCTURE EXPERIMENTALLY --- in this mode algorith seems to start working except but first
		spinor1[0]    = 	cos((dihedral[i])/2);		//	THE INPUT FOR THE DIHEDRAL MUST BE IN RANGE [-PI,PI]
		spinor1[1]    = 	sin((dihedral[i])/2);
		spinor2[0]    = 	cos((M_PI - valence[i])/2);	//	THIS IS CORRECT VALENCE ANGLE FORM: PI - real valence angle (rad)
		spinor2[1]    = 	sin((M_PI - valence[i])/2);	//		*/

		spinor3prod1(spinortot,spinor1,between);
		spinor3prod2(between,spinor2,spinortot);

		//  add to previous atom coordinate
		term1 = 2*bond[i];
		co[i][0] = co[i-1][0] + term1*(spinortot[1]*spinortot[1] + spinortot[0]*spinortot[0] - 0.5);
		co[i][1] = co[i-1][1] + term1*(spinortot[1]*spinortot[2] + spinortot[0]*spinortot[3]);
		co[i][2] = co[i-1][2] + term1*(spinortot[1]*spinortot[3] - spinortot[0]*spinortot[2]);
	}

	message("Calculation finished",switch_debug);


	//	--------------------------------------------------------------------------------------------------------------------------------------



	free(spinor1);
	free(spinor2);
	free(spinortot);
	free(between);
	free(localdihedral);
	free(localval);

	freepointer(copy,10);


}


void rotate_angles3(double **co, double *angles,int *bondnumbers, int bondrange,int number_atoms, double *valence,double *bond,double *dihedral){
	double term1;											//	intermediate term
	int j;
	int switch_debug = 0;									//	switch to show debug information
	message("Start allocation...",0);

	double *spinor1;										//	first spinor for 2-step scheme: dihedral spinor
	spinor1 = (double*) malloc(sizeof(double)*2);
	for(int i = 0; i < 2; i++){
		spinor1[i] = 0.0;
	}

	double *spinor2;										//	second spinor for 2-step scheme: valence angle spinor
	spinor2 = (double*) malloc(sizeof(double)*2);
	for(int i = 0; i < 2; i++){
		spinor2[i] = 0.0;
	}

	double *spinortot;										//	accumulated spinor for 2-step scheme
	spinortot = (double*) malloc(sizeof(double)*4);
	for(int i = 0; i < 4; i++){
		spinortot[i] = 0.0;
	}

	double *between;
	between = (double*) malloc(sizeof(double)*4);
	for(int i = 0; i < 4; i++){
		between[i] = 0.0;
	}

	double *localdihedral;									//	local dihedrals
	localdihedral = (double*) malloc(sizeof(double)*(number_atoms + 2));
	for(int i = 0; i < number_atoms + 2; i++){
		localdihedral[i] = 0.0;
	}

	double *localval;										// local valence angles
	localval = (double*) malloc(sizeof(double)*(number_atoms + 2));
	for(int i = 0; i < number_atoms + 2; i++){
		localval[i] = 0.0;
	}

	double *localbond;										// local bond lengths
	localbond = (double*) malloc(sizeof(double)*(number_atoms + 2));
	for(int i = 0; i < number_atoms + 2; i++){
		localbond[i] = 0.0;
	}

	message("Allocation passed...",switch_debug);

	//	Copy coordinates + add two extra
	double **coor;
	pointeralloc(&coor,number_atoms + 2,3);
	message("Allocation double pointer passed...",switch_debug);

	//	--------------------------------------------------------------------------------------------------------------------------------------


	//	Copying anchor coordinates and keep first two positions free
	for(int i = 0; i < number_atoms; i++){
		//printf("index %i\n",i);
		for(int j = 0; j < 3; j++){
			coor[i + 2][j] = co[i][j];
		}
	}
	message("Copying coordinates passed...",switch_debug);



	//	Creating virtual slot for first two atoms to generate correct internal parameters
	for(int i = 0; i < 3; i++){
		coor[1][i] = co[0][i];									//	copy first atom coordinate in slot 1
		coor[0][i] = co[0][i];									//	copy first atom coordinate in slot 2
	}

	//	XZ SET
	coor[1][0] -= 1.0;
	//coor[1][2] -= 1.0;
	//coor[0][0] -= 2.0;
	//coor[0][2] -= 1.0;

	//	THIS SET WORKS ALSO RELATIVELY FINE
	/*coor[1][2] -= 1.0;
	coor[1][1] -= 1.0;
	coor[0][2] -= 2.0;
	coor[0][1] -= 1.0;*/

	if (switch_debug == 1) {
		show3("Check expanded coordinate","X Y Z coordinates",coor,number_atoms + 2);
	}



	//	Obtain dihedrals and valence angles
	/*finddihedral(coor,localdihedral,number_atoms + 2);
	valencefind(coor,localval,number_atoms + 2);
	findlength(coor,localbond,number_atoms + 2);*/

	finddihedral(coor,localdihedral + 1,number_atoms + 1);
	valencefind(coor,localval + 1,number_atoms + 1);
	findlength(coor,localbond + 1,number_atoms + 1);

	if (switch_debug == 1) {
		//show3("Check loop coordinates","Coordinates",coor,number_atoms + 2);
		show3("Check obtained internal coordinate","Bond - Valence angle - Dihedral angle:",localbond,localval,localdihedral,number_atoms + 2);
	}

	//	CORRECTION DIHEDRALS IN CORRECT [-PI,PI] RANGE
	for (int i = 0; i < bondrange; i++){
		j = bondnumbers[i];
		//localdihedral[j+3] += angles[i];						//	update dihedral angle

		//correctdihedral(localdihedral+j+3);

		localdihedral[j+2] += angles[i];						//	update dihedral angle

		correctdihedral(localdihedral+j+2);

	}

	if (switch_debug == 1){
		//show3("Check loop coordinates","Coordinates",coor,number_atoms + 2);
		show3("Check obtained internal coordinate","Bond - Valence angle - Dihedral angle:",localbond,localval,localdihedral,number_atoms + 2);
		getchar();
	}


	// initialisation start spinor
	spinortot[0] = 1.0;
	spinortot[1] = 0.0;
	spinortot[2] = 0.0;
	spinortot[3] = 0.0;

	message("X-axis Startspinor created",switch_debug);

	//	--------------------------------------------------------------------------------------------------------------------------------------

	//	STARTING ATOM [i-1] for BOND [i] is RELATED TO INTERNAL COORDINATE INDEX [i]
	//for (int i = 3; i < number_atoms + 2; i++) {
	for (int i = 2; i < number_atoms + 1; i++) {
		// THIS IS CRUCIAL PART --- TOGETHER WITH INDEX
		spinor1[0]    = 	cos((localdihedral[i])/2);
		spinor1[1]    = 	sin((localdihedral[i])/2);
		spinor2[0]    = 	cos((M_PI - localval[i])/2);
		spinor2[1]    = 	sin((M_PI - localval[i])/2);

		//	CORRECT
		spinor3prod1(spinortot,spinor1,between);
		spinor3prod2(between,spinor2,spinortot);

		//	CORRECT
		if(i > 4){
			term1 = 2*localbond[i];
			coor[i][0] = coor[i-1][0] + term1*(spinortot[1]*spinortot[1] + spinortot[0]*spinortot[0] - 0.5);
			coor[i][1] = coor[i-1][1] + term1*(spinortot[1]*spinortot[2] + spinortot[0]*spinortot[3]);
			coor[i][2] = coor[i-1][2] + term1*(spinortot[1]*spinortot[3] - spinortot[0]*spinortot[2]);
		}
	}


	//	Copying anchor coordinates and keep first two positions free
	for(int i = 0; i < number_atoms; i++){
		for(int j = 0; j < 3; j++){
			co[i][j] = coor[i + 2 ][j];
		}
	}
	message("Copying coordinates back passed...",switch_debug);
	if (switch_debug == 1){
		show3("Check loop coordinates","Coordinates",co,number_atoms);
		//show3("Check obtained internal coordinate","Bond - Valence angle - Dihedral angle:",localbond,localval,localdihedral,number_atoms + 2);
		getchar();
	}

	//	--------------------------------------------------------------------------------------------------------------------------------------

	free(spinor1);
	free(spinor2);
	free(spinortot);
	free(between);

	free(localdihedral);
	free(localval);
	free(localbond);

	freepointer(coor,number_atoms+2);

}



void rotate_angles4(double **co, double *angles,int *bondnumbers, int bondrange,int number_atoms, double *valence,double *bond,double *dihedral){
	double term1;											//	intermediate term
	int j;
	int switch_debug = 0;									//	switch to show debug information
	message("Start allocation...",0);

	double *spinor1;										//	first spinor for 2-step scheme: dihedral spinor
	spinor1 = (double*) malloc(sizeof(double)*2);
	for(int i = 0; i < 2; i++){
		spinor1[i] = 0.0;
	}

	double *spinor2;										//	second spinor for 2-step scheme: valence angle spinor
	spinor2 = (double*) malloc(sizeof(double)*2);
	for(int i = 0; i < 2; i++){
		spinor2[i] = 0.0;
	}

	double *spinortot;										//	accumulated spinor for 2-step scheme
	spinortot = (double*) malloc(sizeof(double)*4);
	for(int i = 0; i < 4; i++){
		spinortot[i] = 0.0;
	}

	double *between;
	between = (double*) malloc(sizeof(double)*4);
	for(int i = 0; i < 4; i++){
		between[i] = 0.0;
	}

	double *localdihedral;									//	local dihedrals
	localdihedral = (double*) malloc(sizeof(double)*(number_atoms + 1));
	for(int i = 0; i < number_atoms + 2; i++){
		localdihedral[i] = 0.0;
	}

	double *localval;										// local valence angles
	localval = (double*) malloc(sizeof(double)*(number_atoms + 1));
	for(int i = 0; i < number_atoms + 1; i++){
		localval[i] = 0.0;
	}

	double *localbond;										// local bond lengths
	localbond = (double*) malloc(sizeof(double)*(number_atoms + 1));
	for(int i = 0; i < number_atoms + 1; i++){
		localbond[i] = 0.0;
	}

	message("Allocation passed...",switch_debug);

	//	Copy coordinates + add two extra
	double **coor;
	pointeralloc(&coor,number_atoms + 1,3);
	message("Allocation double pointer passed...",switch_debug);

	//	--------------------------------------------------------------------------------------------------------------------------------------


	//	Copying anchor coordinates and keep first two positions free
	for(int i = 0; i < number_atoms; i++){
		//printf("index %i\n",i);
		for(int j = 0; j < 3; j++){
			coor[i + 1][j] = co[i][j];
		}
	}
	message("Copying coordinates passed...",switch_debug);



	//	Creating virtual slot for first two atoms to generate correct internal parameters
	for(int i = 0; i < 3; i++){
		coor[1][i] = co[0][i];									//	copy first atom coordinate in slot 1
		coor[0][i] = co[0][i];									//	copy first atom coordinate in slot 2
	}
	coor[1][0] -= 1.0;


	if (switch_debug == 1) {
		show3("Check expanded coordinate","X Y Z coordinates",coor,number_atoms + 1);
	}


	//	Obtain internal coordinates
	finddihedral(coor,localdihedral + 1,number_atoms + 1);
	valencefind(coor,localval + 1,number_atoms + 1);
	findlength(coor,localbond + 1,number_atoms + 1);

	if (switch_debug == 1) {
		//show3("Check loop coordinates","Coordinates",coor,number_atoms + 2);
		show3("Check obtained internal coordinate","Bond - Valence angle - Dihedral angle:",localbond,localval,localdihedral,number_atoms + 1);
	}

	//	CORRECTION DIHEDRALS IN CORRECT [-PI,PI] RANGE
	for (int i = 0; i < bondrange; i++){
		j = bondnumbers[i];
		localdihedral[j+2] += angles[i];
		correctdihedral(localdihedral+j+2);

	}

	if (switch_debug == 1){
		//show3("Check loop coordinates","Coordinates",coor,number_atoms + 2);
		show3("Check obtained internal coordinate","Bond - Valence angle - Dihedral angle:",localbond,localval,localdihedral,number_atoms + 1);
		getchar();
	}


	// initialisation start spinor
	spinortot[0] = 1.0;
	spinortot[1] = 0.0;
	spinortot[2] = 0.0;
	spinortot[3] = 0.0;

	message("X-axis Startspinor created",switch_debug);

	//	--------------------------------------------------------------------------------------------------------------------------------------

	//	STARTING ATOM [i-1] for BOND [i] is RELATED TO INTERNAL COORDINATE INDEX [i]
	//for (int i = 3; i < number_atoms + 2; i++) {
	for (int i = 2; i < number_atoms + 1; i++) {
		// THIS IS CRUCIAL PART --- TOGETHER WITH INDEX
		spinor1[0]    = 	cos((localdihedral[i])/2);
		spinor1[1]    = 	sin((localdihedral[i])/2);
		spinor2[0]    = 	cos((M_PI - localval[i])/2);
		spinor2[1]    = 	sin((M_PI - localval[i])/2);

		//	CORRECT
		spinor3prod1(spinortot,spinor1,between);
		spinor3prod2(between,spinor2,spinortot);

		//	CORRECT
		if(i > 4){
			term1 = 2*localbond[i];
			coor[i][0] = coor[i-1][0] + term1*(spinortot[1]*spinortot[1] + spinortot[0]*spinortot[0] - 0.5);
			coor[i][1] = coor[i-1][1] + term1*(spinortot[1]*spinortot[2] + spinortot[0]*spinortot[3]);
			coor[i][2] = coor[i-1][2] + term1*(spinortot[1]*spinortot[3] - spinortot[0]*spinortot[2]);
		}
	}


	//	Copying anchor coordinates and keep first two positions free
	for(int i = 0; i < number_atoms; i++){
		for(int j = 0; j < 3; j++){
			co[i][j] = coor[i + 1][j];
		}
	}
	message("Copying coordinates back passed...",switch_debug);
	if (switch_debug == 1){
		show3("Check loop coordinates","Coordinates",co,number_atoms);
		//show3("Check obtained internal coordinate","Bond - Valence angle - Dihedral angle:",localbond,localval,localdihedral,number_atoms + 2);
		getchar();
	}

	//	--------------------------------------------------------------------------------------------------------------------------------------

	free(spinor1);
	free(spinor2);
	free(spinortot);
	free(between);

	free(localdihedral);
	free(localval);
	free(localbond);

	freepointer(coor,number_atoms+1);

}




void rotate_anglesm(double **co, double *angles,int *bondnumbers, int bondrange,int number_atoms, double *valence,double *bond,double *dihedral){
	double term1;											//	intermediate term
	double bond2[3];
	int j;
	int switch_debug = 0;									//	switch to show debug information

	message("Start allocation...",0);

	// MATRIX BACKWARD APPROACH
	double **matrix1;
	double **matrix2;
	double **matrixtot;
	double **matrixrot;
	double **matrixrev;
	pointeralloc(&matrix1,2,2);
	pointeralloc(&matrix2,2,2);
	pointeralloc(&matrixtot,3,3);
	pointeralloc(&matrixrot,3,3);
	pointeralloc(&matrixrev,3,3);

	message("Allocation passed...",switch_debug);

	//	CORRECTION DIHEDRALS IN CORRECT [-PI,PI] RANGE
	for (int i = 0; i < bondrange; i++){						//	POSSIBLE BUGS FOR BIGGER LOOPS WITH INDEXING BONDRANGE
		j = bondnumbers[i];
		dihedral[j+1] += angles[i];
		correctdihedral(dihedral+j+1);
		show1("Bond number chosen",j);
	}

	show1("Bondrange value",bondrange);
	getchar();



	diffvectorseq(co,3,bond2);									//	selecting bond 2
	normalise(bond2);
	rotatematrix(bond2,matrixrot);
	transpose(matrixrot,matrixrev);


	//	Specify first atom
	/*comatrix[0][0] = 0.0;
	comatrix[0][1] = 0.0;
	comatrix[0][2] = 0.0;*/

	// initialisation start matrix
	matrixtot[0][0]  = 1.0;
	matrixtot[0][1]  = 0.0;
	matrixtot[0][2]  = 0.0;
	matrixtot[1][0]  = 0.0;
	matrixtot[1][1]  = 1.0;
	matrixtot[1][2]  = 0.0;
	matrixtot[2][0]  = 0.0;
	matrixtot[2][1]  = 0.0;
	matrixtot[2][2]  = 1.0;


	matrixXmatrix(matrixrev,matrixtot,matrixtot);

	/*for (int i = 2; i < 3; i++) {
			// allocation local matrices: double multiplication
			matrix1[0][0] = cos(dihedral[i]);				//	Index altered to i
			matrix1[1][0] = sin(dihedral[i]);
			matrix1[0][1] = matrix1[1][0];
			matrix1[1][1] = matrix1[0][0];

			matrix2[0][0] = -cos(valence[i]);
			matrix2[1][0] = sin(valence[i]);
			matrix2[0][1] = matrix2[1][0];
			matrix2[1][1] = matrix2[0][0];

			// update accumulated matrix
			matrixprod1(matrixtot,matrix1);
			matrixprod2(matrixtot,matrix2);

	}*/


	for (int i = 3; i < number_atoms; i++) {
		// allocation local matrices: double multiplication
		matrix1[0][0] = cos(dihedral[i]);				//	Index altered to i
		matrix1[1][0] = sin(dihedral[i]);
		matrix1[0][1] = matrix1[1][0];
		matrix1[1][1] = matrix1[0][0];

		matrix2[0][0] = -cos(valence[i]);
		matrix2[1][0] = sin(valence[i]);
		matrix2[0][1] = matrix2[1][0];
		matrix2[1][1] = matrix2[0][0];

		// update accumulated matrix
		matrixprod1(matrixtot,matrix1);
		matrixprod2(matrixtot,matrix2);

		//	rotate + add vector
		co[i][0] = co[i-1][0] +  matrixtot[0][0]*bond[i];
		co[i][1] = co[i-1][1] +  matrixtot[1][0]*bond[i];
		co[i][2] = co[i-1][2] +  matrixtot[2][0]*bond[i];
	}

	free(matrix1);
	free(matrix2);
	free(matrixtot);
	free(matrixrot);
	free(matrixrev);
}


void rotate_bw(double **co,double *bonds,double *valence, double* dihedral, int size){

	// Local spinors
	double spinor1[2];
	double spinor2[2];

	double spinortot[4];
	double between[4];

	double term1;

	//	Specify first atom
	//	co[0][0] = 0.0;						//	not needed	always first coordinate
	//	co[0][1] = 0.0;						//	not needed
	//	co[0][2] = 0.0;						//	not needed

	// initialisation start spinor
	spinortot[0] = 1.0;
	spinortot[1] = 0.0;
	spinortot[2] = 0.0;
	spinortot[3] = 0.0;

	for (int i = 3; i < size + 1; i++) {			//	adapted
		//	Implementation: more operations but less writing to memory----hence little faster
		spinor1[0]    = 	cos(dihedral[i-1]/2);
		spinor1[1]    = 	sin(dihedral[i-1]/2);
		spinor2[0]    = 	cos(valence[i-1]/2);	//	later pass on as direct input
		spinor2[1]    = 	sin(valence[i-1]/2);	//	later pass on as direct input


		spinor3prod1(spinortot,spinor1,between);
		spinor3prod2(between,spinor2,spinortot);

		//  add to previous atom coordinate
		term1 = 2*bonds[i-1];
		co[i-1][0] = co[i-2][0] + term1*(spinortot[1]*spinortot[1] + spinortot[0]*spinortot[0] - 0.5);	//	index atom adapted RH side
		co[i-1][1] = co[i-2][1] + term1*(spinortot[1]*spinortot[2] + spinortot[0]*spinortot[3]);
		co[i-1][2] = co[i-2][2] + term1*(spinortot[1]*spinortot[3] - spinortot[0]*spinortot[2]);
	}

}


void matrixrotate(	double **co,					/*loop coordinates*/
		const double *angles,		/*angle variations from optimization*/
		double **posematrix,		/*pose matrix for anchor*/
		const int *bondnumbers,		/*physical bond numbers*/
		const int bondrange,		/*# of bonds rotating*/
		const int number_atoms,		/*# atoms in loop*/
		const double *valence,		/*valence angles loop*/
		const double *bond,			/*bond lengths loop*/
		double *dihedral,			/*dihedral angles loop*/
		const double *cosarray,		/*preallocated -cosines*/
		const double *sinarray,		/*preallocated sines*/
		double **matrix1,			/*auxiliary matrix 1-preallocated*/
		double **matrix2,			/*auxiliary matrix 2-preallocated*/
		double **matrixtot			/*total matrix-preallocated*/
)
{
	int j;										/*index dihedral correction step*/
	int switch_debug = 0;						/*switch for debug information*/

	message("Start allocation...",0);
	message("Allocation passed...",switch_debug);


	/*show1("dihedral of bond 2",dihedral[3]); getchar();*/

	for (int i = 0; i < bondrange; i++){	/*Correction dihedral angles in [-PI,PI] range*/
		j = bondnumbers[i];
		/*show1("dihedral before",dihedral[j+1]);
		show1("Bond number chosen",j);*/
		dihedral[j + 1] += angles[i];
		/*show1("dihedral after",dihedral[j+1]);*/
		correctdihedral(dihedral + j + 1);
	}


	/*show1("Bond range number",bondrange);	getchar();*/

	matrixtot[0][0]  = 1.0;		/*Initialization start matrix*/
	matrixtot[0][1]  = 0.0;
	matrixtot[0][2]  = 0.0;
	matrixtot[1][0]  = 0.0;
	matrixtot[1][1]  = 1.0;
	matrixtot[1][2]  = 0.0;
	matrixtot[2][0]  = 0.0;
	matrixtot[2][1]  = 0.0;
	matrixtot[2][2]  = 1.0;

	matrixXmatrix(posematrix,matrixtot,matrixtot);	/*Introduce start orientation*/

	for (int i = 3; i < number_atoms; i++) {

		matrix1[0][0] = cos(dihedral[i]);	/*Assign local matrices*/
		matrix1[1][0] = sin(dihedral[i]);
		matrix1[0][1] = matrix1[1][0];
		matrix1[1][1] = matrix1[0][0];

		matrix2[0][0] = cosarray[i];
		matrix2[1][0] = sinarray[i]; //*/
		/*matrix2[0][0] = -cos(valence[i]);
		matrix2[1][0] =  sin(valence[i]); // */
		matrix2[0][1] = matrix2[1][0];
		matrix2[1][1] = matrix2[0][0];

		matrixprod1(matrixtot,matrix1);		/*Update accumulated matrix*/
		matrixprod2(matrixtot,matrix2);

		co[i][0] = co[i-1][0] +  matrixtot[0][0]*bond[i];	/*Rotate and add vector */
		co[i][1] = co[i-1][1] +  matrixtot[1][0]*bond[i];
		co[i][2] = co[i-1][2] +  matrixtot[2][0]*bond[i];
	}
}


void isafrotate(double **co,				/*loop coordinates*/
		const double *angles,		/*angle variations from optimization*/
		double *posespinor,			/*pose spinor for anchor*/
		int chosenbond,
		const int *bondnumbers,		/*physical bond numbers*/
		const int bondrange,		/*# of bonds rotating*/
		const int number_atoms,		/*# atoms in loop*/
		const double *valence,		/*valence angles loop*/
		const double *bond,			/*bond lengths loop*/
		double *dihedral,			/*dihedral angles loop*/
		double *cosarrayd,		/*preallocated -cosines*/
		double *sinarrayd,		/*preallocated sines*/
		const double *cosarray,		/*preallocated -cosines*/
		const double *sinarray,		/*preallocated sines*/
		double *spinor1,			/*auxiliary spinor 1-preallocated*/
		double *spinor2,			/*auxiliary spinor 2-preallocated*/
		double *spinortot			/*total spinor-preallocated*/
)
{

	double term1;	/*intermediate computational term*/
	int j;										/*index dihedral correction step*/
	int switch_debug = 0;						/*switch for debug information*/
	double *between;
	between = (double*) malloc(sizeof(double)*4);
	for(int i = 0; i < 4; i++)
	{
		between[i] = 0.0;
	}

	message("Start allocation...",0);
	message("Allocation passed...",switch_debug);


	//	for (int i = 0; i < bondrange; i++){	/*Correction dihedral angles in [-PI,PI] range*/
	//		j = bondnumbers[i];
	//		dihedral[j + 1] += angles[i];
	//		correctdihedral(dihedral + j + 1);
	//	}

	correctdihedral(dihedral + chosenbond + 1);
	cosarrayd[chosenbond+1]   = cos ( dihedral[chosenbond +1]/2 );
	sinarrayd[chosenbond+1]   = sin ( dihedral[chosenbond +1]/2 );

	spinortot[0] = 1.0;
	spinortot[1] = 0.0;
	spinortot[2] = 0.0;
	spinortot[3] = 0.0;

	quatproduct(posespinor,spinortot,spinortot);

	for (int i = 1; i < number_atoms; i++)
	{

		spinor1[0]    = 	cosarrayd[i];
		spinor1[1]    = 	sinarrayd[i];
		//spinor1[0]    = 	cos(dihedral[i]/2);
		//spinor1[1]    = 	sin(dihedral[i]/2);
		//spinor2[0]    = 	cos((M_PI-valence[i])/2);
		//spinor2[1]    = 	sin((M_PI-valence[i])/2);
		spinor2[0]    = 	cosarray[i];	/*seems to work much better*/
		spinor2[1]    = 	sinarray[i];

		spinor3prod1(spinortot,spinor1,between);	/*update spinortot*/
		spinor3prod2(between,spinor2,spinortot);

		term1 = 2*bond[i];	/*Rotate and add vector */
		co[i][0] = co[i-1][0] + term1*(spinortot[1]*spinortot[1] + spinortot[0]*spinortot[0] - 0.5);
		co[i][1] = co[i-1][1] + term1*(spinortot[1]*spinortot[2] + spinortot[0]*spinortot[3]);
		co[i][2] = co[i-1][2] + term1*(spinortot[1]*spinortot[3] - spinortot[0]*spinortot[2]);
	}



	free(between);
}



void ga_rotate_anglesm2(double **co,				/*loop coordinates*/
		const double *angles,		/*angle variations from optimization*/
		double *posespinor,			/*pose spinor for anchor*/
		const int *bondnumbers,		/*physical bond numbers*/
		const int bondrange,		/*# of bonds rotating*/
		const int number_atoms,		/*# atoms in loop*/
		const double *valence,		/*valence angles loop*/
		const double *bond,			/*bond lengths loop*/
		double *dihedral,			/*dihedral angles loop*/
		const double *cosarray,		/*preallocated -cosines*/
		const double *sinarray,		/*preallocated sines*/
		double *spinor1,			/*auxiliary spinor 1-preallocated*/
		double *spinor2,			/*auxiliary spinor 2-preallocated*/
		double *spinortot			/*total spinor-preallocated*/
)
{

	double term1;	/*intermediate computational term*/
	int j;										/*index dihedral correction step*/
	int switch_debug = 0;						/*switch for debug information*/
	double *between;
	between = (double*) malloc(sizeof(double)*4);
	for(int i = 0; i < 4; i++)
	{
		between[i] = 0.0;
	}

	message("Start allocation...",0);
	message("Allocation passed...",switch_debug);


	for (int i = 0; i < bondrange; i++){	/*Correction dihedral angles in [-PI,PI] range*/
		j = bondnumbers[i];
		dihedral[j + 1] += angles[i];
		correctdihedral(dihedral + j + 1);
	}

	spinortot[0] = 1.0;
	spinortot[1] = 0.0;
	spinortot[2] = 0.0;
	spinortot[3] = 0.0;

	quatproduct(posespinor,spinortot,spinortot);

	for (int i = 1; i < number_atoms; i++)
	{
		spinor1[0]    = 	cos(dihedral[i]/2);
		spinor1[1]    = 	sin(dihedral[i]/2);
		//spinor2[0]    = 	cos((M_PI-valence[i])/2);
		//spinor2[1]    = 	sin((M_PI-valence[i])/2);
		spinor2[0]    = 	cosarray[i];	/*seems to work much better*/
		spinor2[1]    = 	sinarray[i];

		spinor3prod1(spinortot,spinor1,between);	/*update spinortot*/
		spinor3prod2(between,spinor2,spinortot);

		term1 = 2*bond[i];	/*Rotate and add vector */
		co[i][0] = co[i-1][0] + term1*(spinortot[1]*spinortot[1] + spinortot[0]*spinortot[0] - 0.5);
		co[i][1] = co[i-1][1] + term1*(spinortot[1]*spinortot[2] + spinortot[0]*spinortot[3]);
		co[i][2] = co[i-1][2] + term1*(spinortot[1]*spinortot[3] - spinortot[0]*spinortot[2]);
	}



	free(between);
}


void isafrotatetreefs(double **co,				/*loop coordinates*/
		double **p1co,				/*coordinates 1st parallel chain*/
		double **p2co,				/*coordinates 2nd parallel chain*/
		const double *angles,		/*angle variations from optimization*/
		int chosenbond,
		double *posespinor,			/*pose spinor for anchor*/
		const int *bondnumbers,		/*physical bond numbers*/
		const int bondrange,		/*# of bonds rotating*/
		const int number_atoms,		/*# atoms in loop*/
		const double *valence,		/*valence angles loop*/
		const double *bond,			/*bond lengths loop*/
		double *cosarrayd,		/*preallocated -cosines*/
		double *sinarrayd,		/*preallocated sines*/
		const double *cosarray,		/*preallocated -cosines*/
		const double *sinarray,		/*preallocated sines*/
		double *dihedral,			/*dihedral angles loop*/
		double *p1dihedral,			/*dihedral angles loop*/
		double *p2dihedral,			/*dihedral angles loop*/
		int *marker,
		double *spinor1,			/*auxiliary spinor 1-preallocated*/
		double *spinor2,			/*auxiliary spinor 2-preallocated*/
		double *spinortot,			/*total spinor-preallocated*/
		double *s1tot,
		double *s2tot,
		double *between,
		double *full,
		double *s3tot,
		double **s1co
)
{

	//double bondNH = 1.02;
	//double bondCH = 1.10;
	//double bondCO = 1.24;
	double term1;		/*intermediate computational term*/
	double hoek;	/*random dihedral angle for hydrogens on first sidechain carbon*/


	correctdihedral(dihedral + chosenbond + 1);
	cosarrayd[chosenbond+1]   = cos ( dihedral[chosenbond +1]/2 );	/*normally active*/
	sinarrayd[chosenbond+1]   = sin ( dihedral[chosenbond +1]/2 );

	/*Add hydrogen at first nitrogen*/
	//spinortot[0] = 1.0;
	//spinortot[1] = 0.0;
	//spinortot[2] = 0.0;
	//spinortot[3] = 0.0;

	//quatproduct(posespinor,spinortot,spinortot);

	spinor1[0]    = 	cos((dihedral[1] + M_PI)/2);
	spinor1[1]    = 	sin((dihedral[1]+ M_PI)/2);
	spinor2[0]    = 	cos((M_PI-valence[1])/2);
	spinor2[1]    = 	sin((M_PI-valence[1])/2);

	spinor3prod1(full,spinor1,between);
	spinor3prod2(between,spinor2,spinortot);

	//  add to previous atom coordinate
	term1 	 = 2*bondCH;
	p1co[0][0] = co[0][0] + term1*(spinortot[1]*spinortot[1] + spinortot[0]*spinortot[0] - 0.5);
	p1co[0][1] = co[0][1] + term1*(spinortot[1]*spinortot[2] + spinortot[0]*spinortot[3]);
	p1co[0][2] = co[0][2] + term1*(spinortot[1]*spinortot[3] - spinortot[0]*spinortot[2]);

	/*Actual start*/
	spinortot[0] = full[0];
	spinortot[1] = full[1];
	spinortot[2] = full[2];
	spinortot[3] = full[3];

	//quatproduct(posespinor,spinortot,spinortot);

	/*Extension for hydrogens and oxygens of loop*/
	for (int i = 1; i < number_atoms; i++)
	{
		spinor1[0]    = 	cosarrayd[i];
		spinor1[1]    = 	sinarrayd[i];
		spinor2[0]    = 	cosarray[i];	/*seems to work much better*/
		spinor2[1]    = 	sinarray[i];

		//spinor1[0]    = 	cos(dihedral[i]/2);
		//spinor1[1]    = 	sin(dihedral[i]/2);
		//spinor2[0]    = 	cos((M_PI-valence[i])/2);
		//spinor2[1]    = 	sin((M_PI-valence[i])/2);

		spinor3prod1(spinortot,spinor1,between);
		spinor3prod2(between,spinor2,spinortot);

		//  add to previous atom coordinate
		term1 	 = 2*bond[i];
		co[i][0] = co[i-1][0] + term1*(spinortot[1]*spinortot[1] + spinortot[0]*spinortot[0] - 0.5);
		co[i][1] = co[i-1][1] + term1*(spinortot[1]*spinortot[2] + spinortot[0]*spinortot[3]);
		co[i][2] = co[i-1][2] + term1*(spinortot[1]*spinortot[3] - spinortot[0]*spinortot[2]);

		/*Entering tree structure*/
		if (marker[i] == 2)		/*sp3-hybridization or like*/ /*From first CO onward*/
		{
			/*first side chain*/
			spinor1[0]    	= 	cos((dihedral[i+1] + 2*M_PI/3)/2);	/*spinors are fixed -- can be precalculated*/ /*indices dihedral mapping Ok*/
			spinor1[1]    	= 	sin((dihedral[i+1] + 2*M_PI/3)/2);
			spinor2[0]    = 	cosarray[i+1];	/*seems to work much better*/
			spinor2[1]    = 	sinarray[i+1];

			//spinor2[0]   	= 	cos((M_PI - valence[i+1])/2);		/*fixed valences angles ->no variable ones*/
			//spinor2[1]    = 	sin((M_PI - valence[i+1])/2);		/*indices not yet correctly fixed in space*/

			spinor3prod1(spinortot,spinor1,between);
			spinor3prod2(between,spinor2,s1tot);

			term1 			= 	2*bondCH;
			p1co[i][0] 		= 	co[i][0] + term1*(s1tot[1]*s1tot[1] + s1tot[0]*s1tot[0] - 0.5);
			p1co[i][1] 		= 	co[i][1] + term1*(s1tot[1]*s1tot[2] + s1tot[0]*s1tot[3]);
			p1co[i][2] 		= 	co[i][2] + term1*(s1tot[1]*s1tot[3] - s1tot[0]*s1tot[2]);


			/*second side chain*/
			spinor1[0]    	= 	cos((dihedral[i+1] + 4*M_PI/3)/2);					/*sp3-hybridisation*/
			spinor1[1]    	= 	sin((dihedral[i+1] + 4*M_PI/3)/2);
			spinor2[0]    = 	cosarray[i+1];	/*seems to work much better*/
			spinor2[1]    = 	sinarray[i+1];
			//spinor2[0]    	= 	cos((M_PI - valence[i+1])/2);				/*fixed valences angles ->no variable ones*/
			//spinor2[1]    	= 	sin((M_PI - valence[i+1])/2);				/*indices not yet correctly fixed in space*/
			spinor3prod1(spinortot,spinor1,between);
			spinor3prod2(between,spinor2,s2tot);

			//  add to previous atom coordinate
			//term1 			= 	4*bondCH;
			//term1			=   2*bondCH;		/*do some trick!*/
			term1			=	2*bondCACO;		/*altered carbon bond length*/
			p2co[i][0] 		= 	co[i][0] + term1*(s2tot[1]*s2tot[1] + s2tot[0]*s2tot[0] - 0.5);
			p2co[i][1] 		= 	co[i][1] + term1*(s2tot[1]*s2tot[2] + s2tot[0]*s2tot[3]);
			p2co[i][2] 		= 	co[i][2] + term1*(s2tot[1]*s2tot[3] - s2tot[0]*s2tot[2]);


			//NEW CODE ADDITIONAL CODE: MODELING THE SIDECHAIN
			/*1st hydrogen on side chain carbon*/

			hoek =	((double) rand()/(double) RAND_MAX)*M_PI*2 - M_PI;

			spinor1[0]    	= 	cos(hoek/2);	/*spinors are fixed -- can be precalculated*/ /*indices dihedral mapping Ok*/
			spinor1[1]    	= 	sin(hoek/2);
			spinor2[0]   	= 	cos((M_PI - (M_PI*108/180))/2);		/*fixed valences angles ->no variable ones*/
			spinor2[1]    = 	sin((M_PI - (M_PI*108/180))/2);		/*indices not yet correctly fixed in space*/

			spinor3prod1(s2tot,spinor1,between);
			spinor3prod2(between,spinor2,s3tot);

			term1 			= 	2*bondCH;
			s1co[i][0] 		= 	p2co[i][0] + term1*(s3tot[1]*s3tot[1] + s3tot[0]*s3tot[0] - 0.5);
			s1co[i][1] 		= 	p2co[i][1] + term1*(s3tot[1]*s3tot[2] + s3tot[0]*s3tot[3]);
			s1co[i][2] 		= 	p2co[i][2] + term1*(s3tot[1]*s3tot[3] - s3tot[0]*s3tot[2]);


			/*2nd hydrogen on side chain carbon*/
			spinor1[0]    	= 	cos((hoek + 2*M_PI/3)/2);	/*spinors are fixed -- can be precalculated*/ /*indices dihedral mapping Ok*/
			spinor1[1]    	= 	sin((hoek + 2*M_PI/3)/2);
			spinor2[0]   	= 	cos((M_PI - (M_PI*108/180))/2);		/*fixed valences angles ->no variable ones*/
			spinor2[1]    = 	sin((M_PI - (M_PI*108/180))/2);		/*indices not yet correctly fixed in space*/

			spinor3prod1(s2tot,spinor1,between);
			spinor3prod2(between,spinor2,s3tot);

			term1 			= 	2*bondCH;
			s1co[i+1][0] 		= 	p2co[i][0] + term1*(s3tot[1]*s3tot[1] + s3tot[0]*s3tot[0] - 0.5);
			s1co[i+1][1] 		= 	p2co[i][1] + term1*(s3tot[1]*s3tot[2] + s3tot[0]*s3tot[3]);
			s1co[i+1][2] 		= 	p2co[i][2] + term1*(s3tot[1]*s3tot[3] - s3tot[0]*s3tot[2]);



			/*3rd hydrogen on side chain carbon*/
			spinor1[0]    	= 	cos((hoek + 4*M_PI/3)/2);					/*sp3-hybridisation*/
			spinor1[1]    	= 	sin((hoek + 4*M_PI/3)/2);
			spinor2[0]    	= 	cos((M_PI - (M_PI*108/180))/2);				/*fixed valences angles ->no variable ones*/
			spinor2[1]    	= 	sin((M_PI - (M_PI*108/180))/2);				/*indices not yet correctly fixed in space*/

			spinor3prod1(s2tot,spinor1,between);
			spinor3prod2(between,spinor2,s3tot);

			//  add to previous atom coordinate
			term1			=   2*bondCH;
			s1co[i+2][0] 		= 	p2co[i][0] + term1*(s3tot[1]*s3tot[1] + s3tot[0]*s3tot[0] - 0.5);
			s1co[i+2][1] 		= 	p2co[i][1] + term1*(s3tot[1]*s3tot[2] + s3tot[0]*s3tot[3]);
			s1co[i+2][2] 		= 	p2co[i][2] + term1*(s3tot[1]*s3tot[3] - s3tot[0]*s3tot[2]);



			//			printf("Coordinates: %8.3f %8.3f %8.3f\n",p2co[i][0],p2co[i][1],p2co[i][2]);
			//			getchar();

		}
		else if (marker[i] == 1)	/*sp2-hybridization*/
		{
			if (i < number_atoms - 1)
			{
				spinor1[0]   	= 	cos((dihedral[i+1] + M_PI)/2);	/*spinors are fixed -- can be precalculated*/
				spinor1[1]    	= 	sin((dihedral[i+1] + M_PI)/2);
				spinor2[0]    = 	cosarray[i+1];	/*seems to work much better*/
				spinor2[1]    = 	sinarray[i+1];

				//spinor2[0]    	= 	cos((M_PI - valence[i+1])/2);		/*fixed valences angles ->no variable ones*/
				//spinor2[1]    	= 	sin((M_PI - valence[i+1])/2);		/*indices not yet correctly fixed in space*/
			}
			else
			{
				/*but this one must yet be randomized*/
				spinor1[0]   	= 	cos((dihedral[i-2] + M_PI)/2);	/*spinors are fixed -- can be precalculated*/
				spinor1[1]    	= 	sin((dihedral[i-2] + M_PI)/2);
				spinor2[0]    = 	cosarray[i-2];	/*seems to work much better*/
				spinor2[1]    = 	sinarray[i-2];

				//spinor2[0]    	= 	cos((M_PI - valence[i-2])/2);		/*fixed valences angles ->no variable ones*/
				//spinor2[1]    	= 	sin((M_PI - valence[i-2])/2);		/*indices not yet correctly fixed in space*/
			}

			spinor3prod1(spinortot,spinor1,between);
			spinor3prod2(between,spinor2,s1tot);

			term1 		   = 	2*bondCO;
			p1co[i][0] 	   = 	co[i][0] + term1*(s1tot[1]*s1tot[1] + s1tot[0]*s1tot[0] - 0.5);
			p1co[i][1] 	   = 	co[i][1] + term1*(s1tot[1]*s1tot[2] + s1tot[0]*s1tot[3]);
			p1co[i][2] 	   = 	co[i][2] + term1*(s1tot[1]*s1tot[3] - s1tot[0]*s1tot[2]);

			p2co[i][0] = 0.0;	/*fill with zeros for definiteness*/ /*should not be printed!!!*/
			p2co[i][1] = 0.0;
			p2co[i][2] = 0.0;
		}
	}

	/*Print all coordinates*/
	//		for (int j = 0; j < number_atoms; j++)
	//		{
	//			printf("Coordinates: %8.3f %8.3f %8.3f\n",co[j][0],co[j][1],co[j][2]);
	//		}
	//		printf(" \n");
	//
	//		for (int j = 0; j < number_atoms; j++)
	//		{
	//			printf("Coordinates: %8.3f %8.3f %8.3f\n",p1co[j][0],p1co[j][1],p1co[j][2]);
	//		}
	//		printf(" \n");
	//
	//
	//		for (int j = 0; j < number_atoms; j++)
	//		{
	//			printf("Coordinates: %8.3f %8.3f %8.3f\n",p2co[j][0],p2co[j][1],p2co[j][2]);
	//		}
	//getchar();
}

// Set a constant CA-CB distance
void set_cacb_dist2(double **co, double **p1co, int natoms, double dist, char *seq)
{
	double scale;
	double v[3];
	for(int i = 1; i < natoms; i += 3) // Screen CB and CA
	{
		if(seq[i/3] != 'G')
		{
			// Vector CA-->CB
			v[0] = p1co[i][0] - co[i][0];
			v[1] = p1co[i][1] - co[i][1];
			v[2] = p1co[i][2] - co[i][2];
			// Scale factor
			scale = dist / sqrt( pow(v[0],2) + pow(v[1],2) + pow(v[2],2) );
			// Apply scaling factor to CB position
			p1co[i][0] = co[i][0] + v[0] * scale;
			p1co[i][1] = co[i][1] + v[1] * scale;
			p1co[i][2] = co[i][2] + v[2] * scale;
		}
	}
}
void isafrotatetreef_cg(double **co,				// loop coordinates
		double **p1co,				// CB && O coordinates
		int number_atoms,		// atoms in loop
		double *bond,			// bond lengths loop
		double *dihedral,			// dihedral angles loop
		double *cosarrayd,		// preallocated -cosines (BB Dihedral angles, computed inside)
		double *sinarrayd,		// preallocated sines    (BB Dihedral angles, computed inside)
		double *cosarrayd4,		// preallocated -cosines (Dihedral angles related, computed inside)
		double *sinarrayd4,		// preallocated sines    (Dihedral angles related, computed inside)
		double *cosarray,		// preallocated -cosines (Valence angles, computed outside)
		double *sinarray,		// preallocated sines    (Valence angles, computed outside)
		char *seq,				// sequence (one-letter code)
		double *full,
		int mode, 			// mode of operation: 0=(N,CA,C), 1=(N,CA,C,O), 2=(N,CA,C,CB), 3=(N,CA,C,O,CB)
		double factor			// CA-CB common scaling factor (set to <=0 for a default constant distance)
)
{
	double spinor1[2];
	double spinor2[2];
	double spinortot[4];
	double s1tot[4];
	double between[4];
	double dist2 = 0.0; // Current CA-CB distance squared
	int nmax = ((number_atoms/3)-1);
	int cont=0;

	// Pre-computing cos() and sin() from dihedral angles...
	cosarrayd[0]	=	1;
	sinarrayd[0]	=	0;
	cosarrayd[1]	=	1;
	sinarrayd[1]	=	0;
	cosarrayd[2]	=	1;
	sinarrayd[2]	=	0;

	for (int i = 3; i < number_atoms; i++)
	{
		//		sincos( dihedral[i]*0.5,&sinarrayd[i], &cosarrayd[i]);
		cosarrayd[i]	=	cos(dihedral[i]/2.0);
		sinarrayd[i]	=	sin(dihedral[i]/2.0);
	}

	// CB and O precalculations...
	if(mode>0)
	{
		cosarrayd4[1] = cos((dihedral[2] + 4*M_PI/3)/2);
		sinarrayd4[1] = sin((dihedral[2] + 4*M_PI/3)/2);
		// O is M_PI --> 180 deg
		cosarrayd4[2] = cos((dihedral[3] + M_PI)/2); // -sinarrayd[i*3+3]; // cos (D+pi/2)= -sin(pi)  cos((dihedral[i*3+3] + M_PI)/2);
		sinarrayd4[2] = sin((dihedral[3] + M_PI)/2); // co

		for( int i=1; i<nmax; i++) // Screen residues
		{
			// CB is 4*M_PI/3 --> 240 deg
			// The extra +1 considers the angle required to position the next BB-atom (because it defines the CB/O position)
			// cosarrayd4[i*3+1] = cos((dihedral[i*3+2] + 4*M_PI/3)/2);
			// sinarrayd4[i*3+1] = sin((dihedral[i*3+2] + 4*M_PI/3)/2);
			//			cosarrayd4[i*3+1] = cos((dihedral[i*3+2] +  2*valOCON)/2);
			//			sinarrayd4[i*3+1] = sin((dihedral[i*3+2] +  2*valOCON)/2);
			cosarrayd4[i*3+1] =  -cosarrayd[i*3+2]*0.5-SQRT32*sinarrayd[i*3+2]; // cos(d + 2/3 pi) = (-cosD - sqrt(3)sinD)/2 cos((dihedral[i*3+2] + 4*M_PI/3)/2);
			sinarrayd4[i*3+1] =  SQRT32*cosarrayd[i*3+2]-sinarrayd[i*3+2]*0.5; //  sin(d + 2/3 pi) = (sqrt(3)cosD - sinD)/2 sin((dihedral[i*3+2] + 4*M_PI/3)/2);


			// O is M_PI --> 180 deg
			// cosarrayd4[i*3+2] = cos((dihedral[i*3+3] + M_PI)/2);	/*spinors are fixed -- can be precalculated*/
			// sinarrayd4[i*3+2] = sin((dihedral[i*3+3] + M_PI)/2);

			cosarrayd4[i*3+2] =  -sinarrayd[i*3+3]; // cos (D+pi/2)= -sin(pi)  cos((dihedral[i*3+3] + M_PI)/2);
			sinarrayd4[i*3+2] =   cosarrayd[i*3+3]; //  sin (D+pi/2)= -cos(pi)  sin((dihedral[i*3+3] + M_PI)/2);

		}

		// CB is 4*M_PI/3 --> 240 deg
		cosarrayd4[nmax*3+1] = cos((dihedral[nmax*3+2] + 4*M_PI/3)/2);
		sinarrayd4[nmax*3+1] = sin((dihedral[nmax*3+2] + 4*M_PI/3)/2);
		//		cosarrayd4[nmax*3+1] = cos((dihedral[nmax*3+2] + 2*valOCON)/2);
		//		sinarrayd4[nmax*3+1] = sin((dihedral[nmax*3+2] + 2*valOCON)/2);

		// O is M_PI --> 180 deg
		cosarrayd4[nmax*3+2] = cos((dihedral[0] + M_PI)/2);	/*spinors are fixed -- can be precalculated*/
		sinarrayd4[nmax*3+2] = sin((dihedral[0] + M_PI)/2);
	}

	//	for (int i = 0; i < nmax; i++)
	//	{
	//		cosarrayd4[i*3] = cos((dihedral[i*3+2] + 4*M_PI/3)/2);
	//		sinarrayd4[i*3] = sin((dihedral[i*3+2] + 4*M_PI/3)/2);
	//		cosarrayd4[i*3+1] = cos((dihedral[i*3+3] + M_PI)/2);
	//		sinarrayd4[i*3+1] = sin((dihedral[i*3+3] + M_PI)/2);
	//	}
	//	cosarrayd4[nmax*3] = cos((dihedral[nmax*3+2] + 4*M_PI/3)/2);
	//	sinarrayd4[nmax*3] = sin((dihedral[nmax*3+2] + 4*M_PI/3)/2);

	// Initial pose (stored in "full")
	spinortot[0] = full[0];
	spinortot[1] = full[1];
	spinortot[2] = full[2];
	spinortot[3] = full[3];

	for (int i = 1; i < number_atoms; i++) // Screen all atoms but first
	{
		// MOVING BACKBONE (N, CA, C)
		// Dihedrals stuff
		spinor1[0]    = 	cosarrayd[i]; // cos(dihedral[i]/2);
		spinor1[1]    = 	sinarrayd[i]; // sin(dihedral[i]/2);
		// Valences stuff
		spinor2[0]    = 	cosarray[i-1]; // cos((M_PI-valence[i])/2);
		spinor2[1]    = 	sinarray[i-1]; // sin((M_PI-valence[i])/2);

		spinor3prod1(spinortot,spinor1,between);
		spinor3prod2(between,spinor2,spinortot);
		co[i][0] = co[i-1][0] + 2*bond[i]*(spinortot[1]*spinortot[1] + spinortot[0]*spinortot[0] - 0.5);
		co[i][1] = co[i-1][1] + 2*bond[i]*(spinortot[1]*spinortot[2] + spinortot[0]*spinortot[3]);
		co[i][2] = co[i-1][2] + 2*bond[i]*(spinortot[1]*spinortot[3] - spinortot[0]*spinortot[2]);

		// MOVING BRANCHED ATOMS (side-chain atoms)
		if(mode>0)
		{
			if( mode != 1  && i%3 == 1 )		// sp3-hybridization (moving atoms pending from CA)
			{
				// Side chain (-CB)
				if(i==1)
				{
					spinor1[0]    	= 	COS2MPI3; //cos(2*M_PI/3); // 120 degrees (initial reference)
					spinor1[1]    	= 	SIN2MPI3; //sin(2*M_PI/3);
					//					spinor1[0]    	= 	cos(valCBCAC);
					//					spinor1[1]    	= 	sin(valCBCAC);
				}
				else
				{
					//					spinor1[0]    	= 	cos((dihedral[i+1] + 4*M_PI/3)/2);
					//					spinor1[1]    	= 	sin((dihedral[i+1] + 4*M_PI/3)/2);
					spinor1[0]    	= 	cosarrayd4[i];
					spinor1[1]    	= 	sinarrayd4[i];
					//					spinor1[0]    	= 	cos((dihedral[i+1] + 2*valCBCAC)/2);
					//					spinor1[1]    	= 	sin((dihedral[i+1] + 2*valCBCAC)/2);

				}
				spinor2[0]    = 	cosarray[i]; // cos((M_PI - valence[i+1])/2)
				spinor2[1]    = 	sinarray[i]; // sin((M_PI - valence[i+1])/2)

				//				if(seq[i/3]=='P' || seq[i/3]=='p')
				//				{
				//					spinor2[0]    = 	cos( (M_PI - PROvalNCACB)/2.0 ); // NAB-angle  -  cos((M_PI - valence[i+1])/2)
				//					spinor2[1]    = 	sin( (M_PI - PROvalNCACB)/2.0 ); // NAB-angle  -  sin((M_PI - valence[i+1])/2)
				//				}
				//				else
				//				{
				//					spinor2[0]    = 	cos( (M_PI - valNCACB)/2.0 ); // NAB-angle  -  cos((M_PI - valence[i+1])/2)
				//					spinor2[1]    = 	sin( (M_PI - valNCACB)/2.0 ); // NAB-angle  -  sin((M_PI - valence[i+1])/2)
				//				}


				spinor3prod1(spinortot,spinor1,between);
				spinor3prod2(between,spinor2,s1tot);
				//				p1co[i][0] 		= 	co[i][0] + bondCACO2*(s1tot[1]*s1tot[1] + s1tot[0]*s1tot[0] - 0.5);
				//				p1co[i][1] 		= 	co[i][1] + bondCACO2*(s1tot[1]*s1tot[2] + s1tot[0]*s1tot[3]);
				//				p1co[i][2] 		= 	co[i][2] + bondCACO2*(s1tot[1]*s1tot[3] - s1tot[0]*s1tot[2]);

				// Select current CA-CB distance
				if(factor <= 0.0) // constant distance
				{
					dist2 = bondCACB2;
				}
				else // residue dependent distance
				{
					switch (seq[i/3])
					{
					case 'A':     //"ALA" );
						dist2=3.25;
						break;
					case 'C':     //"CYS";
						dist2=4.17;
						break;
					case 'D':		//"ASP";
						dist2=4.25;
						break;
					case 'E':     //"GLU";
						dist2=4.90;
						break;
					case 'F':     //"PHE"
						dist2=4.81;
						break;
						//					  case 'G':		//"GLY"
						//						  // dist2=2.4;
						//						  fprintf(stderr,"rcd> Error, this should never happen...\n");
						//						  exit(1);
						//						  break;
					case 'H':     //"HIS"
						dist2=4.66;
						break;
					case 'I':       //"ILE";
						dist2=4.32;
						break;
					case 'K':       //"LYS" );
						dist2=5.30;
						break;
					case 'L':	      //"LEU" );
						dist2=4.76;
						break;
					case 'M':	      // "MET" );
						dist2=4.95;
						break;
					case 'N':	      // "ASN" );
						dist2=4.29;
						break;
					case 'P':	      // "PRO" );
						dist2=3.58;
						break;
					case 'p':	      // "cis-PRO" );
						dist2=3.58;
						break;
					case 'Q':	      // "GLN" );
						dist2=4.95;
						break;
					case 'R':	      // "ARG" );
						dist2=5.64;
						break;
					case 'S':	      // "SER" );
						dist2=3.56;
						break;
					case 'T':	      //  "THR" );
						dist2=3.79;
						break;
					case 'V':	      // "VAL" );
						dist2=3.90;
						break;
					case 'W':	      //  "TRP" );
						dist2=4.91;
						break;
					case 'Y':	      // "TYR" );
						dist2=4.71;
						break;
					}

					dist2 *= factor;
				}

				// Store new calculated coordinates for CB-atom
				p1co[i][0] 		= 	co[i][0] + dist2 * (s1tot[1]*s1tot[1] + s1tot[0]*s1tot[0] - 0.5);
				p1co[i][1] 		= 	co[i][1] + dist2 * (s1tot[1]*s1tot[2] + s1tot[0]*s1tot[3]);
				p1co[i][2] 		= 	co[i][2] + dist2 * (s1tot[1]*s1tot[3] - s1tot[0]*s1tot[2]);

			}
			if( mode != 2 && i%3 == 2 ) // sp2-hybridization (moving atoms pending from either N or C)
			{
				// Oxygen (-O)
				spinor1[0]    	= 	cosarrayd4[i];
				spinor1[1]    	= 	sinarrayd4[i];
				spinor2[0]    = 	cosarray[i]; // cos((M_PI - valence[i+1])/2)
				spinor2[1]    = 	sinarray[i]; // sin((M_PI - valence[i+1])/2)
				//				spinor2[0]    = 	cos( (M_PI - valCACOO)/2.0 ); // cos((M_PI - valence[i+1])/2)
				//				spinor2[1]    = 	sin( (M_PI - valCACOO)/2.0 ); // sin((M_PI - valence[i+1])/2)

				//				if (i < number_atoms - 1)
				//				{
				////					spinor1[0]   	= 	cos((dihedral[i+1] + M_PI)/2);	/*spinors are fixed -- can be precalculated*/
				////					spinor1[1]    	= 	sin((dihedral[i+1] + M_PI)/2);
				//					spinor1[0]    	= 	cosarrayd4[i];
				//					spinor1[1]    	= 	sinarrayd4[i];
				//					spinor2[0]    = 	cosarray[i+1]; // cos((M_PI - valence[i+1])/2)
				//					spinor2[1]    = 	sinarray[i+1]; // sin((M_PI - valence[i+1])/2)
				//				}
				//				else
				//				{
				////					spinor1[0]   	= 	cos((dihedral[0] + M_PI)/2);	/*spinors are fixed -- can be precalculated*/
				////					spinor1[1]    	= 	sin((dihedral[0] + M_PI)/2);
				//					spinor1[0]    	= 	cosarrayd4[i];
				//					spinor1[1]    	= 	sinarrayd4[i];
				//					// Mon: taking the CA-C-O angle from previous aminoacid ("chapa")
				//					spinor2[0]    = 	cosarray[i-2]; // cos((M_PI - valence[i-2])/2)
				//					spinor2[1]    = 	sinarray[i-2]; // sin((M_PI - valence[i-2])/2)
				//				}

				spinor3prod1(spinortot,spinor1,between);
				spinor3prod2(between,spinor2,s1tot);

				// Store new calculated coordinates for O-atom
				p1co[i][0] 	   = 	co[i][0] + bondCO2*(s1tot[1]*s1tot[1] + s1tot[0]*s1tot[0] - 0.5);
				p1co[i][1] 	   = 	co[i][1] + bondCO2*(s1tot[1]*s1tot[2] + s1tot[0]*s1tot[3]);
				p1co[i][2] 	   = 	co[i][2] + bondCO2*(s1tot[1]*s1tot[3] - s1tot[0]*s1tot[2]);
			}
		}
	}
	//	free(cosarrayd4);
	//	free(sinarrayd4);
}



// "cg" version is what we use now...

// "cg" version is what we use now...
void isafrotatetreef_cg_old(double **co,				// loop coordinates
		double **p1co,				// CB && O coordinates
		int number_atoms,		// atoms in loop
		double *bond,			// bond lengths loop
		double *dihedral,			// dihedral angles loop
		double *cosarrayd,		// preallocated -cosines (BB Dihedral angles, computed inside)
		double *sinarrayd,		// preallocated sines    (BB Dihedral angles, computed inside)
		double *cosarrayd4,		// preallocated -cosines (Dihedral angles related, computed inside)
		double *sinarrayd4,		// preallocated sines    (Dihedral angles related, computed inside)
		double *cosarray,		// preallocated -cosines (Valence angles, computed outside)
		double *sinarray,		// preallocated sines    (Valence angles, computed outside)
		char *seq,				// sequence (one-letter code)
		double *full,
		int mode, 			// mode of operation: 0=(N,CA,C), 1=(N,CA,C,O), 2=(N,CA,C,CB), 3=(N,CA,C,O,CB)
		double factor			// CA-CB common scaling factor (set to <=0 for a default constant distance)
)
{
	double spinor1[2];
	double spinor2[2];
	double spinortot[4];
	double s1tot[4];
	double between[4];
	double dist2 = 0.0; // Current CA-CB distance squared
	int nmax = ((number_atoms/3)-1);
	int cont=0;


	// Pre-computing cos() and sin() from dihedral angles...
	cosarrayd[0]	=	1;
	sinarrayd[0]	=	0;
	cosarrayd[1]	=	1;
	sinarrayd[1]	=	0;
	cosarrayd[2]	=	1;
	sinarrayd[2]	=	0;

	for (int i = 3; i < number_atoms; i++)
	{
		// sincos( dihedral[i]*0.5,&sinarrayd[i], &cosarrayd[i]);
		cosarrayd[i]	=	cos(dihedral[i]/2.0);
		sinarrayd[i]	=	sin(dihedral[i]/2.0);
	}

	// CB and O precalculations...
	if(mode>0)
	{


		for( int i=0; i<nmax*3; i+=3) // Screen residues
		{
			// CB is 4*M_PI/3 --> 240 deg
			// The extra +1 considers the angle required to position the next BB-atom (because it defines the CB/O position)
			cosarrayd4[i+1] =  cos((dihedral[i*3+2] + 4*M_PI/3)/2);
			sinarrayd4[i+1] =  sin((dihedral[i*3+2] + 4*M_PI/3)/2);
			// cos((dihedral[i+1] + 2*valCBCAC)/2);

			// fprintf(stderr,"%f %f\n", cosarrayd4[i*3+1],-cosarrayd[i*3+2]*0.5-SQRT32*sinarrayd[i*3+2]);
			// fprintf(stderr,"%f %f\n", sinarrayd4[i*3+1],SQRT32*cosarrayd[i*3+2]-sinarrayd[i*3+2]*0.5);
			// O is M_PI --> 180 deg
			cosarrayd4[i+2] =  cos((dihedral[i*3+3] + M_PI)/2);
			sinarrayd4[i+2] =  sin((dihedral[i*3+3] + M_PI)/2);
			// fprintf(stderr,"%f %f\n", cosarrayd4[i+2],-cosarrayd[i*3+2]*0.5-SQRT32*sinarrayd[i*3+2]);
			// fprintf(stderr,"%f %f\n", sinarrayd4[i+2],SQRT32*cosarrayd[i*3+2]-sinarrayd[i*3+2]*0.5);
			// cos( (M_PI - valCACOO)/2.0 );

		}

		// getchar();
		// CB is 4*M_PI/3 --> 240 deg
		cosarrayd4[nmax*3+1] = cos((dihedral[nmax*3+2] + 4*M_PI/3)/2);
		sinarrayd4[nmax*3+1] = sin((dihedral[nmax*3+2] + 4*M_PI/3)/2);
		//		cosarrayd4[nmax*3+1] = cos((dihedral[nmax*3+2] + 2*valOCON)/2);
		//		sinarrayd4[nmax*3+1] = sin((dihedral[nmax*3+2] + 2*valOCON)/2);


		// O is M_PI --> 180 deg OJO!! las psi
		cosarrayd4[nmax*3+2] = cos((dihedral[0] + M_PI)/2);	/*spinors are fixed -- can be precalculated*/
		sinarrayd4[nmax*3+2] = sin((dihedral[0] + M_PI)/2);

	}



	// Initial pose (stored in "full")
	spinortot[0] = full[0];
	spinortot[1] = full[1];
	spinortot[2] = full[2];
	spinortot[3] = full[3];

	for (int i = 1; i < number_atoms; i++) // Screen all atoms but first
	{
		// MOVING BACKBONE (N, CA, C)
		// Dihedrals stuff
		spinor1[0]    = 	cosarrayd[i]; // cos(dihedral[i]/2);
		spinor1[1]    = 	sinarrayd[i]; // sin(dihedral[i]/2);
		// Valences stuff
		spinor2[0]    = 	cosarray[i-1]; // cos((M_PI-valence[i])/2);
		spinor2[1]    = 	sinarray[i-1]; // sin((M_PI-valence[i])/2);

		spinor3prod1(spinortot,spinor1,between);
		spinor3prod2(between,spinor2,spinortot);
		co[i][0] = co[i-1][0] + 2*bond[i]*(spinortot[1]*spinortot[1] + spinortot[0]*spinortot[0] - 0.5);
		co[i][1] = co[i-1][1] + 2*bond[i]*(spinortot[1]*spinortot[2] + spinortot[0]*spinortot[3]);
		co[i][2] = co[i-1][2] + 2*bond[i]*(spinortot[1]*spinortot[3] - spinortot[0]*spinortot[2]);

		// MOVING BRANCHED ATOMS (side-chain atoms)
		if(mode>0)
		{
			if( mode != 1  && i%3 == 1 )		// sp3-hybridization (moving atoms pending from CA)
			{
				// Side chain (-CB)
				if(i==1)
				{
					spinor1[0]    	= 	COS2MPI3; // cos(2*M_PI/3); // 120 degrees (initial reference)
					spinor1[1]    	= 	SIN2MPI3; // sin(2*M_PI/3);
					//					spinor1[0]    	= 	cos(valCBCAC);
					//					spinor1[1]    	= 	sin(valCBCAC);
				}
				else
				{
					//					spinor1[0]    	= 	cos((dihedral[i+1] + 4*M_PI/3)/2);
					//					spinor1[1]    	= 	sin((dihedral[i+1] + 4*M_PI/3)/2);
					spinor1[0]    	= 	cosarrayd4[i];
					spinor1[1]    	= 	sinarrayd4[i];
					//					spinor1[0]    	= 	cos((dihedral[i+1] + 2*valCBCAC)/2);
					//					spinor1[1]    	= 	sin((dihedral[i+1] + 2*valCBCAC)/2);

				}
				//				spinor2[0]    = 	cosarray[i]; // cos((M_PI - valence[i+1])/2)
				//				spinor2[1]    = 	sinarray[i]; // sin((M_PI - valence[i+1])/2)

				if(seq[i/3]=='P' || seq[i/3]=='p')
				{
					spinor2[0]    = COSPROvalNCACB;	// cos( (M_PI - PROvalNCACB)/2.0 ); // NAB-angle  -  cos((M_PI - valence[i+1])/2)
					spinor2[1]    = SINPROvalNCACB;	// sin( (M_PI - PROvalNCACB)/2.0 ); // NAB-angle  -  sin((M_PI - valence[i+1])/2)
				}
				else
				{
					spinor2[0]    = COSvalNCACB; //	cos( (M_PI - valNCACB)/2.0 ); // NAB-angle  -  cos((M_PI - valence[i+1])/2)
					spinor2[1]    = SINvalNCACB; //	sin( (M_PI - valNCACB)/2.0 ); // NAB-angle  -  sin((M_PI - valence[i+1])/2)
				}


				spinor3prod1(spinortot,spinor1,between);
				spinor3prod2(between,spinor2,s1tot);
				//				p1co[i][0] 		= 	co[i][0] + bondCACO2*(s1tot[1]*s1tot[1] + s1tot[0]*s1tot[0] - 0.5);
				//				p1co[i][1] 		= 	co[i][1] + bondCACO2*(s1tot[1]*s1tot[2] + s1tot[0]*s1tot[3]);
				//				p1co[i][2] 		= 	co[i][2] + bondCACO2*(s1tot[1]*s1tot[3] - s1tot[0]*s1tot[2]);

				// Select current CA-CB distance
				if(factor <= 0.0) // constant distance
				{
					dist2 = bondCACB2;
				}
				else // residue dependent distance
				{
					switch (seq[i/3])
					{
					case 'A':     //"ALA" );
						dist2=3.25;
						break;
					case 'C':     //"CYS";
						dist2=4.17;
						break;
					case 'D':		//"ASP";
						dist2=4.25;
						break;
					case 'E':     //"GLU";
						dist2=4.90;
						break;
					case 'F':     //"PHE"
						dist2=4.81;
						break;
						//					  case 'G':		//"GLY"
						//						  // dist2=2.4;
						//						  fprintf(stderr,"rcd> Error, this should never happen...\n");
						//						  exit(1);
						//						  break;
					case 'H':     //"HIS"
						dist2=4.66;
						break;
					case 'I':       //"ILE";
						dist2=4.32;
						break;
					case 'K':       //"LYS" );
						dist2=5.30;
						break;
					case 'L':	      //"LEU" );
						dist2=4.76;
						break;
					case 'M':	      // "MET" );
						dist2=4.95;
						break;
					case 'N':	      // "ASN" );
						dist2=4.29;
						break;
					case 'P':	      // "PRO" );
						dist2=3.58;
						break;
					case 'p':	      // "cis-PRO" );
						dist2=3.58;
						break;
					case 'Q':	      // "GLN" );
						dist2=4.95;
						break;
					case 'R':	      // "ARG" );
						dist2=5.64;
						break;
					case 'S':	      // "SER" );
						dist2=3.56;
						break;
					case 'T':	      //  "THR" );
						dist2=3.79;
						break;
					case 'V':	      // "VAL" );
						dist2=3.90;
						break;
					case 'W':	      //  "TRP" );
						dist2=4.91;
						break;
					case 'Y':	      // "TYR" );
						dist2=4.71;
						break;
					}

					dist2 *= factor;
				}

				// Store new calculated coordinates for CB-atom
				p1co[i][0] 		= 	co[i][0] + dist2 * (s1tot[1]*s1tot[1] + s1tot[0]*s1tot[0] - 0.5);
				p1co[i][1] 		= 	co[i][1] + dist2 * (s1tot[1]*s1tot[2] + s1tot[0]*s1tot[3]);
				p1co[i][2] 		= 	co[i][2] + dist2 * (s1tot[1]*s1tot[3] - s1tot[0]*s1tot[2]);

			}
			if( mode != 2 && i%3 == 2 ) // sp2-hybridization (moving atoms pending from either N or C)
			{
				// Oxygen (-O)
				spinor1[0]    	= 	cosarrayd4[i];
				spinor1[1]    	= 	sinarrayd4[i];
				//				spinor2[0]    = 	cosarray[i]; // cos((M_PI - valence[i+1])/2)
				//				spinor2[1]    = 	sinarray[i]; // sin((M_PI - valence[i+1])/2)
				spinor2[0]    = 	COSvalCACOO; // cos( (M_PI - valCACOO)/2.0 ); // cos((M_PI - valence[i+1])/2)
				spinor2[1]    = 	SINvalCACOO; // sin( (M_PI - valCACOO)/2.0 ); // sin((M_PI - valence[i+1])/2)


				spinor3prod1(spinortot,spinor1,between);
				spinor3prod2(between,spinor2,s1tot);

				// Store new calculated coordinates for O-atom
				p1co[i][0] 	   = 	co[i][0] + bondCO2*(s1tot[1]*s1tot[1] + s1tot[0]*s1tot[0] - 0.5);
				p1co[i][1] 	   = 	co[i][1] + bondCO2*(s1tot[1]*s1tot[2] + s1tot[0]*s1tot[3]);
				p1co[i][2] 	   = 	co[i][2] + bondCO2*(s1tot[1]*s1tot[3] - s1tot[0]*s1tot[2]);
			}
		}
	}

}





void isafrotatetreef(double **co,				/*loop coordinates*/
		double **p1co,				/*coordinates 1st parallel chain*/
		double **p2co,				/*coordinates 2nd parallel chain*/
		double *posespinor,			/*pose spinor for anchor*/
		int number_atoms,		/*# atoms in loop*/
		double *valence,		/*valence angles loop*/
		double *bond,			/*bond lengths loop*/
		double *cosarrayd,		/*preallocated -cosines*/
		double *sinarrayd,		/*preallocated sines*/
		double *cosarray,		/*preallocated -cosines*/
		double *sinarray,		/*preallocated sines*/
		double *dihedral,			/*dihedral angles loop*/
		double *p1dihedral,			/*dihedral angles loop*/
		double *p2dihedral,			/*dihedral angles loop*/
		int *marker,
		double *spinor1,			/*auxiliary spinor 1-preallocated*/
		double *spinor2,			/*auxiliary spinor 2-preallocated*/
		double *spinortot,			/*total spinor-preallocated*/
		double *s1tot,
		double *s2tot,
		double *between,
		double *full
)
{

	/*Print all coordinates*/
	//	for (int j = 0; j < number_atoms; j++)
	//	{
	//		printf("Coordinates: %8.3f %8.3f %8.3f\n",co[j][0],co[j][1],co[j][2]);
	//	}
	//	printf(" \n");


	double term1;	/*intermediate computational term*/


	spinor1[0]    = 	cos((dihedral[1] + M_PI)/2);
	spinor1[1]    = 	sin((dihedral[1]+ M_PI)/2);
	spinor2[0]    = 	cos((M_PI-valence[1])/2);
	spinor2[1]    = 	sin((M_PI-valence[1])/2);

	spinor3prod1(full,spinor1,between);
	spinor3prod2(between,spinor2,spinortot);


	// N's H
	term1 	 = 2*bondCH;	/*add to previous atom coordinate*/
	p1co[0][0] = co[0][0] + term1*(spinortot[1]*spinortot[1] + spinortot[0]*spinortot[0] - 0.5);
	p1co[0][1] = co[0][1] + term1*(spinortot[1]*spinortot[2] + spinortot[0]*spinortot[3]);
	p1co[0][2] = co[0][2] + term1*(spinortot[1]*spinortot[3] - spinortot[0]*spinortot[2]);

	// Initial pose (stored in "full")
	spinortot[0] = full[0];	/*Actual start*/
	spinortot[1] = full[1];
	spinortot[2] = full[2];
	spinortot[3] = full[3];

	// Extension for hydrogens and oxygens of loop
	for (int i = 1; i < number_atoms; i++)
	{
		// MOVING BACKBONE (N, CA, C)
		spinor1[0]    = 	cosarrayd[i];
		spinor1[1]    = 	sinarrayd[i];
		spinor2[0]    = 	cosarray[i];	/*seems to work much better*/
		spinor2[1]    = 	sinarray[i];
		//spinor1[0]    = 	cos(dihedral[i]/2);
		//spinor1[1]    = 	sin(dihedral[i]/2);
		//spinor2[0]    = 	cos((M_PI-valence[i])/2);
		//spinor2[1]    = 	sin((M_PI-valence[i])/2);
		spinor3prod1(spinortot,spinor1,between);
		spinor3prod2(between,spinor2,spinortot);
		term1 	 = 2*bond[i];	/*add to previous atom coordinate*/
		co[i][0] = co[i-1][0] + term1*(spinortot[1]*spinortot[1] + spinortot[0]*spinortot[0] - 0.5);
		co[i][1] = co[i-1][1] + term1*(spinortot[1]*spinortot[2] + spinortot[0]*spinortot[3]);
		co[i][2] = co[i-1][2] + term1*(spinortot[1]*spinortot[3] - spinortot[0]*spinortot[2]);

		//		fprintf(stderr,"i= %2d  marker[i]= %d\n",i,marker[i]);
		// MOVING BRANCHED ATOMS (side-chain atoms)
		if (marker[i] == 2)		// sp3-hybridization (moving atoms pending from CA)
		{
			// first side chain (-H)
			spinor1[0]    	= 	cos((dihedral[i+1] + 2*M_PI/3)/2);	/*spinors are fixed -- can be precalculated*/ /*indices dihedral mapping Ok*/
			spinor1[1]    	= 	sin((dihedral[i+1] + 2*M_PI/3)/2);
			spinor2[0]    = 	cosarray[i+1];	/*seems to work much better*/
			spinor2[1]    = 	sinarray[i+1];
			//			spinor2[0]   	= 	cos((M_PI - valence[i+1])/2);		/*fixed valences angles ->no variable ones*/
			//			spinor2[1]    = 	sin((M_PI - valence[i+1])/2);		/*indices not yet correctly fixed in space*/
			spinor3prod1(spinortot,spinor1,between);
			spinor3prod2(between,spinor2,s1tot);
			term1 			= 	2*bondCH;
			p1co[i][0] 		= 	co[i][0] + term1*(s1tot[1]*s1tot[1] + s1tot[0]*s1tot[0] - 0.5);
			p1co[i][1] 		= 	co[i][1] + term1*(s1tot[1]*s1tot[2] + s1tot[0]*s1tot[3]);
			p1co[i][2] 		= 	co[i][2] + term1*(s1tot[1]*s1tot[3] - s1tot[0]*s1tot[2]);

			// second side chain (-R)
			spinor1[0]    	= 	cos((dihedral[i+1] + 4*M_PI/3)/2);					/*sp3-hybridisation*/
			spinor1[1]    	= 	sin((dihedral[i+1] + 4*M_PI/3)/2);
			spinor2[0]    = 	cosarray[i+1];	/*seems to work much better*/
			spinor2[1]    = 	sinarray[i+1];
			//			spinor2[0]    	= 	cos((M_PI - valence[i+1])/2);				/*fixed valences angles ->no variable ones*/
			//			spinor2[1]    	= 	sin((M_PI - valence[i+1])/2);				/*indices not yet correctly fixed in space*/
			spinor3prod1(spinortot,spinor1,between);
			spinor3prod2(between,spinor2,s2tot);
			term1			=   2*bondCH;				/*add to previous atom coordinate*/	//THIS IS ALKYL GROUP OF PEPTIDE!!!
			p2co[i][0] 		= 	co[i][0] + term1*(s2tot[1]*s2tot[1] + s2tot[0]*s2tot[0] - 0.5);
			p2co[i][1] 		= 	co[i][1] + term1*(s2tot[1]*s2tot[2] + s2tot[0]*s2tot[3]);
			p2co[i][2] 		= 	co[i][2] + term1*(s2tot[1]*s2tot[3] - s2tot[0]*s2tot[2]);
		}
		else if (marker[i] == 1)	// sp2-hybridization (moving atoms pending from either N or C)
		{
			if (i < number_atoms - 1)
			{
				spinor1[0]   	= 	cos((dihedral[i+1] + M_PI)/2);	/*spinors are fixed -- can be precalculated*/
				spinor1[1]    	= 	sin((dihedral[i+1] + M_PI)/2);
				spinor2[0]    = 	cosarray[i+1];	/*seems to work much better*/
				spinor2[1]    = 	sinarray[i+1];
				//spinor2[0]    	= 	cos((M_PI - valence[i+1])/2);		/*fixed valences angles ->no variable ones*/
				//spinor2[1]    	= 	sin((M_PI - valence[i+1])/2);		/*indices not yet correctly fixed in space*/
			}
			else // Last atom (Last dihedral is different --> dihedral[0])
			{
				spinor1[0]   	= 	cos((dihedral[0] + M_PI)/2);	/*spinors are fixed -- can be precalculated*/
				spinor1[1]    	= 	sin((dihedral[0] + M_PI)/2);
				spinor2[0]    = 	cosarray[i-2];	// Mon: taking the CA-C-O angle from previous aminoacid ("chapa")
				spinor2[1]    = 	sinarray[i-2];
				//spinor2[0]    	= 	cos((M_PI - valence[i-2])/2);		/*fixed valences angles ->no variable ones*/
				//spinor2[1]    	= 	sin((M_PI - valence[i-2])/2);		/*indices not yet correctly fixed in space*/
			}

			spinor3prod1(spinortot,spinor1,between);
			spinor3prod2(between,spinor2,s1tot);

			term1 		   = 	2*bondCO;
			p1co[i][0] 	   = 	co[i][0] + term1*(s1tot[1]*s1tot[1] + s1tot[0]*s1tot[0] - 0.5);
			p1co[i][1] 	   = 	co[i][1] + term1*(s1tot[1]*s1tot[2] + s1tot[0]*s1tot[3]);
			p1co[i][2] 	   = 	co[i][2] + term1*(s1tot[1]*s1tot[3] - s1tot[0]*s1tot[2]);

			p2co[i][0] = 0.0;	/*fill with zeros for definiteness*/ /*should not be printed!!!*/
			p2co[i][1] = 0.0;
			p2co[i][2] = 0.0;
		}
	}

	//		/*Print all coordinates*/
	//		for (int j = 0; j < number_atoms; j++)
	//		{
	//			printf("Coordinates: %8.3f %8.3f %8.3f\n",co[j][0],co[j][1],co[j][2]);
	//		}
	//		printf(" \n");
	//
	//		for (int j = 0; j < number_atoms; j++)
	//		{
	//			printf("Coordinates: %8.3f %8.3f %8.3f\n",p1co[j][0],p1co[j][1],p1co[j][2]);
	//		}
	//		printf(" \n");
	//
	//
	//		for (int j = 0; j < number_atoms; j++)
	//		{
	//			printf("Coordinates: %8.3f %8.3f %8.3f\n",p2co[j][0],p2co[j][1],p2co[j][2]);
	//		}
	//		getchar();
}


void isafrotatetreef2(double **co,				/*loop coordinates*/
		double **p1co,				/*coordinates 1st parallel chain*/
		double **p2co,				/*coordinates 2nd parallel chain*/
		double *posespinor,			/*pose spinor for anchor*/
		int number_atoms,		/*# atoms in loop*/
		double *valence,		/*valence angles loop*/
		double *bond,			/*bond lengths loop*/
		double *cosarrayd,		/*preallocated -cosines*/
		double *sinarrayd,		/*preallocated sines*/
		double *cosarray,		/*preallocated -cosines*/
		double *sinarray,		/*preallocated sines*/
		double *dihedral,			/*dihedral angles loop*/
		double *p1dihedral,			/*dihedral angles loop*/
		double *p2dihedral,			/*dihedral angles loop*/
		int *marker,
		double *spinor1,			/*auxiliary spinor 1-preallocated*/
		double *spinor2,			/*auxiliary spinor 2-preallocated*/
		double *spinortot,			/*total spinor-preallocated*/
		double *s1tot,
		double *s2tot,
		double *between,
		double *full
)
{

	double term1;	/*intermediate computational term*/


	spinor1[0]    = 	cos((dihedral[1] + M_PI)/2);
	spinor1[1]    = 	sin((dihedral[1]+ M_PI)/2);
	spinor2[0]    = 	cos((M_PI-valence[1])/2);
	spinor2[1]    = 	sin((M_PI-valence[1])/2);

	spinor3prod1(full,spinor1,between);
	spinor3prod2(between,spinor2,spinortot);


	// N's H
	term1 	 = 2*bondCH;	/*add to previous atom coordinate*/
	p1co[0][0] = co[0][0] + term1*(spinortot[1]*spinortot[1] + spinortot[0]*spinortot[0] - 0.5);
	p1co[0][1] = co[0][1] + term1*(spinortot[1]*spinortot[2] + spinortot[0]*spinortot[3]);
	p1co[0][2] = co[0][2] + term1*(spinortot[1]*spinortot[3] - spinortot[0]*spinortot[2]);

	// Initial pose (stored in "full")
	spinortot[0] = full[0];	/*Actual start*/
	spinortot[1] = full[1];
	spinortot[2] = full[2];
	spinortot[3] = full[3];

	// Extension for hydrogens and oxygens of loop
	for (int i = 1; i < number_atoms; i++)
	{
		// MOVING BACKBONE (N, CA, C)
		spinor1[0]    = 	cosarrayd[i];
		spinor1[1]    = 	sinarrayd[i];
		spinor2[0]    = 	cosarray[i];	/*seems to work much better*/
		spinor2[1]    = 	sinarray[i];
		//spinor1[0]    = 	cos(dihedral[i]/2);
		//spinor1[1]    = 	sin(dihedral[i]/2);
		//spinor2[0]    = 	cos((M_PI-valence[i])/2);
		//spinor2[1]    = 	sin((M_PI-valence[i])/2);
		spinor3prod1(spinortot,spinor1,between);
		spinor3prod2(between,spinor2,spinortot);
		term1 	 = 2*bond[i];	/*add to previous atom coordinate*/
		co[i][0] = co[i-1][0] + term1*(spinortot[1]*spinortot[1] + spinortot[0]*spinortot[0] - 0.5);
		co[i][1] = co[i-1][1] + term1*(spinortot[1]*spinortot[2] + spinortot[0]*spinortot[3]);
		co[i][2] = co[i-1][2] + term1*(spinortot[1]*spinortot[3] - spinortot[0]*spinortot[2]);

		//		fprintf(stderr,"i= %2d  marker[i]= %d\n",i,marker[i]);
		// MOVING BRANCHED ATOMS (side-chain atoms)
		if (marker[i] == 2)		// sp3-hybridization (moving atoms pending from CA)
		{
			// first side chain (-H)
			spinor1[0]    	= 	cos((dihedral[i+1] + 2*M_PI/3)/2);	/*spinors are fixed -- can be precalculated*/ /*indices dihedral mapping Ok*/
			spinor1[1]    	= 	sin((dihedral[i+1] + 2*M_PI/3)/2);
			//			spinor2[0]    = 	cosarray[i+1];	/*seems to work much better*/
			//			spinor2[1]    = 	sinarray[i+1];
			spinor2[0]   	= 	cos((M_PI - valence[i+1])/2);		/*fixed valences angles ->no variable ones*/
			spinor2[1]    = 	sin((M_PI - valence[i+1])/2);		/*indices not yet correctly fixed in space*/
			spinor3prod1(spinortot,spinor1,between);
			spinor3prod2(between,spinor2,s1tot);
			term1 			= 	2*bondCH;
			p1co[i][0] 		= 	co[i][0] + term1*(s1tot[1]*s1tot[1] + s1tot[0]*s1tot[0] - 0.5);
			p1co[i][1] 		= 	co[i][1] + term1*(s1tot[1]*s1tot[2] + s1tot[0]*s1tot[3]);
			p1co[i][2] 		= 	co[i][2] + term1*(s1tot[1]*s1tot[3] - s1tot[0]*s1tot[2]);

			// second side chain (-R)
			spinor1[0]    	= 	cos((dihedral[i+1] + 4*M_PI/3)/2);					/*sp3-hybridisation*/
			spinor1[1]    	= 	sin((dihedral[i+1] + 4*M_PI/3)/2);
			spinor2[0]    = 	cosarray[i+1];	/*seems to work much better*/
			spinor2[1]    = 	sinarray[i+1];
			//			spinor2[0]    	= 	cos((M_PI - valence[i+1])/2);				/*fixed valences angles ->no variable ones*/
			//			spinor2[1]    	= 	sin((M_PI - valence[i+1])/2);				/*indices not yet correctly fixed in space*/
			spinor3prod1(spinortot,spinor1,between);
			spinor3prod2(between,spinor2,s2tot);
			term1			=   2*bondCH;				/*add to previous atom coordinate*/	//THIS IS ALKYL GROUP OF PEPTIDE!!!
			p2co[i][0] 		= 	co[i][0] + term1*(s2tot[1]*s2tot[1] + s2tot[0]*s2tot[0] - 0.5);
			p2co[i][1] 		= 	co[i][1] + term1*(s2tot[1]*s2tot[2] + s2tot[0]*s2tot[3]);
			p2co[i][2] 		= 	co[i][2] + term1*(s2tot[1]*s2tot[3] - s2tot[0]*s2tot[2]);
		}
		else if (marker[i] == 1)	// sp2-hybridization (moving atoms pending from either N or C)
		{
			if (i < number_atoms - 1)
			{
				spinor1[0]   	= 	cos((dihedral[i+1] + M_PI)/2);	/*spinors are fixed -- can be precalculated*/
				spinor1[1]    	= 	sin((dihedral[i+1] + M_PI)/2);
				spinor2[0]    = 	cosarray[i+1];	/*seems to work much better*/
				spinor2[1]    = 	sinarray[i+1];
				//spinor2[0]    	= 	cos((M_PI - valence[i+1])/2);		/*fixed valences angles ->no variable ones*/
				//spinor2[1]    	= 	sin((M_PI - valence[i+1])/2);		/*indices not yet correctly fixed in space*/
			}
			else // Last atom (Last dihedral is different --> dihedral[0])
			{
				spinor1[0]   	= 	cos((dihedral[0] + M_PI)/2);	/*spinors are fixed -- can be precalculated*/
				spinor1[1]    	= 	sin((dihedral[0] + M_PI)/2);
				spinor2[0]    = 	cosarray[i-2];	// Mon: taking the CA-C-O angle from previous aminoacid ("chapa")
				spinor2[1]    = 	sinarray[i-2];
				//spinor2[0]    	= 	cos((M_PI - valence[i-2])/2);		/*fixed valences angles ->no variable ones*/
				//spinor2[1]    	= 	sin((M_PI - valence[i-2])/2);		/*indices not yet correctly fixed in space*/
			}

			spinor3prod1(spinortot,spinor1,between);
			spinor3prod2(between,spinor2,s1tot);

			term1 		   = 	2*bondCO;
			p1co[i][0] 	   = 	co[i][0] + term1*(s1tot[1]*s1tot[1] + s1tot[0]*s1tot[0] - 0.5);
			p1co[i][1] 	   = 	co[i][1] + term1*(s1tot[1]*s1tot[2] + s1tot[0]*s1tot[3]);
			p1co[i][2] 	   = 	co[i][2] + term1*(s1tot[1]*s1tot[3] - s1tot[0]*s1tot[2]);

			p2co[i][0] = 0.0;	/*fill with zeros for definiteness*/ /*should not be printed!!!*/
			p2co[i][1] = 0.0;
			p2co[i][2] = 0.0;
		}
	}
}


void isafrotatetree(double **co,				/*loop coordinates*/
		double **p1co,				/*coordinates 1st parallel chain*/
		double **p2co,				/*coordinates 2nd parallel chain*/
		int chosenbond,
		double *posespinor,			/*pose spinor for anchor*/
		const int *bondnumbers,		/*physical bond numbers*/
		const int bondrange,		/*# of bonds rotating*/
		const int number_atoms,		/*# atoms in loop*/
		const double *valence,		/*valence angles loop*/
		const double *bond,			/*bond lengths loop*/
		double *cosarrayd,		/*preallocated -cosines*/
		double *sinarrayd,		/*preallocated sines*/
		const double *cosarray,		/*preallocated -cosines*/
		const double *sinarray,		/*preallocated sines*/
		double *dihedral,			/*dihedral angles loop*/
		double *p1dihedral,			/*dihedral angles loop*/
		double *p2dihedral,			/*dihedral angles loop*/
		int *marker,
		double *spinor1,			/*auxiliary spinor 1-preallocated*/
		double *spinor2,			/*auxiliary spinor 2-preallocated*/
		double *spinortot,			/*total spinor-preallocated*/
		double *s1tot,
		double *s2tot,
		double *between,
		double *full
)
{
	double term1;	/*intermediate computational term*/

	correctdihedral(dihedral + chosenbond + 1);

	/*Add hydrogen at first nitrogen*/
	//spinortot[0] = 1.0;
	//spinortot[1] = 0.0;
	//spinortot[2] = 0.0;
	//spinortot[3] = 0.0;

	//quatproduct(posespinor,spinortot,spinortot);

	spinor1[0]    = 	cos((dihedral[1] + M_PI)/2);
	spinor1[1]    = 	sin((dihedral[1]+ M_PI)/2);
	spinor2[0]    = 	cos((M_PI-valence[1])/2);
	spinor2[1]    = 	sin((M_PI-valence[1])/2);

	spinor3prod1(full,spinor1,between);
	spinor3prod2(between,spinor2,spinortot);

	//  add to previous atom coordinate
	term1 	 = 2*bondCH;
	p1co[0][0] = co[0][0] + term1*(spinortot[1]*spinortot[1] + spinortot[0]*spinortot[0] - 0.5);
	p1co[0][1] = co[0][1] + term1*(spinortot[1]*spinortot[2] + spinortot[0]*spinortot[3]);
	p1co[0][2] = co[0][2] + term1*(spinortot[1]*spinortot[3] - spinortot[0]*spinortot[2]);

	/*Actual start*/
	spinortot[0] = full[0];
	spinortot[1] = full[1];
	spinortot[2] = full[2];
	spinortot[3] = full[3];

	//quatproduct(posespinor,spinortot,spinortot);

	/*Extension for hydrogens and oxygens of loop*/
	for (int i = 1; i < number_atoms; i++)
	{
		spinor1[0]    = 	cos(dihedral[i]/2);
		spinor1[1]    = 	sin(dihedral[i]/2);
		spinor2[0]    = 	cos((M_PI-valence[i])/2);
		spinor2[1]    = 	sin((M_PI-valence[i])/2);

		spinor3prod1(spinortot,spinor1,between);
		spinor3prod2(between,spinor2,spinortot);

		//  add to previous atom coordinate
		term1 	 = 2*bond[i];
		co[i][0] = co[i-1][0] + term1*(spinortot[1]*spinortot[1] + spinortot[0]*spinortot[0] - 0.5);
		co[i][1] = co[i-1][1] + term1*(spinortot[1]*spinortot[2] + spinortot[0]*spinortot[3]);
		co[i][2] = co[i-1][2] + term1*(spinortot[1]*spinortot[3] - spinortot[0]*spinortot[2]);

		/*Entering tree structure*/
		if (marker[i] == 2)		/*sp3-hybridization or like*/ /*From first CO onward*/
		{
			/*first side chain*/
			spinor1[0]    	= 	cos((dihedral[i+1] + 2*M_PI/3)/2);	/*spinors are fixed -- can be precalculated*/ /*indices dihedral mapping Ok*/
			spinor1[1]    	= 	sin((dihedral[i+1] + 2*M_PI/3)/2);
			spinor2[0]   	= 	cos((M_PI - valence[i+1])/2);		/*fixed valences angles ->no variable ones*/
			spinor2[1]    = 	sin((M_PI - valence[i+1])/2);		/*indices not yet correctly fixed in space*/

			spinor3prod1(spinortot,spinor1,between);
			spinor3prod2(between,spinor2,s1tot);

			term1 			= 	2*bondCH;
			p1co[i][0] 		= 	co[i][0] + term1*(s1tot[1]*s1tot[1] + s1tot[0]*s1tot[0] - 0.5);
			p1co[i][1] 		= 	co[i][1] + term1*(s1tot[1]*s1tot[2] + s1tot[0]*s1tot[3]);
			p1co[i][2] 		= 	co[i][2] + term1*(s1tot[1]*s1tot[3] - s1tot[0]*s1tot[2]);


			/*second side chain*/
			spinor1[0]    	= 	cos((dihedral[i+1] + 4*M_PI/3)/2);					/*sp3-hybridisation*/
			spinor1[1]    	= 	sin((dihedral[i+1] + 4*M_PI/3)/2);
			spinor2[0]    	= 	cos((M_PI - valence[i+1])/2);				/*fixed valences angles ->no variable ones*/
			spinor2[1]    	= 	sin((M_PI - valence[i+1])/2);				/*indices not yet correctly fixed in space*/
			spinor3prod1(spinortot,spinor1,between);
			spinor3prod2(between,spinor2,s2tot);

			//  add to previous atom coordinate
			//term1 			= 	4*bondCH;
			term1			=   2*bondCH;		/*do some trick!*/
			p2co[i][0] 		= 	co[i][0] + term1*(s2tot[1]*s2tot[1] + s2tot[0]*s2tot[0] - 0.5);
			p2co[i][1] 		= 	co[i][1] + term1*(s2tot[1]*s2tot[2] + s2tot[0]*s2tot[3]);
			p2co[i][2] 		= 	co[i][2] + term1*(s2tot[1]*s2tot[3] - s2tot[0]*s2tot[2]);


			//			printf("Coordinates: %8.3f %8.3f %8.3f\n",p2co[i][0],p2co[i][1],p2co[i][2]);
			//			getchar();

		}
		else if (marker[i] == 1)	/*sp2-hybridization*/
		{
			if (i < number_atoms - 1)
			{
				spinor1[0]   	= 	cos((dihedral[i+1] + M_PI)/2);		/*spinors are fixed -- can be precalculated*/
				spinor1[1]    	= 	sin((dihedral[i+1] + M_PI)/2);
				spinor2[0]    	= 	cos((M_PI - valence[i+1])/2);		/*fixed valences angles ->no variable ones*/
				spinor2[1]    	= 	sin((M_PI - valence[i+1])/2);		/*indices not yet correctly fixed in space*/
			}
			else
			{
				/*but this one must yet be randomized*/
				spinor1[0]   	= 	cos((dihedral[i-2] + M_PI)/2);	/*spinors are fixed -- can be precalculated*/
				spinor1[1]    	= 	sin((dihedral[i-2] + M_PI)/2);
				spinor2[0]    	= 	cos((M_PI - valence[i-2])/2);		/*fixed valences angles ->no variable ones*/
				spinor2[1]    	= 	sin((M_PI - valence[i-2])/2);		/*indices not yet correctly fixed in space*/
			}

			spinor3prod1(spinortot,spinor1,between);
			spinor3prod2(between,spinor2,s1tot);

			term1 		   = 	2*bondCO;
			p1co[i][0] 	   = 	co[i][0] + term1*(s1tot[1]*s1tot[1] + s1tot[0]*s1tot[0] - 0.5);
			p1co[i][1] 	   = 	co[i][1] + term1*(s1tot[1]*s1tot[2] + s1tot[0]*s1tot[3]);
			p1co[i][2] 	   = 	co[i][2] + term1*(s1tot[1]*s1tot[3] - s1tot[0]*s1tot[2]);

			p2co[i][0] = 0.0;	/*fill with zeros for definiteness*/ /*should not be printed!!!*/
			p2co[i][1] = 0.0;
			p2co[i][2] = 0.0;
		}
	}

	/*Print all coordinates*/
	//		for (int j = 0; j < number_atoms; j++)
	//		{
	//			printf("Coordinates: %8.3f %8.3f %8.3f\n",co[j][0],co[j][1],co[j][2]);
	//		}
	//		printf(" \n");
	//
	//		for (int j = 0; j < number_atoms; j++)
	//		{
	//			printf("Coordinates: %8.3f %8.3f %8.3f\n",p1co[j][0],p1co[j][1],p1co[j][2]);
	//		}
	//		printf(" \n");
	//
	//
	//		for (int j = 0; j < number_atoms; j++)
	//		{
	//			printf("Coordinates: %8.3f %8.3f %8.3f\n",p2co[j][0],p2co[j][1],p2co[j][2]);
	//		}
	//getchar();
}

void ga_rotate_anglesm2tree(double **co,				/*loop coordinates*/
		double **p1co,				/*coordinates 1st parallel chain*/
		double **p2co,				/*coordinates 2nd parallel chain*/
		const double *angles,		/*angle variations from optimization*/
		double *posespinor,			/*pose spinor for anchor*/
		const int *bondnumbers,		/*physical bond numbers*/
		const int bondrange,		/*# of bonds rotating*/
		const int number_atoms,		/*# atoms in loop*/
		const double *valence,		/*valence angles loop*/
		const double *bond,			/*bond lengths loop*/
		double *dihedral,			/*dihedral angles loop*/
		double *p1dihedral,			/*dihedral angles loop*/
		double *p2dihedral,			/*dihedral angles loop*/
		int *marker,
		double *spinor1,			/*auxiliary spinor 1-preallocated*/
		double *spinor2,			/*auxiliary spinor 2-preallocated*/
		double *spinortot			/*total spinor-preallocated*/
)
{

	/*Auxiliary spinors*/
	double *between;
	double *s1tot;
	double *s2tot;
	between = (double*) malloc(sizeof(double)*4);
	s1tot   = (double*) malloc(sizeof(double)*4); /*in principle only size 2 integrated approach*/
	s2tot   = (double*) malloc(sizeof(double)*4);  /*in principle only size 2 integrated approach*/
	for(int i = 0; i < 4; i++)
	{
		between[i] = 0.0;
		s1tot[i]   = 0.0;
		s2tot[i]   = 0.0;
	}

	//double bondNH = 1.02;
	//double bondCH = 1.10;
	//double bondCO = 1.24;
	double term1;	/*intermediate computational term*/


	/*Add hydrogen at first nitrogen*/
	spinortot[0] = 1.0;
	spinortot[1] = 0.0;
	spinortot[2] = 0.0;
	spinortot[3] = 0.0;

	quatproduct(posespinor,spinortot,spinortot);

	spinor1[0]    = 	cos((dihedral[1] + M_PI)/2);
	spinor1[1]    = 	sin((dihedral[1]+ M_PI)/2);
	spinor2[0]    = 	cos((M_PI-valence[1])/2);
	spinor2[1]    = 	sin((M_PI-valence[1])/2);

	spinor3prod1(spinortot,spinor1,between);
	spinor3prod2(between,spinor2,spinortot);

	//  add to previous atom coordinate
	term1 	 = 2*bondCH;
	p1co[0][0] = co[0][0] + term1*(spinortot[1]*spinortot[1] + spinortot[0]*spinortot[0] - 0.5);
	p1co[0][1] = co[0][1] + term1*(spinortot[1]*spinortot[2] + spinortot[0]*spinortot[3]);
	p1co[0][2] = co[0][2] + term1*(spinortot[1]*spinortot[3] - spinortot[0]*spinortot[2]);

	/*Actual start*/
	spinortot[0] = 1.0;
	spinortot[1] = 0.0;
	spinortot[2] = 0.0;
	spinortot[3] = 0.0;

	quatproduct(posespinor,spinortot,spinortot);

	/*Extension for hydrogens and oxygens of loop*/
	for (int i = 1; i < number_atoms; i++)
	{
		spinor1[0]    = 	cos(dihedral[i]/2);
		spinor1[1]    = 	sin(dihedral[i]/2);
		spinor2[0]    = 	cos((M_PI-valence[i])/2);
		spinor2[1]    = 	sin((M_PI-valence[i])/2);

		spinor3prod1(spinortot,spinor1,between);
		spinor3prod2(between,spinor2,spinortot);

		//  add to previous atom coordinate
		term1 	 = 2*bond[i];
		co[i][0] = co[i-1][0] + term1*(spinortot[1]*spinortot[1] + spinortot[0]*spinortot[0] - 0.5);
		co[i][1] = co[i-1][1] + term1*(spinortot[1]*spinortot[2] + spinortot[0]*spinortot[3]);
		co[i][2] = co[i-1][2] + term1*(spinortot[1]*spinortot[3] - spinortot[0]*spinortot[2]);

		/*Entering tree structure*/
		if (marker[i] == 2)		/*sp3-hybridization or like*/ /*From first CO onward*/
		{
			/*first side chain*/
			spinor1[0]    	= 	cos((dihedral[i+1] + 2*M_PI/3)/2);	/*spinors are fixed -- can be precalculated*/ /*indices dihedral mapping Ok*/
			spinor1[1]    	= 	sin((dihedral[i+1] + 2*M_PI/3)/2);
			spinor2[0]   	= 	cos((M_PI - valence[i+1])/2);		/*fixed valences angles ->no variable ones*/
			spinor2[1]    	= 	sin((M_PI - valence[i+1])/2);		/*indices not yet correctly fixed in space*/

			spinor3prod1(spinortot,spinor1,between);
			spinor3prod2(between,spinor2,s1tot);

			term1 			= 	2*bondCH;
			p1co[i][0] 		= 	co[i][0] + term1*(s1tot[1]*s1tot[1] + s1tot[0]*s1tot[0] - 0.5);
			p1co[i][1] 		= 	co[i][1] + term1*(s1tot[1]*s1tot[2] + s1tot[0]*s1tot[3]);
			p1co[i][2] 		= 	co[i][2] + term1*(s1tot[1]*s1tot[3] - s1tot[0]*s1tot[2]);


			/*second side chain*/
			spinor1[0]    	= 	cos((dihedral[i+1] + 4*M_PI/3)/2);					/*sp3-hybridisation*/
			spinor1[1]    	= 	sin((dihedral[i+1] + 4*M_PI/3)/2);
			spinor2[0]    	= 	cos((M_PI - valence[i+1])/2);				/*fixed valences angles ->no variable ones*/
			spinor2[1]    	= 	sin((M_PI - valence[i+1])/2);				/*indices not yet correctly fixed in space*/
			spinor3prod1(spinortot,spinor1,between);
			spinor3prod2(between,spinor2,s2tot);

			//  add to previous atom coordinate
			//term1 			= 	2*bondCH;
			term1			=   4*bondCH;		/*do some trick!*/
			p2co[i][0] 		= 	co[i][0] + term1*(s2tot[1]*s2tot[1] + s2tot[0]*s2tot[0] - 0.5);
			p2co[i][1] 		= 	co[i][1] + term1*(s2tot[1]*s2tot[2] + s2tot[0]*s2tot[3]);
			p2co[i][2] 		= 	co[i][2] + term1*(s2tot[1]*s2tot[3] - s2tot[0]*s2tot[2]);


			//			printf("Coordinates: %8.3f %8.3f %8.3f\n",p2co[i][0],p2co[i][1],p2co[i][2]);
			//			getchar();

		}
		else if (marker[i] == 1)	/*sp2-hybridization*/
		{
			if (i < number_atoms - 1)
			{
				spinor1[0]   	= 	cos((dihedral[i+1] + M_PI)/2);	/*spinors are fixed -- can be precalculated*/
				spinor1[1]    	= 	sin((dihedral[i+1] + M_PI)/2);
				spinor2[0]    	= 	cos((M_PI - valence[i+1])/2);		/*fixed valences angles ->no variable ones*/
				spinor2[1]    	= 	sin((M_PI - valence[i+1])/2);		/*indices not yet correctly fixed in space*/
			}
			else
			{
				/*but this one must yet be randomized*/
				spinor1[0]   	= 	cos((dihedral[i-2] + M_PI)/2);	/*spinors are fixed -- can be precalculated*/
				spinor1[1]    	= 	sin((dihedral[i-2] + M_PI)/2);
				spinor2[0]    	= 	cos((M_PI - valence[i-2])/2);		/*fixed valences angles ->no variable ones*/
				spinor2[1]    	= 	sin((M_PI - valence[i-2])/2);		/*indices not yet correctly fixed in space*/
			}

			spinor3prod1(spinortot,spinor1,between);
			spinor3prod2(between,spinor2,s1tot);

			term1 		   = 	bondCO2;
			p1co[i][0] 	   = 	co[i][0] + term1*(s1tot[1]*s1tot[1] + s1tot[0]*s1tot[0] - 0.5);
			p1co[i][1] 	   = 	co[i][1] + term1*(s1tot[1]*s1tot[2] + s1tot[0]*s1tot[3]);
			p1co[i][2] 	   = 	co[i][2] + term1*(s1tot[1]*s1tot[3] - s1tot[0]*s1tot[2]);

			p2co[i][0] = 0.0;	/*fill with zeros for definiteness*/ /*should not be printed!!!*/
			p2co[i][1] = 0.0;
			p2co[i][2] = 0.0;
		}
	}

	free(between);
	free(s1tot);
	free(s2tot);

	/*Print all coordinates*/
	//		for (int j = 0; j < number_atoms; j++)
	//		{
	//			printf("Coordinates: %8.3f %8.3f %8.3f\n",co[j][0],co[j][1],co[j][2]);
	//		}
	//		printf(" \n");
	//
	//		for (int j = 0; j < number_atoms; j++)
	//		{
	//			printf("Coordinates: %8.3f %8.3f %8.3f\n",p1co[j][0],p1co[j][1],p1co[j][2]);
	//		}
	//		printf(" \n");
	//
	//
	//		for (int j = 0; j < number_atoms; j++)
	//		{
	//			printf("Coordinates: %8.3f %8.3f %8.3f\n",p2co[j][0],p2co[j][1],p2co[j][2]);
	//		}
	//getchar();
}



void spinorrotate(double **co,				/*loop coordinates*/
		const double *angles,		/*angle variations from optimization*/
		double *pspinor,			/*pose spinor for anchor*/
		const int *bondnumbers,		/*physical bond numbers*/
		const int bondrange,		/*# of bonds rotating*/
		const int number_atoms,		/*# atoms in loop*/
		const double *valence,		/*valence angles loop*/
		const double *bond,			/*bond lengths loop*/
		double *dihedralangle,			/*dihedral angles loop*/
		const double *cosarray,		/*preallocated -cosines*/
		const double *sinarray,		/*preallocated sines*/
		double *spinor1,			/*auxiliary spinor 1-preallocated*/
		double *spinor2,			/*auxiliary spinor 2-preallocated*/
		double *spinortot			/*total spinor-preallocated*/
)
{

	double term1;	/*intermediate computational term*/
	int j;										/*index dihedral correction step*/
	int switch_debug = 0;						/*switch for debug information*/
	double *between;
	between = (double*) malloc(sizeof(double)*4);
	for(int i = 0; i < 4; i++)
	{
		between[i] = 0.0;
	}

	message("Start allocation...",0);
	message("Allocation passed...",switch_debug);

	for (int i = 0; i < bondrange; i++){	/*Correction dihedral angles in [-PI,PI] range*/
		j = bondnumbers[i];
		dihedralangle[j + 1] += angles[i];
		correctdihedral(dihedralangle + j + 1);
	}


	//	spinortot[0] = 1.0;
	//	spinortot[1] = 0.0;
	//	spinortot[2] = 0.0;
	//	spinortot[3] = 0.0;
	//
	//	spinorproduct(posespinor,spinortot,spinortot);

	spinortot[0] = pspinor[0];
	spinortot[1] = pspinor[1];
	spinortot[2] = pspinor[2];
	spinortot[3] = pspinor[3];



	for (int i = 3; i < number_atoms; i++)
	{
		spinor1[0]    = 	cos(dihedralangle[i]/2);
		spinor1[1]    = 	sin(dihedralangle[i]/2);
		//spinor2[0]    = 	cosarray[i];	/*seems to work much better*/
		//spinor2[1]    = 	sinarray[i];
		spinor2[0]    = 	cos((M_PI - valence[i])/2);
		spinor2[1]    = 	sin((M_PI - valence[i])/2);

		spinor3prod1(spinortot,spinor1,between);	/*update spinortot*/
		spinor3prod2(between,spinor2,spinortot);

		term1 = 2*bond[i];	/*Rotate and add vector */
		co[i][0] = co[i-1][0] + term1*(spinortot[1]*spinortot[1] + spinortot[0]*spinortot[0] - 0.5);
		co[i][1] = co[i-1][1] + term1*(spinortot[1]*spinortot[2] + spinortot[0]*spinortot[3]);
		co[i][2] = co[i-1][2] + term1*(spinortot[1]*spinortot[3] - spinortot[0]*spinortot[2]);
	}
	free(between);
}


void garotate(double **co,				/*loop coordinates*/
		double **bond,
		double **bondn,
		const double *length,		/*bond length array*/
		const double *angles,		/*angle variations from optimization*/
		const int *marker,
		const int number_atoms,		/*# atoms in loop*/
		double *between,			/*auxiliary spinor*/
		double *spinor1,			/*auxiliary spinor 1-preallocated*/
		double *spinor2,			/*auxiliary spinor 2-preallocated*/
		double *spinortot			/*total spinor-preallocated*/
)
{
	for (int i = 3; i <  number_atoms; i++)	/*Correction dihedral angles in [-PI,PI] range*/
	{	if (i != 3 && marker[i-1] == 1 )
	{
		spinor(bondn[i-1],angles[i-1],spinor1);	/*still not tuned correctly*/
		spinorproduct(spinortot,spinor1,spinortot);
		revspinor(spinortot,spinor2);
	}
	else if (i == 3)
	{
		spinor(bondn[i-1],angles[i-1],spinortot);	/*still not tuned correctly*/
		revspinor(spinortot,spinor2);
	}
	geom32opt(bond[i],spinortot,between);
	geom23opt(spinor2,between,bond[i]);
	double norm;
	norm = 1/sqrt(bond[i][0]*bond[i][0] + bond[i][1]*bond[i][1] + bond[i][2]*bond[i][2]);
	bondn[i][0] = bond[i][0]*norm;
	bondn[i][1] = bond[i][1]*norm;
	bondn[i][2] = bond[i][2]*norm;
	bond[i][0]  = length[i]*bondn[i][0];
	bond[i][1]  = length[i]*bondn[i][1];
	bond[i][2]  = length[i]*bondn[i][2];
	//normalize(bond[i],bondn[i]);

	for(int j = 0; j < 3; j++)	/*can be put outside this loop -> maybe accuracy issue ...*/
	{
		co[i][j] = co[i-1][j] + bond[i][j]; 	/*bond[0] is not written*/
	}
	}
}



void ga_rotate_anglesm4(double **co,				/*loop coordinates*/
		const double *angles,		/*angle variations from optimization*/
		const int *bondnumbers,		/*physical bond numbers*/
		const int bondrange,		/*# of bonds rotating*/
		const int number_atoms,		/*# atoms in loop*/
		double *spinor1,			/*auxiliary spinor 1-preallocated*/
		double *spinor2,			/*auxiliary spinor 2-preallocated*/
		double *spinortot			/*total spinor-preallocated*/
)
{



	double term1;	/*intermediate computational term*/
	int j;										/*index dihedral correction step*/
	int switch_debug = 0;						/*switch for debug information*/
	double *between;
	between = (double*) malloc(sizeof(double)*4);
	for(int i = 0; i < 4; i++)
	{
		between[i] = 0.0;
	}



	message("Start allocation...",0);
	message("Allocation passed...",switch_debug);



	double **colocal;
	pointeralloc(&colocal,number_atoms + 2,3);

	for(int i = 0; i < number_atoms; i++)
	{	for(int j = 0; j < 3; j++)
	{
		colocal[i + 2][j] = co[i][j];
	}
	}
	colocal[1][0] = colocal[2][0] - 1.0;
	colocal[1][1] = colocal[2][1] - 0.0;
	colocal[1][2] = colocal[2][2] - 0.0;
	colocal[0][0] = colocal[1][0] - 1.0;	/*alternative approach also possible*/
	colocal[0][1] = colocal[1][1] - 1.0;
	colocal[0][2] = colocal[1][2] + 0.0;



	double *lengthset 		= new double[number_atoms + 2];		//	Array bond lengths (Angstrom) for internal coordinate set
	double *valenceset   	= new double[number_atoms + 2];		//	Array valence angles (rad) for internal coordinate set
	double *dihedralset  	= new double[number_atoms + 2];		//	Array dihedral angles for internal coordinate set
	findlength(colocal,lengthset,number_atoms + 2);
	valencefind(colocal,valenceset,number_atoms + 2);
	finddihedral(colocal,dihedralset,number_atoms + 2);


	//show3("Check loop coordinates","Coordinates",colocal,number_atoms + 2);
	//show3("Check obtained internal coordinate","Bond - Valence angle - Dihedral angle:",lengthset,valenceset,dihedralset,number_atoms + 2);




	for (int i = 0; i < bondrange; i++){	/*Correction dihedral angles in [-PI,PI] range*/
		j = bondnumbers[i];
		//printf("Bond number given %d\n",j);
		dihedralset[j + 1 + 2] += angles[i];		/*correct 2 higher than before*/
		correctdihedral(dihedralset + j + 1 + 2);	/*correct 2 higher than before*/
	}

	//show3("Check obtained internal coordinate","Bond - Valence angle - Dihedral angle:",lengthset,valenceset,dihedralset,number_atoms + 2);
	//getchar();

	spinortot[0] = 1.0;
	spinortot[1] = 0.0;
	spinortot[2] = 0.0;
	spinortot[3] = 0.0;

	for (int i = 3; i < number_atoms; i++)
	{
		spinor1[0]    = 	cos(dihedralset[i + 2]/2);
		spinor1[1]    = 	sin(dihedralset[i + 2]/2);
		spinor2[0]    = 	cos((M_PI - valenceset[i + 2])/2);
		spinor2[1]    = 	sin((M_PI - valenceset[i + 2])/2);

		spinor3prod1(spinortot,spinor1,between);	/*update spinortot*/
		spinor3prod2(between,spinor2,spinortot);

		term1 = 2*lengthset[i+2];	/*Rotate and add vector */
		colocal[i+2][0] = colocal[i+1][0] + term1*(spinortot[1]*spinortot[1] + spinortot[0]*spinortot[0] - 0.5);
		colocal[i+2][1] = colocal[i+1][1] + term1*(spinortot[1]*spinortot[2] + spinortot[0]*spinortot[3]);
		colocal[i+2][2] = colocal[i+1][2] + term1*(spinortot[1]*spinortot[3] - spinortot[0]*spinortot[2]);
	}



	for (int i = 2; i < number_atoms +  2; i++)
	{
		co[i-2][0] = colocal[i][0];
		co[i-2][1] = colocal[i][1];
		co[i-2][2] = colocal[i][2];
	}




	free(between);
	freepointer(colocal,number_atoms + 2);
}



void posespinor( double **co, int number_atoms, double *spinortot, double *spinor1,	double *spinor2)
{	/*auxiliary spinor 1-preallocated*/
	/*auxiliary spinor 2-preallocated*/
	double *between;
	between = (double*) malloc(sizeof(double)*4);
	for(int i = 0; i < 4; i++)
	{
		between[i] = 0.0;
	}

	double **colocal;
	pointeralloc(&colocal,number_atoms + 2,3);
	for(int i = 0; i < number_atoms; i++)
	{	for(int j = 0; j < 3; j++)
	{
		colocal[i + 2][j] = co[i][j];
	}
	}
	colocal[1][0] = colocal[2][0] - 1.0;
	colocal[1][1] = colocal[2][1] - 0.0;
	colocal[1][2] = colocal[2][2] - 0.0;
	colocal[0][0] = colocal[1][0] - 1.0;	/*alternative approach also possible*/
	colocal[0][1] = colocal[1][1] - 1.0;
	colocal[0][2] = colocal[1][2] + 0.0;

	double *lengthset 		= new double[number_atoms + 2];		//	Array bond lengths (Angstrom) for internal coordinate set
	double *valenceset   	= new double[number_atoms + 2];		//	Array valence angles (rad) for internal coordinate set
	double *dihedralset  	= new double[number_atoms + 2];		//	Array dihedral angles for internal coordinate set
	findlength(colocal,lengthset,number_atoms + 2);
	valencefind(colocal,valenceset,number_atoms + 2);
	finddihedral(colocal,dihedralset,number_atoms + 2);



	//show3("Check loop coordinates","Coordinates",colocal,number_atoms);
	//show3("Check obtained internal coordinate","Bond - Valence angle - Dihedral angle:",lengthset,valenceset,dihedralset,number_atoms);
	//	for (int i = 0; i < number_atoms; i++)
	//	{
	//		printf("Obtained dihedrals in posespinor %3i %8.3f\n",i,dihedralset[i]);
	//	}
	//	getchar();

	spinortot[0] = 1.0;
	spinortot[1] = 0.0;
	spinortot[2] = 0.0;
	spinortot[3] = 0.0;

	spinor1[0]    = 	cos(dihedralset[3]/2);
	spinor1[1]    = 	sin(dihedralset[3]/2);
	spinor2[0]    = 	cos((M_PI - valenceset[3])/2);
	spinor2[1]    = 	sin((M_PI - valenceset[3])/2);

	spinor3prod1(spinortot,spinor1,between);	/*update spinortot*/
	spinor3prod2(between,spinor2,spinortot);

	spinor1[0]    = 	cos(dihedralset[4]/2);
	spinor1[1]    = 	sin(dihedralset[4]/2);
	spinor2[0]    = 	cos((M_PI - valenceset[4])/2);
	spinor2[1]    = 	sin((M_PI - valenceset[4])/2);

	spinor3prod1(spinortot,spinor1,between);	/*update spinortot*/
	spinor3prod2(between,spinor2,spinortot);

	freepointer(colocal,number_atoms + 2);
}

void startspinor(double **co, double *posespinor)
{
	double *spinor1;
	double *spinor2;
	double *between;
	spinor1 = (double*) malloc(sizeof(double)*2);
	spinor2 = (double*) malloc(sizeof(double)*2);
	between = (double*) malloc(sizeof(double)*4);
	for(int i = 0; i < 2; i++)
	{
		spinor1[i] = 0.0;
		spinor2[i] = 0.0;
	}
	for(int i = 0; i < 2; i++)
	{
		between[i] = 0.0;
	}


	double *xaxis;
	double *yaxis;
	double *zaxis;
	xaxis = (double*) malloc(sizeof(double)*3);
	yaxis = (double*) malloc(sizeof(double)*3);
	zaxis = (double*) malloc(sizeof(double)*3);
	for(int i = 0; i < 3; i++)
	{
		xaxis[i] = 0.0;
		yaxis[i] = 0.0;
		zaxis[i] = 0.0;
	}
	xaxis[0]	= 	1.0;
	yaxis[1]	=	1.0;
	zaxis[2]	=	1.0;

	double *bond2;
	bond2 = (double*) malloc(sizeof(double)*3);
	for(int i = 0; i < 3; i++)
	{
		bond2[i] = co[2][i] - co[1][i];
	}

	double valangle;
	double dihangle;

	valangle = acos(dotprod(xaxis,bond2)/(norm(xaxis)*norm(bond2)));
	if (valangle > M_PI) {
		valangle = 2*M_PI - valangle;
	}
	if (valangle > M_PI/2) {
		valangle = M_PI - valangle;
	}

	bond2[0] = 0.0;
	dihangle = acos(dotprod(zaxis,bond2)/(norm(zaxis)*norm(bond2)));
	if (dihangle > M_PI) {
		dihangle = 2*M_PI - dihangle;
	}
	if (dotprod(yaxis,bond2) > 0) {
		dihangle = M_PI - dihangle;
	} else
	{
		dihangle *= -1;
	}

	/*creating start spinor*/
	posespinor[0] = 1.0;
	posespinor[1] = 0.0;
	posespinor[2] = 0.0;
	posespinor[3] = 0.0;

	spinor1[0]    = 	cos(dihangle/2);
	spinor1[1]    = 	sin(dihangle/2);
	spinor2[0]    = 	cos(valangle/2);	/*seems to work much better*/
	spinor2[1]    = 	sin(valangle/2);

	spinor3prod1(posespinor,spinor1,between);	/*update spinortot*/
	spinor3prod2(posespinor,spinor2,posespinor);

	free(spinor1);
	free(spinor2);
	free(between);
}

void rotate_bw2(double **co,double *bonds,double *cos_valence,double *sin_valence,double *cos_dihedral,double *sin_dihedral, int size){

	// Local spinors
	double spinor1[2];
	double spinor2[2];

	double spinortot[4];
	double between[4];

	double term1;

	//	Specify first atom
	//	co[0][0] = 0.0;						//	not needed	always first coordinate
	//	co[0][1] = 0.0;						//	not needed
	//	co[0][2] = 0.0;						//	not needed

	// initialisation start spinor	pass as argument
	spinortot[0] = 1.0;
	spinortot[1] = 0.0;
	spinortot[2] = 0.0;
	spinortot[3] = 0.0;

	for (int i = 1; i < size  ; i++) {			//	adapted
		//	Implementation: more operations but less writing to memory----hence little faster
		spinor1[0]    = 	cos_dihedral[i];
		spinor1[1]    = 	sin_dihedral[i];
		spinor2[0]    =    	cos_valence[i];		//	cos(valence[i-1]/2);	//	later pass on as direct input
		spinor2[1]    = 	sin_valence[i]; 		//	sin(valence[i-1]/2);	//	later pass on as direct input

		spinor3prod1(spinortot,spinor1,between);
		spinor3prod2(between,spinor2,spinortot);

		//  add to previous atom coordinate
		if (i > 3) {								//	but i = 3 must be used to do calculation ... starter
			term1 = 2*bonds[i];
			co[i][0] = co[i-1][0] + term1*(spinortot[1]*spinortot[1] + spinortot[0]*spinortot[0] - 0.5);	//	index atom adapted RH side
			co[i][1] = co[i-1][1] + term1*(spinortot[1]*spinortot[2] + spinortot[0]*spinortot[3]);
			co[i][2] = co[i-1][2] + term1*(spinortot[1]*spinortot[3] - spinortot[0]*spinortot[2]);
		}
	}

}


/*===============================================================================
 *
 *  BESTANGLE.CPP
 *
 *	Created on: Dec 13, 2010
 *      Author: Pieter Chys
 *
 *===============================================================================
 * BESTANGLE.CPP finds the best angle to rotate for the polymer
 * end point towards a chosen fixed end point. It is a optimization
 * algorithm in the plane that is used and makes use of matrix
 * algebra to do the job.  The coordinate frame is varied
 * as to make only an optimization in Y and Z coordinates necessary.
 * The rotation before analysis makes the x planes fixed so that
 * it can be excluded in the analysis to make.
 *
 *--------------------------------------------------------------------------------
 *
 *  INPUT: - int referring to the number of atoms in the chain
 *  	   - int referring to the bond number
 *  	   - double** pointer to the coordinates of the whole chain
 *  	   - double array[3] which has the fixed end point coordinates
 *
 *-------------------------------------------------------------------------------
 *
 *  OUTPUT: - double is returned: the value of rotation angle to
 *            go through
 *
 *
 *-------------------------------------------------------------------------------
 *
 *  DEPENDENCIES:
 *
 *===============================================================================
 */

void startmatrix(	double **co,		/*loop coordinates*/
		double **matrix,	/*output matrix*/
		double **matrixrot)	/*temporary matrix-preallocated*/
{
	rotatematrix2(co + 1,matrixrot);
	transpose(matrixrot,matrix);
}

void startspinor(	double **co,		/*loop coordinates*/
		double **matrix,	/*output matrix*/
		double **matrixrot,	/*temporary matrix-preallocated*/
		double *spinor)	/*still implement if trace is wrong*/
{
	rotatematrix2(co + 1,matrixrot);
	transpose(matrixrot,matrix);

	double w;
	double w4;
	double x;
	double y;
	double z;

	w = sqrt(1 + matrix[0][0] + matrix[1][1] + matrix[2][2])/2;
	w4 = 4.0 * w;
	x = (matrix[2][1] - matrix[1][2]) / w4;
	y = (matrix[0][2] - matrix[2][0]) / w4;
	z = (matrix[1][0] - matrix[0][1]) / w4;

	double tracesum;							// is correct
	tracesum = matrix[0][0] + matrix[1][1] + matrix[2][2] + 1;
	printf("Trace sum %8.3f\n",tracesum);

	//double length;							// is correct
	//length = sqrt(x*x + y*y + z*z + w*w);
	//printf("Length quaternion %8.3f\n",length);

	spinor[0] =   w;
	spinor[1] = + z;
	spinor[2] = + x;
	spinor[3] = - y;	/*or either with different sign*/
}


double bestangle(double** co,int bond_number,int number_atoms,double* point)
{
	double polpoint[3],endpoint[3],bond2[3];
	double** matrix;

	//	Initialization of rotation matrix
	pointeralloc(&matrix,3,3);

	//	Actual algorithm
	subtract(co[number_atoms-1],co[bond_number],polpoint);
	subtract(point,co[bond_number],endpoint);
	diffvectorseq(co,bond_number+1,bond2);						//	Modified at 27/01/2011


	rotatematrix(bond2,matrix); 			//	But mind bond1 ABOVE (not dealt with) -> bond_number 1 should be avoided
	double result1[3],result2[3];
	double teta;
	matrixXvector(matrix,polpoint,result1);
	matrixXvector(matrix,endpoint,result2);

	//	Free allocated matrix
	freepointer(matrix,3);

	teta = acos(result1[1]/(sqrt((result1[1]*result1[1]+result1[2]*result1[2]))));
	if (result1[2] < 0) {
		teta = -teta;
	}

	return atan2(result2[2],result2[1])-teta;
}

/*===============================================================================
 *
 *  BESTANGLE_ANCHOR.CPP
 *
 *	Created on: Jan 27, 2010
 *      Author: Pieter Chys
 *
 *===============================================================================
 * BESTANGLE_ANCHOR.CPP finds the best angle to rotate for the polymer
 * end point towards a chosen fixed anchor. An anchor is a set of endpoints
 * which are fixed 3-dimensional coordinates. It is similar optimization
 * algorithm as BESTANGLE.CPP in the plane but extended to three points.
 * The coordinate frame is first varied as in BESTANGLE.CPP but then the
 * computational procedure is different to account for three points, albeit
 * similar in concept. The reference from literature is :
 *
 * Canutescu, A.A. and Dunback Jr., R.L., Cyclic coordinate descent: a
 * robotics algorithm for protein loop closure, Protein Science (2003),
 * 12:963-972
 *
 * The algorithm is slightly changed and the formula -proposed to be
 * not used- we use here with some tricks to reduce computations.
 *
 *--------------------------------------------------------------------------------
 *
 *  INPUT: - int referring to the number of atoms in the chain
 *  	   - int referring to the bond number
 *  	   - double** pointer to the coordinates of the whole chain
 *  	   - double** pointer which has the 3 fixed end point coordinates
 *  		 of the anchor
 *
 *-------------------------------------------------------------------------------
 *
 *  OUTPUT: - double is returned: the value of rotation angle to
 *            go through
 *
 *
 *-------------------------------------------------------------------------------
 *
 *  DEPENDENCIES:
 *
 *===============================================================================
 */

double bestangle_anchor(double** co,int bond_number,int number_atoms,double** anchor) // 	Code seems to work but --sometimes wrong angle
{

	double** matrix;																			// 	rotation matrix used
	double** polpoint;																		// 	coordinates of the loop end atoms
	double** anchorpoint;																	// 	coordinates of the fixed anchor to go
	double** result1;																		//  rotated loop end coordinates
	double** result2;																		//	rotated anchor coordinates
	double bond2[3];																			//	bond around which we rotate
	double** result3;

	//printf("Debugging mark 1 ...\n");

	//	Initialization of rotation matrix
	pointeralloc(&matrix,3,3);
	pointeralloc(&polpoint,3,3);
	pointeralloc(&anchorpoint,3,3);
	pointeralloc(&result1,3,3);
	pointeralloc(&result2,3,3);
	pointeralloc(&result3,3,3);

	//printf("%s %d\n","Bond_number+1 value used",bond_number-1);	//	bond_number is fixed -> problem!!!!



	//	Actual algorithm
	for (int i = 0; i < 3; i++) {																//	make correct active anchor coordinates
		subtract(co[bond_number-1],co[number_atoms - 3 + i],polpoint[i]);								//	crucial_____XXXXX OLD CODE
		//	subtract(co[bond_number],co[number_atoms-3+i],polpoint[i]);
		// segmentation fault can happen
	}
	//printf("%-30s %8.3f %8.3f %8.3f\n","polymer anchor coordinates",polpoint[1][0], polpoint[1][1],polpoint[2][2]);



	for (int i = 0; i < 3; i++) {																//	make correct fixed anchor coordinates
		//	subtract(co[bond_number-1],anchor[i],anchorpoint[i]);									//	in translated frame ... OLD ONE
		//printf("Bond number - 1 %i\n",bond_number-1);
		subtract(co[bond_number-1],anchor[i],anchorpoint[i]);									//	in translated frame ... OLD ONE
	}
	//	MISTAKE IS HERE IN SUBTRACT

	//printf("Debugging mark 2 ...\n");
	//printf("%s %d\n","Bond_number+1 value used",bond_number-1);

	diffvectorseq(co,bond_number+1,bond2);

	//printf("%-30s %8.3f %8.3f %8.3f\n","bond 2 co",bond2[0], bond2[1],bond2[2]);



	rotatematrix(bond2,matrix); 										//	But mind bond1 ABOVE (not dealt with) -> bond_number 1 should be avoided
	for (int i = 0; i < 3; i++) {										//	apply to each of the anchor points the previous procedure
		matrixXvector(matrix,polpoint[i],result1[i]);				//	problem here somewhere....
		matrixXvector(matrix,anchorpoint[i],result2[i]);
	}

	//	Free allocated matrix
	freepointer(matrix,3);


	//	MATRIX ALGEBRA APPROACH

	for (int i = 0; i < 3; i++) {
		//double a;
		result2[i][0] = 0.0;												//	setting x coordinate at zero
		result1[i][0] = 0.0;
	}

	double denominator = 0;
	for (int i = 0; i < 3; i++) {
		denominator += dotprod(result1[i],result2[i]);
	}


	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			result3[i][0] = 0.0;
			result3[i][1] = -result1[i][2];
			result3[i][2] = result1[i][1];
		}
	}

	double numerator = 0;
	for (int i = 0; i < 3; i++) {
		numerator += dotprod(result2[i],result3[i]);
	}

	freepointer(polpoint,3);
	freepointer(anchorpoint,3);
	freepointer(result1,3);
	freepointer(result2,3);
	freepointer(result3,3);

	return atan2(numerator,denominator);

}


/*===============================================================================
 *
 *  BESTANGLE_ANCHOR.CPP
 *
 *	Created on: 1 September, 2011
 *      Author: Pieter Chys
 *
 *===============================================================================
 * BESTANGLE_ANCHOR_NEW.CPP finds the best angle to for the polymer end anchor
 * (of the last three points) to rotate towards a chosen fixed anchor. The optimization
 * minimizes the distance between the two sets of points.  The fixed anchor is the set
 * of three end points which are fixed 3-dimensional coordinates. It is similar optimization
 * algorithm as BESTANGLE.CPP (with planar optimization) but extended for three points.
 * The coordinate frame is first changed as in BESTANGLE.CPP but hereafter the
 * computational procedure is different to account for the three points, albeit
 * similar in concept. The algorithm has many **pointer input  arguments as to
 * increase speed of execution. The function BESTANGLE_ANCHOR_NEW.CPP is a direct
 * descendant which additionally implements optimization for the last-minus-one
 * bond to be rotated.
 *
 *
 * The reference from literature that was initially the goal is:
 *
 * Canutescu, A.A. and Dunback Jr., R.L., Cyclic coordinate descent: a
 * robotics algorithm for protein loop closure, Protein Science (2003),
 * 12:963-972
 *
 * However, the final algorithm here is radically different but the
 * derivation in this reference is anyway interesting to learn from.
 *
 *--------------------------------------------------------------------------------
 *
 *  INPUT: - int referring to the number of atoms in the chain
 *  	   - int referring to the bond number -> switch for which path is chosen
 *  	   - double** pointer to the coordinates of the whole chain
 *  	   - double** pointer which has the 3 fixed end point coordinates
 *  		 of the anchor
 *  	   - double** pointers to POLPOINT, ANCHORPOINT specifying the loop and
 *  	     anchor coordinates
 *         - double** pointers to result1, result2 and result3, used in the
 *           minimization procedure
 *
 *
 *-------------------------------------------------------------------------------
 *
 *  OUTPUT: - double is returned: the value of rotation angle to
 *            go through
 *
 *
 *-------------------------------------------------------------------------------
 *
 *  DEPENDENCIES: - many functions from my_library
 *                - the function is the predecessor of BESTANGLE_ANCHOR_NEW.CPP
 *
 *===============================================================================
 */


double bestangle_anchor(double** co,int bond_number,int number_atoms,double** anchor,						//	Preallocated version
		double** matrix,double** polpoint,double** anchorpoint,double** result1,double** result2, double** result3) {
	double bond2[3];																			//	bond around which we rotate

	//	Actual algorithm
	for (int i = 0; i < 3; i++) {																//	make correct active anchor coordinates
		subtract(co[bond_number-1],co[number_atoms-3+i],polpoint[i]);							//	crucial_____XXXXX OLD CODE
		// segmentation fault can happen
	}

	for (int i = 0; i < 3; i++) {																//	make correct fixed anchor coordinates
		subtract(co[bond_number-1],anchor[i],anchorpoint[i]);									//	in translated frame ... OLD ONE
	}

	diffvectorseq(co,bond_number+1,bond2);

	rotatematrix(bond2,matrix);

	rotatematrix(bond2,matrix); 										//	But mind bond1 ABOVE (not dealt with) -> bond_number 1 should be avoided
	for (int i = 0; i < 3; i++) {										//	apply to each of the anchor points the previous procedure
		matrixXvector(matrix,polpoint[i],result1[i]);				//	problem here somewhere....
		matrixXvector(matrix,anchorpoint[i],result2[i]);
	}

	//	MATRIX ALGEBRA APPROACH
	for (int i = 0; i < 3; i++) {
		result2[i][0] = 0.0;												//	setting x coordinate at zero
		result1[i][0] = 0.0;
	}

	double denominator = 0;
	for (int i = 0; i < 3; i++) {
		denominator += dotprod(result1[i],result2[i]);
	}


	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			result3[i][0] = 0.0;
			result3[i][1] = -result1[i][2];
			result3[i][2] = result1[i][1];
		}
	}

	double numerator = 0;
	for (int i = 0; i < 3; i++) {
		numerator += dotprod(result2[i],result3[i]);
	}

	return atan2(numerator,denominator);

}





void intersec(double *a,double *b,double c[3],double d[3])
{
	double t;
	double f[3];
	int i;

	for (i = 0; i < 3; i++)
	{
		f[i] = b[i]-a[i];
	}

	t = ((c[0]-a[0])*f[0] + (c[1]-a[1])*f[1] + (c[2]-a[2])*f[2])/(f[0]*f[0]+f[1]*f[1]+f[2]*f[2]);

	for (i = 0; i < 3; i++)
	{
		d[i] = a[i] + (b[i]-a[i])*t;
	}
}


void orient(double v0[],double a[],double b[],double v[])
{
	int i;
	double v1[3];
	double v2[3];
	double f;

	for (i = 0; i < 3; i++)
	{
		v1[i] = a[i];
		v2[i] = b[i];
	}


	for (i = 0; i < 3; i++)
	{
		v1[i] = v1[i] - v0[i];
		v2[i] = v2[i] - v0[i];
	}
	v[0] =  v1[1]*v2[2] -v1[2]*v2[1];
	v[1] = -v1[0]*v2[2] +v1[2]*v2[0];
	v[2] =  v1[0]*v2[1] -v1[1]*v2[0];


	f = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);

	for (i = 0 ;i < 3; i++)
	{
		v[i] /= f;
	}

}


double bestangle_ana(double **co, int bondnumber,double **co_i, int nr_atoms)
{
	double x;
	double y;
	double dr;
	double r[3];
	double ro[3];
	double fo[3];
	double o[3];
	double ab[3];
	double s[3];
	double angle;

	x = 0.0;
	y = 0.0;

	for (int m = 0; m < 3 ;m++)
	{
		for (int l = 0; l < 3 ; l++)
		{
			//r[l]=residue[back-1].main[m+1][l];
			r[l] = co[nr_atoms - 3 + m][l];		//actual anchor coordinates
		}

		intersec(co[bondnumber],co[bondnumber-1],r,o); 		// a and b are successive coordinates in main chain
		// a,b,c are input --- o is output
		// o is a transformed coordinate from r


		for (int l = 0; l < 3; l++)
		{
			ro[l]   = r[l]        - o[l];
			//fo[l] = target[m][l]- o[l];
			fo[l]	= co_i[nr_atoms - 3 +m][l] - o[l];	// target coordinates minus created o-coordinates
		}

		dr = sqrt(ro[0]*ro[0]+ro[1]*ro[1]+ro[2]*ro[2]);

		if (dr < TINY)
			continue;

		for (int l = 0; l < 3; l++)
		{
			ab[l] = o[l] + co[bondnumber-1][l] - co[bondnumber][l];
		}

		orient(o,ab,r,s);

		y	+= dr*dotprod(fo,s);
		x	+= dotprod(fo,ro);
	}

	angle 	= -atan2(y,x); 	// optimized angle ...
	//angle 	*= RADIAN;

	return angle;

}



/*===============================================================================
 *
 *  BESTANGLE_ANCHOR_NEW.CPP
 *
 *	Created on: 1 September, 2011
 *      Author: Pieter Chys
 *
 *===============================================================================
 * BESTANGLE_ANCHOR_NEW.CPP finds the best angle to for the polymer end anchor
 * (of the last three points) to rotate towards a chosen fixed anchor. The optimization
 * minimizes the distance between the two sets of points.  The fixed anchor is the set
 * of three end points which are fixed 3-dimensional coordinates. It is similar optimization
 * algorithm as BESTANGLE.CPP (with planar optimization) but extended for three points.
 * The coordinate frame is first changed as in BESTANGLE.CPP but hereafter the
 * computational procedure is different to account for the three points, albeit
 * similar in concept. This function is the same as BESTANGLE_ANCHOR but contains
 * an additional procedure to include motion of the last point of the protein loop
 * by rotation around the last-minus-one bond. Since this is only optimization
 * for 1 point an alternative procedure was implemented for this. Thus two paths
 * exists in this version depending on which rotation bond we are dealing with.
 * If the rotation bond is the last-minus-one the standard optimization is switched
 * for the alternative minimization. The algorithm has many **pointer input
 * arguments as to increase speed of execution.
 *
 *
 * The reference from literature that was initially the goal is:
 *
 * Canutescu, A.A. and Dunback Jr., R.L., Cyclic coordinate descent: a
 * robotics algorithm for protein loop closure, Protein Science (2003),
 * 12:963-972
 *
 * However, the final algorithm here is radically different but the
 * derivation in this reference is anyway interesting to learn from.
 *
 *--------------------------------------------------------------------------------
 *
 *  INPUT: - int referring to the number of atoms in the chain
 *  	   - int referring to the bond number -> switch for which path is chosen
 *  	   - double** pointer to the coordinates of the whole chain
 *  	   - double** pointer which has the 3 fixed end point coordinates
 *  		 of the anchor
 *  	   - double** pointers to POLPOINT, ANCHORPOINT specifying the loop and
 *  	     anchor coordinates
 *         - double** pointers to result1, result2 and result3, used in the
 *           minimization procedure
 *
 *
 *-------------------------------------------------------------------------------
 *
 *  OUTPUT: - double is returned: the value of rotation angle to
 *            go through
 *
 *
 *-------------------------------------------------------------------------------
 *
 *  DEPENDENCIES: - many functions from my_library
 *                - the function is an immediate descendant of BEST_ANGLE_ANCHOR.CPP
 *
 *===============================================================================
 */


double bestangle_anchor_new(double** co,int bond_number,int number_atoms,double** anchor,						//	Preallocated version
		double** matrix,double** polpoint,double** anchorpoint,double** result1,double** result2, double** result3)
{

	double angle;	 //	optimal angle we want
	double bond2[3]; // 	bond around which we rotate
	double denominator;
	double numerator;

	for (int i = 0; i < 3; i++)
	{		//	make correct active anchor coordinates
		subtract(co[bond_number-1],co[number_atoms-3+i],polpoint[i]);
	}

	for (int i = 0; i < 3; i++)
	{		//	make correct fixed anchor coordinates
		subtract(co[bond_number-1],anchor[i],anchorpoint[i]);	 //	in translated frame ... OLD ONE
	}

	diffvectorseq(co,bond_number + 1,bond2);

	rotatematrix(bond2,matrix);


	for (int i = 0; i < 3; i++)
	{
		matrixXvector(matrix,polpoint[i],result1[i]);
		matrixXvector(matrix,anchorpoint[i],result2[i]);
	}

	//	MATRIX ALGEBRA APPROACH
	for (int i = 0; i < 3; i++)
	{
		result2[i][0] = 0.0;		//	setting x coordinate at zero
		result1[i][0] = 0.0;
	}


	denominator = 0;
	for (int i = 0; i < 3; i++)
	{
		denominator += dotprod(result1[i],result2[i]);
	}


	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			result3[i][0] 	= 0.0;	//	Construct orthonormal with cross product of result1 and [1 0 0] and
			result3[i][1] 	= -result1[i][2];	//	work out analytically
			result3[i][2] 	= result1[i][1];
		}
	}

	numerator = 0;
	for (int i = 0; i < 3; i++)
	{
		numerator += dotprod(result2[i],result3[i]);
	}

	angle = atan2(numerator,denominator);



	//		if(angle > M_PI || angle < -M_PI || isnan(angle))
	if( isnan(angle) )
	{
		fprintf(stderr,"\n\nbest angle= %f  num= %f   den= %f bond %d number_atoms %d\n",angle,numerator,denominator, bond_number, number_atoms);
		fprintf(stderr," polpoint %f %f %f\n",polpoint[0],polpoint[1],polpoint[2]);
		fprintf(stderr," anchorpoint %f %f %f\n",anchorpoint[0],anchorpoint[1],anchorpoint[2]);
		fprintf(stderr," bond2 %f %f %f\n\n",bond2[0],bond2[1],bond2[2]);
		exit(1);
	}

	//		if (isnan(angle))
	//			angle=0.0;

	//		else if(angle > M_PI)  angle=M_PI;
	//		 else if (angle < -M_PI ) angle=-M_PI;

	return angle;
}


void bestangles(double** co,double *angles,int *bondsequence,int bondrange,int number_atoms,double** anchor,						//	Preallocated version
		double** matrix,double** polpoint,double** anchorpoint,double** result1,double** result2, double** result3) {

	double angle;							//	optimal angle we want

	//	THIS BRANCH HAS NOT YET BEEN MODIFIED
	/*if (bond_number == number_atoms - 2){
		double t;								//	normal increment == length normal
		double radius;							//	radius last bond rotation in circle plane
		double distance;						//	distance  projection circle plane and point on circle
		double length3;							//	length last bond
		double angle3;							//	valence angle last bond
		double dotpr;							//	dot product from angle checking

		double cotogo[3];						//	ideal point to go
		double cross[3];						//	normal vector to check dihedral angle sign
		double normal[3];						//	normal vector to circle plane
		double fixed[3];						//	last point of anchor
		double tocenter[3];						//	vector that goes to center
		double center[3];						//	center of circle last bond
		double bond2[3];						//	bond to rotate around
		double bond3[3];						//  bond that will be rotated
		double projection[3];					//  point projected on circle plane (fixed)
		double difference[3];					//	difference projected and circle point
		double x[3];							//  coordinate optimal point on circle
		double binormal[3];						//	coordinate circle coordinate actual bond
		double localsum;						//	temporary variable

		//	ACTIVE ANCHOR COORDINATES
		diffvectorseq(co,bond_number + 1, bond2);			// indices fine	// first bond active anchor
		diffvectorseq(co,bond_number + 2, bond3);			// indices fine	// second bond active anchor
		normalize(bond2,normal);

		length3 	= 	norm(bond3);
		angle3 		= 	acos(dotprod(bond2,bond3)/(norm(bond2)*norm(bond3)));	//	correct angle ... smaller/bigger than 90 degrees	// valence angle between bond 1 and 2 in anchor --- POSSIBLE BUG???
		radius 		=	sin(angle3)*length3;
		// DEBUGGING FUNCTION
		//printf("%-50s %8.3f\n","Valence angle in anchor is",angle3);			//	seems Ok
		//printf("-50s %8.3f\n", "Radius in anchor is",radius);					//	seems Ok
		//fflush(stdout);

		t 	= 0.0;
		localsum = 0.0;
		for(int i = 0; i < 3; i++){
			// CAN BE SPEED UP LATER IN SINGLE LINES
			tocenter[i]	=	length3*cos(angle3)*normal[i];
			center[i] 	=	co[bond_number][i] + tocenter[i];		//	changes // seems Ok until here ...
			binormal[i] =   co[bond_number + 1][i] - center[i];		//	changes

			fixed[i]  	= 	anchor[2][i];							//	seems Ok ...
			//t     		+=	normal[i]*(center[i] - fixed[i]);
			t     		+=	normal[i]*(fixed[i]-center[i]);			//	BUG: SHOULD BE REVERSED
			localsum	+= normal[i];								//	NEEDS TO BE ADDED TO LET IT WORK!!!!! BUG BEFORE ....
		}
		t				/= localsum;								//	NEEDS TO BE ADDED TO LET IT WORK!!!!! BUG BEFORE ...
		//printf("%-50s %8.3f\n","Parameter t",t);
		//fflush(stdout);


		for(int i = 0; i < 3; i++){
			//projection[i] = fixed[i] + t*normal[i];
			projection[i]	= center[i] + t*normal[i];				//	CORRECTED!!!!!! BUG BEFORE....
			//difference[i] = projection[i] - center[i];
			difference[i] = fixed[i] - projection[i];
		}

		distance = norm(difference);

		for(int i = 0; i < 3; i++){
			x[i] 	  =  radius*difference[i]/distance;	// projection vector in A parallel plane !!!!! DEBUGGING ...
			//cotogo[i] =  center[i]	+  x[i];
		}

		//	DEBUGGING
		printf("\n");
		printf("%50s %8.3f %8.3f %8.3f\n","Original end coordinate",co[number_atoms - 1][0],co[number_atoms - 1][1],co[number_atoms - 1][2]);
		printf("%50s %8.3f %8.3f %8.3f\n","Ideal end coordinate",cotogo[0],cotogo[1],cotogo[2]);
		printf("%50s %8.3f %8.3f %8.3f\n","Actual radius vector",binormal[0],binormal[1],binormal[2]);
		//printf("%50s %8.3f %8.3f %8.3f\n","Ideal radius vector",x[0],x[1],x[2]);
		//fflush(stdout);


		crossprod(bond3,bond2,cross);								//	calculating cross product
		dotpr = dotprod(cross,x);									//  check value dot product

		angle = acos(dotprod(x,binormal)/(norm(x)*norm(binormal)));	//	how about correction for angle????

		if (dotpr > 0){
			angle = - fabs(angle);
		} else if (dotpr < 0) {
			angle = + fabs(angle);
		}


		//	DEBUGGING
		printf("%50s %8.3f\n","angle to move",angle);
		printf("\n");
		fflush(stdout);

	} else {
	 */
	double bond2[3];																			//	bond around which we rotate

	//	FIXED POLYMER ANCHOR COORDINATES
	//	NOT NEEDED: YOU HAVE CO->EL TO SPECIFY THEM
	/*for (int i = 0; i < 3; i++) {																//	make correct active anchor coordinates
			subtract(co[bond_number-1],co[number_atoms-3+i],polpoint[i]);							//	crucial_____XXXXX OLD CODE
																									// segmentation fault can happen
		}*/

	//	CREATE EXACT BOND NUMBERS FOR REAL CHAIN ABOUT WHICH MUST BE ROTATED
	//	CAN ALSO BE ALLOCATED OUTSIDE OF FUNCTION AS FIXED 1-DIMENSIONAL ARRAY
	/*int bond_number;
		int *bondsequence;
		bondsequence = new int[bondrange];
		int bond_preimage = 1;
		for (int i = 0; i < bondrange; i++){																//	check for still within range
				bondsequence[i]	=  	(int) ( floor(3.0*((double) (bond_preimage))/2.0) + 1.0 );				//	mapping to real protein chain
				bond_preimage++;
		}*/

	//	COPY ANCHOR COORDINATES LOCALLY
	double **copyanchor;
	pointeralloc(&copyanchor,3,3);
	for (int i = 0; i < 3; i++ ){
		for(int j = 0; j < 3; j++)
			copyanchor[i][j] = anchor[i][j];
	}

	//printf("ANCHOR COORDINATES\n");										//	COORDINATES ARE OK
	//for(int i = 0; i < 3; i++){
	//printf("Anchor coordinates %8.3f %8.3f %8.3f\n",copyanchor[i][0],copyanchor[i][1],copyanchor[i][2]);
	//fflush(stdout);							//	LOOKS OK
	//}

	//printf("Anchor coordinates copied to local variable\n");
	//fflush(stdout);

	int bond_number;

	//	ITERATE OVER THE DIFFERENT BOND NUMBERS AND OPTIMIZE EACH ANGLE
	for(int k = 0; k < bondrange; k++){
		//printf("Running index in BESTANGLES %i\n",k);
		//fflush(stdout);

		bond_number = bondsequence[k];

		//printf("Bond number  %i\n",bond_number);
		//fflush(stdout);				-----> Bond numbers Ok

		//	PRINT INPUT COORDINATES
		/*printf("LOOP COORDINATES\n");										//	COORDINATES ARE OK
			  	  	  for(int i = 0; i < number_atoms; i++){
							printf("Coordinates %8.3f %8.3f %8.3f\n",co[i][0],co[i][1],co[i][2]);
							fflush(stdout);
						}*/

		//	Actual algorithm
		for (int i = 0; i < 3; i++) {																//	make correct flexible anchor coordinates
			//printf("Loop index %i\n",i);
			//fflush(stdout);			//	LOOKS OK
			subtract(co[bond_number-1],copyanchor[i],anchorpoint[i]);									//	in translated frame ... OLD ONE
		}

		//printf("Subtraction passed\n");
		//fflush(stdout);


		diffvectorseq(co,bond_number + 1,bond2);

		rotatematrix(bond2,matrix);

		//printf("Construction matrix passed\n");
		//fflush(stdout);



		//rotatematrix(bond2,matrix); 										//	But mind bond1 ABOVE (not dealt with) -> bond_number 1 should be avoided
		for (int i = 0; i < 3; i++) {										//	apply to each of the anchor points the previous procedure
			//matrixXvector(matrix,polpoint[i],result1[i]);					//	problem here somewhere....
			matrixXvector(matrix,anchorpoint[i],result1[i]);
			matrixXvector(matrix,co[number_atoms-3+i],result2[i]);			// Index in co --- possible source of mistakes // indices 1 and 2 of result switched to deal with fixed anchor rotation

		}

		//	MATRIX ALGEBRA APPROACH
		for (int i = 0; i < 3; i++) {
			result2[i][0] = 0.0;												//	setting x coordinate at zero
			result1[i][0] = 0.0;
		}

		double denominator = 0;
		for (int i = 0; i < 3; i++) {
			denominator += dotprod(result1[i],result2[i]);
		}


		for (int i = 0; i < 3; i++) {						//	What is this protocol exactly ? -> Does not seem to fit with theory  -> is OK
			for (int j = 0; j < 3; j++) {
				result3[i][0] = 0.0;						//	Construct orthonormal with cross product of result1 and [1 0 0] and
				result3[i][1] = -result1[i][2];				//	work out analytically
				result3[i][2] = result1[i][1];
			}
		}

		double numerator = 0;
		for (int i = 0; i < 3; i++) {
			numerator += dotprod(result2[i],result3[i]);
		}

		angle = atan2(numerator,denominator);

		//	Correction to get values in [-PI,PI] range
		if (angle > M_PI) {
			angle = angle - 2*M_PI;
		} else if (angle < -M_PI) {
			angle = angle + 2*M_PI;
		}



		angles[k] = -angle;

		//	ROTATE THE ANCHOR COORDINATES
		//	faster quaternion/matrix=spinor/matrix code ----- seems to work also---slightly faster
		double rotor[4];
		quaternion(bond2,angle/2,rotor);
		quatmatrix(rotor,matrix);
		for(int i = 0; i < 3;i++) {
			matrixXvector(matrix,anchorpoint[i],copyanchor[i]);
		}

	}

	//printf("Inside loop finished in BESTANGLES\n");
	//fflush(stdout);

	//printf("\n");
	//for(int i = 0; i < bondrange; i++){
	//printf("Angle values obtained in rad %8.3f\n",angles[i]);
	//fflush(stdout);
	//}
	//getchar();

	//delete[] bondsequence;	//	You may not do this ---> CODE HAS CHANGED AND IS NOW ALLOCATED OUTSIDE
	freepointer(copyanchor,3);
	//}

	//return 0;
}



void bestanglesm(	double** co,		/*loop coordinates*/
		double *angles,		/*output dihedral angle variations*/
		double *angles2,	/*physically representing bond angles*/
		int *bondsequence,	/*array with rotation bonds*/
		int bondrange,		/*# bonds to rotate*/
		int number_atoms,	/*# atoms in loop*/
		double** anchor,	/*fixed end anchor coordinates*/
		double** matrix,	/*auxiliary matrix*/
		double** result1,	/*matrix for spatial points 1*/
		double** result2,	/*matrix for spatial points 2*/
		double** result3,	/*matrix for spatial points 3*/
		double*  bond,		/*rotation bond in optimization*/
		double *rotor,		/*rotor for spatial transformation*/
		double *result,		/*auxiliary array*/
		double **copyanchor, /*matrix for copy anchor*/
		double **polanchor,  /*matrix for polymer anchor????*/
		double **anchorpoint2	/*other extra matrix*/
)
{

	double 	angle;															//	Optimal angle to rotate at given bond
	//double 	bond[3];													//	Coordinates of bond around which we rotate
	int 	bond_number;													//	Index of local bond number
	double 	numerator;
	double 	denominator;
	int 	switch_debug = 0;												//	Switch for debugging code

	for (int i = 0; i < 3; i++ ){
		for(int j = 0; j < 3; j++){
			copyanchor[i][j] 	= 	anchor[i][j];
			//anchorpoint2[i][j]	=	anchorpoint[i][j];
		}
	}

	if(switch_debug == 1){													//	The anchor coordinates seem to be copied always correctly
		message("NEW ITERATION",1);
		show3("Fixed anchor coordinates","X Y Z coordinates",anchor,3);
		show3("Copied anchor coordinates","X Y Z coordinates",copyanchor,3);
	}

	message("Anchor coordinates copied ...",0);

	if(switch_debug == 1){													//	Copied anchor of loop is Ok --- last three atoms
		message("NEW ITERATION",1);
		show3("Loop coordinates","X Y Z coordinates",co,number_atoms);
		show3("Loop anchor coordinates","X Y Z coordinates",co+number_atoms - 3 ,3);
		getchar();
	}



	//	ITERATE OVER THE DIFFERENT BOND ROTATIONS AND OPTIMIZE EACH ANGLE
	for(int k = 0; k < bondrange; k++){

		//show1("Index for loop",k);
		bond_number = bondsequence[k];										//	Correct numbers for small protein loop

		if(switch_debug == 1){												//	Checking correct bond numbers
			show1("Bond number to rotate",bond_number);
			getchar();
		}

		//	Actual algorithm
		for (int i = 0; i < 3; i++) {										//	make correct flexible anchor coordinates
			subtract(co[bond_number-1],copyanchor[i],anchorpoint2[i]);		//	in translated frame ... OLD ONE
			subtract(co[bond_number-1],co[number_atoms - 3 + i],polanchor[i]);		//	in translated frame ... OLD ONE
		}																		//	----->  Index [bond_number] also possible

		if(switch_debug == 1){
			show3("Loop coordinates","X Y Z coordinates",co,number_atoms);
			show3("Copied anchor coordinates","X Y Z coordinates",copyanchor,3);
			show3("Subtracted fixed anchor","X Y Z coordinates",anchorpoint2,3);
			show3("Subtracted loop anchor","X Y Z coordinates",polanchor,3);
			getchar();
		}


		//	OPTIMIZATION PROTOCOL
		//-----------------------------------------------------------------------------------------------------------------------------------
		diffvectorseq(co,bond_number + 1,bond);
		if(switch_debug == 1){
			show3("Loop coordinates","X Y Z coordinates",co,number_atoms);
			show1("X rotate bond",*bond);
			show1("X rotate bond",*(bond+1));
			show1("X rotate bond",*(bond+2));
		}
		normalize(bond,bond);												//	No errors inside ---- obligatory for later
		if(switch_debug == 1){
			show1("X rotate bond",*bond);
			show1("X rotate bond",*(bond+1));
			show1("X rotate bond",*(bond+2));
			getchar();
		}
		//=========================================================================================================================================
		rotatematrix(bond,matrix);		//	NEEDS DEBUGGING? OR NOT?



		for (int i = 0; i < 3; i++) {										//	apply to each of the anchor points the previous procedure
			matrixXvector(matrix,anchorpoint2[i],result1[i]);				//	transformed anchor coordinates
			matrixXvector(matrix,polanchor[i],result2[i]);					// Can be put out of loop ...
		}

		//	MATRIX ALGEBRA APPROACH
		for (int i = 0; i < 3; i++) {
			result2[i][0] = 0.0;												//	setting x coordinate at zero
			result1[i][0] = 0.0;
		}

		denominator = 0;
		for (int i = 0; i < 3; i++) {
			denominator += dotprod(result1[i],result2[i]);
		}


		for (int i = 0; i < 3; i++) {												//	What is this protocol exactly ? -> Does not seem to fit with theory  -> is OK
			for (int j = 0; j < 3; j++) {
				result3[i][0] = 0.0;												//	Construct orthonormal with cross product of result1 and [1 0 0] and
				result3[i][1] = -result1[i][2];										//	work out analytically
				result3[i][2] = result1[i][1];
			}
		}

		numerator = 0;
		for (int i = 0; i < 3; i++) {
			numerator += dotprod(result2[i],result3[i]);
		}

		angle = atan2(numerator,denominator);
		correctdihedral(&angle);													//	Correction to get values in [-PI,PI] range --- is this needed to do????


		//-----------------------------------------------------------------------------------------------------------------------------------

		//	Debugging code for optimization
		/*srand( (unsigned) time(NULL) );
		angle 	= 	((double) rand()/(double) RAND_MAX)*M_PI*2 - M_PI;*/


		//	pure spinor code ---- code that works
		/*double 	*revrotor;														//	Coordinates of bond around which we rotate
		revrotor = (double*) malloc(sizeof(double)*4);
		for(int i = 0; i < 4; i++)
		{
			revrotor[i] = 0.0;
		}
		spinor(bond,angle,rotor);
		revspinor(rotor,revrotor);
		for(int i = 0; i < 3 ;i++)
		{
				geom32(anchorpoint2[i],rotor,result);
				geom23(revrotor,result,copyanchor[i]);
		}
		free(revrotor);
		 */


		//double rotor[4];															//	Rotate effectively the anchor
		quaternion(bond,angle/2,rotor);
		quatmatrix(rotor,matrix);
		for(int i = 0; i < 3;i++) {
			matrixXvector(matrix,anchorpoint2[i],copyanchor[i]);
		}

		if(switch_debug == 1){
			message("AFTER OPTIMIZATION",1);
			show3("moved fixed anchor","X Y Z coordinates",copyanchor,3);
		}


		//  Reconstruct full active anchor coordinates
		for (int i = 0; i < 3; i++) {
			add(co[bond_number-1],copyanchor[i],copyanchor[i]);
		}

		if(switch_debug == 1){
			show3("Loop coordinates","X Y Z coordinates",co,number_atoms);
			show3("Added anchor coordinates","X Y Z coordinates",copyanchor,3);
			getchar();
		}

		if(switch_debug == 1){
			//message("After calculation",1);
			show3("Added and moved fixed anchor","X Y Z coordinates",copyanchor,3);

		}


		angles[k] 			=  -angle; /*working for spinor scheme*/
		angles2[bond_number] = -angle;	/*new scheme to be used --- sign correct? */


	}


	//getchar();

	if(switch_debug == 1){
		message("AFTER FUNCTION",1);
		show3("Fixed anchor coordinates","X Y Z coordinates",anchor,3);
	}
	//getchar();

	message("Loop finished in obtaining angles ...",0);
	message("Arrays deallocated ...",0);

}

void bestanglesm3(double** co,double *angles,int *bondsequence,int bondrange,int number_atoms,double** anchor,						//	Preallocated version
		double** matrix,double** polpoint,double** anchorpoint,double** result1,double** result2, double** result3) {

	double 	angle;															//	Optimal angle to rotate at given bond
	double 	bond[3];														//	Coordinates of bond around which we rotate
	int 	bond_number;													//	Index of local bond number
	double 	numerator;
	double 	denominator;

	//	COPY ANCHOR COORDINATES ----- LOCALLY AND INITIALIZE  //			//	Is always the same each time called--correct atoms invoked
	double **copyanchor;													//	Real anchor ---- copied to enter algorithm
	double **polanchor;														//	Flexible anchor ---- used in transformations
	pointeralloc(&copyanchor,3,3);
	pointeralloc(&polanchor,3,3);

	for (int i = 0; i < 3; i++ ){
		for(int j = 0; j < 3; j++){
			copyanchor[i][j] = anchor[i][j];
		}
	}

	message("Anchor coordinates copied ...",0);

	//	ITERATE OVER THE DIFFERENT BOND ROTATIONS AND OPTIMIZE EACH ANGLE
	for(int k = 0; k < 3; k++){

		//show1("Index for loop",k);

		bond_number = bondsequence[k];										//	Correct numbers for small protein loop

		//	Actual algorithm
		for (int i = 0; i < 3; i++) {										//	make correct flexible anchor coordinates
			subtract(co[bond_number-1],copyanchor[i],anchorpoint[i]);		//	in translated frame ... OLD ONE
			subtract(co[bond_number-1],co[number_atoms - 3 + i],polanchor[i]);		//	in translated frame ... OLD ONE
		}																	//	----->  Index [bond_number] also possible

		diffvectorseq(co,bond_number + 1,bond);
		normalize(bond,bond);												//	No errors inside ---- obligatory for later
		rotatematrix(bond,matrix);

		for (int i = 0; i < 3; i++) {										//	apply to each of the anchor points the previous procedure
			matrixXvector(matrix,anchorpoint[i],result1[i]);				//	transformed anchor coordinates
			matrixXvector(matrix,polanchor[i],result2[i]);					// Can be put out of loop ...
		}

		//	MATRIX ALGEBRA APPROACH
		for (int i = 0; i < 3; i++) {
			result2[i][0] = 0.0;												//	setting x coordinate at zero
			result1[i][0] = 0.0;
		}

		denominator = 0;
		for (int i = 0; i < 3; i++) {
			denominator += dotprod(result1[i],result2[i]);
		}


		for (int i = 0; i < 3; i++) {												//	What is this protocol exactly ? -> Does not seem to fit with theory  -> is OK
			for (int j = 0; j < 3; j++) {
				result3[i][0] = 0.0;												//	Construct orthonormal with cross product of result1 and [1 0 0] and
				result3[i][1] = -result1[i][2];										//	work out analytically
				result3[i][2] = result1[i][1];
			}
		}

		numerator = 0;
		for (int i = 0; i < 3; i++) {
			numerator += dotprod(result2[i],result3[i]);
		}

		angle = atan2(numerator,denominator);
		//correctdihedral(&angle);													//	Correction to get values in [-PI,PI] range --- is this needed to do????



		//	Debugging code for optimization
		/*srand( (unsigned) time(NULL) );
		angle 	= 	((double) rand()/(double) RAND_MAX)*M_PI*2 - M_PI;*/


		double rotor[4];															//	Rotate effectively the anchor
		quaternion(bond,angle/2,rotor);
		quatmatrix(rotor,matrix);
		for(int i = 0; i < 3;i++) {
			matrixXvector(matrix,anchorpoint[i],copyanchor[i]);
		}

		//  Reconstruct full active anchor coordinates
		for (int i = 0; i < 3; i++) {
			add(co[bond_number-1],copyanchor[i],copyanchor[i]);
		}

		angles[k] = - angle;


	}

	message("Loop finished in obtaining angles ...",0);

	freepointer(copyanchor,3);
	freepointer(polanchor,3);

	message("Arrays deallocated ...",0);

}



void fix_anchor(double **co,double **matrix,double **anchor,double **anchor_fix_centroid,int nr_atoms_loop) {
	double centroid[3];			//	Algorithm does not work properly
	double a1;
	double a2;
	double a3;
	double a4;
	double b1;
	double b2;
	double x;
	double y;
	double z;

	printf("%-30s %8.3f %8.3f %8.3f\n","Anchors fixed 1",anchor_fix_centroid[0][0],anchor_fix_centroid[0][1],anchor_fix_centroid[0][2]);
	printf("%-30s %8.3f %8.3f %8.3f\n","Anchors fixed 2",anchor_fix_centroid[1][0],anchor_fix_centroid[1][1],anchor_fix_centroid[1][2]);
	printf("%-30s %8.3f %8.3f %8.3f\n","Anchors fixed 3",anchor_fix_centroid[2][0],anchor_fix_centroid[2][1],anchor_fix_centroid[2][2]);


	printf("%-30s %8.3f %8.3f %8.3f\n","Anchors 1",anchor[0][0],anchor[0][1],anchor[0][2]);
	printf("%-30s %8.3f %8.3f %8.3f\n","Anchors 2",anchor[1][0],anchor[1][1],anchor[1][2]);
	printf("%-30s %8.3f %8.3f %8.3f\n","Anchors 3",anchor[2][0],anchor[2][1],anchor[2][2]);

	//	Calculate centroid variable anchor
	for(int j = 0; j < 3; j++) {
		centroid[j] = 0.0;
		for(int i = 0; i < 3; i++) {
			centroid[j] += anchor[i][j];
		}
		centroid[j] /= 3;

	}

	printf("%-30s %8.3f %8.3f %8.3f\n","centroid 1",centroid[0],centroid[1],centroid[2]);


	// Translate loop coordinates
	for (int i = 0; i < nr_atoms_loop; i++) {
		for (int j = 0; j < 3; j++) {
			co[i][j] -= centroid[j];

		}
	}

	// Translate anchor coordinates						//	Must be manipulated as well!!!
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			anchor[i][j] -= centroid[j];

		}
	}


	printf("%-30s %8.3f %8.3f %8.3f\n","Anchors 1",anchor[0][0],anchor[0][1],anchor[0][2]);
	printf("%-30s %8.3f %8.3f %8.3f\n","Anchors 2",anchor[1][0],anchor[1][1],anchor[1][2]);
	printf("%-30s %8.3f %8.3f %8.3f\n","Anchors 3",anchor[2][0],anchor[2][1],anchor[2][2]);

	//	Find rotation matrix									//	Important that this code works----when it does---work on speed
	for (int i = 0; i < 3; i++) {
		a1	= anchor[2][1]*anchor[0][2]      			- 	anchor[2][2]*anchor[0][1];
		a2	= anchor[1][2]*anchor_fix_centroid[0][i] 	- 	anchor[0][2]*anchor_fix_centroid[1][i];
		a3	= anchor[0][2]*anchor[1][1] 				- 	anchor[1][2]*anchor[0][1];
		a4	= anchor[2][2]*anchor_fix_centroid[0][i] 	- 	anchor[0][2]*anchor_fix_centroid[1][i];

		b1	= anchor[0][2]*anchor[1][0] 				- 	anchor[0][0]*anchor[1][2];
		b2	= anchor[0][2]*anchor[2][0] 				- 	anchor[2][2]*anchor[0][0];

		x	= ( a1*a2 - a3*a4 )/ ( a1*b1 - a3*b2 );
		y	= ( a2 -  b1*x )/a3;
		z   = (anchor_fix_centroid[0][i] - (anchor[0][0]*x) - (anchor[0][1]*y))/anchor[0][2];

		matrix[i][0] = x;
		matrix[i][1] = y;
		matrix[i][2] = z;

		printf("division terms %8.3f %8.3f %8.3f\n", a1*b1 - a3*b2, a3,anchor[0][2]);
	}


	for(int i = 0; i< nr_atoms_loop; i++) {
		matrixXvector(matrix,co[i],co[i]);
	}

	//	Translated loop coordinates back
	for (int i = 0; i < nr_atoms_loop; i++) {
		for (int j = 0; j < 3; j++) {
			co[i][j] += centroid[j];

		}
	}
}




void gafixanchor(double **co,double **anchor,double **anchor_fix_centroid,int nr_atoms_loop,double **normals1,double **normals2, double **normals3) {
	double centroid[3];
	double phi1;
	double phi2;
	double spinor1[4];
	double spinor2[4];
	double spinor3[4];
	double rspinor1[4];
	double rspinor2[4];
	double rspinor3[4];
	double result[4];

	// /*
	printf("Initial coordinates\n");
	printf("%-30s %8.3f %8.3f %8.3f\n","Anchors fixed 1",anchor_fix_centroid[0][0],anchor_fix_centroid[0][1],anchor_fix_centroid[0][2]);
	printf("%-30s %8.3f %8.3f %8.3f\n","Anchors fixed 2",anchor_fix_centroid[1][0],anchor_fix_centroid[1][1],anchor_fix_centroid[1][2]);
	printf("%-30s %8.3f %8.3f %8.3f\n","Anchors fixed 3",anchor_fix_centroid[2][0],anchor_fix_centroid[2][1],anchor_fix_centroid[2][2]);
	// */

	//	Calculate centroid variable anchor
	for(int j = 0; j < 3; j++) {
		centroid[j] = 0.0;
		for(int i = 0; i < 3; i++) {
			centroid[j] += anchor[i][j];
		}
		centroid[j] /= 3;

	}

	// Translate loop coordinates
	for (int i = 0; i < nr_atoms_loop; i++) {
		for (int j = 0; j < 3; j++) {
			co[i][j] -= centroid[j];

		}
	}

	// Translate anchor coordinates						//	Must be manipulated as well!!!
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			anchor[i][j] -= centroid[j];

		}
	}

	//	Basic algorithm
	//	Normals to triangle planes
	geom11out(anchor[0],anchor[2],normals1[0]);
	geomi2(normals1[0],normals2[0]);
	normalize(normals2[0],normals3[0]);	// active normal
	geom11out(anchor_fix_centroid[0],anchor_fix_centroid[2],normals1[1]);
	geomi2(normals1[1],normals2[1]);
	normalize(normals2[1],normals3[1]);	// fixed normal

	//	Double normal
	geom11out(normals3[0],normals3[1],normals1[2]); //	double normal
	geomi2(normals1[2],normals2[2]);
	normalize(normals2[2],normals3[2]);

	//	Spinor 1
	phi1 = acos(dotprod(normals3[0],normals3[1]));

	spinor(normals3[2],phi1,spinor1);
	revspinor(spinor1,rspinor1);

	//	Spinor 2
	for (int i = 0; i < 3; i++) {
		geom32(anchor[i],spinor1,result);
		geom23(rspinor1,result,anchor[i]);
	}

	phi2 = acos(dotprod(anchor[0],anchor_fix_centroid[0])/(norm(anchor[0])*norm(anchor_fix_centroid[0])));

	//	To test correct angle
	geom11out(anchor[0],anchor_fix_centroid[0],normals1[0]);
	geomi2(normals1[0],normals2[0]);
	normalize(normals2[0],normals3[0]);

	if (dotprod(normals3[0],normals3[1]) > 0.0) {

	} else {
		phi2 = -phi2;
	}

	spinor(normals3[1],phi2,spinor2);
	revspinor(spinor2,rspinor2);

	/* DOES NOT COMPLETELY WORK
	// Complete spinor
	spinorproduct(spinor1,spinor2,spinor3); //	mistake in spinorproduct---debug
	revspinor(spinor3,rspinor3);
	// */


	/*
	for (int i = 0; i < nr_atoms_loop; i++) {
		geom32(co[i],spinor3,result);
		geom23(rspinor3,result,co[i]);
	}
	// */

	// /*
	for (int i = 0; i < nr_atoms_loop; i++) {
		geom32(co[i],spinor1,result);
		geom23(rspinor1,result,co[i]);
		geom32(co[i],spinor2,result);
		geom23(rspinor2,result,co[i]);
	}

	// */
	// /*
	printf("End coordinates\n");
	printf("%-30s %8.3f %8.3f %8.3f\n","co anchor 1",co[nr_atoms_loop-3+0][0],co[nr_atoms_loop-3+0][1],co[nr_atoms_loop-3+0][2]);
	printf("%-30s %8.3f %8.3f %8.3f\n","co anchor 2",co[nr_atoms_loop-3+1][0],co[nr_atoms_loop-3+1][1],co[nr_atoms_loop-3+1][2]);
	printf("%-30s %8.3f %8.3f %8.3f\n","co anchor 3",co[nr_atoms_loop-3+2][0],co[nr_atoms_loop-3+2][1],co[nr_atoms_loop-3+2][2]);
	// */

	//	Translated loop coordinates back
	for (int i = 0; i < nr_atoms_loop; i++) {
		for (int j = 0; j < 3; j++) {
			co[i][j] += centroid[j];
		}
	}

}



/*===============================================================================
 *
 *	BONDSELECT3.CPP
 *
 *	Created on: February 21, 2010
 *      Author: Pieter Chys
 *
 *===============================================================================
 *
 * BONDSELECT4.CPP selects randomly a bond from a given & continuous range (1-range) and
 * maps it back to the real position in the protein chain. This mapping is in agreement
 * with the model of the protein chain in the GA_algorthm families and maps thus back
 * to a discrete range with CA-CO the first remapped rotation bond. The input
 * is the range of available bonds. This is simply the total number of bonds that
 * can be used without consideration of the exact spread over the real chain.
 *
 *--------------------------------------------------------------------------------
 *
 *  INPUT: - int range = total number of bonds available
 *
 *-------------------------------------------------------------------------------
 *
 *  OUTPUT: - explicitly int
 *  		- int which is the real bond in the protein loop to be rotated
 *
 *-------------------------------------------------------------------------------
 *
 *  DEPENDENCIES:
 *
 *===============================================================================
 */

int	bond_select3(int range) {
	double 	number 	= 	(double) ceil( ((double)rand()/(double)RAND_MAX) * (double) range ) ;
	int 	bond	=  	(int) ( floor(3.0*number/2.0) + 1.0 );
	return bond;
}

/*===============================================================================
 *
 *	BONDSELECT4.CPP
 *
 *	Created on: February 21, 2010
 *      Author: Pieter Chys
 *
 *===============================================================================
 *
 * BONDSELECT4.CPP selects randomly a bond from a given & continuous range (1-range) and
 * maps it back to the real position in the protein chain. This mapping is in agreement
 * with the model of the protein chain in the GA_algorthm families and maps thus back
 * to a discrete range with CA-CO the first remapped rotation bond. The input
 * is the range of available bonds. This is simply the total number of bonds that
 * can be used without consideration of the exact spread over the real chain. In
 * addition this bondselect algorithm avoids that the previous bond is chosen. As
 * such the algorithm is equal to algorithm BONDSELECT3.CPP but with inclusion
 * of avoiding the previously chosen bond.
 *
 *--------------------------------------------------------------------------------
 *
 *  INPUT: - int range = total number of bonds available
 *  	   - int *pointer which allows to change the value of the previous bond
 *
 *-------------------------------------------------------------------------------
 *
 *  OUTPUT: - explicitly int
 *  		- int which is the real bond in the protein loop to be rotated
 *  		- implicitly int *pointer to the previous bond value
 *
 *
 *-------------------------------------------------------------------------------
 *
 *  DEPENDENCIES:
 *
 *===============================================================================
 */

int	bond_select4(int range,int *bond_previous) {
	double 	number 	= 	(double) ceil( ((double)rand()/(double)RAND_MAX) * (double) range ) ;
	int 	bond	=  	(int) ( floor(3.0*number/2.0) + 1.0 );
	while (bond==*bond_previous) {
		double 	number 	= 	(double) ceil( ((double)rand()/(double)RAND_MAX) * (double) range ) ;
		bond			=  	(int) ( floor(3.0*number/2.0) + 1.0 );
	}
	*bond_previous = bond;
	return bond;
}

int	bond_select5(int range,int *bond_previous) {
	double 	number 	= 	(double) ceil( ((double)rand()/(double)RAND_MAX) * (double) range ) ;
	int 	bond	=  	(int) ( floor(3.0*number/2.0) + 1.0 );
	while (bond == *bond_previous || bond == 2) {
		double 	number 	= 	(double) ceil( ((double)rand()/(double)RAND_MAX) * (double) range ) ;
		bond			=  	(int) ( floor(3.0*number/2.0) + 1.0 );
	}
	*bond_previous = bond;
	return bond;
}

//	6 is same as 4 but with more output arguments

/*===============================================================================
 *
 *	BONDSELECT6.CPP
 *
 *	Created on: September 1, 2011
 *      Author: Pieter Chys
 *
 *===============================================================================
 *
 * BONDSELECT6.CPP selects a bond from a given and continuous range (1 -> range) and
 * outputs it to use as choice for the bond around which rotation and optimization
 * will happen. This version picks randomly bonds from the full range (1 -> range)
 * and makes a mapping.  The mapping outputs to the real position of the bond in
 * the protein chain. This mapping is in agreement with the model of the protein
 * chain in the GA_algorthm families and maps thus back to a discrete range with
 * CA-CO the first rotation bond. This selection algorithm contains additionally
 * a tweak that prevents that the previous bond is skipped always. The input is
 * the range of available bonds. This s simply the total number of bonds that
 * can be used without consideration of the exact spread over the real chain. This
 * version is used for the sequential variant BONDSELECT7.CPP which was derived
 * from it.
 *
 *--------------------------------------------------------------------------------
 *
 *  INPUT: - int range 			= 	total number of bonds available
 *  	   - int *bond_preimage =    allows to pass the value of the chosen bond
 *  	   - int *bond_previous =    allows to pass the value of the previous bond
 *  								 and prevent it from being used
 *
 *
 *-------------------------------------------------------------------------------
 *
 *  OUTPUT: - explicitly int which is the bond number (mapped to protein chain)
 *  		- implicitly int *pointer to the bond pre-real number chosen
 *  		- implicitly int *pointer to the previously chosen bond in protein chain
 *
 *
 *-------------------------------------------------------------------------------
 *
 *  DEPENDENCIES:
 *
 *===============================================================================
 */


int	bond_select6(int range,int *bond_previous,int *bond_preimage)
{
	double 	number 	= 	(double) ceil( ((double)rand()/(double)RAND_MAX) * (double) range );
	int 	bond	=  	(int) ( floor(3.0*number/2.0) + 1.0 );
	while (bond == *bond_previous || bond == 1)
	{//	correction for extreme case of bond = 1 ... how possible is this?
		//while (bond == *bond_previous || bond == 0) {
		//	correction for extreme case of bond = 1 ... BEFORE ... spinor FW
		double 	number 	= 	(double) ceil( ((double)rand()/(double)RAND_MAX) * (double) range ) ;
		bond			=  	(int) ( floor(3.0*number/2.0) + 1.0 );
	}
	*bond_preimage = (int) number;
	*bond_previous = bond;
	return bond;
}


/*===============================================================================
 *
 *	BONDSELECT.CPP
 *
 *	Created on: May 29th, 2012
 *      Author: Pieter Chys
 *
 *===============================================================================
 *
 * BONDSELECT.CPP selects a bond from a given and continuous range (1 -> range) and
 * outputs it to use as choice for the bond around which rotation and optimization
 * will happen. This version picks randomly a number from the full range (1 -> range)
 * and makes the bond mapping. The mapping outputs to the real position of the bond in
 * the protein chain. This mapping is in agreement with the model of the protein
 * chain in the GA_algorthm families and maps thus back to a discrete range with
 * CA-CO the first rotation bond. This selection algorithm contains additionally
 * a tweak that makes that the previous bond is always skipped. The input is
 * the range of available bonds.  Here a marker array is used which also determines
 * whether the computed bond can be used. If the marker for the bond has value zero
 * it cannot be used and further cycling must occur. However, if the marker value is
 * equal to 1 the bond can be taken.  The marker array is passed as pointer and
 * contains all markers for the different bonds. The range is simply the total
 * number of bonds that can be used without consideration of the exact spread
 * over the real chain. This version is a straight derivation from BOND_SELECT6.CPP
 * to account for unallowed bond rotations under the normal PHI and PSI angles.
 *
 *--------------------------------------------------------------------------------
 *
 *  INPUT: - int range 			= 	total number of bonds available
 *  	   - int *bond_preimage =   allows to pass the value of the chosen bond
 *  	   - int *bond_previous =   allows to pass the value of the previous bond
 *  								and prevent it from being used
 *         - int *marker		=	array with length equal to the number of atoms
 *                                  in the loop containing 1 for PHI and PSI angles
 *                                  and 0 for the omega angles
 *
 *
 *-------------------------------------------------------------------------------
 *
 *  OUTPUT: - explicitly int which is the bond number (mapped to protein chain)
 *  		- implicitly int *pointer to the bond pre-real number chosen
 *  		- implicitly int *pointer to the previously chosen bond in protein chain
 *
 *
 *-------------------------------------------------------------------------------
 *
 *  DEPENDENCIES:	- an int *marker array must be passed and have been constructed
 *                    to use in conjunction
 *
 *===============================================================================
 */

int	bondselect(int range,int *marker,int *previous,int *preimage,int *infinite)
{
	double 	number;		/*picking a number from an integer range: e.g. [1-N] and pick number 1 -> N from it*/
	int bond;			/*actual rotation bond that is constructed from NUMBER*/
	int count;

	number 	= 	double( ceil(( double(rand())/ double(RAND_MAX)) * double(range)) );
	bond	=  	int( floor(3.0*number/2.0) + 1.0 );

	count 	= 0;
	while (bond == *previous || bond == 1 || marker[bond] == 0)
	{
		//printf("Not finding the correct angle\n");
		number 	= 	double ( ceil((double (rand())/ double (RAND_MAX)) * double (range)) );
		bond	=  	int ( floor(3.0*number/2.0) + 1.0 );
		count++;

		if (count > 110)
		{
			break;
		}
	}

	*infinite = count;
	*preimage = int (number);
	*previous = bond;

	return bond;
}

/*===============================================================================
 *
 *	BONDSELECT7.CPP
 *
 *	Created on: September 1, 2011
 *      Author: Pieter Chys
 *
 *===============================================================================
 *
 * BONDSELECT7.CPP selects a bond from a given and continuous range (1->range) and
 * outputs it to use as choice for the bond around which rotation and optimization
 * will happen. This version progresses sequentially through the range: 1,2, ... ,
 * range and resets the counter to one hereafter. The number from the range is taken
 * and mapped to the real position in the protein chain. This mapping is in agreement
 * with the model of the protein chain in the GA_algorthm families and maps thus back
 * to a discrete range with CA-CO the first mapped rotation bond. The input
 * is the range of available bonds. This is simply the total number of bonds that
 * can be used without consideration of the exact spread over the real chain. This
 * version is developed from BONDSELECT6.CPP
 *
 *--------------------------------------------------------------------------------
 *
 *  INPUT: - int range 			= 	total number of bonds available
 *  	   - int *bond_preimage =    allows to change the value of the previous bond
 *  								 in the sequence
 *
 *-------------------------------------------------------------------------------
 *
 *  OUTPUT: - explicitly int which is the bond number (mapped to protein chain)
 *  		- implicitly int *pointer to the previous bond number value
 *
 *
 *-------------------------------------------------------------------------------
 *
 *  DEPENDENCIES:
 *
 *===============================================================================
 */


int	bond_select7(int range,int *bond_preimage)
{

	//COPY FROM BELOW

	//	double 	number;		/*picking a number from an integer range: e.g. [1-N] and pick number 1 -> N from it*/
	int bond;			/*actual rotation bond that is constructed from NUMBER*/
	//
	//	number 	= 	double( ceil(( double(rand())/ double(RAND_MAX)) * double(range)) );
	//		bond	=  	int( floor(3.0*number/2.0) + 1.0 );
	//
	//
	//		while (bond == *previous || bond == 1 || marker[bond] == 0)
	//		{
	//			//printf("Not finding the correct angle\n");
	//			number 	= 	double ( ceil((double (rand())/ double (RAND_MAX)) * double (range)) );
	//			bond	=  	int ( floor(3.0*number/2.0) + 1.0 );
	//		}
	//
	//		*preimage = int (number);
	//		*previous = bond;
	//
	//		return bond;

	if (*bond_preimage < range)
	{															//	check for still within range
		*bond_preimage += 1;																//	update from sequence number
		bond	=  	(int) ( floor(3.0*((double) (*bond_preimage))/2.0) + 1.0 );				//	mapping to real protein chain

	} else {
		*bond_preimage = 1;																	//	reset to start sequence
		bond	=  	(int) ( floor(3.0*((double) (*bond_preimage))/2.0) + 1.0 );				//	mapping  for this case
	}
	//	printf("%d\n",bond);

	return bond;
}


int	bondselect(int range,int *marker,int *previous,int *preimage)
{
	double 	number;		/*picking a number from an integer range: e.g. [1-N] and pick number 1 -> N from it*/
	int bond;			/*actual rotation bond that is constructed from NUMBER*/


	number 	= 	double( ceil(( double(rand())/ double(RAND_MAX)) * double(range)) );
	bond	=  	int( floor(3.0*number/2.0) + 1.0 );


	while (bond == *previous || bond == 1 || marker[bond] == 0)
	{
		//printf("Not finding the correct angle\n");
		number 	= 	double ( ceil((double (rand())/ double (RAND_MAX)) * double (range)) );
		bond	=  	int ( floor(3.0*number/2.0) + 1.0 );
	}

	*preimage = int (number);
	*previous = bond;

	return bond;
}


/*===============================================================================
 *
 *	COSNTRUCTMARKER.CPP
 *
 *	Created on: May 29th, 2012
 *      Author: Pieter Chys
 *
 *===============================================================================
 *
 * CONSTRUCTMARKER.CPP constructs a marker array which is filled with boolean
 * values 1 or 0.  Filled values of 1 correspond with a rotation bond which is
 * active whereas filled values of 0 mean the rotation bond is inactive.  The index
 * of the bond always corresponds to the index of the atom at which the bond ends.
 * Bond 2 e.g. ends at atom 3 which has an index in arrays of [2]. Hence the value
 * in the int *marker for this bond is marker[2].  In general all PHI and PSI bonds
 * are active and thus correspond with bool values 1. On the other hand OMEGA angles
 * correspond to bool values 0.  The array may be modified appropriately but always
 * with the convention that 1 refers to an active bond and 0 to an inactive rotation
 * bond. The first and last residue are correctly modified since they contain inactive
 * rotation bonds.
 *
 *--------------------------------------------------------------------------------
 *
 *  INPUT: - int *marker =	array with length equal to the number of atoms
 *                          in the loop containing 1 for PHI and PSI angles
 *                          and 0 for the OMEGA angles
 *         - int nratoms = number of atoms in the loop and this is also the size
 *                         of the marker array which has to have that dimension
 *
 *
 *-------------------------------------------------------------------------------
 *
 *  OUTPUT: - explicitly void
 *  		- implicitly int *pointer to the marker array

 *
 *
 *-------------------------------------------------------------------------------
 *
 *  DEPENDENCIES:
 *
 *===============================================================================
 */

void constructmarker(int *marker, int nratoms)
{
	for (int i = 0; i < nratoms/3; i++)
	{
		marker[3*i] 	=	0;						/*Assignment of OMEGA angle-inactive*/
		marker[3*i + 1]	=	1;						/*Assignment of PHI angle-active*/
		marker[3*i + 2]	=   1;						/*Assignment of PSI angle-active*/
	}

	/*1st residue*/
	marker[1]	=	0;								/*first PHI bond-inactive*/

	/*last residue*/
	marker[nratoms-1] = 0;							/*last PSI bond-inactive*/
}

void constructpsi(int *marker, int nratoms)
{
	for (int i = 0; i < nratoms/3; i++)
	{
		marker[3*i] 	=	0;	/*assignment of OMEGA angle-inactive*/
		marker[3*i + 1]	=	0;	/*in-activate PHI angle*/
		marker[3*i + 2]	=   1;	/*activate PSI angle*/
	}

	/*1st residue*/
	marker[1]	=	0;		/*first PHI bond-inactive*/

	/*last residue*/
	marker[nratoms-1] = 0;	/*last PSI bond-inactive*/
	marker[nratoms-2] = 0;  /*last PHI bond-inactive*/
}




int	bondbreaksel(int start,int range,int *bond_previous,int *bond_preimage,int counter,int flag)
{
	//	START is the n-th and last rotation bond  which cannot rotate
	//  RANGE is the full range of rotation bonds
	//  FLAG is a switch for different loop lengths!
	//  COUNTER is the counter for th current individual closure
	//int step;
	double number;
	int bond;

	//step	=	    (int) ( floor(3.0*start/2.0) + 1.0 );
	number 	= 	    (double) ceil( ((double)rand()/(double)RAND_MAX) * (double) (range-start) ) ;
	bond	=	  	(int) ( floor(3.0*(number+double(start))/2.0) + 1.0 );

	while (bond == *bond_previous || bond == 1)
	{				//	correction for extreme case of bond = 1 ... how possible is this?
		//while (bond == *bond_previous || bond == 0) {				//	correction for extreme case of bond = 1 ... BEFORE ... spinor FW
		number 			=  	(double) ceil( ((double)rand()/(double)RAND_MAX) * (double) (range-start) ) ;
		bond			=  	(int) ( floor(3.0*(number+double(start))/2.0) + 1.0 );
	}

	*bond_preimage = (int) number;
	*bond_previous = bond;

	return bond;
}




int	bond_select_oddeven(int *angles,int range,int *bond_previous) {
	int 	number 	= 	(int) (ceil( ((double)rand()/(double)RAND_MAX) * (double) range ) - 1.0) ;
	int 	bond	=  	(int) ( floor(3.0*angles[number]/2.0) + 1.0 );
	while (bond==*bond_previous) {
		number 		= 	(int) (ceil( ((double)rand()/(double)RAND_MAX) * (double) range ) - 1.0) ;
		bond		=  	(int) ( floor(3.0*angles[number]/2.0) + 1.0 );
	}
	*bond_previous = bond;
	return bond;
}


/*===============================================================================
 *
 *	CENTROID.CPP
 *
 *	Created on: April 26, 2011
 *      Author: Pieter Chys
 *
 *===============================================================================
 *
 * CENTROID.CPP calculates the centroid of an atom system.  This is basically
 * the center of mass that is calculated without invoking the mass terms.
 *
 *--------------------------------------------------------------------------------
 *
 *  INPUT: - double **coordinates points to the atom coordinates which contain the x,y,z
 *  		 coordinates
 *  	   - int number_atoms = number atoms of the atom system
 *
 *
 *-------------------------------------------------------------------------------
 *
 *  OUTPUT: - double *center which contains the centroids' x,y and z coordinates
 *
 *-------------------------------------------------------------------------------
 *
 *  DEPENDENCIES:	- mind you that double **coordinates is often (size x 4) in
 *  				  format
 *
 *===============================================================================
 */

void	centroid(double **coordinates, double *center,int number_atoms){
	center[0] = 0.0;
	center[1] = 0.0;
	center[2] = 0.0;
	for(int i = 1; i < number_atoms; i++) {
		center[0] += coordinates[i][0];
		center[1] += coordinates[i][1];
		center[2] += coordinates[i][2];
	}
	center[0] /= number_atoms;
	center[1] /= number_atoms;
	center[2] /= number_atoms;
}


/*===============================================================================
 *
 *	RELATIVECO.CPP
 *
 *	Created on: April 26, 2011
 *      Author: Pieter Chys
 *
 *===============================================================================
 *
 * RELATIVECO.CPP calculates the relative coordinates of an atom system with respect
 * to the centroid.  This is the absolute coordinates minus the centroid. This
 * could be modified for the center of mass if you want.
 *
 *--------------------------------------------------------------------------------
 *
 *  INPUT: - double **coordinates are the absolute atom coordinates in (number_atoms x 3/4)
 *  		 format
 *	 	   - int number_atoms = number atoms of the atom system
 *
 *
 *-------------------------------------------------------------------------------
 *
 *  OUTPUT: - double **new_coordinates are the new and relative atom coordinates
 *            in (number_atoms x 3/4) format
 *
 *-------------------------------------------------------------------------------
 *
 *  DEPENDENCIES:	- mind you that double **coordinates is often (size x 4) in
 *  				  format
 *                  - the center is locally created
 *
 *===============================================================================
 */

void 	relativeco(double **coordinates,double **new_coordinates,int number_atoms){
	double center[3];
	centroid(coordinates,center,number_atoms);
	for(int i = 0; i < number_atoms; i++){
		for(int j = 0; j < 3; j++){
			new_coordinates[i][j] = coordinates[i][j] - center[j];
		}
	}
}



// DEBUGGING CODE TRICKS

//f//printf(stderr,"bond_number=%d number_atoms-bond_number-1=%d number_atoms=%d\n",bond_number,number_atoms-bond_number-1,number_atoms);

/*//RUNNING COMPLETE
//===========================================================================
		for(int i=0; i<(number_atoms-bond_number-1); i++)
		{
				f//printf(stderr,"%d\n",i);
				free(part[i]);
		}
		return;
//============================================================================
 */

/*===============================================================================
 *
 *  MOVECHAIN.CPP
 *
 *	Created on: Dec 13, 2010
 *      Author: Pieter Chys
 *
 *===============================================================================
 * MOVECHAIN.CPP takes the coordinates of a polymer chain and
 * translates them with a given vector.
 *
 *
 *--------------------------------------------------------------------------------
 *
 *  INPUT: - int referring to the number of bonds in the chain
 *  	   - double** pointer to the coordinates of the whole chain
 *  	   - double* pointer to the translation to be made
 *
 *-------------------------------------------------------------------------------
 *
 *  OUTPUT: - function return is explicitly VOID
 *  		- double** pointer which refers to the new polymer
 *  		  coordinates
 *
 *-------------------------------------------------------------------------------
 *
 *  DEPENDENCIES:
 *
 *===============================================================================
 */

void movechain(double** co,double** step,int atom_number,int number_atoms,int choice)
{
	/*
	for(int i=0;i<number_atoms;i++)
	{
		for(int j=0;j<3;j++)
		{
			if (choice==1){
				co[i][j] = co[i][j]+step[atom_number-1][j];
			}
			else {co[i][j] = co[i][j]-step[atom_number-1][j];};
		}
	}
	 */
	switch (choice) {
	case 1:
		for(int i = 0; i < number_atoms; i++) {
			for(int j = 0; j < 3; j++) {
				co[i][j] = co[i][j]+step[atom_number-1][j];
				//co[i][j] +=  step[atom_number-1][j];
			}
		}
		break;
	case 2:
		for(int i = 0; i < number_atoms; i++) {
			for(int j = 0; j < 3; j++) {
				co[i][j] = co[i][j]-step[atom_number-1][j];
				//co[i][j] -= step[atom_number-1][j];
			}
		}
		break;
	}
}

/*===============================================================================
 *
 *  SELECTPART.CPP
 *
 *	Created on: Dec 13, 2010
 *      Author: Pieter Chys
 *
 *===============================================================================
 * SELECTPART.CPP takes the coordinates of a specified
 * segment of the polymer chain. For the present the part
 * is selected at the terminal end-always running till the
 * end of the chain.  The segment starts at atom START_ATOM
 * and finishes at NUMBER_ATOMS. Everything is done by
 * using pointers.
 *
 *--------------------------------------------------------------------------------
 *
 *  INPUT: - int referring to the number of bonds in the chain
 *  	   - int referring to the first atom of the segment
 *  	   - double** pointer to the coordinates of the whole chain
 *  	   - double** pointer to the newe part
 *
 *-------------------------------------------------------------------------------
 *
 *  OUTPUT: - function return is explicitly VOID
 *  		- double** pointer which refers to the selected portion
 *            of the chain
 *
 *-------------------------------------------------------------------------------
 *
 *  DEPENDENCIES:
 *
 *===============================================================================
 */

void selectpart(double** co, int start_atom, int number_atoms, double** part) {
	for(int i = 0; i < number_atoms-start_atom+1; i++) {
		for(int j = 0; j < 3; j++) {
			part[i][j]=co[i+start_atom-1][j];
		}
	}

}

double rmsd_anchor3(double** anchor_flexible,double** anchor_fixed)	//:TRICKY:	variant using all three anchor atoms
{
	double sum = 0.0;
	for (int i = 0; i < 3; i++)
	{	for (int j = 0; j < 3; j++)
	{
		sum += pow(anchor_fixed[i][j] - anchor_flexible[i][j],2);
	}
	}

	return sqrt(sum/3);
}


double rmsd_anchor2(double** anchor_flexible,double** anchor_fixed)	//:TRICKY:	variant using last two anchor atoms (CA and C -> forward; N and CA -> reverse)
{
	double sum = 0.0;
	for (int i = 1; i < 3; i++)
	{	for (int j = 0; j < 3; j++)
	{
		sum += pow(anchor_fixed[i][j] - anchor_flexible[i][j],2);
	}
	}

	return sqrt(sum/2);
}


void rand_shift(double **co,int nr_atoms,double rmsd) {
	double translate[3];

	// Calculate translate coordinates
	for (int i = 0; i < 3; i++) {
		translate[i] = ((double) rand()/(double) RAND_MAX) * rmsd/sqrt(3.0);
	}

	// New coordinates co
	for(int i = 0; i < nr_atoms; i++) {
		for(int j = 0; j < 3; j++) {
			co[i][j] 	+= translate[j];

		}
	}

}

/*===============================================================================
 *
 *	OPT_SHIFT.CPP
 *
 *	Created on: February 10, 2010
 *      Author: Pieter Chys
 *
 *===============================================================================
 *
 * OPT_SHIFT.CPP calculates the optimal rigid body translation
 * for a protein loop and executes it on the active coordinates.
 * The new coordinates are returned in the original matrix. The
 * algorithm is along simple analytical lines of coordinate
 * minimization with a variable translation coordinate. Basically,
 * the difference coordinates between the fixed anchor and
 * real anchor coordinates must be made and be divided by 6.
 * This yields the translation coordinates in one dimension. This
 * must be done in 3 dimensions after which the translation coordinates
 * are added to the real coordinates.
 *
 * PREALLOCATION OF difference[6][3] is still possible to speed up code
 * if wanted!
 *
 *--------------------------------------------------------------------------------
 *
 *  INPUT: - double **pointer to the fixed begin anchor coordinates (3 x 3 matrix) ANCHOR_BEGIN
 *  	   - double **pointer to the fixed end anchor coordinates (3 x 3 matrix) ANCHOR_END
 *  	   - double **pointer to the loop coordinates (N x 3 matrix) CO
 *  	   - int nr_atoms_loop which represents the number of atoms in the loop
 *
 *-------------------------------------------------------------------------------
 *
 *  OUTPUT: - explicitly void
 *  		- implicitly double **pointer to the new loop coordinates CO
 *
 *
 *-------------------------------------------------------------------------------
 *
 *  DEPENDENCIES:
 *
 *===============================================================================
 */

void opt_shift(double **anchor_begin,double **anchor_end, double **co,int nr_atoms_loop) {
	double translate[3];
	double difference[6][3];

	//	Initialize difference vectors
	for(int i = 0; i < 3; i++) {
		for(int j = 0; j < 3; j++) {
			difference[i][j] 	= anchor_begin[i][j]-	co[i][j];
			difference[3+i][j]	= anchor_end[i][j]	-	co[nr_atoms_loop-3 +i][j];
		}
	}

	// Calculate translate coordinates
	for (int i = 0; i < 3; i++) {
		translate[i] = 0.0;
		for (int j = 0; j < 6; j++) {
			translate[i] += difference[j][i];
		}
		translate[i] /= 6.0;
	}

	// New coordinates co
	for(int i = 0; i < nr_atoms_loop; i++) {
		for(int j = 0; j < 3; j++) {
			co[i][j] 	+= translate[j];

		}
	}

}





//OLD CODES: SOME WORK PROPERLY

//  CODE TO ROTATE TWO QUATERNIONS

//	quaternion[0] =	qt2[0]*qt1[0] - qt2[1]*qt1[1] - qt2[2]*qt1[2] - qt2[3]*qt1[3];
//	quaternion[1] =	qt2[0]*qt1[1] + qt2[1]*qt1[0] + qt2[2]*qt1[3] - qt2[3]*qt1[2];
//	quaternion[2] =	qt2[0]*qt1[2] - qt2[1]*qt1[3] + qt2[2]*qt1[0] + qt2[3]*qt1[1];
//	quaternion[3] =	qt2[0]*qt1[3] + qt2[1]*qt1[2] - qt2[2]*qt1[1] + qt2[3]*qt1[0];





//WORKING INTEGRATED BACKWARD SPINORCODE
//void ga_rotate_anglesm5(double **co,				/*loop coordinates*/
//						double **bond,
//						double **bondn,
//						const double *angles,		/*angle variations from optimization*/
//						int *marker,
//						const int *bondnumbers,		/*physical bond numbers*/
//						const int bondrange,		/*# of bonds rotating*/
//						const int number_atoms,		/*# atoms in loop*/
//						double *between,			/*auxiliary spinor*/
//						double *spinor1,			/*auxiliary spinor 1-preallocated*/
//						double *spinor2,			/*auxiliary spinor 2-preallocated*/
//						double *spinortot			/*total spinor-preallocated*/
//						)
//{
//
//	int switch_debug = 0;						/*switch for debug information*/
//	message("Start allocation...",0);
//	message("Allocation passed...",switch_debug);
//
//
////	for(int j = 0; j < 3; j++)
////	{
////		bond[0][j] = 0.0;
////		bondn[0][j]= 0.0;
////	}
////
////	for(int i = 1; i < number_atoms; i++)
////	{	for(int j = 0; j < 3; j++)
////		{
////			bond[i][j] = co[i][j] - co[i-1][j]; 	/*bond[0] is not written*/
////			bond[i][3] = 0.0;
////			bondn[i][3]= 0.0;
////			normalize(bond[i],bondn[i]);
////		}
////	}
//
//
////	show3("Initial loop coordinates","Coordinates",co,number_atoms);
////	show3("Bond coordinates","Coordinates",bond,number_atoms);
////	show3("Normalized bond coordinates","Coordinates",bondn,number_atoms);
//
//
//
//	for (int i = 3; i <  number_atoms; i++){	/*Correction dihedral angles in [-PI,PI] range*/
//		if (i != 3 && marker[i-1] == 1 )
//		{
//			spinor(bondn[i-1],angles[i-1],spinor1);	/*still not tuned correctly*/
//			spinorproduct(spinortot,spinor1,spinortot);
//			revspinor(spinortot,spinor2);
//		}
//		else if (i == 3)
//		{
//			spinor(bondn[i-1],angles[i-1],spinortot);	/*still not tuned correctly*/
//			revspinor(spinortot,spinor2);
//		}
//		geom32opt(bond[i],spinortot,between);
//		geom23opt(spinor2,between,bond[i]);
//		normalize(bond[i],bondn[i]);
//	}
//
//	//show3("Updated bond coordinates","Coordinates",bond,number_atoms);
//
//
//	for(int i = 3; i < number_atoms; i++)
//	{	for(int j = 0; j < 3; j++)
//		{
//			co[i][j] = co[i-1][j] + bond[i][j]; 	/*bond[0] is not written*/
//		}
//	}
//
//	//show3("Returned loop coordinates","Coordinates",co,number_atoms);
//	//getchar();
//}



//BACKUP ORIGINAL CODE BESTANGLE_ANCHOR_NEW	21/09/2012
//WITH DIRTY OLD COMMENTING AND OLD CODE LINES

//double bestangle_anchor_new(double** co,int bond_number,int number_atoms,double** anchor,						//	Preallocated version
//			double** matrix,double** polpoint,double** anchorpoint,double** result1,double** result2, double** result3) {
//
//	double angle;							//	optimal angle we want
//
//	if (bond_number == number_atoms - 2)						//PROTOCOL FOR IN-ANCHOR DEGREE OF FREEDOM
//	{
//		double t;								//	normal increment == length normal
//		double radius;							//	radius last bond rotation in circle plane
//		double distance;						//	distance  projection circle plane and point on circle
//		double length3;							//	length last bond
//		double angle3;							//	valence angle last bond
//		double dotpr;							//	dot product from angle checking
//
//		double cotogo[3];						//	ideal point to go
//		double cross[3];						//	normal vector to check dihedral angle sign
//		double normal[3];						//	normal vector to circle plane
//		double fixed[3];						//	last point of anchor
//		double tocenter[3];						//	vector that goes to center
//		double center[3];						//	center of circle last bond
//		double bond2[3];						//	bond to rotate around
//		double bond3[3];						//  bond that will be rotated
//		double projection[3];					//  point projected on circle plane (fixed)
//		double difference[3];					//	difference projected and circle point
//		double x[3];							//  coordinate optimal point on circle
//		double binormal[3];						//	coordinate circle coordinate actual bond
//		double localsum;						//	temporary variable
//
//		//	ACTIVE ANCHOR COORDINATES
//		diffvectorseq(co,bond_number + 1, bond2);			// indices fine	// first bond active anchor
//		diffvectorseq(co,bond_number + 2, bond3);			// indices fine	// second bond active anchor
//		normalize(bond2,normal);
//
//		length3 	= 	norm(bond3);
//		angle3 		= 	acos(dotprod(bond2,bond3)/(norm(bond2)*norm(bond3)));	//	correct angle ... smaller/bigger than 90 degrees	// valence angle between bond 1 and 2 in anchor --- POSSIBLE BUG???
//		radius 		=	sin(angle3)*length3;
//		// DEBUGGING FUNCTION
//		//printf("%-50s %8.3f\n","Valence angle in anchor is",angle3);			//	seems Ok
//		//printf("-50s %8.3f\n", "Radius in anchor is",radius);					//	seems Ok
//		//fflush(stdout);
//
//		t 	= 0.0;
//		localsum = 0.0;
//		for(int i = 0; i < 3; i++)
//		{
//			// CAN BE SPEED UP LATER IN SINGLE LINES
//			tocenter[i]	=	length3*cos(angle3)*normal[i];
//			center[i] 	=	co[bond_number][i] + tocenter[i];		//	changes // seems Ok until here ...
//			binormal[i] =   co[bond_number + 1][i] - center[i];		//	changes
//
//			fixed[i]  	= 	anchor[2][i];							//	seems Ok ...
//			//t     		+=	normal[i]*(center[i] - fixed[i]);
//			t     		+=	normal[i]*(fixed[i]-center[i]);			//	BUG: SHOULD BE REVERSED
//			localsum	+= 	normal[i];								//	NEEDS TO BE ADDED TO LET IT WORK!!!!! BUG BEFORE ....
//		}
//		t				/= localsum;								//	NEEDS TO BE ADDED TO LET IT WORK!!!!! BUG BEFORE ...
//		//printf("%-50s %8.3f\n","Parameter t",t);
//		//fflush(stdout);
//
//
//		for(int i = 0; i < 3; i++){
//			//projection[i] = fixed[i] + t*normal[i];
//			projection[i]	= center[i] + t*normal[i];				//	CORRECTED!!!!!! BUG BEFORE....
//			//difference[i] = projection[i] - center[i];
//			difference[i] = fixed[i] - projection[i];
//		}
//
//		distance = norm(difference);
//
//		for(int i = 0; i < 3; i++)
//		{
//			x[i] 	  =  radius*difference[i]/distance;	// projection vector in A parallel plane !!!!! DEBUGGING ...
//			//cotogo[i] =  center[i]	+  x[i];
//		}
//
//		//	DEBUGGING
//		/*printf("\n");
//		printf("%50s %8.3f %8.3f %8.3f\n","Original end coordinate",co[number_atoms - 1][0],co[number_atoms - 1][1],co[number_atoms - 1][2]);
//		printf("%50s %8.3f %8.3f %8.3f\n","Ideal end coordinate",cotogo[0],cotogo[1],cotogo[2]);
//		printf("%50s %8.3f %8.3f %8.3f\n","Actual radius vector",binormal[0],binormal[1],binormal[2]);*/
//		//printf("%50s %8.3f %8.3f %8.3f\n","Ideal radius vector",x[0],x[1],x[2]);
//		//fflush(stdout);
//
//
//		crossprod(bond3,bond2,cross);								//	calculating cross product
//		dotpr = dotprod(cross,x);									//  check value dot product
//
//		angle = acos(dotprod(x,binormal)/(norm(x)*norm(binormal)));	//	how about correction for angle????
//
//		if (dotpr > 0){
//			angle = - fabs(angle);
//		} else if (dotpr < 0) {
//			angle = + fabs(angle);
//		}
//
//
//		//	DEBUGGING
//		/*printf("%50s %8.3f\n","angle to move",angle);
//		printf("\n");
//		fflush(stdout);*/
//
//	} else {												//PROTOCOL FOR  REGULAR BONDS
//
//		double bond2[3];																			//	bond around which we rotate
//
//		//	Actual algorithm
//		for (int i = 0; i < 3; i++)
//		{																//	make correct active anchor coordinates
//			subtract(co[bond_number-1],co[number_atoms-3+i],polpoint[i]);							//	crucial_____XXXXX OLD CODE
//																									// segmentation fault can happen
//		}
//
//		for (int i = 0; i < 3; i++)
//		{																//	make correct fixed anchor coordinates
//			subtract(co[bond_number-1],anchor[i],anchorpoint[i]);									//	in translated frame ... OLD ONE
//		}
//
//		diffvectorseq(co,bond_number + 1,bond2);
//
//		rotatematrix(bond2,matrix);
//
//		//rotatematrix(bond2,matrix); 										//	But mind bond1 ABOVE (not dealt with) -> bond_number 1 should be avoided
//		for (int i = 0; i < 3; i++) {										//	apply to each of the anchor points the previous procedure
//			matrixXvector(matrix,polpoint[i],result1[i]);				//	problem here somewhere....
//			matrixXvector(matrix,anchorpoint[i],result2[i]);
//		}
//
//		//	MATRIX ALGEBRA APPROACH
//		for (int i = 0; i < 3; i++)
//		{
//			result2[i][0] = 0.0;												//	setting x coordinate at zero
//			result1[i][0] = 0.0;
//		}
//
//		double denominator;
//		denominator = 0;
//		for (int i = 0; i < 3; i++)
//		{
//			denominator += dotprod(result1[i],result2[i]);
//		}
//
//
//		for (int i = 0; i < 3; i++) {						//	What is this protocol exactly ? -> Does not seem to fit with theory  -> is Ok
//			for (int j = 0; j < 3; j++) {
//				result3[i][0] = 0.0;						//	Construct orthonormal with cross product of result1 and [1 0 0] and
//				result3[i][1] = -result1[i][2];				//	work out analytically
//				result3[i][2] = result1[i][1];
//			}
//		}
//
//		double numerator = 0;
//		for (int i = 0; i < 3; i++) {
//			numerator += dotprod(result2[i],result3[i]);
//		}
//
//		angle = atan2(numerator,denominator);
//	}
//
//	return angle;
//}

