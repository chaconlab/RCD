/*
 * geometry.h
 *
 *  Created on: Oct 26, 2013
 *      Author: pablo
 */
//
//  geometry.h
//
//
//  Created by Pablo on 11/9/12.
//
//

#ifndef _geometry_h
#define _geometry_h


//INTERNAL COORDINATE CONSTANTS
/********************************************************************************
 *
 *			1st set of internal coordinates based on reference:
 *
 *
 *		 Fersht Alan, Enzyme Structure and Mechanism, 2nd edition,
 *		 1985, W.H. Freeman and Company, New York
 *
 *
 *
 ********************************************************************************/



//BOND LENGTHS (Angstrom)
//#define bondNCA2 	 		1.46			/*N-CA bond  */
//#define bondCACO2    		1.51			/*CA-CO bond */
//#define bondCON2 	 		1.33	  		/*CO-N bond  */

//#define bondCC				1.54
//#define bondCO				1.24
//#define bondNH				1.02
//#define bondCH				1.10
//#define bondCN				1.488

//VALENCE ANGLES (rad)
//#define valCONCA2   			122*M_PI/180		/*valence angle CO-N-CA*/
//#define valNCACO2   			111*M_PI/180		/*valence angle N-CA-CO*/
//#define valCACON2   			116*M_PI/180		/*valence angle CA-CO-N*/

//#define valNCH	 			116*M_PI/180		/*average value chosen: 116 degrees*/
//#define valNCC	 			116*M_PI/180
//#define valCCO	 			116*M_PI/180
//#define valCNH	 			116*M_PI/180
//#define valCCC	 			0.667*M_PI
//#define valCCH	 			0.667*M_PI
//#define valCCN	 			0.667*M_PI	        /*based on CCH and CCC; must be around 0.65 rad*/

//DIHEDRAL ANGLES (rad)
//#define dihCON				M_PI			/*dihedral for peptide plane = 180 degrees*/




/********************************************************************************
 *
 *	2nd set of internal coordinates based on reference:
 *	Engh Richar A. and Huber Robert, Accurate bond and angle
 * 	for X-ray protein structure refinement, Acta Cryst, 1991, A47,
 *       392-400
 *	Remark:	- Used as result of need to compare with Soto paper
 *	- Sigma's also used as perturbations later on
 *	- Check also Dunbrack values and variations
 *
 ********************************************************************************/


/*Ramachandran valid angle domains*/
#define LOWERPHI -175*M_PI/180
#define UPPERPHI  -40*M_PI/180
#define LOWERPSI -60*M_PI/180
#define UPPERPSI  175*M_PI/180



// REGULAR RESIDUES (EXCEPT GLYCINE AND PROLINE)
//
// Values taken from the most populated values (manually) of Dumbrack's maps.
// http://dunbrack.fccc.edu/nmhrcm/KernRegr_CDL_v1.2.1_Oct14-2009.txt
//
// mon@yolem:~/PROCESS/rcd/CDL$ more summary.txt
//
// --------
// GENERAL:
// --------
//NonPGIV_nonxpro.txt
//# ResTypeGroup    Phi  Psi S   Num   mCNA sCNA   mNAB sNAB   mNAC sNAC   mBAC sBAC   mACO sACO   mACN sACN   mOCN sOCN    mCN    sCN    mNA    sNA    mAB    sAB    mAC    sAC    mCO    sCO
//#Avg(1296):        -5   -5 0  8360 121.69 1.85 110.57 1.67 110.73 1.98 110.57 1.93 120.53 1.38 116.83 1.63 122.57 1.32 1.3322 0.0140 1.4576 0.0129 1.5302 0.0168 1.5235 0.0133 1.2352 0.0127
//#Sig(1296):       104  104 0  4608   0.54 0.17   0.42 0.13   0.98 0.30   0.77 0.17   0.49 0.15   0.63 0.18   0.27 0.12 0.0013 0.0010 0.0016 0.0009 0.0020 0.0020 0.0016 0.0012 0.0011 0.0011
//
// BOND LENGTHS (Angstrom)
// CO-N bond  #  CN    C(-1) - N(0)
#define bondCON 	 	1.3322
#define sigma_bondCON 	0.0140
// N-CA bond  #  NA    N(0)  - Ca(0)
#define bondNCA 	 	1.4576
#define sigma_bondNCA 	0.0129
// CA-CO bond #  AC    Ca(0) - C(0)
#define bondCACO    	1.5235
#define sigma_bondCACO  0.0133
#define bondCACO2  		3.0470
// CA-CB bond #  AB    Ca(0) - Cb(0)
#define bondCACB    	1.5302
#define sigma_bondCACB  0.0168
#define bondCACB2    	3.0604
// CO=O bond  #  CO    C(0) - O(0)
#define bondCO  		1.2352
#define sigma_bondCO2  	0.0127
#define bondCO2  		2.4704
//
// VALENCE ANGLES (rad)
// valence angle CO-N-CA   # CNA   C(-1) - N(0)  - Ca(0)
#define valCONCA   		121.69*M_PI/180
#define sigma_valCONCA  1.85*M_PI/180
// valence angle N-CA-CO   # NAC   N(0)  - Ca(0) - C(0)
#define valNCACO   		110.73*M_PI/180
#define sigma_valNCACO  1.98*M_PI/180
// valence angle CA-CO-N   # ACN   Ca(0) - C(0)  - N(+1)
#define valCACON   		116.83*M_PI/180
#define sigma_valCACON  1.63*M_PI/180
// valence angle O-CO-N    # OCN   O(0)  - C(0)  - N(+1)
#define valOCON   		122.57*M_PI/180
#define sigma_valOCON   1.32*M_PI/180
// valence NAB   N(0)  - Ca(0) - Cb(0)
#define valNCACB   		110.57*M_PI/180
#define sigma_valNCACB    1.67*M_PI/180
// valence ACO   Ca(0) - C(0)  - O(0)
#define valCACOO   		120.53*M_PI/180
#define sigma_valCACOO    1.38*M_PI/180


#define COSvalNCACB     cos( (M_PI - valNCACB)/2.0 );
#define COSvalCACOO     cos( (M_PI - valCACOO)/2.0 );
#define COSPROvalNCACB  cos( (M_PI - PROvalNCACB)/2.0 );
#define COSvalNCACB     cos( (M_PI - valNCACB)/2.0 );

#define SINvalNCACB     sin( (M_PI - valNCACB)/2.0 );
#define SINvalCACOO     sin( (M_PI - valCACOO)/2.0 );
#define SINPROvalNCACB  sin( (M_PI - PROvalNCACB)/2.0 );
#define SINvalNCACB     sin( (M_PI - valNCACB)/2.0 );

#define COS2MPI3 cos(2.0*M_PI/3.0);
#define SIN2MPI3 sin(2.0*M_PI/3.0);
#define SQRT32 sqrt(3.0)/2.0
//
// --------
// ILE/VAL:
// --------
//IleVal_nonxpro.txt
//# ResTypeGroup   Phi  Psi S Num    mCNA sCNA   mNAB sNAB   mNAC sNAC   mBAC sBAC   mACO sACO   mACN sACN   mOCN sOCN    mCN    sCN    mNA    sNA    mAB    sAB    mAC    sAC    mCO    sCO
//#Avg(1296)        -5   -5 0 1634 121.97 1.76 111.25 1.64 109.39 2.00 111.29 1.62 120.77 1.23 116.62 1.43 122.56 1.24 1.3319 0.0135 1.4590 0.0125 1.5401 0.0135 1.5230 0.0127 1.2361 0.0119
//#Sig(1296)       104  104 0  550   0.30 0.14   0.24 0.13   0.67 0.25   0.18 0.16   0.25 0.06   0.34 0.09   0.15 0.07 0.0008 0.0010 0.0009 0.0005 0.0015 0.0006 0.0007 0.0005 0.0009 0.0005
//
// BOND LENGTHS (Angstrom)
// CO-N bond  #  CN    C(-1) - N(0)
#define IVbondCON 	 			1.3319
#define sigma_IVbondCON 	 	0.0135
// N-CA bond  #  NA    N(0)  - Ca(0)
#define IVbondNCA 	 			1.4590
#define sigma_IVbondNCA 	 	0.0125
// CA-CO bond #  AC    Ca(0) - C(0)
#define IVbondCACO    			1.5230
#define sigma_IVbondCACO    	0.0127
// CA-CB bond #  AB    Ca(0) - Cb(0)
#define IVbondCACB    			1.5401
#define sigma_IVbondCACB  		0.0135
#define IVbondCACB2    			3.0802
// CO=O bond  #  CO    C(0)  - O(0)
#define IVbondCO  				1.2361
#define sigma_IVbondCO2  		0.0119
#define IVbondCO2  				2.4722
//
// VALENCE ANGLES (rad)
// valence angle CO-N-CA   # CNA   C(-1) - N(0)  - Ca(0)
#define IVvalCONCA   			121.97*M_PI/180
#define sigma_IVvalCONCA   		  1.76*M_PI/180
// valence angle N-CA-CO   # NAC   N(0)  - Ca(0) - C(0)
#define IVvalNCACO   			109.39*M_PI/180
#define sigma_IVvalNCACO   		  2.00*M_PI/180
// valence angle CA-CO-N   # ACN   Ca(0) - C(0)  - N(+1)
#define IVvalCACON   			116.62*M_PI/180
#define sigma_IVvalCACON  	 	  1.43*M_PI/180
// valence angle O-CO-N    # OCN   O(0)  - C(0)  - N(+1)
#define IVvalOCON   			122.56*M_PI/180
#define sigma_IVvalOCON  		  1.24*M_PI/180
//
// --------
// GLYCINE:
// --------
//Gly_nonxpro.txt
//#ResTypeGroup Phi  Psi S  Num   mCNA sCNA  -  -   mNAC sNAC  -  -   mACO sACO   mACN sACN  mOCN  sOCN    mCN    sCN    mNA    sNA  -  -    mAC    sAC    mCO    sCO
//#Avg(1296)     -5   -5 0 1131 121.45 1.89 -1 -1 113.06 2.22 -1 -1 120.66 1.63 116.56 1.90 122.74 1.27 1.3303 0.0144 1.4492 0.0143 -1 -1 1.5145 0.0139 1.2345 0.0134
//#Sig(1296):   104  104 0  527   0.39 0.27  0  0   0.74 0.35  0  0   0.53 0.26   0.67 0.31   0.23 0.12 0.0017 0.0013 0.0016 0.0013  0  0 0.0016 0.0013 0.0013 0.0011
//
// BOND LENGTHS (Angstrom)
// CO-N bond  #  CN    C(-1) - N(0)
#define GLYbondCON 	 			1.3303
#define sigma_GLYbondCON 	 	0.0144
// N-CA bond  #  NA    N(0)  - Ca(0)
#define GLYbondNCA 	 			1.4492
#define sigma_GLYbondNCA 	 	0.0143
// CA-CO bond #  AC    Ca(0) - C(0)
#define GLYbondCACO    			1.5145
#define sigma_GLYbondCACO    	0.0139
// CO=O bond  #  CO    C(0) - O(0)
#define GLYbondCO  				1.2345
#define sigma_GLYbondCO2  		0.0134
#define GLYbondCO2  			2.4690
//
// VALENCE ANGLES (rad)
// valence angle CO-N-CA   # CNA   C(-1) - N(0)  - Ca(0)
#define GLYvalCONCA   			121.45*M_PI/180
#define sigma_GLYvalCONCA   	1.89*M_PI/180
// valence angle N-CA-CO   # NAC   N(0)  - Ca(0) - C(0)
#define GLYvalNCACO   			113.06*M_PI/180
#define sigma_GLYvalNCACO   	2.22*M_PI/180
// valence angle CA-CO-N   # ACN   Ca(0) - C(0)  - N(+1)
#define GLYvalCACON   			116.56*M_PI/180
#define sigma_GLYvalCACON   	1.90*M_PI/180
// valence angle O-CO-N    # OCN   O(0)  - C(0)  - N(+1)
#define GLYvalOCON   			122.74*M_PI/180
#define sigma_GLYvalOCON  		1.27*M_PI/180
//
// --------
// PROLINE:
// --------
//Pro_nonxpro.txt
//#ResTypeGroup Phi  Psi S Num   mCNA sCNA   mNAB sNAB   mNAC sNAC   mBAC sBAC   mACO sACO   mACN sACN   mOCN sOCN    mCN    sCN    mNA    sNA    mAB    sAB    mAC    sAC    mCO    sCO
//#Avg(1296):    -5   -5 0 612 119.85 1.25 103.24 1.05 112.47 2.04 111.55 1.64 120.61 1.80 116.70 2.04 122.64 1.35 1.3341 0.0234 1.4687 0.0128 1.5332 0.0142 1.5196 0.0141 1.2351 0.0130
//#Sig(1296):   104  104 0 127   0.15 0.08   0.10 0.04   0.30 0.12   0.18 0.06   0.29 0.12   0.32 0.15   0.09 0.05 0.0031 0.0076 0.0008 0.0004 0.0003 0.0005 0.0008 0.0006 0.0005 0.0005
//
// BOND LENGTHS (Angstrom)
// CO-N bond  #  CN    C(-1) - N(0)
#define PRObondCON 	 			1.3341
#define sigma_PRObondCON 	 	0.0234
// N-CA bond  #  NA    N(0)  - Ca(0)
#define PRObondNCA 	 			1.4687
#define sigma_PRObondNCA 	 	0.0128
// CA-CO bond #  AC    Ca(0) - C(0)
#define PRObondCACO    			1.5196
#define sigma_PRObondCACO    	0.0141
// CA-CB bond #  AB    Ca(0) - Cb(0)
#define PRObondCACB    			1.5332
#define sigma_PRObondCACB  		0.0142
#define PRObondCACB2    		3.0664
// CO=O bond  #  CO    C(0) - O(0)
#define PRObondCO  				1.2351
#define sigma_PRObondCO2  		0.0130
#define PRObondCO2  			2.4702
//
// VALENCE ANGLES (rad)
// valence angle CO-N-CA   # CNA   C(-1) - N(0)  - Ca(0)
#define PROvalCONCA   			119.85*M_PI/180
#define sigma_PROvalCONCA   	1.25*M_PI/180
// valence angle N-CA-CO   # NAC   N(0)  - Ca(0) - C(0)
#define PROvalNCACO   			112.47*M_PI/180
#define sigma_PROvalNCACO   	2.04*M_PI/180
// valence angle CA-CO-N   # ACN   Ca(0) - C(0)  - N(+1)
#define PROvalCACON   			116.70*M_PI/180
#define sigma_PROvalCACON   	2.04*M_PI/180
// valence angle O-CO-N    # OCN   O(0)  - C(0)  - N(+1)
#define PROvalOCON   			122.64*M_PI/180
#define sigma_PROvalOCON  		1.35*M_PI/180
// valence NAB   N(0)  - Ca(0) - Cb(0)
#define PROvalNCACB   			103.24*M_PI/180
#define sigma_PROvalNCACB   	 1.05*M_PI/180
//
//
// ------------------
// GENERAL "pre-PRO":
// ------------------
//NonPGIV_xpro.txt
//#ResTypeGroup  Phi  Psi S Num   mCNA sCNA   mNAB sNAB   mNAC sNAC   mBAC sBAC   mACO sACO   mACN sACN   mOCN sOCN    mCN    sCN    mNA    sNA    mAB    sAB    mAC    sAC    mCO    sCO
//#Avg(1296):     -5   -5 0 470 121.82 2.40 110.38 1.77 109.80 2.17 110.17 1.93 120.16 1.34 118.44 1.55 121.32 1.14 1.3310 0.0204 1.4569 0.0140 1.5298 0.0157 1.5237 0.0125 1.2378 0.0127
//#Sig(1296):    104  104 0 137   0.29 0.33   0.17 0.11   0.39 0.19   0.30 0.15   0.22 0.10   0.27 0.15   0.11 0.10 0.0011 0.0031 0.0012 0.0011 0.0012 0.0009 0.0009 0.0005 0.0010 0.0007
//
// BOND LENGTHS (Angstrom)
// CO-N bond  #  CN    C(-1) - N(0)
#define pPRO_bondCON 	 	1.3310
#define pPRO_sigma_bondCON 	0.0204
// N-CA bond  #  NA    N(0)  - Ca(0)
#define pPRO_bondNCA 	 	1.4569
#define pPRO_sigma_bondNCA 	0.0140
// CA-CO bond #  AC    Ca(0) - C(0)
#define pPRO_bondCACO    	1.5237
#define pPRO_sigma_bondCACO 0.0125
#define pPRO_bondCACO2  	3.0474
// CA-CB bond #  AB    Ca(0) - Cb(0)
#define pPRO_bondCACB    	1.5298
#define pPRO_sigma_bondCACB 0.0157
#define pPRO_bondCACB2    	3.0596
// CO=O bond  #  CO    C(0) - O(0)
#define pPRO_bondCO  		1.2378
#define pPRO_sigma_bondCO2  0.0127
#define pPRO_bondCO2  		2.4756
//
// VALENCE ANGLES (rad)
// valence angle CO-N-CA   # CNA   C(-1) - N(0)  - Ca(0)
#define pPRO_valCONCA   	121.82*M_PI/180
#define pPRO_sigma_valCONCA   2.40*M_PI/180
// valence angle N-CA-CO   # NAC   N(0)  - Ca(0) - C(0)
#define pPRO_valNCACO   	109.80*M_PI/180
#define pPRO_sigma_valNCACO   2.17*M_PI/180
// valence angle CA-CO-N   # ACN   Ca(0) - C(0)  - N(+1)
#define pPRO_valCACON   	118.44*M_PI/180
#define pPRO_sigma_valCACON   1.55*M_PI/180
// valence angle O-CO-N    # OCN   O(0)  - C(0)  - N(+1)
#define pPRO_valOCON   		121.32*M_PI/180
#define pPRO_sigma_valOCON    1.14*M_PI/180
//
// ------------------
// ILE/VAL "pre-PRO":
// ------------------
//IleVal_xpro.txt
//# ResType     Phi  Psi S Num  mCNA sCNA   mNAB sNAB   mNAC sNAC   mBAC sBAC   mACO sACO   mACN sACN   mOCN sOCN    mCN    sCN    mNA    sNA    mAB    sAB    mAC    sAC    mCO    sCO
//#Avg(1296):    -5   -5 0 80 122.13 1.84 111.21 1.39 108.88 2.14 111.36 1.81 119.95 1.33 118.88 1.53 121.10 1.14 1.3292 0.0113 1.4612 0.0123 1.5372 0.0129 1.5245 0.0105 1.2397 0.0126
//#Sig(1296):   104  104 0 10   0.17 0.10   0.09 0.05   0.21 0.14   0.17 0.10   0.11 0.07   0.15 0.08   0.06 0.04 0.0004 0.0003 0.0005 0.0004 0.0007 0.0005 0.0004 0.0003 0.0006 0.0003
//
// BOND LENGTHS (Angstrom)
// CO-N bond  #  CN    C(-1) - N(0)
#define pPRO_IVbondCON 	 			1.3292
#define pPRO_sigma_IVbondCON 	 	0.0113
// N-CA bond  #  NA    N(0)  - Ca(0)
#define pPRO_IVbondNCA 	 			1.4612
#define pPRO_sigma_IVbondNCA 	 	0.0125
// CA-CO bond #  AC    Ca(0) - C(0)
#define pPRO_IVbondCACO    			1.5230
#define pPRO_sigma_IVbondCACO    	0.0123
// CA-CB bond #  AB    Ca(0) - Cb(0)
#define pPRO_IVbondCACB    			1.5372
#define pPRO_sigma_IVbondCACB  		0.0129
#define pPRO_IVbondCACB2    		3.0744
// CO=O bond  #  CO    C(0)  - O(0)
#define pPRO_IVbondCO  				1.2397
#define pPRO_sigma_IVbondCO2  		0.0126
#define pPRO_IVbondCO2  			2.4794
//
// VALENCE ANGLES (rad)
// valence angle CO-N-CA   # CNA   C(-1) - N(0)  - Ca(0)
#define pPRO_IVvalCONCA   			122.13*M_PI/180
#define pPRO_sigma_IVvalCONCA   	  1.84*M_PI/180
// valence angle N-CA-CO   # NAC   N(0)  - Ca(0) - C(0)
#define pPRO_IVvalNCACO   			108.88*M_PI/180
#define pPRO_sigma_IVvalNCACO   	  2.14*M_PI/180
// valence angle CA-CO-N   # ACN   Ca(0) - C(0)  - N(+1)
#define pPRO_IVvalCACON   			118.88*M_PI/180
#define pPRO_sigma_IVvalCACON  	 	  1.53*M_PI/180
// valence angle O-CO-N    # OCN   O(0)  - C(0)  - N(+1)
#define pPRO_IVvalOCON   			121.10*M_PI/180
#define pPRO_sigma_IVvalOCON  		  1.14*M_PI/180
//
// ------------------
// GLYCINE "pre-PRO":
// ------------------
//Gly_xpro.txt
//#ResTypeup Phi  Psi S Num  mCNA sCNA  -  -   mNAC sNAC  -  -   mACO sACO   mACN sACN  mOCN  sOCN    mCN    sCN    mNA    sNA  -  -    mAC    sAC    mCO    sCO
//#Avg(1296)  -5   -5 0 44 121.87 1.57 -1 -1 112.34 2.04 -1 -1 121.52 1.43 116.69 1.71 121.77 1.00 1.3304 0.0149 1.4426 0.0137 -1 -1 1.5146 0.0143 1.2351 0.0135
//#Sig(1296) 104  104 0  3   0.04 0.04  0  0   0.11 0.05  0  0   0.08 0.04   0.11 0.04   0.03 0.03 0.0004 0.0004 0.0003 0.0002  0  0 0.0002 0.0003 0.0002 0.0002
//
// BOND LENGTHS (Angstrom)
// CO-N bond  #  CN    C(-1) - N(0)
#define pPRO_GLYbondCON 	 	1.3304
#define pPRO_sigma_GLYbondCON 	0.0149
// N-CA bond  #  NA    N(0)  - Ca(0)
#define pPRO_GLYbondNCA 	 	1.4426
#define pPRO_sigma_GLYbondNCA 	0.0137
// CA-CO bond #  AC    Ca(0) - C(0)
#define pPRO_GLYbondCACO    	1.5146
#define pPRO_sigma_GLYbondCACO    	0.0143
// CO=O bond  #  CO    C(0) - O(0)
#define pPRO_GLYbondCO  		1.2351
#define pPRO_sigma_GLYbondCO2  	0.0135
#define pPRO_GLYbondCO2  		2.4702
//
// VALENCE ANGLES (rad)
// valence angle CO-N-CA   # CNA   C(-1) - N(0)  - Ca(0)
#define pPRO_GLYvalCONCA   			121.87*M_PI/180
#define pPRO_sigma_GLYvalCONCA        1.57*M_PI/180
// valence angle N-CA-CO   # NAC   N(0)  - Ca(0) - C(0)
#define pPRO_GLYvalNCACO   			112.34*M_PI/180
#define pPRO_sigma_GLYvalNCACO   	  2.04*M_PI/180
// valence angle CA-CO-N   # ACN   Ca(0) - C(0)  - N(+1)
#define pPRO_GLYvalCACON   			116.69*M_PI/180
#define pPRO_sigma_GLYvalCACON   	  1.71*M_PI/180
// valence angle O-CO-N    # OCN   O(0)  - C(0)  - N(+1)
#define pPRO_GLYvalOCON   			121.77*M_PI/180
#define pPRO_sigma_GLYvalOCON  		  1.00*M_PI/180
//
// ------------------
// PROLINE "pre-PRO":
// ------------------
//Pro_xpro.txt
//#ResType   Phi  Psi S Num  mCNA sCNA   mNAB sNAB   mNAC sNAC   mBAC sBAC   mACO sACO   mACN sACN   mOCN sOCN    mCN    sCN    mNA    sNA    mAB    sAB    mAC    sAC    mCO    sCO
//#Avg(1296)  -5   -5 0 12 120.38 1.03 103.08 0.97 110.70 1.22 110.92 1.22 120.56 1.45 117.93 1.20 121.46 1.18 1.3292 0.0118 1.4649 0.0203 1.5357 0.0137 1.5171 0.0093 1.2404 0.0112
//#Sig(1296) 104  104 0  0   0.03 0.01   0.00 0.02   0.01 0.01   0.01 0.01   0.01 0.03   0.01 0.01   0.01 0.03 0.0002 0.0001 0.0003 0.0003 0.0001 0.0002 0.0001 0.0001 0.0002 0.0001
//
// BOND LENGTHS (Angstrom)
// CO-N bond  #  CN    C(-1) - N(0)
#define pPRO_PRObondCON 	 		1.3292
#define pPRO_sigma_PRObondCON 	 	0.0118
// N-CA bond  #  NA    N(0)  - Ca(0)
#define pPRO_PRObondNCA 	 		1.4649
#define pPRO_sigma_PRObondNCA 	 	0.0203
// CA-CO bond #  AC    Ca(0) - C(0)
#define pPRO_PRObondCACO    		1.5171
#define pPRO_sigma_PRObondCACO    	0.0093
// CA-CB bond #  AB    Ca(0) - Cb(0)
#define pPRO_PRObondCACB    		1.5357
#define pPRO_sigma_PRObondCACB  	0.0137
#define pPRO_PRObondCACB2    		3.0714
// CO=O bond  #  CO    C(0) - O(0)
#define pPRO_PRObondCO  			1.2404
#define pPRO_sigma_PRObondCO2  		0.0112
#define pPRO_PRObondCO2  			2.4808
//
// VALENCE ANGLES (rad)
// valence angle CO-N-CA   # CNA   C(-1) - N(0)  - Ca(0)
#define pPRO_PROvalCONCA   			120.38*M_PI/180
#define pPRO_sigma_PROvalCONCA   	  1.03*M_PI/180
// valence angle N-CA-CO   # NAC   N(0)  - Ca(0) - C(0)
#define pPRO_PROvalNCACO   			110.70*M_PI/180
#define pPRO_sigma_PROvalNCACO   	  1.22*M_PI/180
// valence angle CA-CO-N   # ACN   Ca(0) - C(0)  - N(+1)
#define pPRO_PROvalCACON   			117.93*M_PI/180
#define pPRO_sigma_PROvalCACON   	  1.20*M_PI/180
// valence angle O-CO-N    # OCN   O(0)  - C(0)  - N(+1)
#define pPRO_PROvalOCON   			121.46*M_PI/180
#define pPRO_sigma_PROvalOCON  		  1.18*M_PI/180




//BOND LENGTHS (Angstrom)
//#define bondNCA 	 	1.457	// N-CA bond   sigma +- 0.013   # NA    N(0)  - Ca(0)
//#define bondCACO    	1.523	// CA-CO bond  sigma +-         # AC    Ca(0) - C(0)
//#define bondCACO2  		3.046
//#define bondCACB    	1.535	// CA-CB bond  sigma +-         # AB    Ca(0) - Cb(0)   <--- Mon (I've added this...)
//#define bondCACB2    	3.070
//#define bondCON 	 	1.332	// CO-N bond   sigma +-         # CN    C(-1) - N(0)
//#define bondCO  		1.235	// CO=O bond   sigma +-         # CO    C(0) - O(0)
//#define bondCO2  		2.470
//#define sigma_bondNCA 	0.013
//#define sigma_bondCACO  0.017
//#define sigma_bondCACB  0.017
//#define sigma_bondCON 	0.014
// Other H-based (old) bond lengths...
#define bondNH  		1.02
#define bondCH  		1.10
#define bondCH2  		2.20

//VALENCE ANGLES (rad)
//#define valCONCA   		121.7*M_PI/180 // valence angle CO-N-CA   # CNA   C(-1) - N(0)  - Ca(0)
//#define valNCACO   		111.0*M_PI/180 // valence angle N-CA-CO   # NAC   N(0)  - Ca(0) - C(0)
//#define valCACON   		117.2*M_PI/180 // valence angle CA-CO-N   # ACN   Ca(0) - C(0)  - N(+1)
//#define valOCON   		122.6*M_PI/180 // valence angle O-CO-N    # OCN   O(0)  - C(0)  - N(+1)
//#define sigma_valCONCA   	1.8*M_PI/180
//#define sigma_valNCACO   	2.1*M_PI/180
//#define sigma_valCACON   	1.5*M_PI/180
//#define sigma_valOCON   	1.3*M_PI/180

//
//// GLYCINE:
//// -------
////BOND LENGTHS (Angstrom)
//#define GLYbondNCA 	 			1.449  // N-CA bond
//#define GLYbondCACO    			1.5143 // CA-CO bond
//#define GLYbondCON 	 			1.3305 // CO-N bond
//#define sigma_GLYbondNCA 	 	0.0145
//#define sigma_GLYbondCACO    	0.0141
//#define sigma_GLYbondCON 	 	0.0146
////VALENCE ANGLES (rad) REGULAR PEPTIDES
//#define GLYvalCONCA   			121.4*M_PI/180 // valence angle CO-N-CA   # CNA   C(-1) - N(0)  - Ca(0)
//#define GLYvalNCACO   			113.2*M_PI/180 // valence angle N-CA-CO   # NAC   N(0)  - Ca(0) - C(0)
//#define GLYvalCACON   			116.7*M_PI/180 // valence angle CA-CO-N   # ACN   Ca(0) - C(0)  - N(+1)
//#define sigma_GLYvalCONCA   	1.7*M_PI/180
//#define sigma_GLYvalNCACO   	2.2*M_PI/180
//#define sigma_GLYvalCACON   	2.0*M_PI/180
//
// PROLINE:
// --------
//BOND LENGTHS (Angstrom)
//#define PRObondNCA 	 		1.4687 // N-CA bond
//#define PRObondCACO    		1.5195 // CA-CO bond
//#define PRObondCON 	 		1.3339 // CO-N bond
//#define sigma_PRObondNCA 	0.0128
//#define sigma_PRObondCACO   0.0142
//#define sigma_PRObondCON 	0.0234
////VALENCE ANGLES (rad) REGULAR PEPTIDES
//#define PROvalCONCA   		119.8*M_PI/180		/*valence angle CO-N-CA	    */
//#define PROvalNCACO   		112.5*M_PI/180		/*valence angle N-CA-CO	    */
//#define PROvalCACON   		116.7*M_PI/180		/*valence angle CA-CO-N	    */
//#define sigma_PROvalCONCA   1.3*M_PI/180
//#define sigma_PROvalNCACO   2.1*M_PI/180
//#define sigma_PROvalCACON   2.1*M_PI/180
//
// PRO-PrePROLINE:
// --------
//BOND LENGTHS (Angstrom)
//#define PRObondNCA 	 		1.4687 // N-CA bond
//#define PRObondCACO    		1.5195 // CA-CO bond
//#define PRObondCON 	 		1.3339 // CO-N bond
//#define sigma_PRObondNCA 	0.0128
//#define sigma_PRObondCACO   0.0142
//#define sigma_PRObondCON 	0.0234
////VALENCE ANGLES (rad) REGULAR PEPTIDES
//#define PROvalCONCA   		120.4*M_PI/180		/*valence angle CO-N-CA	    */
//#define PROvalNCACO   		110.7*M_PI/180		/*valence angle N-CA-CO	    */
//#define PROvalCACON   		117.9*M_PI/180		/*valence angle CA-CO-N	    */
//#define sigma_PROvalCONCA   1.0*M_PI/180
//#define sigma_PROvalNCACO   2.2*M_PI/180
//#define sigma_PROvalCACON   1.2*M_PI/180


//// Rosetta's standard parameters for "centroid_pre_talaris2013", "fa_standard_05.2009_icoor", and "fa_standard_05.2009_icoor"
//// ~/prgs/rosetta_2014.35.57232_bundle/main/database/chemical/residue_type_sets/centroid_pre_talaris2013/residue_types
//#define bondNCA 	 	1.458	// N-CA bond   sigma +-
//#define bondCACO    	1.52326	// CA-CO bond  sigma +-
//#define bondCON 	 	1.3287	// CO-N bond   sigma +-


//BOND LENGTHS (Angstrom)
//#define bondNCA 	 	1.458	/*N-CA bond      sigma	+- 0.019	*/
//#define bondCACO    	1.525	/*CA-CO bond	 sigma  +- 0.021	*/
//#define bondCON 	 	1.329	/*CO-N bond      sigma	+- 0.014	*/

//#define sigma_bondNCA 	 	0.012 // 0.019
//#define sigma_bondCACO    	0.010 // 0.021
//#define sigma_bondCON 	 	0.010 // 0.014

////VALENCE ANGLES (rad)
//#define valCONCA   		121.7*M_PI/180		/*valence angle CO-N-CA	       */
//#define valNCACO   		111.2*M_PI/180		/*valence angle N-CA-CO	       */
//#define valCACON   		116.2*M_PI/180		/*valence angle CA-CO-N	       */
//
//#define sigma_valCONCA   	1.8*M_PI/180
//#define sigma_valNCACO   	2.8*M_PI/180
//#define sigma_valCACON   	2.0*M_PI/180
//
//#define sigma_valCONCA_B   		1.7*M_PI/180	//	VALUES FROM BERKHOLZ RELATED TO ENGH AND HUBER
//#define sigma_valNCACO_B   		1.5*M_PI/180    //  (SMALLER SIGMA AND SMALLER RANGE, OPTIMIZED VALUES)
//#define sigma_valCACON_B   		1.3*M_PI/180
//
//
////GLYCINE
//
////BOND LENGTHS (Angstrom)
//#define GLYbondNCA 	 			1.41			/*N-CA bond      	*/
//#define GLYbondCACO    			1.516			/*CA-CO bond	 	*/
//#define GLYbondCON 	 			1.329	  		/*CO-N bond		*/
//
//#define sigma_GLYbondNCA 	 	0.016
//#define sigma_GLYbondCACO    	0.018
//#define sigma_GLYbondCON 	 	0.014
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
////PROLINE
//
////BOND LENGTHS (Angstrom)
//#define PRObondNCA 	 		1.466			/*N-CA bond  	*/
//#define PRObondCACO    		1.525			/*CA-CO bond	*/
//#define PRObondCON 	 		1.341	  		/*CO-N bond	*/
//
//#define sigma_PRObondNCA 	 	0.015
//#define sigma_PRObondCACO    	0.021
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

//#define LOWERPHI -175*M_PI/180
//#define UPPERPHI  1-40*M_PI/180
//#define LOWERPSI -1-60*M_PI/180
//#define UPPERPSI  175*M_PI/180

#endif
