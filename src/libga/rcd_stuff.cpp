/*
 * rcd_stuff.cpp
 *
 *  Created on: Jun 16, 2016
 *      Author: mon
 */

#include "rcd_stuff.h"

/*===============================================================================
 *
 *	BONDSELECT, by Pieter
 *
 *===============================================================================
 *
 * to account for unallowed bond rotations under the normal PHI and PSI angles.
 *
 *--------------------------------------------------------------------------------
 *
 *  INPUT: - int range 			= 	total number of bonds available
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

// CRandomMersenne *rg; // Mersenne Twister global object

// Mon: Remove REVERSE here!!! (not used)
int	bondselect_marker(int nr_atoms_loop, int *marker,int *previous, int reverse)
{
	int marker_count;
	int marker_list[100];

	marker_count=0; // number of rotable bonds
	for (int z = 2; z < nr_atoms_loop-1; z++)
	{
		if (marker[z]==1) // if rotable...
		{
			if (z!=*previous) // if not previous
			{
				marker_list[marker_count]=z;
				marker_count++;
			}
		}
	}

	if (marker_count <= 0 )
	{
		fprintf(stderr,"Msg(bondselect_marker): Huston we've got a problem!\n marker_count=%d  *previous=%d",marker_count,*previous);
		exit(1);
	}

	if (marker_count <= 1 )
		return  *previous=marker_list[0];
	else
		*previous=marker_list[rg->IRandom(0,marker_count-1)];

	return *previous;
}

// Mon: Remove REVERSE here!!! (not used)
int	bondselect_marker2(int nr_atoms_loop, int *marker,int *previous, int reverse)
{
	int marker_count;
	int marker_list[100];

	marker_count=0;
	for (int z = 2; z < nr_atoms_loop-1; z++)	 {
		if (marker[z]==1) {
			{marker_list[marker_count]=z;
			marker_count++;}
		}
	}

	if (marker_count == 1 ) { return  *previous=marker_list[0]; }
	else {*previous=marker_list[rg->IRandom(0,marker_count-1)];}

	return *previous;
}

void randomizeIC(double *bond_length,double *valence_angle,double *bond_length2,double *valence_angle2,int nratoms,double lengthptb,double valenceptb,char *seq, char seq_ct1, bool flat)
{
	int i;
	for (i = 1; i < (nratoms/3)-1; i++) // Screen all loop residues but Nt and Ct anchors.
	{
		if( !(seq[i+1] == 'P' || seq[i+1] == 'p') ) // NOT pre-PRO
		{
			if (seq[i] == 'I' || seq[i] == 'V') // Iso-Leucine or Valine
			{
				// Bond lengths
				bond_length[3*i] 		= 	bond_length2[3*i]     +   lengthptb * sigma_IVbondCON * normRand(flat);  // CN   IVbondCON
				bond_length[3*i+1] 	   	= 	bond_length2[3*i+1]   +   lengthptb * sigma_IVbondNCA * normRand(flat);  // NA   IVbondNCA
				bond_length[3*i+2]     	= 	bond_length2[3*i+2]   +   lengthptb * sigma_IVbondCACO * normRand(flat); // AC   IVbondCACO
				// Bond angles
				valence_angle[3*i+3]    = 	valence_angle2[3*i+3]	  +  valenceptb * sigma_IVvalCACON * normRand(flat);   // ACN (of the next residue)   IVvalCACON
				valence_angle[3*i+1]    = 	valence_angle2[3*i+1]	  +  valenceptb * sigma_IVvalCONCA * normRand(flat);   // CNA      IVvalCONCA
				valence_angle[3*i+2]    = 	valence_angle2[3*i+2]	  +  valenceptb * sigma_IVvalNCACO * normRand(flat);   // NAC      IVvalNCACO
			}
			else if (seq[i] == 'G') // Glycine
			{
				// Bond lengths
				bond_length[3*i] 		= 	bond_length2[3*i]     +   lengthptb * sigma_GLYbondCON * normRand(flat);  // CN     GLYbondCON
				bond_length[3*i+1] 	   	= 	bond_length2[3*i+1]   +   lengthptb * sigma_GLYbondNCA * normRand(flat);  // NA     GLYbondNCA
				bond_length[3*i+2]     	= 	bond_length2[3*i+2]   +   lengthptb * sigma_GLYbondCACO * normRand(flat); // AC     GLYbondCACO
				// Bond angles
				valence_angle[3*i+3]    = 	valence_angle2[3*i+3]	  +  valenceptb * sigma_GLYvalCACON * normRand(flat);   // ACN (of the next residue)   GLYvalCACON
				valence_angle[3*i+1]    = 	valence_angle2[3*i+1]	  +  valenceptb * sigma_GLYvalCONCA * normRand(flat);   // CNA     GLYvalCONCA
				valence_angle[3*i+2]    = 	valence_angle2[3*i+2]	  +  valenceptb * sigma_GLYvalNCACO * normRand(flat);   // NAC     GLYvalNCACO
			}
			else if (seq[i] == 'P' || seq[i] == 'p') // Proline
			{
				// Bond lengths
				bond_length[3*i] 		= 	bond_length2[3*i]     +   lengthptb * sigma_PRObondCON * normRand(flat);  // CN     PRObondCON
				bond_length[3*i+1] 	   	= 	bond_length2[3*i+1]   +   lengthptb * sigma_PRObondNCA * normRand(flat);  // NA     PRObondNCA
				bond_length[3*i+2]     	= 	bond_length2[3*i+2]   +   lengthptb * sigma_PRObondCACO * normRand(flat); // AC     PRObondCACO
				// Bond angles
				valence_angle[3*i+3]    = 	valence_angle2[3*i+3]	  +  valenceptb * sigma_PROvalCACON * normRand(flat);   // ACN (of the next residue)   PROvalCACON
				valence_angle[3*i+1]    = 	valence_angle2[3*i+1]	  +  valenceptb * sigma_PROvalCONCA * normRand(flat);   // CNA     PROvalCONCA
				valence_angle[3*i+2]    = 	valence_angle2[3*i+2]	  +  valenceptb * sigma_PROvalNCACO * normRand(flat);   // NAC     PROvalNCACO
			}
			else // NOT I, V, G or P (by default)
			{
				// Bond lengths
				bond_length[3*i] 		= 	bond_length2[3*i]     +   lengthptb * sigma_bondCON * normRand(flat);  // CN     bondCON
				bond_length[3*i+1] 	   	= 	bond_length2[3*i+1]   +   lengthptb * sigma_bondNCA * normRand(flat);  // NA     bondNCA
				bond_length[3*i+2]     	= 	bond_length2[3*i+2]   +   lengthptb * sigma_bondCACO * normRand(flat); // AC     bondCACO
				// Bond angles
				valence_angle[3*i+3]    = 	valence_angle2[3*i+3]	  +  valenceptb * sigma_valCACON * normRand(flat);   // ACN (of the next residue)   valCACON
				valence_angle[3*i+1]    = 	valence_angle2[3*i+1]	  +  valenceptb * sigma_valCONCA * normRand(flat);   // CNA     valCONCA
				valence_angle[3*i+2]    = 	valence_angle2[3*i+2]	  +  valenceptb * sigma_valNCACO * normRand(flat);   // NAC     valNCACO
			}
		}
		else // pre-PRO
		{
			if (seq[i] == 'I' || seq[i] == 'V') // Iso-Leucine or Valine
			{
				// Bond lengths
				bond_length[3*i] 		= 	bond_length2[3*i]     +   lengthptb * pPRO_sigma_IVbondCON * normRand(flat);  // CN     pPRO_IVbondCON
				bond_length[3*i+1] 	   	= 	bond_length2[3*i+1]   +   lengthptb * pPRO_sigma_IVbondNCA * normRand(flat);  // NA     pPRO_IVbondNCA
				bond_length[3*i+2]     	= 	bond_length2[3*i+2]   +   lengthptb * pPRO_sigma_IVbondCACO * normRand(flat); // AC     pPRO_IVbondCACO
				// Bond angles
				valence_angle[3*i+3]    = 	valence_angle2[3*i+3]	  +  valenceptb * pPRO_sigma_IVvalCACON * normRand(flat);   // ACN (of the next residue)   pPRO_IVvalCACON
				valence_angle[3*i+1]    = 	valence_angle2[3*i+1]	  +  valenceptb * pPRO_sigma_IVvalCONCA * normRand(flat);   // CNA     pPRO_IVvalCONCA
				valence_angle[3*i+2]    = 	valence_angle2[3*i+2]	  +  valenceptb * pPRO_sigma_IVvalNCACO * normRand(flat);   // NAC     pPRO_IVvalNCACO
			}
			else if (seq[i] == 'G') // Glycine
			{
				// Bond lengths
				bond_length[3*i] 		= 	bond_length2[3*i]     +   lengthptb * pPRO_sigma_GLYbondCON * normRand(flat);  // CN     pPRO_GLYbondCON
				bond_length[3*i+1] 	   	= 	bond_length2[3*i+1]   +   lengthptb * pPRO_sigma_GLYbondNCA * normRand(flat);  // NA     pPRO_GLYbondNCA
				bond_length[3*i+2]     	= 	bond_length2[3*i+2]   +   lengthptb * pPRO_sigma_GLYbondCACO * normRand(flat); // AC     pPRO_GLYbondCACO
				// Bond angles
				valence_angle[3*i+3]    = 	valence_angle2[3*i+3]	  +  valenceptb * pPRO_sigma_GLYvalCACON * normRand(flat);   // ACN (of the next residue)   pPRO_GLYvalCACON
				valence_angle[3*i+1]    = 	valence_angle2[3*i+1]	  +  valenceptb * pPRO_sigma_GLYvalCONCA * normRand(flat);   // CNA     pPRO_GLYvalCONCA
				valence_angle[3*i+2]    = 	valence_angle2[3*i+2]	  +  valenceptb * pPRO_sigma_GLYvalNCACO * normRand(flat);   // NAC     pPRO_GLYvalNCACO
			}
			else if (seq[i] == 'P' || seq[i] == 'p') // Proline
			{
				// Bond lengths
				bond_length[3*i] 		= 	bond_length2[3*i]     +   lengthptb * pPRO_sigma_PRObondCON * normRand(flat);  // CN     pPRO_PRObondCON
				bond_length[3*i+1] 	   	= 	bond_length2[3*i+1]   +   lengthptb * pPRO_sigma_PRObondNCA * normRand(flat);  // NA     pPRO_PRObondNCA
				bond_length[3*i+2]     	= 	bond_length2[3*i+2]   +   lengthptb * pPRO_sigma_PRObondCACO * normRand(flat); // AC     pPRO_PRObondCACO
				// Bond angles
				valence_angle[3*i+3]    = 	valence_angle2[3*i+3]	  +  valenceptb * pPRO_sigma_PROvalCACON * normRand(flat);   // ACN (of the next residue)   pPRO_PROvalCACON
				valence_angle[3*i+1]    = 	valence_angle2[3*i+1]	  +  valenceptb * pPRO_sigma_PROvalCONCA * normRand(flat);   // CNA     pPRO_PROvalCONCA
				valence_angle[3*i+2]    = 	valence_angle2[3*i+2]	  +  valenceptb * pPRO_sigma_PROvalNCACO * normRand(flat);   // NAC     pPRO_PROvalNCACO
			}
			else // NOT I, V, G or P (by default)
			{
				// Bond lengths
				bond_length[3*i] 		= 	bond_length2[3*i]     +   lengthptb * pPRO_sigma_bondCON * normRand(flat);  // CN     pPRO_bondCON
				bond_length[3*i+1] 	   	= 	bond_length2[3*i+1]   +   lengthptb * pPRO_sigma_bondNCA * normRand(flat);  // NA     pPRO_bondNCA
				bond_length[3*i+2]     	= 	bond_length2[3*i+2]   +   lengthptb * pPRO_sigma_bondCACO * normRand(flat); // AC     pPRO_bondCACO
				// Bond angles
				valence_angle[3*i+3]    = 	valence_angle2[3*i+3]	  +  valenceptb * pPRO_sigma_valCACON * normRand(flat);   // ACN (of the next residue)   pPRO_valCACON
				valence_angle[3*i+1]    = 	valence_angle2[3*i+1]	  +  valenceptb * pPRO_sigma_valCONCA * normRand(flat);   // CNA     pPRO_valCONCA
				valence_angle[3*i+2]    = 	valence_angle2[3*i+2]	  +  valenceptb * pPRO_sigma_valNCACO * normRand(flat);   // NAC     pPRO_valNCACO
			}
		}
	}

	// Randomizing the remaining Ct anchor CN bond length and CNA bond angle
	i=(nratoms/3)-1;
	if( !(seq_ct1 == 'P' || seq_ct1 == 'p') ) // NOT pre-PRO
	{
		if (seq[i] == 'I' || seq[i] == 'V') // Iso-Leucine or Valine
		{
			bond_length[3*i] 		= 	bond_length2[3*i]     +   lengthptb * sigma_IVbondCON * normRand(flat);  // CN     IVbondCON
			valence_angle[3*i+1]    = 	valence_angle2[3*i+1]	  +  valenceptb * sigma_IVvalCONCA * normRand(flat);   // CNA        IVvalCONCA
		}
		else if (seq[i] == 'G') // Glycine
		{
			bond_length[3*i] 		= 	bond_length2[3*i]     +   lengthptb * sigma_GLYbondCON * normRand(flat);  // CN     GLYbondCON
			valence_angle[3*i+1]    = 	valence_angle2[3*i+1]	  +  valenceptb * sigma_GLYvalCONCA * normRand(flat);   // CNA        GLYvalCONCA
		}
		else if (seq[i] == 'P' || seq[i] == 'p') // Proline
		{
			bond_length[3*i] 		= 	bond_length2[3*i]     +   lengthptb * sigma_PRObondCON * normRand(flat);  // CN     PRObondCON
			valence_angle[3*i+1]    = 	valence_angle2[3*i+1]	  +  valenceptb * sigma_PROvalCONCA * normRand(flat);   // CNA        PROvalCONCA
		}
		else // NOT I, V, G or P (by default)
		{
			bond_length[3*i] 		= 	bond_length2[3*i]     +   lengthptb * sigma_bondCON * normRand(flat);  // CN     bondCON
			valence_angle[3*i+1]    = 	valence_angle2[3*i+1]	  +  valenceptb * sigma_valCONCA * normRand(flat);   // CNA        valCONCA
		}
	}
	else // pre-PRO
	{
		if (seq[i] == 'I' || seq[i] == 'V') // Iso-Leucine or Valine
		{
			bond_length[3*i] 		= 	bond_length2[3*i]     +   lengthptb * pPRO_sigma_IVbondCON * normRand(flat);  // CN     pPRO_IVbondCON
			valence_angle[3*i+1]    = 	valence_angle2[3*i+1]	  +  valenceptb * pPRO_sigma_IVvalCONCA * normRand(flat);   // CNA        pPRO_IVvalCONCA
		}
		else if (seq[i] == 'G') // Glycine
		{
			bond_length[3*i] 		= 	bond_length2[3*i]     +   lengthptb * pPRO_sigma_GLYbondCON * normRand(flat);  // CN     pPRO_GLYbondCON
			valence_angle[3*i+1]    = 	valence_angle2[3*i+1]	  +  valenceptb * pPRO_sigma_GLYvalCONCA * normRand(flat);   // CNA        pPRO_GLYvalCONCA
		}
		else if (seq[i] == 'P' || seq[i] == 'p') // Proline
		{
			bond_length[3*i] 		= 	bond_length2[3*i]     +   lengthptb * pPRO_sigma_PRObondCON * normRand(flat);  // CN     pPRO_PRObondCON
			valence_angle[3*i+1]    = 	valence_angle2[3*i+1]	  +  valenceptb * pPRO_sigma_PROvalCONCA * normRand(flat);   // CNA        pPRO_PROvalCONCA
		}
		else // NOT I, V, G or P (by default)
		{
			bond_length[3*i] 		= 	bond_length2[3*i]     +   lengthptb * pPRO_sigma_bondCON * normRand(flat);  // CN     pPRO_bondCON
			valence_angle[3*i+1]    = 	valence_angle2[3*i+1]	  +  valenceptb * pPRO_sigma_valCONCA * normRand(flat);   // CNA        pPRO_valCONCA
		}
	}
	// CONSIDER REMOVAL...
	// Backup Ct-anchor (NOT randomizing)
	bond_length[nratoms-2]	= bond_length2[nratoms -2]; // N-CA    		(NA)
	bond_length[nratoms-1]	= bond_length2[nratoms -1]; // CA-C    		(AC)
	valence_angle[nratoms -1] = valence_angle2[nratoms -1]; // N-CA-C 	(NAC)


	// Randomizing the remaining Nt anchor ACN bond angle
	if( !(seq[1] == 'P' || seq[1] == 'p') ) // NOT pre-PRO
	{
		if (seq[0] == 'I' || seq[0] == 'V') // Iso-Leucine or Valine
			valence_angle[3]    = 	valence_angle2[3]	  +  valenceptb * sigma_IVvalCACON * normRand(flat);   // ACN        IVvalCACON
		else if (seq[0] == 'G') // Glycine
			valence_angle[3]    = 	valence_angle2[3]	  +  valenceptb * sigma_GLYvalCACON * normRand(flat);   // ACN        GLYvalCACON
		else if (seq[0] == 'P' || seq[0] == 'p') // Proline
			valence_angle[3]    = 	valence_angle2[3]	  +  valenceptb * sigma_PROvalCACON * normRand(flat);   // ACN        PROvalCACON
		else // NOT I, V, G or P (by default)
			valence_angle[3]    = 	valence_angle2[3]	  +  valenceptb * sigma_valCACON * normRand(flat);   // ACN        valCACON
	}
	else // pre-PRO
	{
		if (seq[0] == 'I' || seq[0] == 'V') // Iso-Leucine or Valine
			valence_angle[3]    = 	valence_angle2[3]	  +  valenceptb * pPRO_sigma_IVvalCACON * normRand(flat);   // ACN        pPRO_IVvalCACON
		else if (seq[0] == 'G') // Glycine
			valence_angle[3]    = 	valence_angle2[3]	  +  valenceptb * pPRO_sigma_GLYvalCACON * normRand(flat);   // ACN        pPRO_GLYvalCACON
		else if (seq[0] == 'P' || seq[0] == 'p') // Proline
			valence_angle[3]    = 	valence_angle2[3]	  +  valenceptb * pPRO_sigma_PROvalCACON * normRand(flat);   // ACN        pPRO_PROvalCACON
		else // NOT I, V, G or P (by default)
			valence_angle[3]    = 	valence_angle2[3]	  +  valenceptb * pPRO_sigma_valCACON * normRand(flat);   // ACN        pPRO_valCACON
	}
	// CONSIDER REMOVAL...
	// Backup Nt-anchor (NOT randomizing)
	bond_length[1] = bond_length2[1]; // N-CA
	bond_length[2] = bond_length2[2]; // CA-C
	valence_angle[2] = valence_angle2[2]; // N-CA-C

	// CONSIDER REMOVAL...
	// We don't know these...so zero.
	bond_length[0] = 0.0; // consider removal
	valence_angle[0] = 0.0; // consider removal
	valence_angle[1] = 0.0; // consider removal																								//ROW 2
}

// A random number generator that picks numbers from a normal distribution by default (Marsaglia's Polar Method)
// Or it can return a random numbers according to a flat [-1:1[ random distribution if "flat=true"
float normRand(bool flat)
{
	if(flat) // Flat distribution
		return 2.0 * rg->Random() - 1.0;
	else // Normal distribution
	{
		float rsq, v2;
		do
		{
			rsq = 2.0 * rg->Random() - 1.0;
			v2 = 2.0 * rg->Random() - 1.0;
			rsq = rsq * rsq + v2 * v2;
		}
		while ( rsq >= 1.0 || rsq == 0.0 );
		return sqrt( -2.0 * log( rsq ) / rsq ) * v2;
	}
}

// Loading coordinates from PDB-file
void loop_co( Macromolecule *loop,  double** co, double** p1co, int *residuemarker, int * seq, bool rmsd_O)
{
	Fragment * resi;
	int j,resn, res_count;
	char * at_name;
	Segment * seg;
	Atom * atomi;
	//	pdbIter *iter1, *iter2;
	pdbIter *iter2;
	int atoms_for_checking;
	Tcoor pos;
	bool  misCA, misC, misN, misO, misCB;

	//Iterador para recorrer segmentos
	//	iter1 = new pdbIter( loop );
	iter2 = new pdbIter( loop );

	// Mon 4/11/2016: Not-sorted atoms bug fixed!
	//	//Bucle para recorrer segmentos
	//	for ( iter1->pos_segment = 0; !iter1->gend_segment(); iter1->next_segment() )
	//	{
	//		seg = ( Segment * ) iter1->get_segment();
	//		iter2 = new pdbIter( seg );

	//Bucle para recorrer residuos
	for ( iter2->pos_fragment = 0; !iter2->gend_fragment(); iter2->next_fragment() )
	{
		resi = ( Fragment * ) iter2->get_fragment();
		resn = resnum_from_resname( resi->getName() );
		seq[iter2->pos_fragment]=resn;
		// fprintf( stderr, "Processing aa %s resn= %d  %d\n", resi->getName(), resn, iter2->pos_fragment );

		if ( resn > N_AMINO )
		{
			fprintf( stderr, "Warning %s residue is not an aa\n", resi->getName());
			//getchar();
		}

		//MARK NON-REGULAR AMINO ACIDS GLYCINE AND PROLINE
		if ( resn==5)  /*Mark a GLYCINE residue	*/
		{
			residuemarker[iter2->pos_fragment] = 1;
			// fprintf( stderr, " residue \"%s\" ", resi->getName());

		} else if ( resn == 12 ) /*Mark a PROLINE residue	*/
		{
			residuemarker[iter2->pos_fragment] = 2;
			//	      				fprintf( stderr, " residue \"%s\" ", resi->getName());
		}

		//Iterador sobre los atomos del residuo
		pdbIter * iter3 = new pdbIter( resi );

		res_count=3*(iter2->pos_fragment );
		misCA=true;
		misC=true;
		misN=true;
		misO=true;
		misCB=true;

		for ( iter3->pos_atom = 0; !iter3->gend_atom(); iter3->next_atom() )
		{
			atomi = ( Atom * ) iter3->get_atom();
			at_name = atomi->getName();
			atomi->getPosition(pos);
			//fprintf( stderr, " Atom \"%s\" of residue %s %d at_check= %d %d it3= %d %f %f %f\n",
			//		  at_name, resi->getName(), resn, AA[resn].nheavyatoms, atoms_for_checking, iter3->pos_fragment *3, pos[0], pos[1], pos[2]);

			if ( strcmp(at_name," N  ") == 0  )
			{
				co[res_count][0] 	= pos[0];
				co[res_count][1] 	= pos[1];
				co[res_count][2] 	= pos[2];
				misN=false;
			}
			else  if ( strcmp(at_name," CA ") == 0  )
			{
				co[res_count+1][0] 	= pos[0];
				co[res_count +1][1] 	= pos[1];
				co[res_count +1][2] 	= pos[2];
				misCA=false;
			}
			else if ( strcmp(at_name," C  ") == 0  )
			{
				co[res_count +2][0] 	= pos[0];
				co[res_count +2][1] 	= pos[1];
				co[res_count +2][2] 	= pos[2];
				misC=false;
			}
			else if ( strcmp(at_name," O  ") == 0  )
			{
				p1co[res_count +2][0] 	= pos[0];
				p1co[res_count +2][1] 	= pos[1];
				p1co[res_count +2][2] 	= pos[2];
				misO=false;
			}
			else if ( strcmp(at_name," CB ") == 0  )
			{
				p1co[res_count +1][0] 	= pos[0];
				p1co[res_count +1][1] 	= pos[1];
				p1co[res_count +1][2] 	= pos[2];
				misCB=false;
			}

		}
		delete iter3;

		// check missing atoms
		if (misCA ||  misC || misN )
		{
			if (misCA)
				fprintf( stderr, "Error: Atom CA of residue %s is missing (loop residue %d). Check input coordinates.\n",resi->getName(),iter2->pos_fragment);
			if (misC)
				fprintf( stderr, "Error: Atom C of residue %s is missing (loop residue %d. Check input coordinates.\n",resi->getName(),iter2->pos_fragment);
			if (misN)
				fprintf( stderr, "Error: Atom N of residue %s is missing (loop residue %d). Check input coordinates.\n",resi->getName(),iter2->pos_fragment);
			exit(0);
		}

		// check missing atoms
		if (misO && rmsd_O)
			fprintf( stderr, "Warning: Atom O of residue %s is missing (loop residue %d). Check input coordinates.\n",resi->getName(),iter2->pos_fragment);

		//			if ((misCB)&& ( resn!=5))
		//				fprintf( stderr, "Warning: Atom CB of residue %s is missing (loop residue %d). Check input coordinates.\n",resi->getName(),iter2->pos_fragment);
	}
	delete iter2;
	//	}
	//	delete iter1;
}

// Create residuemarker and seqnum from Macromolecule
void residuemarker_seq( Macromolecule *loop, int *residuemarker, int *seqnum)
{
	Fragment * resi;
	int resn;
	Segment * seg;
	pdbIter *iter1, *iter2;

	//Iterador para recorrer segmentos
	iter1 = new pdbIter( loop );

	//Bucle para recorrer segmentos
	for ( iter1->pos_segment = 0; !iter1->gend_segment(); iter1->next_segment() )
	{
		seg = ( Segment * ) iter1->get_segment();
		iter2 = new pdbIter( seg );

		//Bucle para recorrer residuos
		for ( iter2->pos_fragment = 0; !iter2->gend_fragment(); iter2->next_fragment() )
		{
			resi = ( Fragment * ) iter2->get_fragment();
			resn = resnum_from_resname( resi->getName() );
			seqnum[iter2->pos_fragment]=resn;
			// fprintf( stderr, "Processing aa %s resn= %d  %d\n", resi->getName(), resn, iter2->pos_fragment );

			if ( resn > N_AMINO )
			{
				fprintf( stderr, "Warning %s residue is not an aa\n", resi->getName());
				//getchar();
			}

			//MARK NON-REGULAR AMINO ACIDS GLYCINE AND PROLINE
			if ( resn==5)  /*Mark a GLYCINE residue	*/
			{
				residuemarker[iter2->pos_fragment] = 1;
				// fprintf( stderr, " residue \"%s\" ", resi->getName());

			} else if ( resn == 12 ) /*Mark a PROLINE residue	*/
			{
				residuemarker[iter2->pos_fragment] = 2;
				//	      				fprintf( stderr, " residue \"%s\" ", resi->getName());
			}
		}
		delete iter2;
	}
	delete iter1;
}

int clash(vlVolume *vlm, double** co, double**  p1co,  int nr_atoms_loop, int reverse, bool check_p1co)
{
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

	for (int index = 3; index < nr_atoms_loop - 2; index++)
	{
		x=(int)roundf( ((co[index][0]- pos.x()) / unitsx) );
		y=(int)roundf( ((co[index][1]- pos.y()) / unitsy) );
		z=(int)roundf( ((co[index][2]- pos.z()) / unitsz) );

		if( (x < limitx && y < limity && z < limitz && x >= 0 && y >= 0 && z >= 0) && (*(m_pData + (z*stepz+y*stepy+x))>0.00001) )
			return 1;

		if(check_p1co) // Check "O"
		{
			x=(int)roundf( ((p1co[index][0]- pos.x()) / unitsx) );
			y=(int)roundf( ((p1co[index][1]- pos.y()) / unitsy) );
			z=(int)roundf( ((p1co[index][2]- pos.z()) / unitsz) );

			if( (x < limitx && y < limity && z < limitz && x >= 0 && y >= 0 && z >= 0) && (*(m_pData + (z*stepz+y*stepy+x))>0.00001) )
				return 1;
		}
	}

	return 0;
}

// Check clashes between loop and pocket-map
//  mode  --> mode of operation: 0=(N,CA,C), 1=(N,CA,C,O), 2=(N,CA,C,CB), 3=(N,CA,C,O,CB)
//  count --> set true to enable counting all clashes, otherwise (false) it will return upon first clash detection
int clash_cg(vlVolume *vlm, double **co, double **p1co, int nr_atoms_loop, int *marker2, int reverse, int mode, bool count)
{
	vlPoint3f pos;
	float *m_pData = (float*)vlm->getVoxelVoidPtr(vlPoint3ui(0,0,0));
	vlm->getPosition(&pos);
	int x,y,z;
	int clashes = 0; // counting clashes

	if(reverse==0) // forward
	{
		if( mode == 3)
		{
			for (int i = nr_atoms_loop - 3; i>=3 ; i--)
				if (i%3 != 0) // Checking O & CB
				{
					if (!((i%3 == 1)&&(marker2[i] == 2)))
					{
						x=(int)roundf( ((p1co[i][0]- pos.x()) / unitsx) );
						y=(int)roundf( ((p1co[i][1]- pos.y()) / unitsy) );
						z=(int)roundf( ((p1co[i][2]- pos.z()) / unitsz) );
						if( (x < limitx && y < limity && z < limitz && x >= 0 && y >= 0 && z >= 0) && (*(m_pData + (z*stepz+y*stepy+x))>0.00001) )
							if(count)
								clashes++; // count clashes
							else
								return 1; // return on first clash
					}
				}
		}
		else if( mode == 2)
		{
			for (int i = nr_atoms_loop - 3; i>=3 ; i--)
			{
				if ((i%3 == 1)&&(marker2[i] != 2)) // checking CB and non GLY
				{
					x=(int)roundf( ((p1co[i][0]- pos.x()) / unitsx) );
					y=(int)roundf( ((p1co[i][1]- pos.y()) / unitsy) );
					z=(int)roundf( ((p1co[i][2]- pos.z()) / unitsz) );
					if( (x < limitx && y < limity && z < limitz && x >= 0 && y >= 0 && z >= 0) && (*(m_pData + (z*stepz+y*stepy+x))>0.00001) )
						if(count)
							clashes++; // count clashes
						else
							return 1; // return on first clash
				}
			}
		}
		else if( mode == 1) {
			for (int i = nr_atoms_loop - 3; i>=3 ; i--)
			{
				if (i%3 == 2) // checking 0
				{
					x=(int)roundf( ((p1co[i][0]- pos.x()) / unitsx) );
					y=(int)roundf( ((p1co[i][1]- pos.y()) / unitsy) );
					z=(int)roundf( ((p1co[i][2]- pos.z()) / unitsz) );
					if( (x < limitx && y < limity && z < limitz && x >= 0 && y >= 0 && z >= 0) && (*(m_pData + (z*stepz+y*stepy+x))>0.00001) )
						if(count)
							clashes++; // count clashes
						else
							return 1; // return on first clash
				}
			}
		}
		// Backbone
		for (int i = nr_atoms_loop - 3; i>=3 ; i--)
		{
			x=(int)roundf( ((co[i][0]- pos.x()) / unitsx) );
			y=(int)roundf( ((co[i][1]- pos.y()) / unitsy) );
			z=(int)roundf( ((co[i][2]- pos.z()) / unitsz) );

			if( (x < limitx && y < limity && z < limitz && x >= 0 && y >= 0 && z >= 0) && (*(m_pData + (z*stepz+y*stepy+x))>0.00001) )
				if(count)
					clashes++; // count clashes
				else
					return 1; // return on first clash
		}
	}
	else // reverse
	{
		if( mode == 3)
		{
			for (int i = 3; i < nr_atoms_loop - 2; i++)  //  divisible by 3
				if (i%3 != 2) // Checking O & CB
				{
					if (!((i%3 == 0)&&(marker2[i] == 2)))
					{
						x=(int)roundf( ((p1co[i][0]- pos.x()) / unitsx) );
						y=(int)roundf( ((p1co[i][1]- pos.y()) / unitsy) );
						z=(int)roundf( ((p1co[i][2]- pos.z()) / unitsz) );
						if( (x < limitx && y < limity && z < limitz && x >= 0 && y >= 0 && z >= 0) && (*(m_pData + (z*stepz+y*stepy+x))>0.00001) )
							if(count)
								clashes++; // count clashes
							else
								return 1; // return on first clash
					}
				}
		}
		else if( mode == 2)
		{
			for (int i = 3; i < nr_atoms_loop - 2; i++)
			{
				if ((i%3 == 0)&&(marker2[i] != 2)) // checking CB no gly
				{
					x=(int)roundf( ((p1co[i][0]- pos.x()) / unitsx) );
					y=(int)roundf( ((p1co[i][1]- pos.y()) / unitsy) );
					z=(int)roundf( ((p1co[i][2]- pos.z()) / unitsz) );
					if( (x < limitx && y < limity && z < limitz && x >= 0 && y >= 0 && z >= 0) && (*(m_pData + (z*stepz+y*stepy+x))>0.00001) )
						if(count)
							clashes++; // count clashes
						else
							return 1; // return on first clash
				}
			}
		}
		else if( mode == 1)
		{
			for (int i = 3; i < nr_atoms_loop - 2; i++)
			{
				if (i%3 == 1) // checking 0
				{
					x=(int)roundf( ((p1co[i][0]- pos.x()) / unitsx) );
					y=(int)roundf( ((p1co[i][1]- pos.y()) / unitsy) );
					z=(int)roundf( ((p1co[i][2]- pos.z()) / unitsz) );
					if( (x < limitx && y < limity && z < limitz && x >= 0 && y >= 0 && z >= 0) && (*(m_pData + (z*stepz+y*stepy+x))>0.00001) )
						if(count)
							clashes++; // count clashes
						else
							return 1; // return on first clash
				}
			}
		}
		// Backbone
		for (int i = 3; i < nr_atoms_loop - 2; i++)
		{
			x=(int)roundf( ((co[i][0]- pos.x()) / unitsx) );
			y=(int)roundf( ((co[i][1]- pos.y()) / unitsy) );
			z=(int)roundf( ((co[i][2]- pos.z()) / unitsz) );

			if( (x < limitx && y < limity && z < limitz && x >= 0 && y >= 0 && z >= 0) && (*(m_pData + (z*stepz+y*stepy+x))>0.00001) )
				if(count)
					clashes++; // count clashes
				else
					return 1; // return on first clash

		}
	}

	return clashes;
}

int clash_cg_ibond(vlVolume *vlm, double **co, double **p1co, int bond_number, int nr_atoms_loop, int *marker2, int reverse, int mode, bool count)
{
	vlPoint3f pos;
	float *m_pData = (float*)vlm->getVoxelVoidPtr(vlPoint3ui(0,0,0));
	vlm->getPosition(&pos);
	int x,y,z;
	int clashes = 0; // counting clashes



	if(reverse==0) // forward
	{
		// if (bond_number<2) bond_number=2;

		if( mode == 3)
		{
			 for (int i = nr_atoms_loop - 3; i>=bond_number ; i--)
				if (i%3 != 0) // Checking O & CB
				{
					if (!((i%3 == 1)&&(marker2[i] == 2)))
					{
						x=(int)roundf( ((p1co[i][0]- pos.x()) / unitsx) );
						y=(int)roundf( ((p1co[i][1]- pos.y()) / unitsy) );
						z=(int)roundf( ((p1co[i][2]- pos.z()) / unitsz) );
						if( (x < limitx && y < limity && z < limitz && x >= 0 && y >= 0 && z >= 0) && (*(m_pData + (z*stepz+y*stepy+x))>0.00001) )
							if(count)
								clashes++; // count clashes
							else
								return 1; // return on first clash
					}
				}
		}
		else if( mode == 2)
		{
			for (int i = nr_atoms_loop - 3; i>=bond_number ; i--)
			{
				if ((i%3 == 1)&&(marker2[i] != 2)) // checking CB and non GLY
				{
					x=(int)roundf( ((p1co[i][0]- pos.x()) / unitsx) );
					y=(int)roundf( ((p1co[i][1]- pos.y()) / unitsy) );
					z=(int)roundf( ((p1co[i][2]- pos.z()) / unitsz) );
					if( (x < limitx && y < limity && z < limitz && x >= 0 && y >= 0 && z >= 0) && (*(m_pData + (z*stepz+y*stepy+x))>0.00001) )
						if(count)
							clashes++; // count clashes
						else
							return 1; // return on first clash
				}
			}
		}
		else if( mode == 1) {
			for (int i = nr_atoms_loop - 3; i>=bond_number ; i--)
			{
				if (i%3 == 2) // checking 0
				{
					x=(int)roundf( ((p1co[i][0]- pos.x()) / unitsx) );
					y=(int)roundf( ((p1co[i][1]- pos.y()) / unitsy) );
					z=(int)roundf( ((p1co[i][2]- pos.z()) / unitsz) );
					if( (x < limitx && y < limity && z < limitz && x >= 0 && y >= 0 && z >= 0) && (*(m_pData + (z*stepz+y*stepy+x))>0.00001) )
						if(count)
							clashes++; // count clashes
						else
							return 1; // return on first clash
				}
			}
		}
		// Backbone forward
		for (int i = bond_number; i < nr_atoms_loop - 2; i++)
		{
			x=(int)roundf( ((co[i][0]- pos.x()) / unitsx) );
			y=(int)roundf( ((co[i][1]- pos.y()) / unitsy) );
			z=(int)roundf( ((co[i][2]- pos.z()) / unitsz) );

			if( (x < limitx && y < limity && z < limitz && x >= 0 && y >= 0 && z >= 0) && (*(m_pData + (z*stepz+y*stepy+x))>0.00001) )
				if(count)
					clashes++; // count clashes
				else
					return 1; // return on first clash
		}
	}
	else // reverse
	{
		// if (bond_number<3) bond_number=3;

		if( mode == 3)
		{
			for (int i = bond_number; i < nr_atoms_loop - 2; i++)  //  divisible by 3
				if (i%3 != 2) // Checking O & CB
				{
					if (!((i%3 == 0)&&(marker2[i] == 2)))
					{
						x=(int)roundf( ((p1co[i][0]- pos.x()) / unitsx) );
						y=(int)roundf( ((p1co[i][1]- pos.y()) / unitsy) );
						z=(int)roundf( ((p1co[i][2]- pos.z()) / unitsz) );
						if( (x < limitx && y < limity && z < limitz && x >= 0 && y >= 0 && z >= 0) && (*(m_pData + (z*stepz+y*stepy+x))>0.00001) )
							if(count)
								clashes++; // count clashes
							else
								return 1; // return on first clash
					}
				}
		}
		else if( mode == 2)
		{
			for (int i = bond_number; i < nr_atoms_loop - 2; i++)
			{
				if ((i%3 == 0)&&(marker2[i] != 2)) // checking CB no gly
				{
					x=(int)roundf( ((p1co[i][0]- pos.x()) / unitsx) );
					y=(int)roundf( ((p1co[i][1]- pos.y()) / unitsy) );
					z=(int)roundf( ((p1co[i][2]- pos.z()) / unitsz) );
					if( (x < limitx && y < limity && z < limitz && x >= 0 && y >= 0 && z >= 0) && (*(m_pData + (z*stepz+y*stepy+x))>0.00001) )
						if(count)
							clashes++; // count clashes
						else
							return 1; // return on first clash
				}
			}
		}
		else if( mode == 1)
		{
			for (int i = bond_number; i < nr_atoms_loop - 2; i++)
			{
				if (i%3 == 1) // checking 0
				{
					x=(int)roundf( ((p1co[i][0]- pos.x()) / unitsx) );
					y=(int)roundf( ((p1co[i][1]- pos.y()) / unitsy) );
					z=(int)roundf( ((p1co[i][2]- pos.z()) / unitsz) );
					if( (x < limitx && y < limity && z < limitz && x >= 0 && y >= 0 && z >= 0) && (*(m_pData + (z*stepz+y*stepy+x))>0.00001) )
						if(count)
							clashes++; // count clashes
						else
							return 1; // return on first clash
				}
			}
		}
		// Backbone reverse
		for (int i = nr_atoms_loop - 3; i>=bond_number ; i--)
		{
			x=(int)roundf( ((co[i][0]- pos.x()) / unitsx) );
			y=(int)roundf( ((co[i][1]- pos.y()) / unitsy) );
			z=(int)roundf( ((co[i][2]- pos.z()) / unitsz) );

			if( (x < limitx && y < limity && z < limitz && x >= 0 && y >= 0 && z >= 0) && (*(m_pData + (z*stepz+y*stepy+x))>0.00001) )
				if(count)
					clashes++; // count clashes
				else
					return 1; // return on first clash

		}
	}

	return clashes;
}


// Clash test of backbone without anchors.
int clash_co(vlVolume *vlm, double** co, int nr_atoms_loop)
{
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

	for (int index = 3; index < nr_atoms_loop - 2; index++)
	{
		x=(int)roundf( ((co[index][0]- pos.x()) / unitsx) );
		y=(int)roundf( ((co[index][1]- pos.y()) / unitsy) );
		z=(int)roundf( ((co[index][2]- pos.z()) / unitsz) );

		if( (x < limitx && y < limity && z < limitz && x >= 0 && y >= 0 && z >= 0) && (*(m_pData + (z*stepz+y*stepy+x))>0.00001) )
			return 1;
	}

	return 0;
}

// Check clashes between loop and pocket-map
int clash2(vlVolume *vlm, double** co, double**  p1co,  int nr_atoms_loop, int reverse, bool check_p1co)
{
	char *format = "%-6s%5d %4s %3s %1s%4d%4s%8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n";
	vlPoint3f pos;
	float *m_pData = (float*)vlm->getVoxelVoidPtr(vlPoint3ui(0,0,0));
	vlm->getPosition(&pos);
	float gx, gy, gz;
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
	int index_t0, index_tend;
	if (reverse == 1)	//:TODO: FOR REVERSE PROTOCOL
	{
		index_t0	= 2;
		index_tend = 1;
	}
	else
	{
		index_t0	= 3;
		index_tend = 1;
	}

	for (int index = index_t0; index < nr_atoms_loop - index_tend; index++)
	{

		gx = ((co[index][0]- pos.x()) / unitsx);
		gy = ((co[index][1]- pos.y()) / unitsy);
		gz = ((co[index][2]- pos.z()) / unitsz);

		x=(int)roundf(gx);
		y=(int)roundf(gy);
		z=(int)roundf(gz);

		if( (x < limitx && y < limity && z < limitz && x >= 0 && y >= 0 && z >= 0) )
		{
			if (*(m_pData + (z*stepz+y*stepy+x))>0.00001)
				return 1;
		}

		if(check_p1co) // Check "O"
		{
			gx = ((p1co[index][0]- pos.x()) / unitsx);
			gy = ((p1co[index][1]- pos.y()) / unitsy);
			gz = ((p1co[index][2]- pos.z()) / unitsz);

			x=(int)roundf(gx);
			y=(int)roundf(gy);
			z=(int)roundf(gz);

			if( (x < limitx && y < limity && z < limitz && x >= 0 && y >= 0 && z >= 0) )
			{
				if (*(m_pData + (z*stepz+y*stepy+x))>0.00001)
					return 1;
			}
		}
	}

	return 0;
}


int intra_clash_cg(double **co,double **p1co, double fact_co, double fact_ocb, int nr_atoms_loop,int atomshift, int *marker2, int reverse, int cg_mode)
{
	double d; // required by ISNEAR
	int offset; // Oxigen atom offset
	int shiftp1co = 6; // Sequential shift for O and CB atoms (p1co's)

	if(atomshift < 6)
		shiftp1co = 3;

	if(reverse==0)
		offset = 2;
	else
		offset = 0;

	if(cg_mode > 0)
	{
		if(cg_mode != 2) // O
		{
			for (int i = offset; i < nr_atoms_loop; i+=3) // Just screen "O"s
			{
				for (int j = i + shiftp1co ; j < nr_atoms_loop; j+=3) // Just screen "O"s
					ISNEAR(p1co[i],p1co[j],fact_ocb)
					for (int j = 0; j <= i-atomshift; j++)  // Oi vs. BBj-
						ISNEAR(p1co[i],co[j],fact_co)
						for (int j = i + atomshift; j < nr_atoms_loop; j++)   // Oi vs. BBj+
							ISNEAR(p1co[i],co[j],fact_co)
			}
		}

		// Mon: Consider including sequence information in CB size here...
		//
		if(cg_mode > 1) // CB (the same for forward and backward...)
		{
			for (int i = 1; i < nr_atoms_loop; i+=3) // Just screen CB's
			{
				if(marker2[i] != 2) // NON-GLY
				{
					for (int j = i + shiftp1co ; j < nr_atoms_loop; j+=3) // Just screen CB-CB
						if(marker2[j] != 2)
							ISNEAR(p1co[i],p1co[j],fact_ocb * 1.2)

							for (int j = i + offset + shiftp1co + 1; j < nr_atoms_loop; j+=3) // CB-O
								ISNEAR(p1co[i],p1co[j],fact_ocb)

								for (int j = 0; j <= i-atomshift; j++)   // CB-BB
									ISNEAR(p1co[i],co[j],fact_ocb)

									for (int j = i + atomshift; j < nr_atoms_loop; j++)
										ISNEAR(p1co[i],co[j],fact_ocb)
				}
			}
		}
	}

	for (int i = 0; i < nr_atoms_loop; i++)
		for (int j = i + atomshift; j < nr_atoms_loop; j++) // BBi vs BBj
			ISNEAR(co[i],co[j],fact_co)

			return 0;
}

int intra_clash_cg_ibond(double **co,double **p1co, double fact_co, double fact_ocb,  int bond_number, int nr_atoms_loop,int atomshift, int *marker2, int reverse, int cg_mode)
{
	double d; // required by ISNEAR
	int offset; // Oxigen atom offset
	int shiftp1co = 6; // Sequential shift for O and CB atoms (p1co's)
    int resi=bond_number/3; // CHECK THIS TOO LATE TO THINK


	for (int i = bond_number; i < nr_atoms_loop; i++) {

		for (int j = i - atomshift; j = 0; j--) // BBi vs BBj
					ISNEAR(co[i],co[j],fact_co)

		for (int j = i + atomshift; j < nr_atoms_loop; j++) // BBi vs BBj
			ISNEAR(co[i],co[j],fact_co)
	}

			return 0;
}


/**
 * Project a "co" minimal backbone representation (just N, CA, and C atoms) into a pre-existent map.
 * (Mon made 27/8/2018)
 * vol  = Input/Output volume
 * res  = Resolution [Angstroms]
 * nat  = Number of loop atoms (nr_atoms_loop)
 * co   = N, CA, C, coordinates array
 * mask = Mask of indices (array of integers) for fast cross-correlation. WARNING: A negative integer must indicate the array end!
 */
void co2map_real( vlVolume *vol, float res, int nat, double **co)
{
	vlPoint3ui *position;
	float eweight, value;
	float g_cte, kg;
	int g_width, x[2], y[2], z[2];
	float gx, gy, gz, dist;
	int i, j, k;

	float unit = vol->units().x(); // Warning: only valid for isotropic sampling: units.x=units.y=units.z

	g_cte = res / unit; // Warning: only valid for isotropic sampling: units.x=units.y=units.z
	g_width = ( int )( roundFloat_( g_cte * 3.0 ) );

	if ( g_cte / 2.0 * sqrt( 2.0 * log( 2.0 ) ) / sqrt( ( float )3 ) < 1.0 )
	{
		printf( "Increase the sampling and try again\n" );
		exit( 0 );
	}

	g_cte = ( PDB_PI / g_cte ) * ( PDB_PI / g_cte ); // constant for a gaussian cutoff

	// for res defined at y=1/2 Gausian width  g_cte=log(2)*(M_PI/g_cte)*(M_PI/g_cte);
	kg = pow( g_cte / PDB_PI, 1.5 ); // normalization constant

	// Get the geometric center ("gc") of the atomic coordinates ("co"), excluding anchors!!! (absolute reference frame, in Angstroms)
	float gc[3];
	gc[0] = 0.0;
	gc[1] = 0.0;
	gc[2] = 0.0;
	for (int i = 0; i < nat ; i++) // including anchors
	{
		gc[0] += co[i][0];
		gc[1] += co[i][1];
		gc[2] += co[i][2];
	}
	gc[0] /= nat;
	gc[1] /= nat;
	gc[2] /= nat;

	vlPoint3f ori; // Origin (minimum corner)
	vol->getPosition(&ori); // Origins

	// Gaussian expansion of PDB atoms

	// Get the geometric center ("gc") of the atomic coordinates ("co"), excluding anchors!!!
	int rest; // this indicates atom identity (0=N, 1=CA, 2=C)
	position = new vlPoint3ui( 0, 0, 0 );

	for (int n = 0; n < nat ; n++)
	{
		rest = n%3; // compute rest to know atom type

		// Number of electrons ???
		if(rest == 0) // N atom (NH pseudo-atom)
			eweight = 7;
		else if(rest == 1) // CA atom
			eweight = 6;
		else // C atom (C+O pseudo-atom)
			eweight = 6;

		// fprintf(stderr,"i= %2d  Atom %d (w= %f) --> co %f %f %f\n",n,rest,eweight,co[n][0],co[n][1],co[n][2]);

		// Atom relative position in number of voxels
		gx = ((co[n][0]- ori.x()) / unit) ;
		gy = ((co[n][1]- ori.y()) / unit) ;
		gz = ((co[n][2]- ori.z()) / unit) ;

		// Calculate kernel gaussian size in grid units
		x[0] = ( int )roundFloat_( gx - g_width );
		x[1] = ( int )roundFloat_( gx + g_width );
		y[0] = ( int )roundFloat_( gy - g_width );
		y[1] = ( int )roundFloat_( gy + g_width );
		z[0] = ( int )roundFloat_( gz - g_width );
		z[1] = ( int )roundFloat_( gz + g_width );

		// Gaussian expansion
		for ( k = z[0]; k < z[1]; k++ )
			for ( j = y[0]; j < y[1]; j++ )
				for ( i = x[0]; i < x[1]; i++ )
				{
					// compute distance
					dist = ( ( i - gx ) * ( i - gx ) + ( j - gy ) * ( j - gy ) + ( k - gz ) * ( k - gz ) );

					// get position
					position->x( i );
					position->y( j );
					position->z( k );
					vol->getVoxel( *position, value );

					// Recipe to get direct access "index" (from getVoxel in vlvoldata_linear.h)
					// Return the voxel data
					// return *(this->m_pData + (position.z()*this->m_step.z()+position.y()*this->m_step.y()+position.x()));

					// set position
					value += kg * eweight * exp( -dist * g_cte );
					vol->setVoxel( *position, value );
				}
	}
}


/**
 * Project a "co" and (optionally) "p1co" backbone representation (just N, CA, C, O, and CB atoms) into a pre-existent map.
 * (Mon made 29/8/2018)
 * vol           = Input/Output volume
 * res           = Resolution [Angstroms]
 * nat           = Number of loop atoms (nr_atoms_loop)
 * co            = N, CA, and C coordinates array
 * pco           = (optional) O and CB coordinates array
 * cg            = (optional) Coarse graining level
 * residuemarker = (optional) Array with the residue types. GLY=1 (size= loop number of atoms, N, CA, C only!)
 */
void co2map_real( vlVolume *vol, float res, int nat, double **co, double **pco, int cg, int *residuemarker)
{
	vlPoint3ui *position;
	float eweight, value;
	float g_cte, kg;
	int g_width, x[2], y[2], z[2];
	float gx, gy, gz, dist;
	int i, j, k;

	float unit = vol->units().x(); // Warning: only valid for isotropic sampling: units.x=units.y=units.z

	g_cte = res / unit; // Warning: only valid for isotropic sampling: units.x=units.y=units.z
	g_width = ( int )( roundFloat_( g_cte * 3.0 ) );

	if ( g_cte / 2.0 * sqrt( 2.0 * log( 2.0 ) ) / sqrt( ( float )3 ) < 1.0 )
	{
		printf( "Increase the sampling and try again\n" );
		exit( 0 );
	}

	g_cte = ( PDB_PI / g_cte ) * ( PDB_PI / g_cte ); // constant for a gaussian cutoff

	// for res defined at y=1/2 Gausian width  g_cte=log(2)*(M_PI/g_cte)*(M_PI/g_cte);
	kg = pow( g_cte / PDB_PI, 1.5 ); // normalization constant

	// Get the geometric center ("gc") of the atomic coordinates ("co"), excluding anchors!!! (absolute reference frame, in Angstroms)
	float gc[3];
	gc[0] = 0.0;
	gc[1] = 0.0;
	gc[2] = 0.0;
	for (int i = 0; i < nat ; i++) // including anchors
	{
		gc[0] += co[i][0];
		gc[1] += co[i][1];
		gc[2] += co[i][2];
	}
	gc[0] /= nat;
	gc[1] /= nat;
	gc[2] /= nat;

	vlPoint3f ori; // Origin (minimum corner)
	vol->getPosition(&ori); // Origins

	// Get the geometric center ("gc") of the atomic coordinates ("co"), excluding anchors!!!
	int rest; // this indicates atom identity (0=N, 1=CA, 2=C)
	position = new vlPoint3ui( 0, 0, 0 );

	// Gaussian expansion of PDB atoms
	for (int n = 0; n < nat ; n++) // Screen minimal backbone atoms (N, CA, C)
	{
		rest = n%3; // compute rest to know atom type

		// Number of electrons ???
		if(rest == 0) // N atom (NH pseudo-atom)
			eweight = 7;
		else if(rest == 1) // CA atom
			eweight = 6;
		else // C atom (C+O pseudo-atom)
			eweight = 6;

		// fprintf(stderr,"i= %2d  Atom %d (w= %f) --> co %f %f %f\n",n,rest,eweight,co[n][0],co[n][1],co[n][2]);

		// Atom relative position in number of voxels
		gx = ((co[n][0]- ori.x()) / unit) ;
		gy = ((co[n][1]- ori.y()) / unit) ;
		gz = ((co[n][2]- ori.z()) / unit) ;

		// Calculate kernel gaussian size in grid units
		x[0] = ( int )roundFloat_( gx - g_width );
		x[1] = ( int )roundFloat_( gx + g_width );
		y[0] = ( int )roundFloat_( gy - g_width );
		y[1] = ( int )roundFloat_( gy + g_width );
		z[0] = ( int )roundFloat_( gz - g_width );
		z[1] = ( int )roundFloat_( gz + g_width );

		// Gaussian expansion
		for ( k = z[0]; k < z[1]; k++ )
			for ( j = y[0]; j < y[1]; j++ )
				for ( i = x[0]; i < x[1]; i++ )
				{
					// compute distance
					dist = ( ( i - gx ) * ( i - gx ) + ( j - gy ) * ( j - gy ) + ( k - gz ) * ( k - gz ) );

					// get position
					position->x( i );
					position->y( j );
					position->z( k );
					vol->getVoxel( *position, value );

					// Recipe to get direct access "index" (from getVoxel in vlvoldata_linear.h)
					// Return the voxel data
					// return *(this->m_pData + (position.z()*this->m_step.z()+position.y()*this->m_step.y()+position.x()));

					// set position
					value += kg * eweight * exp( -dist * g_cte );
					vol->setVoxel( *position, value );
				}
	}

	// "pco"'s CB
	if(cg >= 2) //  2: N,CA,C,CB    3: N,CA,C,O,CB
	{
		eweight = 6; // CB density

		for (int n = 1; n < nat ; n += 3) // Screen CB's
			if(residuemarker[n/3] != 1) // If non-GLY, it has CB
			{
				// fprintf(stderr,"i= %2d  Atom %d (w= %f) --> co %f %f %f\n",n,rest,eweight,co[n][0],co[n][1],co[n][2]);

				// Atom relative position in number of voxels
				gx = ((pco[n][0]- ori.x()) / unit) ;
				gy = ((pco[n][1]- ori.y()) / unit) ;
				gz = ((pco[n][2]- ori.z()) / unit) ;

				// Calculate kernel gaussian size in grid units
				x[0] = ( int )roundFloat_( gx - g_width );
				x[1] = ( int )roundFloat_( gx + g_width );
				y[0] = ( int )roundFloat_( gy - g_width );
				y[1] = ( int )roundFloat_( gy + g_width );
				z[0] = ( int )roundFloat_( gz - g_width );
				z[1] = ( int )roundFloat_( gz + g_width );

				// Gaussian expansion
				for ( k = z[0]; k < z[1]; k++ )
					for ( j = y[0]; j < y[1]; j++ )
						for ( i = x[0]; i < x[1]; i++ )
						{
							// compute distance
							dist = ( ( i - gx ) * ( i - gx ) + ( j - gy ) * ( j - gy ) + ( k - gz ) * ( k - gz ) );

							// get position
							position->x( i );
							position->y( j );
							position->z( k );
							vol->getVoxel( *position, value );

							// set position
							value += kg * eweight * exp( -dist * g_cte );
							vol->setVoxel( *position, value );
						}
			}
	}

	// "pco"'s O
	if(cg == 1 || cg == 3) //  1: N,CA,C,O    3: N,CA,C,O,CB
	{
		eweight = 8; // O density

		for (int n = 2; n < nat ; n += 3) // Screen O's (all residues have O)
		{

			// fprintf(stderr,"i= %2d  Atom %d (w= %f) --> co %f %f %f\n",n,rest,eweight,co[n][0],co[n][1],co[n][2]);

			// Atom relative position in number of voxels
			gx = ((pco[n][0]- ori.x()) / unit) ;
			gy = ((pco[n][1]- ori.y()) / unit) ;
			gz = ((pco[n][2]- ori.z()) / unit) ;

			// Calculate kernel gaussian size in grid units
			x[0] = ( int )roundFloat_( gx - g_width );
			x[1] = ( int )roundFloat_( gx + g_width );
			y[0] = ( int )roundFloat_( gy - g_width );
			y[1] = ( int )roundFloat_( gy + g_width );
			z[0] = ( int )roundFloat_( gz - g_width );
			z[1] = ( int )roundFloat_( gz + g_width );

			// Gaussian expansion
			for ( k = z[0]; k < z[1]; k++ )
				for ( j = y[0]; j < y[1]; j++ )
					for ( i = x[0]; i < x[1]; i++ )
					{
						// compute distance
						dist = ( ( i - gx ) * ( i - gx ) + ( j - gy ) * ( j - gy ) + ( k - gz ) * ( k - gz ) );

						// get position
						position->x( i );
						position->y( j );
						position->z( k );
						vol->getVoxel( *position, value );

						// set position
						value += kg * eweight * exp( -dist * g_cte );
						vol->setVoxel( *position, value );
					}
		}
	}
}

// Required by "co2map_real"
int roundFloat_( float x )
{

	int aux1, aux2;
	float dist1, dist2;
	aux1 = ( int )ceil( x );
	aux2 = ( int )floor( x );
	dist1 = x - aux1;
	if ( dist1 < 0 ) dist1 *= -1;
	dist2 = x - aux2;
	if ( dist2 < 0 ) dist2 *= -1;

	if ( dist2 < dist1 ) return aux2;
	else
		return aux1;
}


// Write Multi-PDB file (N, CA, C, O, H, HA2)
void multiPDB(double **co,double **p1co,double **p2co,int number_loops,int number_residues, char *name, char **secu, int start)
{	//:INFO: debugged once ---> looks Ok ---> only nomenclature H atoms*/

	FILE* p_file;
	p_file = fopen(name,"w");

	if (p_file == NULL)
		printf("%s opening failed\n",name);

	char *format = "%-6s%5d %4s %3s %1s%4d%4s%8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n";		//:REMARK:	PDB format---exact format

	int index;
	for (int i = 0; i < number_loops; i++)
	{
		fprintf(p_file,"%-6s%8d\n","MODEL ",i+1);

		for (int j = 0; j < number_residues; j++)	/*iterate # residues*/ /*local j used and override global j*/
		{
			index = (i*number_residues*3)+(j*3);	/*atoms physically different from total shift of 4*/
			fprintf(p_file,format,"ATOM",j+ 1,"N  ",secu[j],"A",start+j," ",co[index][0],co[index][1],co[index][2],0.00,0.00);					/*backbone*/
			fprintf(p_file,format,"ATOM",j+ 2,"CA ",secu[j],"A",start+j," ",co[index+1][0],co[index+1][1],co[index+1][2],0.00,0.00);
			fprintf(p_file,format,"ATOM",j+ 3,"C  ",secu[j],"A",start+j," ",co[index+2][0],co[index+2][1],co[index+2][2],0.00,0.00);
			fprintf(p_file,format,"ATOM",j+ 4,"O  ",secu[j],"A",start+j," ",p1co[index+2][0],p1co[index+2][1],p1co[index+2][2],0.00,0.00);		/*O atom*/
			fprintf(p_file,format,"ATOM",j+ 5,"H  ",secu[j],"A",start+j," ",p1co[index][0],p1co[index][1],p1co[index][2],0.00,0.00);			/*H atom nitrogen*/
			fprintf(p_file,format,"ATOM",j+ 6,"HA2",secu[j],"A",start+j," ",p1co[index+1][0],p1co[index+1][1],p1co[index+1][2],0.00,0.00);		/*H atom 1 nitrogen CA*/
		}

		fprintf(p_file,"%-6s\n","ENDMDL");
	}

	fclose(p_file);
}

// Write Multi-PDB file (N, CA, C, O, H, HA2)
void multiPDB_cg(double **co,double **p1co,int number_loops,int number_residues, char *name, char **secu, int start, int mode, char chain)
{	//:INFO: debugged once ---> looks Ok ---> only nomenclature H atoms*/

	FILE* p_file;
	p_file = fopen(name,"w");

	if (p_file == NULL)
		printf("%s opening failed\n",name);

	char *format = "%-6s%5d %4s %3s %c%4d%4s%8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n";		//:REMARK:	PDB format---exact format

	int index;
	for (int i = 0; i < number_loops; i++)
	{
		int atom=1; // Same numbering for all loops
		fprintf(p_file,"%-6s%8d\n","MODEL ",i+1);
		for (int j = 0; j < number_residues; j++)	/*iterate # residues*/ /*local j used and override global j*/
		{
			index = (i*number_residues*3)+(j*3);	/*atoms physically different from total shift of 4*/
			fprintf(p_file,format,"ATOM",atom,"N  ",secu[j],chain,start+j," ",co[index][0],co[index][1],co[index][2],1.00,0.00);					/*backbone*/
			atom++;
			fprintf(p_file,format,"ATOM",atom,"CA ",secu[j],chain,start+j," ",co[index+1][0],co[index+1][1],co[index+1][2],1.00,0.00);
			atom++;
			fprintf(p_file,format,"ATOM",atom,"C  ",secu[j],chain,start+j," ",co[index+2][0],co[index+2][1],co[index+2][2],1.00,0.00);
			atom++;
			if(mode>0)
			{
				if(mode==1 || mode==3)
				{
					fprintf(p_file,format,"ATOM",atom,"O  ",secu[j],chain,start+j," ",p1co[index+2][0],p1co[index+2][1],p1co[index+2][2],1.00,0.00);		/*O atom*/
					atom++;
				}
				if(mode==2 || mode==3)
				{
					if (strcmp(secu[j],"GLY")!=0)
						fprintf(p_file,format,"ATOM",atom,"CB ",secu[j],chain,start+j," ",p1co[index+1][0],p1co[index+1][1],p1co[index+1][2],1.00,0.00);		/*H atom 1 nitrogen CA*/
					atom++;
				}
			}
		}
		fprintf(p_file,"%-6s\n","ENDMDL");
	}
	fclose(p_file);
}
// Write Multi-PDB file (N, CA, C)
void multiPDB_write(double **co,int number_loops,int number_residues_loop, char *name)
{
	FILE* p_file;
	p_file = fopen(name,"w");

	if (p_file == NULL)
	{
		printf("%s opening failed\n",name);
		printf("Program execution aborted ...\n");
		exit (1);
	}

	int index;
	char *format = "%-6s%5d %4s %3s %1s%4d%4s%8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n";	  //:REMARK:	PDB format---exact format critical

	for (int i = 0; i < number_loops-1; i++)
	{
		fprintf(p_file,"%-6s%8d\n","MODEL ",i+1);

		for (int j = 0; j < number_residues_loop; j++)
		{
			index = (i*number_residues_loop*3) + (j*3);
			fprintf(p_file,format,"ATOM",(j*3)+1,"N  ","ALA","A",j+1," ",co[index][0],co[index][1],co[index][2],1.00,0.00);
			fprintf(p_file,format,"ATOM",(j*3)+2,"CA ","ALA","A",j+1," ",co[index+1][0],co[index+1][1],co[index+1][2],1.00,0.00);
			fprintf(p_file,format,"ATOM",(j*3)+3,"C  ","ALA","A",j+1," ",co[index+2][0],co[index+2][1],co[index+2][2],1.00,0.00);
		}
		fprintf(p_file,"%-6s\n","ENDMDL");
	}

	fclose(p_file);
}

//************************************************************
// NEW functions (ALEX)
//************************************************************

// Computes the accumulated probability for a 2D probability density function (PDF) map.
float *accunorm(float **pdf, int size)
{
	float *acc;
	acc = (float *) malloc( sizeof(float) * size * size);
	double sum=0.0;

	// Accumulate
	for( int i=0; i<size; i++) // rows
		for( int j=0; j<size; j++) // cols
		{
			sum += pdf[i][j];
			acc[i*size + j] = sum;
		}

	// Normalize
	for ( int i=0; i<size*size; i++)
		acc[i] /= sum;
	return acc;
}

// Now read_dunbrack() is in rama.cpp
// Reads all Dunbrack's PDFs from file (created with "rang" program)
//float *****read_dunbrack(char *file, int *p_size, float *p_step)
//{
//	FILE *f_map; // File handler
//	if( !(f_map = fopen(file,"rb")) ) // Read binary file
//	{
//		printf("Reading input error! file %s\n",file);
//		exit(1);
//	}
//
//	// Write header with "size" and "step"
//	fread(p_size,sizeof(int),1,f_map);
//	fread(p_step,sizeof(float),1,f_map);
//
//	// Memory allocation
//	float *****pdfs; // 5*pointer containing all pdfs
//	int mem=0;
//	pdfs = (float *****) malloc( sizeof(float ****) * NAAS); // Central aminoacid level
//	mem += sizeof(float ****) * NAAS;
//	for( int c=0; c<NAAS; c++)
//	{
//		pdfs[c] = (float ****) malloc( sizeof(float ***) * 2); // Whether it is a Left-side (0) or Right-side (1) aminoacid
//		mem += sizeof(float ***) * 2;
//		for( int lr=0; lr<2; lr++)
//		{
//			pdfs[c][lr] = (float ***) malloc( sizeof(float **) * (NAAS+1)); // Left/Right-side aminoacid (+1 includes "ALL")
//			mem += sizeof(float **) * NAAS+1;
//			for( int s=0; s<NAAS+1; s++)
//			{
//				pdfs[c][lr][s] = (float **) malloc( sizeof(float *) * (*p_size)); // Pointer array
//				mem += sizeof(float *) * (*p_size);
//				for( int i=0; i<*p_size; i++)
//				{
//					pdfs[c][lr][s][i] = (float *) malloc( sizeof(float) * (*p_size)); // data
//					mem += sizeof(float) * (*p_size);
//					for( int j=0; j<*p_size; j++)
//						pdfs[c][lr][s][i][j] = 0.0; // initialize pdf
//				}
//			}
//		}
//	}
//	// fprintf(stdout,"Msg(read_dunbrack): Allocated %d bytes!\n",mem);
//
//	for( int c=0; c<NAAS; c++)
//		for( int lr=0; lr<2; lr++)
//			for( int s=0; s<NAAS+1; s++)
//				for( int i=0; i<*p_size; i++)
//					fread(pdfs[c][lr][s][i], sizeof(float), *p_size, f_map); // write data
//
//	fclose(f_map);
//	return pdfs;
//}

// Generate Phi/Psi dihedral angle pair "in Ramachandran".
void rang_index2D_dun(float *pdf_accum, int size, int step, float *phi, float *psi) // aleator --> rg
{
	float random = (float) rg->Random(); // Output random float number in the interval 0 <= x < 1
	int i=0;
	while(pdf_accum[i] < random && i<size*size)
		i++;

	*phi = ((i/size)*step-180+((float) step)/2)*M_PI/180.0; // in radians
	*psi = ((i%size)*step-180+((float) step)/2)*M_PI/180.0; // in radians

}

// Now MatrixProduct() is in rama.cpp
// Element-wise multiplication of two squared matrices (size x size)
//float **MatrixProduct(float **M1, float **M2, int size)
//{
//	float **prod;
//	prod = (float **) malloc( sizeof(float *) * size);
//
//	for( int i=0; i<size; i++)
//	{
//		prod[i] = (float *) malloc( sizeof(float) * size);
//		for( int j=0; j<size; j++)
//			prod[i][j] = 0.0; // initialize
//	}
//
//	for( int i=0; i<size; i++)
//		for( int j=0; j<size; j++)
//			prod[i][j] = M1[i][j]*M2[i][j];
//
//	return prod;
//}

// Element-wise division of two squared matrices (size x size)
//float **MatrixDivision(float **M1, float **M2, int size)
//{
//	float **divgc;
//	divgc = (float **) malloc( sizeof(float *) * size);
//	for( int i=0; i<size; i++)
//	{
//		divgc[i] = (float *) malloc( sizeof(float) * size);
//		for( int j=0; j<size; j++)
//			divgc[i][j] = 0.0; // initialize
//	}
//
//	for( int i=0; i<size; i++)
//		for( int j=0; j<size; j++)
//			divgc[i][j] = M1[i][j]/M2[i][j];
//
//	return divgc;
//}

// map_gen() is now in rama.cpp
// Generates the final PDF from Left (L) and Right (R) PDFs for a given Central (C) aminoacid and Dunbrack's tabulated PDFs.
// (it allocates output map memory)
//float **map_gen(int C, int L, int R, float *****pdfs, int size)
//{
//	float **prod;
//	float **map;
//
//	prod = MatrixProduct(pdfs[C][0][L],pdfs[C][1][R],size);
//	// Bug (4/1/2019)
//	// map = MatrixDivision(prod,pdfs[C][1][8],size); // ALL-R [21][1]
//	map = MatrixDivision(prod,pdfs[C][1][21],size); // L or R = ALL --> [21] (both seem equivalent..., see "rang" program)
//
//	for(int i=0; i<size; i++)
//		free(prod[i]);
//	free(prod);
//
//	return map;
//}

// Generates a binary map given an absolute threshold "thr"
// (it allocates output map memory)
float **map_bin(float **map, float thr , int size)
{
	float **mapb;
	mapb = (float **) malloc( sizeof(float *) * size);
	for( int i=0; i<size; i++)
		mapb[i] = (float *) malloc( sizeof(float) * size);

	for(int i=0;i<size;i++)
		for(int j=0;j<size;j++)
			if(map[i][j] > thr)
				mapb[i][j] = 1;
			else
				mapb[i][j] = 0;

	return mapb;
}

// Generates a binary map given a relative threshold ( thr=]0,1[ )
// Warning, map values are sorted first (not using quick-sort algorithm, thus not very efficient...)
// (it allocates output map memory)
float **map_bin2(float **map, float thr, int size)
{
	double acum=0.0;
	float d_max;
	int i_max,j_max;
	float **mapb;

	mapb = (float **) malloc( sizeof(float *) * size);
	for( int i=0; i<size; i++)
	{
		mapb[i] = (float *) malloc( sizeof(float) * size);
		for(int j=0; j<size; j++)
			mapb[i][j]=0.0; // initialize binary map
	}

	while (acum < thr)
	{
		d_max = 0.0;
		for(int i=0;i<size;i++)
			for(int j=0;j<size;j++)
				if(mapb[i][j] < 1.0)
					if(map[i][j] > d_max)
					{
						d_max = map[i][j];
						i_max = i;
						j_max = j;
					}
		acum += d_max;
		mapb[i_max][j_max] = 1.0;
	}

	return mapb;
}

// Computes the cutoff value that accounts for the selected probability using a relative threshold ( thr=]0,1[ )
// Warning, map values are sorted first (not using quick-sort algorithm, thus not very efficient...)
// (it allocates output map memory)
float map_cutoff(float **map, float thr, int size)
{
	double acum=0.0;
	float d_max;
	int i_max,j_max;
	float **mapb;

	mapb = (float **) malloc( sizeof(float *) * size);
	for( int i=0; i<size; i++)
	{
		mapb[i] = (float *) malloc( sizeof(float) * size);
		for(int j=0; j<size; j++)
			mapb[i][j]=0.0; // initialize binary map
	}

	while (acum <= thr)
	{
		d_max = 0.0;
		for(int i=0;i<size;i++)
			for(int j=0;j<size;j++)
				if(mapb[i][j] < 1.0)
					if(map[i][j] > d_max)
					{
						d_max = map[i][j];
						i_max = i;
						j_max = j;
					}
		acum += d_max;
		mapb[i_max][j_max] = 1.0;
	}

	for( int i=0; i<size; i++)
		free(mapb[i]);
	free(mapb);

	return d_max;
}

// Generates a binary map given a relative threshold ( thr=]0,1[ )
// Warning, map values are sorted first (not using quick-sort algorithm, thus not very efficient...)
// (it allocates output map memory)
float **map_bin3(float **map, float thr, int size)
{
	float *ord_array;
	float **out_map;
	double aux = 0.0, acum = 0.0;
	bool cambio=true;
	float probthr;
	int i=0;
	int j=0;

	out_map = (float **) malloc( sizeof(float *) *size);
	for( int i=0; i<size; i++)
		out_map[i] = (float *) malloc( sizeof(float) * size);

	ord_array = (float *)malloc( sizeof(float) * size * size);
	for( int i=0; i<size*size; i++)
		ord_array[i]=0.0;
	for( int i=0; i<size; i++)
		for( int j=0; j<size; j++)
			ord_array[i*size + j] = map[i][j];

	while(cambio)
	{
		cambio = false;
		for( int i=0;i<size*size;i++)
			if (ord_array[i] < ord_array[i+1])
			{
				aux = ord_array[i];
				ord_array[i] = ord_array[i+1];
				ord_array[i+1] = aux;
				cambio = true;
			}
	}

	for( int i=0;i<size*size;i++)
		if (acum < thr)
			acum += ord_array[i];
		else
		{
			probthr = ord_array[i-1]; // Get the probability threshold for the current distribution
			break;
		}

	for( int i=0; i<size; i++)
		for ( int j=0; j<size; j++)
			if (map[i][j] > probthr)
				out_map[i][j] = 1.0; // Binarize the current map
			else
				out_map[i][j] = 0.0; // Binarize the current map

	return out_map;
}

// Normalizes "in situ" a map (size x size) to sum=1.0. It does not allocate memory.
//void norm_map(float **map, int size)
//{
//	double acum=0.0;
//
//	// Get sum
//	for(int i=0;i<size;i++)
//		for(int j=0;j<size;j++)
//			acum += map[i][j];
//
//	// Normalize
//	for(int i=0;i<size;i++)
//		for(int j=0;j<size;j++)
//			map[i][j] /= acum;
//}

// Set to zero those elements with values higher or lower than some threshold. (not allocates memory)
// less = true --> lower values than "thr" will be set to zero, otherwise higher values than "thr" will be set to zero
void map_setzero(float **map, int size, float thr, bool less)
{
	if(less)
	{
		for(int i=0;i<size;i++)
			for(int j=0;j<size;j++)
				if(map[i][j] < thr)
					map[i][j] = 0.0;
	}
	else
	{
		for(int i=0;i<size;i++)
			for(int j=0;j<size;j++)
				if(map[i][j] > thr)
					map[i][j] = 0.0;
	}
}

// Copy a 2D map of "size" side elements. (it allocates output map memory)
float **map_copy(float **map, int size)
{
	float **mapb;

	mapb = (float **) malloc( sizeof(float *) * size);
	for( int i=0; i<size; i++)
	{
		mapb[i] = (float *) malloc( sizeof(float) * size);
		for(int j=0; j<size; j++)
			mapb[i][j]=map[i][j]; // copy map
	}

	return mapb;
}

// Copy a 2D map of "size" side elements. (it allocates output map memory)
float **map_copyh3(float ***map, int k, int size)
{
	float **mapb;
	double suma;
	suma=0.0;
	mapb = (float **) malloc( sizeof(float *) * size);
	for( int i=0; i<size; i++)
	{
		mapb[i] = (float *) malloc( sizeof(float) * size);
		for(int j=0; j<size; j++)
		{
			mapb[i][j]=map[k][i][j]; // copy map
			suma+=mapb[i][j];
		}
	}

	for( int i=0; i<size; i++)
		for(int j=0; j<size; j++)
			mapb[i][j]/=suma;



	return mapb;
}

// Copy a 2D map of "size" side elements. (it allocates output map memory)
float **map_mergeh3(float ***map, float **map2, int k, int size, float factor)
{
	float **mapb;
	double suma;
	suma=0.0;
	mapb = (float **) malloc( sizeof(float *) * size);
	for( int i=0; i<size; i++)
	{
		mapb[i] = (float *) malloc( sizeof(float) * size);
		for(int j=0; j<size; j++)
		{
			mapb[i][j]=factor*map[k][i][j]+(1-factor)*map2[i][j]; // copy map
			suma+=mapb[i][j];
		}
	}

	for( int i=0; i<size; i++)
		for(int j=0; j<size; j++)
			mapb[i][j]/=suma;


	return mapb;
}



// Writes a pdf map in some file (size x size)
//void show_map(float **map,int size, char *file)
//{
//	FILE *f_file;
//	if( !(f_file = fopen(file,"w")) )
//	{
//		printf("Reading input error! file %s\n",f_file);
//		exit(1);
//	}
//
//	for(int i=0;i<size;i++)
//		for(int j=0;j<size;j++)
//			fprintf(f_file, "%5d %5d %12e \n", i, j, map[i][j]);
//
//	fclose(f_file);
//}

// Plots a 2D map and (optionally) phi and psi values.
void paint_map_point(float **map,int size, char *file, int indexphi, int indexpsi)
{
	float ang_i,ang_j, ang_indexphi, ang_indexpsi;
	FILE *f_file;
	if( !(f_file = fopen(file,"w")) )
	{
		printf("Reading input error! file %s\n",f_file);
		exit(1);
	}

	for(int i=0;i<size;i++)
	{
		for(int j=0;j<size;j++)
		{
			ang_i = (180.0/M_PI)*(index2angle(i,size));
			ang_j = (180.0/M_PI)*(index2angle(j,size));

			if((i == indexphi) && (j == indexpsi))
			{
				ang_indexphi = (180.0/M_PI)*(index2angle(indexphi,size));
				ang_indexpsi = (180.0/M_PI)*index2angle(indexpsi,size);
				fprintf(f_file, "%10.5f %10.5f %10.5f \n", ang_i, ang_j, 5.0);
				// fprintf(stdout, "Map: %10.5f %10.5f %10.5f \n", ang_i, ang_j, map[i][j]);
			}
			else
			{
				fprintf(f_file, "%10.5f %10.5f %10.5f \n", ang_i, ang_j, map[i][j]);
			}
		}
	}
	fclose(f_file);
}

// Plot the Ramachandran maps for all aminoacids with all dihedral angles.
// INPUT:  map --> Ramachandran map
//         size --> Map size
//         dihedral_angle --> Dihedral angles array
//         nres --> Number of residues including anchors (nr_loop+2)
//         reverse --> 0= forward, 1= backward
//         base --> Output files base name
// OUTPUT: One Ramachandran map per aminoacid with current dihedral angles represented therein.
void paint_maps_points(float ***maps, int size, int nres, char *base, double *dihedral_angle, int reverse )
{
	bool debug = false;
	int indexphi, indexpsi, index;
	double anglephi, anglepsi;
	char file[100];
	char num[10];
	int j;

	for(int i = 0; i < nres; i++ )
	{
		if(dihedral_angle != NULL)
		{
			if(reverse == 0)
			{
				anglephi = dihedral_angle[i*3+2];
				if(i == nres-1)
					anglepsi = dihedral_angle[0];
				else
					anglepsi = dihedral_angle[i*3+3];
				index = i;
			}
			else
			{
				anglepsi = dihedral_angle[i*3+2];
				if(i == nres-1)
					anglephi = dihedral_angle[0];
				else
					anglephi = dihedral_angle[i*3+3];
				index = (nres - 1) - i;
			}

			indexphi = angle2index(anglephi,size);
			indexpsi = angle2index(anglepsi,size);

			if(debug)
				fprintf(stderr,"AA %d (absolute): anglephi= %f   anglepsi= %f   indexphi= %d indexpsi = %d \n", index, anglephi*180.0/M_PI, anglepsi*180.0/M_PI,indexphi,indexpsi);
		}
		else
			index = i;

		sprintf(num,"%02d",index+1);
		strcpy(file,base);
		strcat(file,num);
		strcat(file,".txt");

		if(dihedral_angle != NULL)
			paint_map_point(maps[index], size, file, indexphi,indexpsi);
		else
			paint_map_point(maps[index], size, file);
	}
}

// Plot the Ramachandran maps for all aminoacids with all dihedral angles. G
// INPUT:  dihedral_angle --> Dihedral angles array of the "nsol" solutions (nsol*nres*3)
//         size --> Map size
//         nres --> Number of residues including anchors (nr_loop+2)
//         reverse --> 0= forward, 1= backward
//         base --> Output files base name
//         nsol --> Number of solutions
// OUTPUT: One Ramachandran map per aminoacid with current dihedral angles represented therein.
void paint_all_dihedrals(double *dihedral_angle, int size, int nres, int reverse, char *base, int nsol )
{
	bool debug = false;
	int indexphi, indexpsi, index;
	double anglephi, anglepsi;
	char file[500];
	char num[10];
	int j;
	int natoms = nres*3;
	float ang_i,ang_j, ang_indexphi, ang_indexpsi;
	FILE *f_file;

	for(int i = 0; i < nres; i++ ) // screen aminoacids
	{
		if(reverse == 0)
			index = i;
		else
			index = (nres - 1) - i;

		sprintf(num,"%02d",index+1);
		strcpy(file,base);
		strcat(file,num);
		strcat(file,".txt");

		if( !(f_file = fopen(file,"w")) )
		{
			printf("Reading input error! file %s\n",f_file);
			exit(1);
		}

		for(int k=0; k<nsol; k++) // screen solutions
		{
			if(reverse == 0) // Forward
			{
				anglephi = dihedral_angle[k*natoms+i*3+2];
				if(i == nres-1)
					anglepsi = dihedral_angle[k*natoms];
				else
					anglepsi = dihedral_angle[k*natoms+i*3+3];
			}
			else // Backward
			{
				anglepsi = dihedral_angle[k*natoms+i*3+2];
				if(i == nres-1)
					anglephi = dihedral_angle[k*natoms];
				else
					anglephi = dihedral_angle[k*natoms+i*3+3];
			}
			fprintf(f_file, "%10.5f %10.5f\n", anglephi*180/M_PI, anglepsi*180/M_PI);
		}
		fclose(f_file);
	}
}

// CHECK (and count) ALL DIHEDRAL ANGLES
// count_bad --> true= return the total number of aminoacids with bad dihedrals, false= return 1 as soon as any aminoacid has bad dihedrals (default)
int check_angles(float ***binmap, double *dihedral, int size, int nr_atoms_loop, int reverse, bool count_bad, int *marker, float *cutoffs)
{
	int indexpsi,indexphi;
	int bad_angles = 0;
	int nres = nr_atoms_loop/3;
	int j;
	float cutoff = 0.5;

	if (reverse == 0) // Forward
	{// FORWARD --> PHI-PSI-PHI-PSI...
		for(int i = 0; i < nres; i++)
		{
			if(cutoffs != NULL)
				cutoff = cutoffs[i];

			indexphi = angle2index(dihedral[i*3+2], size);

			if(i == nres-1)
				indexpsi = angle2index(dihedral[0], size);
			else
				indexpsi = angle2index(dihedral[i*3+3], size);

			if (binmap[i][indexphi][indexpsi] < cutoffs[i])
			{
				// fprintf(stderr,"Rama failed at i=%d indexphi=%d indexpsi=%d   binmap=%f < cutoff=%f  rev=%d \n",i,indexphi,indexpsi,binmap[i][indexphi][indexpsi],cutoff,reverse);
				if(!count_bad)
					return 1;
				bad_angles++; // count bad pair of angles
			}
		}
	}
	else // Backward
	{
		for(int i = 0; i < nres; i++)
		{// REVERSE --> PSI-PHI-PSI-PHI...
			j = (nres-1) - i;

			if(cutoffs != NULL)
				cutoff = cutoffs[j];

			indexpsi = angle2index(dihedral[i*3+2], size);

			if(i == nres-1)
				indexphi = angle2index(dihedral[0], size);
			else
				indexphi = angle2index(dihedral[i*3+3], size);

			if (binmap[j][indexphi][indexpsi] < cutoff)
			{
				//	fprintf(stderr,"Rama failed at j=%d indexphi=%d indexpsi=%d   binmap=%f < cutoff=%f  rev=%d \n",j,indexphi,indexpsi,binmap[j][indexphi][indexpsi],cutoff,reverse);
				if(!count_bad)
					return 1;
				bad_angles++; // count bad pair of angles
			}
		}
	}

	return bad_angles;
}

// CHECK JUST A GIVEN DIHEDRAL ANGLES PAIR (returns 1 if bad angles or 0 if OK)
int check_angles_final(float **binmap, double ang_phi, double ang_psi, int size, int nr_atoms_loop, float cutoff)
{
	int indexphi,indexpsi;

	indexphi = angle2index(ang_phi,size);
	indexpsi = angle2index(ang_psi,size);

	if (binmap[indexphi][indexpsi] < cutoff)
		return 1; // If fails

	return 0; // If OK
}

// Given Phi (or Psi), it tells whether there is any valid Psi (or Phi) "in Ramachandran".
// angle --> Given Phi or Psi (fixed)
// isphi --> =true if "angle" is Phi, =false otherwise.
// count_bad --> =true, counts the number of valid Ramachandran map positions (default = false)
int check1angle(float **binmap, double angle, bool isphi, int size, int reverse, bool count_bad, float cutoff)
{
	int nbad = 0; // counts bads
	int index = angle2index(angle,size);

	if(isphi)
	{
		for(int i=0; i<size; i++)
			if(binmap[index][i] >= cutoff)
			{
				if(!count_bad)
					return 0;
			}
			else
				nbad ++;
	}
	else
	{
		for(int i=0; i<size; i++)
			if(binmap[i][index] >= cutoff)
			{
				if(!count_bad)
					return 0;
			}
			else
				nbad ++;
	}

	return nbad;
}

// Computes the dihedral angle around bond between at2 and at3
double get_dihedral(double *at1, double *at2, double *at3, double *at4)
{
	double dihedral;
	double s;
	double ct;
	double bond1[3];
	double bond2[3];
	double bond3[3];
	double normal1[3];
	double normal2[3];

	// Get the three bond vectors
	for (int i = 0; i < 3; i++)
	{
		bond1[i] = at1[i] - at2[i];
		bond2[i] = at3[i] - at2[i]; // Bond that contains current dihedral angle
		bond3[i] = at3[i] - at4[i];
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

	//	In case normalization fails
	if ( ct > 1.0 )
		ct = 1.0;
	else if ( ct < -1.0 )
		ct = -1.0;

	dihedral = acos( ct );

	s = bond2[0] * ( normal1[2] * normal2[1] - normal1[1] * normal2[2] ) +
			bond2[1] * ( normal1[0] * normal2[2] - normal1[2] * normal2[0]  ) +
			bond2[2] * ( normal1[1] * normal2[0] - normal1[0] * normal2[1] );

	if ( s < 0.0 )
		dihedral = -dihedral;

	dihedral = ( dihedral > 0.0 ) ? M_PI - dihedral : -( M_PI + dihedral );

	//	Correction algorithm to get value in [-Pi,+Pi] range --- BUT IS THIS REALLY NEEDED?
	if (dihedral > M_PI)
		dihedral -= 2*M_PI;
	else if (dihedral < -M_PI)
		dihedral += 2*M_PI;

	return dihedral;
}

//int angle2index(double angle, int size, int site)
//// Avoid quantization problems when "angle" is -180 or +180 (-M_PI or +M_PI)
//// 	angle = ]-M_PI,+M_PI[
//// 	returned index is bounded to [0,size-1]
//{
//	// double angle0 = angle; // Print luego
//	if (angle  >=  M_PI)
//		angle += -2*M_PI + 1.0E-6; // Avoid (angle = -M_PI)
//	else
//		if (angle <= -M_PI)
//			angle += 2*M_PI - 1.0E-6; // Avoid (angle = +M_PI)
//
//	int index = angle * size/(2*M_PI) + size/2; // Convert from continuous angle to discrete bin
//
////	// Remove this later
////	if(index >= size || index < 0)
////	{
////		fprintf(stderr,"Hei, angle2index() crashed! angle=%lf  index= %d  site=%d\n",angle,index,site);
////		exit(0);
////	}
//
//	return index; // Index is bounded to [0,size-1]
//}

double index2angle(int index, int size)
// returned angle is at the middle of the respective index bin ( ((step/2)*(M_PI/180)) )
{
	double ang = ((index - size/2) * (2*M_PI) / size ) + (((360.0/size)/2)*(M_PI/180.0)); // Convert from discrete bin to continuous angle
	return ang;
}

double index2angle_beta(int index, int size)
// returned angle is between the current index and the next one
{
	float random = (float) rg->Random();
	if(random == 5.0)
	{
		random -= 0.00001;
	}
	else if (random == 0.0)
	{
		random += 0.00001;
	}
	double ang = ((index - size/2) * (2*M_PI) / size ) + ((360.0/size)*random*(M_PI/180.0)); // Convert from discrete bin to continuous angle
	return ang;
}

double getoptphi(int indexphi, int indexpsi, int aux_index, float ***binmaps, double *dihedral_angle, int ind, int ibond, int size, int *opt_index, float cutoff)
// Function that returns the optimal value of phi angle
{
	double opt_angle;
	double angle = 0.0;

	if(indexphi < aux_index)
	{
		for(int i = indexphi+1; i <= aux_index; i++)
		{
			if(binmaps[ind][i][indexpsi] > cutoff)
			{
				*opt_index = i;
				opt_angle = index2angle(i,size);
				angle = opt_angle - dihedral_angle[ibond+1];
				break;
			}
		}
	}
	else if (indexphi > aux_index)
	{
		for(int i = indexphi-1; i >= aux_index; i--)
		{
			if(binmaps[ind][i][indexpsi] > cutoff)
			{
				*opt_index = i;
				opt_angle = index2angle(i,size);
				angle= opt_angle - dihedral_angle[ibond+1];
				break;
			}
		}
	}
	else
		*opt_index = aux_index;

	return angle;
}

double getoptpsi(int indexphi, int indexpsi, int aux_index, float ***binmaps, double *dihedral_angle, int ind, int ibond, int size, int *opt_index, float cutoff)
// Function that returns the optimal value of psi angle
{
	double opt_angle;
	double angle = 0.0;

	if(indexpsi < aux_index)
	{
		for(int i = indexpsi+1; i <= aux_index; i++)
		{
			if(binmaps[ind][indexphi][i] > cutoff)
			{
				*opt_index = i;
				opt_angle = index2angle(i,size);
				angle = opt_angle - dihedral_angle[ibond+1];
				break;
			}
		}
	}
	else if (indexpsi > aux_index)
	{
		for(int i = indexpsi-1; i >= aux_index; i--)
		{
			if(binmaps[ind][indexphi][i] > cutoff)
			{
				*opt_index = i;
				opt_angle = index2angle(i,size);
				angle = opt_angle - dihedral_angle[ibond+1];
				break;
			}
		}
	}
	else
		*opt_index = aux_index;
	return angle;
}

void show_points_angle(float anglephi, float anglepsi, int size, char *file)
{
	FILE *f_file;
	if( !(f_file = fopen(file,"w")) )
	{
		printf("Reading input error! file %s\n",f_file);
		exit(1);
	}

	int indexphi,indexpsi;
	indexphi = angle2index(anglephi,size);
	indexpsi = angle2index(anglepsi,size);
	//	fprintf(f_file, "#   x     y        point value\n");

	for(int i=0;i<size;i++)
	{
		for(int j=0;j<size;j++)
		{
			if((i == indexphi) && (j == indexpsi))
			{
				fprintf(f_file, "%5d %5d %12e \n", i, j, 1.0);
			}
			else
			{
				fprintf(f_file, "%5d %5d %12e \n", i, j, 0.0);
			}
		}
	}
	fclose(f_file);
}

void accu_1D(float **pdf, int size, double angle, bool isphi, float **p_acc)
{
	if(*p_acc == NULL)
		*p_acc = (float *) malloc( sizeof(float) * size); // 1D acc

	double sum = 0.0;
	int index = angle2index(angle,size);

	for( int i=0; i<size; i++) // rows or cols
	{
		if(isphi) // We fix the phi (Nt anchor) angle.
		{
			sum += pdf[index][i];
			(*p_acc)[i] = sum;
		}
		else // We fix the psi (Ct anchor) angle.
		{
			sum += pdf[i][index];
			(*p_acc)[i] = sum;
		}
	}

	// Normalize
	for (int i=0; i<size; i++)
		(*p_acc)[i] /= sum;
}

// Input "index" instead of "angle" (this one does not allocate memory)
// isphi = true, if the fixed angle is Phi
void accu_1D(float **pdf, int size, int index, bool isphi, float *acc)
{
	double sum = 0.0;

	for( int i=0; i<size; i++) // rows or cols
	{
		if(isphi) // We fix the phi (Nt anchor) angle.
		{
			sum += pdf[index][i];
			acc[i] = sum;
		}
		else // We fix the psi (Ct anchor) angle.
		{
			sum += pdf[i][index];
			acc[i] = sum;
		}
	}

	// Normalize
	for (int i=0; i<size; i++)
		acc[i] /= sum;
}


float rang_index1D_dun(float *pdf_accum, int size) // aleator --> rg
{
	float random = (float) rg->Random(); // Output random float number in the interval 0 <= x < 1
	int i=0;

	while(i<size && pdf_accum[i] < random)
		i++;

	return( index2angle(i,size) );
}

// Quick screen dump of dihedral angles
void dump_dihedrals(char *name, double *dihedral_angle, int nr_atoms_loop)
{
	// Showing dihedral_angle at the end of "valid loop"
	fprintf(stderr,"%s:\n",name);
	for(int i=0; i<nr_atoms_loop; i++)
		fprintf(stderr," %6.3f", dihedral_angle[i]);
	fprintf(stderr,"\n");

}

// Compute RMSD between two conformations of "nr_atoms_loop" atoms and "cg_mode" coarse graining
double rmsd_co_p1co(double **co1,double **p1co1, double **co2,double **p1co2, int nr_atoms_loop, int cg_mode)
{
	double simil = 0.0;
	int natoms = nr_atoms_loop;
	for(int k=0; k<nr_atoms_loop; k++)
		simil += pow(co1[k][0] - co2[k][0] ,2) + pow(co1[k][1] - co2[k][1] ,2) + pow(co1[k][2] - co2[k][2] ,2);
	if(cg_mode>0)
	{
		if(cg_mode !=2)
		{
			natoms += nr_atoms_loop/3;
			for(int k=2; k<nr_atoms_loop; k+=3) // Screen "O"s
				simil += pow(p1co1[k][0] - p1co2[k][0] ,2) + pow(p1co1[k][1] - p1co2[k][1] ,2) + pow(p1co1[k][2] - p1co2[k][2] ,2);
		}

		if(cg_mode !=1)
		{
			natoms += nr_atoms_loop/3;
			for(int k=1; k<nr_atoms_loop; k+=3) // Screen "CB"s
				simil += pow(p1co1[k][0] - p1co2[k][0] ,2) + pow(p1co1[k][1] - p1co2[k][1] ,2) + pow(p1co1[k][2] - p1co2[k][2] ,2);
		}
	}

	return sqrt(simil/natoms); // Square root is slow!
}

double bond_length2( double * a, double * b )
{
	return sqrt( pow( a[0] - b[0], 2 ) + pow( a[1] - b[1], 2 ) + pow( a[2] - b[2], 2 ) );
}


// Compute RMSD (N,CA,C,O) from two "co" and "p1co" arrays
double total_rmsd_NCACO(double** co_i,double** co,double** pco_i,double** pco, int number_residues)
{
	double sum;						/*sub sum floating point number																			*/
	int index;

	sum=0;
	// RMSD is computed using "moving" residues only (not anchors)
	for (int j = 1; j < number_residues+1; j++)	/*iterate # residues*/ /*local j used and override global j*/
	{
		index = (j*3);

		for(int i = 0; i < 3; i++)
		{
			sum += pow( co[index][i] - co_i[index][i], 2 ); // N
			sum += pow( co[index+1][i] - co_i[index+1][i], 2 ); // CA
			sum += pow( co[index+2][i] - co_i[index+2][i], 2 ); // C
			sum += pow( pco[index+2][i] - pco_i[index+2][i], 2 ); // O
			//  sum += pow( pco[index+1][i] - pco_i[index+1][i], 2 ); // CB
		}
	}

	return sqrt((sum/(number_residues*4.0)));
}

// Compute RMSD (N,CA,C) from two "co" arrays
double total_rmsd_NCAC(double** co_i,double** co, int number_residues)
{
	double sum;						/*sub sum floating point number																			*/
	int index;

	sum=0;
	// RMSD is computed using "moving" residues only (not anchors)
	for (int j = 1; j < number_residues+1; j++)	/*iterate # residues*/ /*local j used and override global j*/
	{
		index = (j*3);

		for(int i = 0; i < 3; i++)
		{
			sum += pow( co[index][i] - co_i[index][i], 2 ); // N
			sum += pow( co[index+1][i] - co_i[index+1][i], 2 ); // CA
			sum += pow( co[index+2][i] - co_i[index+2][i], 2 ); // C
		}
	}

	return sqrt((sum/(number_residues*3.0)));
}

// Quick screen dump of dihedral angles
void dump_IC(double *bondlength, double *valence,double *dihedral_angle, int nr_atoms_loop, char **secu,float FirstPhi,float LastPsi)
{

	// Showing all angles at the end of "valid loop"
	float c;
	c=180/M_PI;
	dihedral_angle[2]=FirstPhi;
	for(int i=0; i<nr_atoms_loop/3; i++)
	{
		fprintf(stdout,"pdbloop> %4s ", secu[i]);
		if (i==(nr_atoms_loop/3-1))
			fprintf(stdout,"%9.3f %9.3f %9.3f ", dihedral_angle[i*3+2]*c,LastPsi*c,0.0);
		else
			fprintf(stdout,"%9.3f %9.3f %9.3f ", dihedral_angle[i*3+2]*c,dihedral_angle[i*3+3]*c,dihedral_angle[i*3+4]*c);
		fprintf(stdout,"%7.3f %7.3f %7.3f ", valence[i*3]*c, valence[i*3+1]*c,valence[i*3+2]*c);
		if (i==(nr_atoms_loop/3-1))
			fprintf(stdout,"%4.3f %4.3f %4.3f  ", bondlength[i*3+1], bondlength[i*3+2], 0.0);
		else
			fprintf(stdout,"%4.3f %4.3f %4.3f ", bondlength[i*3+1], bondlength[i*3+2],bondlength[i*3+3]);
		fprintf(stdout,"\n");
	}

	// Standard
	fprintf(stdout,"pdbloop>   X                                ");
	fprintf(stdout,"%7.3f %7.3f %7.3f ", valCACON*c, valCONCA*c, valNCACO*c  );
	fprintf(stdout,"%4.3f %4.3f %4.3f ",  bondNCA, bondCACO, bondCON);
	fprintf(stdout,"\n");
	fprintf(stdout,"pdbloop>                                    ");
	fprintf(stdout,"%7.3f %7.3f %7.3f ", sigma_valCACON*c, sigma_valCONCA*c, sigma_valNCACO*c );
	fprintf(stdout,"%4.3f %4.3f %4.3f ",  sigma_bondNCA, sigma_bondCACO, sigma_bondCON);
	fprintf(stdout,"\n");
	// GLY
	fprintf(stdout,"pdbloop>   G                                ");
	fprintf(stdout,"%7.3f %7.3f %7.3f ", GLYvalCACON*c, GLYvalCONCA*c, GLYvalNCACO*c  );
	fprintf(stdout,"%4.3f %4.3f %4.3f ", GLYbondNCA, GLYbondCACO, GLYbondCON);
	fprintf(stdout,"\n");
	fprintf(stdout,"pdbloop>                                    ");
	fprintf(stdout,"%7.3f %7.3f %7.3f ", sigma_GLYvalCACON*c, sigma_GLYvalCONCA*c, sigma_GLYvalNCACO*c );
	fprintf(stdout,"%4.3f %4.3f %4.3f ", sigma_GLYbondNCA, sigma_GLYbondCACO, sigma_GLYbondCON);
	fprintf(stdout,"\n");
	// PRO
	fprintf(stdout,"pdbloop>   P                                ");
	fprintf(stdout,"%7.3f %7.3f %7.3f ", PROvalCACON*c, PROvalCONCA*c, PROvalNCACO*c  );
	fprintf(stdout,"%4.3f %4.3f %4.3f ", PRObondNCA, PRObondCACO, PRObondCON);
	fprintf(stdout,"\n");
	fprintf(stdout,"pdbloop>                                    ");
	fprintf(stdout,"%7.3f %7.3f %7.3f ", sigma_PROvalCACON*c, sigma_PROvalCONCA*c, sigma_PROvalNCACO*c );
	fprintf(stdout,"%4.3f %4.3f %4.3f ", sigma_PRObondNCA, sigma_PRObondCACO, sigma_PRObondCON);
	fprintf(stdout,"\n");
}

// Save into text file all dihedral angles (Phi, Psi and Omega, in degrees)
void save_dihedrals(double *dihedral_angle, int nr_atoms_loop, char *file, char *seq)
{
	// Showing all angles at the end of "valid loop"
	float c = 180/M_PI;
	FILE *f_file;
	if( !(f_file = fopen(file,"w")) )
	{
		fprintf(stderr,"Unable to write file %s\n Forcing exit!\n", file);
		exit(1);
	}

	if(seq == NULL)
	{
		for(int i=0; i < (nr_atoms_loop/3)-1; i++)
			fprintf(f_file,"%9.3f %9.3f %9.3f\n", dihedral_angle[i*3+2]*c,dihedral_angle[i*3+3]*c,dihedral_angle[i*3+4]*c);
		fprintf(f_file,"%9.3f %9.3f %9.3f\n", dihedral_angle[nr_atoms_loop-1]*c, dihedral_angle[0]*c, 0.0);
	}
	else
	{
		for(int i=0; i < (nr_atoms_loop/3)-1; i++)
			fprintf(f_file,"%c: %9.3f %9.3f %9.3f\n", seq[i], dihedral_angle[i*3+2]*c,dihedral_angle[i*3+3]*c,dihedral_angle[i*3+4]*c);
		fprintf(f_file,"%c: %9.3f %9.3f %9.3f\n", seq[(nr_atoms_loop/3)-1], dihedral_angle[nr_atoms_loop-1]*c, dihedral_angle[0]*c, 0.0);
	}

	fclose(f_file);
}

//// Ramachandran energy derived from neighbor-dependent PDFs from Dunbrack's paper.
//float rama_energy(double *dihedral_angle,int nr_atoms_loop,float ***maps,int size)
//{
//	float c = 180/M_PI;
//	float energy = 0.0;
//
//	for(int i=0; i < (nr_atoms_loop/3)-1; i++) // screen residues
//	{
////		fprintf(stderr,"Phi%d= %6.2f  Psi%d= %6.2f  p= %f\n",i,dihedral_angle[i*3+2]*c,i,dihedral_angle[i*3+3]*c,
////				maps[ i ][ angle2index(dihedral_angle[i*3+2], size) ][ angle2index(dihedral_angle[i*3+3], size) ] );
//		energy -= log10(maps[ i ][ angle2index(dihedral_angle[i*3+2], size) ][ angle2index(dihedral_angle[i*3+3], size) ]);
//	}
////	fprintf(stderr,"Phi%d= %6.2f, Psi%d= %6.2f  p= %f\n",(nr_atoms_loop/3)-1, dihedral_angle[nr_atoms_loop-1]*c, (nr_atoms_loop/3)-1, dihedral_angle[0]*c,
////			maps[ (nr_atoms_loop/3)-1 ][ angle2index(dihedral_angle[nr_atoms_loop-1], size) ][ angle2index(dihedral_angle[0], size) ] );
//	energy -= log10(maps[ (nr_atoms_loop/3)-1 ][ angle2index(dihedral_angle[nr_atoms_loop-1], size) ][ angle2index(dihedral_angle[0], size) ]);
//
//	return energy;
//}

// Get the most likely index of Phi, if bestphi=true, (or Psi if bestphi=false) for a given Psi (or Phi) "index" from a PDF distribution "map" of side "size".
int best_index(float **map, bool bestphi, int index, int size)
{
	// FILE *myf = fopen("profile.txt","w");

	float best = 0.0;
	int ibest = -1;
	if(bestphi) // fix Psi, search for best Phi
		for(int i=0; i<size; i++)
		{
			// fprintf(myf,"%d %f %f\n",i,-177.5+5*i,map[i][index]);
			if(map[i][index] > best)
			{
				best = map[i][index];
				ibest = i;
			}
		}
	else // fix Phi, search for best Psi
		for(int i=0; i<size; i++)
		{
			// fprintf(myf,"%d %f %f\n",i,-177.5+5*i,map[index][i]);
			if(map[index][i] > best)
			{
				best = map[index][i];
				ibest = i;
			}
		}
	// fprintf(myf,"# ibest= %d\n",ibest);
	// fclose(myf);

	if(ibest<0)
	{
		fprintf(stderr,"Error, plain Ramachandran region for a given Phi or Psi! Please, check your thresholds...\n");
		exit(1);
	}

	return ibest;
}

void quicksort(double **arr, int low, int high)
{
	double *y; // dummy variable for swap
	int i = low; // left
	int j = high; // right
	double z = *arr[(low + high) / 2]; // pivot

	do {
		while(*arr[i] < z) // find member above
			i++;

		while(*arr[j] > z) // // find member below
			j--;

		if(i <= j) // Swap two elements
		{
			swap(arr[i],arr[j],y);
			i++;
			j--;
		}
	} while(i <= j);

	// recurse
	if(low < j)
		quicksort(arr, low, j);

	if(i < high)
		quicksort(arr, i, high);
}


// Quick-sort recursive routine to everything in between `low' <-> `high'
void quicksort(double *arr, int low, int high, int *index)
{
	double y; // dummy variable for swap
	int i = low; // left
	int j = high; // right
	double z = arr[(low + high) / 2]; // pivot




	do {
		while(arr[i] < z) // find member above
			i++;

		while(arr[j] > z) // // find member below
			j--;

		if(i <= j) // Swap two elements
		{
			swap(arr[i],arr[j],y);
			swap(index[i],index[j],y);
			i++;
			j--;
		}
	} while(i <= j);

	// recurse
	if(low < j)
		quicksort(arr, low, j, index);

	if(i < high)
		quicksort(arr, i, high, index);
}


