/*
 * parameters.h
 *
 *  Created on: Oct 26, 2013
 *      Author: pablo
 */
//
//  parameters.h
//
//
//  Created by Pablo on 11/9/12.
//
//

#ifndef _parameters_h
#define _parameters_h


  /*
    ------------------------------------------------

    PROGRAM VARIABLES

    ------------------------------------------------
  */


  //:INFO:	Different types of counter
  int save_counter;	/*Counter to keep track of saving the exact loop closure solutions*/
  int k;		/*Solutions counter loop closure procedure*/
  int totalk;		/*Total counter loop counter loop closure*/
  int self_marker;	/*marker to indicate that a self-loop clash has occurred; if equal to 1 = yes, if equal to no*/
  int atom_marker;	/*mark which atom a clash occurred*/
  int counter_perturbed;/*to count how many perturebed loops must be rejected*/
  int counter_optimized;/*to count how many optimized loops must be rejected*/
  int counter_passed;   /*to count how many optimized loops have really passed*/



  int marker_input;		/*counter number input arguments*/
  int switch_stop;		/*switch to stop at several places during program execution for debugging/checking*/ /*NOT YET FULLY DEBUGGED*/
  int switch_verbose;		/*Switch on printing debugging information to screen during program execution*/
  //:TODO:	Implement a different level for end users	-  Keep at minimum - Mind you in iterations
  //			if only once in a while (each 100 or more)

  int switch_path;		/*Switch to save one loop closure trajectory (1 = on, 0 = off); similar to molecular dynamics	*/
  int switch_file;		/*Switch on saving of output files with results (1 = on, 0 = off)*/
  //:TODO:	work a bit on these options
  int switch_situs;		/*switch on saving of protein grid in situs file (1 = on, 0 = off)*/
  int switch_solutions;	/*switch on to save all found loop closure solutions (1 = on, 0 = off)*/
  int grid_backbone;		/*switch to only select backbone atoms for filling protein grid*/
  int seq_initgrid;		/*switch on procedure to sequentially put loop clashfree on the grid*/
  int switch_dof;		 //	Switch on extra degree of freedom to be used (1 = not used, 0 = extra DoF)
  int switch_bond;	 //	Switch which bond selection algorithm is used
  //	= 1 : random bond choice minus previous: BONDSELECT6.CPP
  //	= 0 : sequential cycling through range:  BONDSELECT7.CPP
  int switch_perturbation; //	Switch to introduce perturbation step at start
  int no_loop;		 /*switch to let program know that no residue coordinates are present and that only anchor positions are used*/
  int switching_side;      /*to switch direction of protocol in 50%-50% manner for all loop closures */
  int switching_side2;	 /*to switch direction of protocol when failure occurs of loop closure */
  int switching_mode;

  //	Parameters
  // ------------------------------------------------------------------------------------------
  int randomize;		/*randomization for each loop in forward mode with native loop, only working under: REVERSE = 0, NO_LOOP = 0*/
  int random_omega;
				//:TODO:	only working under the above condition at present
  int randbond;	        /*switch to indicate that bond lengths are randomized according to literature values (including standard deviation sigma's)*/
  int randval;		/*switch to indicate that valence angles are randomized according to literature values (including standarr deviation sigma's)*/
  int inbetween;		/*integer stating whether the atom coordinates in between the two loop anchors are present or not
				  in the only-anchor runs; in principle this is so and thus the state is 1 (present)				  */
  int ramachandran_check;		/*randomize start loop in Ramachandran valid range*/
  int ramachandran_inside;	/*Ramachandran zone protocol inside loop closure routine*/
 // int selfloop_check;			/*switch for checking collisions in loop itself*/
  int sidechain_check;		/*switch for side chain collision check with protein: == 1 is on, == 0 is off*/
  int mainchain_check;		/*switch for main chain collision check with protein: == 1 is on, == 0 is off*/
  int post_check;				/*switch for main chain collision check with protein after loop but inside closure: == 1 is on, == 0 is off*/
  int include_hydrogens;		/*switch to select hydrogens in loop environment or not*/

  double degrees_perturbation;/*degrees that angle is perturbed*/
  double angle;				/*Optimized dihedral angle (rad)*/
  int randcount;				/*number of times (multiplied by number of rotatable bonds) native loop is randomized for each closure*/
  int max_path;				// Maximum number of completed RMSD iterations in algorithm that is allowed
  float path_slope; // Slope of the maximum number of RCD iterations
  int cycles;					/*Number of loop solutions that must be computed with RMSD < RMSD_CRIT
						  Default: 1000																	*/
  int number_path; 			/*Number of the index for which a loop closure trajectory is saved; is function of max_path*/
  int	profile_counter;		//	counter for saving rmsd-profile -> fully debugged -> I think so: FOR MULTIPLE PDB'S?
  int max_perturbation;		/*Maximum number allowed interations perturbation part of algorithm*/
  double rmsd_crit;   		/*Minimum RMSD value within which loop anchor must come to achieve actual loop closure*/
  double stepsize;			/*Distance between two grid points (correspons to a voxel size)*/
  double softness;  // Softness factor grid [0-1]: fraction pure hard sphere approach for atoms on grid
  // Value for pure hard sphere approach 	= 	1
  // Universal radius used for all atoms 	= 	1.5 Angstrom
  // A softness factor 0.8 is normally a lower limit for atoms; beyond this value energy
  // goes  asymptotically to infinity






  //CLASH DETECTION VARIABLES
  int infinite;			/*counter in escape protocol of Ramachandran-inside filter 			*/
  int check;				/*check for clashes 1 in post-processing; equal to 0 if no clash	*/
  int check2;	    		/*check for clashes 2 					; equal to 0 if no clash	*/
  int kindex;				/*running atom index 1 for clash detection  						*/
  int kindex2;			/*running atom index 2 for clash detection				    */
  int counter_clash ;		/*counter number of clashes in post-processing*/
  int cutoff;				/*number of atoms at the end of complete loop not used for clash detection:
					 *this is relevant for the end anchor itself since the position of this upon real
					 *loop closure is known and fixed */
  int cutoffp;	  		/*post-filter cutoff in number of atoms from the anchor end; see cutoff             */
  int atomshift;			/*number of intervening atoms not screened for self-loop collision detection:
					 *normal values are skips of 3 and 4 atoms				   */
  double clash_factor;	/*Multiplication factor soft sphere approach in self-loop clash detection [0-1]
***************************************************************************************
*This factor is equivalent to the softness factor for the grid but applies to the
*self-loop	clash detection protocol; the meaning is identical
*Thus, full hard sphere approach 		  = 1.0	   (1.0 if defined by squares in Eq.)
*      soft sphere approach (Soto paper) = 0.50   (0.7071 )
* 		soft sphere default			      = 0.80   (0.894  )
****************************************************************************************/


  /*perturbation variables*/
  int native;	  // switch for using the native or default peptide coordinates
  double bondptb; //factor which x 95% confidence interval for perturbation from ideal peptide bond lengths
  double valptb; // factor which x 95% confidence interval for perturbation from ideal valence angles
  double dihptb; // factor which X 95% confidence interval for perturbation from ideal dihedral angles
  	         // forgoing factors allow to artificially increase the range of extreme values outside typical


int gridlimit;				/*Maximum allowed size of the cubic grid*/
int sizegrid;				/*number of grid points in one dimension on 3D grid*/
int numberlist;				/*number of atoms found inside the cubic grid*/
int presizelist;			/*dimension coordinate array: fixed size ... not dynamically allocated*/

int sizegrid2;				/*number of grid points in one dimension on 3D grid*/
int numberlist2;			/*number of atoms found inside the cubic grid*/
int presizelist2;			/*dimension coordinate array: fixed size ... not dynamically allocated*/








#endif
