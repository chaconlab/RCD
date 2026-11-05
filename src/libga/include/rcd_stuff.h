/*
 * rcd_stuff.h
 *
 *  Created on: Jun 16, 2016
 *      Author: mon
 */

#ifndef INCLUDE_RCD_STUFF_H_
#define INCLUDE_RCD_STUFF_H_


#include <stdlib.h>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <cstring>
#include <vector>
#include <cmath>
#include <ctime>
#include <sys/time.h>	//:TRICKY:  rare header
#include <sys/stat.h>
#include <sstream>	//:TRICKY:  rare header  	needed for int-string conversion
#include <string.h>

// #include "parameters.h"
// #include "strum.h"
// #include "loop.h"

//#include "libga/strum.h"
//#include "libga/loop.h"

#include <libvolume/include/floatops.h>
#include <libpdb/include/Macromolecule.h>
#include <libpdb/include/pdbIter.h>
#include "libtools/include/Mersenne.h"
#include "libenergy/include/pd2.h"
#include <cmdl/CmdLine.h>

//#include "rasps.h"
//#include "libga/timer.h"	     /*no more linked to library Nacio*/
#include "geometry.h"		 // IC's stuff
#include "ga.h"		       /*external library my_library*/
#include "own_math.h"
#include "chain_manip.h"
#include "statistics.h"
#include "stat.h"
#include "matrix.h"
#include "loopdata.h"
#include "libenergy/include/icosa.h" // ICOSA energy
#include "libenergy/include/korpe.h" // KORP energy
//#include "librasp/include/rasp.h" // RASP side chain repacking
#include "libenergy/include/rama.h" // Rama energy


//int MAXIT, MAX_ITER_SECANT;


//CRandomMersenne *rg; // Mersenne Twister global object
//float unitsx;
//float unitsy;
//float unitsz;
//int limitx;
//int limity;
//int limitz;
//int stepy;
//int stepz;

extern CRandomMersenne *rg; // Mersenne Twister global object
extern float unitsx;
extern float unitsy;
extern float unitsz;
extern int limitx;
extern int limity;
extern int limitz;
extern int stepy;
extern int stepz;

#define NAAS 21 // Number of aminoacids (20 + Cis-PRO) // TEST
#define MAX_PDBS 200 // Maximum number of PDBs in input file (Mon: we should use dynamic memory allocation...)

#define rad2deg 57.29577951308232087721		//CONVERT RAD TO DEGREE


#define ISNEAR(I,J,C) {                                    \
							d = pow(I[0]-J[0],2);          \
							if( d < C )                    \
							{                              \
								d += pow(I[1]-J[1],2);     \
								if( d < C )                \
								{                          \
									d += pow(I[2]-J[2],2); \
									if( d < C )            \
										return 1;          \
								}                          \
							}                              \
						}
// Required by "quicksort"
#define swap(a,b,t) ((t)=(a),(a)=(b),(b)=(t))

//using namespace std;


float LoopClosure( double **co, double** anchor_fixed, double *bond_length, double *valence_angle, double * dihedral_angles, int nr_atoms_loop, int reverse);
double bond_length2( double * a, double * b );

int	bondselect_marker(int nratoms, int *marker,int *previous, int sense);
int	bondselect_marker2(int nr_atoms_loop, int *marker,int *previous, int sense);
void randomizeIC(double *bond_length,double *valence_angle,double *bond_length2,double *valence_angle2,int nratoms,double lengthptb,double valenceptb,char *seq, char seq_ct1, bool flat=false);
float normRand(bool flat=false); // flat= true --> flat random distribution in [-1:1] range.
void loop_co( Macromolecule *loop,  double** co, double** p1co,int *residuemarker, int *seq, bool rmsd_O=false);
// Create residuemarker and seqnum from Macromolecule
void residuemarker_seq( Macromolecule *loop, int *residuemarker, int *seqnum);
int clash(vlVolume *vlm, double** co,  double **p1co, int nr_atoms_loop, int sense, bool check_p1co=true);
// Check clashes between loop and pocket-map
//  mode  --> mode of operation: 0=(N,CA,C), 1=(N,CA,C,O), 2=(N,CA,C,CB), 3=(N,CA,C,O,CB)
//  count --> set true to enable counting all clashes, otherwise (false) it will return upon first clash detection
int clash_cg(vlVolume *vlm, double **co, double **p1co,  int nr_atoms_loop, int *marker2, int reverse, int mode, bool count = false);
int clash_cg_ibond(vlVolume *vlm, double **co, double **p1co, int bond_number, int nr_atoms_loop, int *marker2, int reverse, int mode, bool count = false);


// Clash test of backbone without anchors.
int clash_co(vlVolume *vlm, double** co, int nr_atoms_loop);
int intra_clash_cg(double** co,double **p1co, double fact_co, double fact_ocb, int nr_atoms_loop,int atomshift, int *marker2, int reverse, int cg_mode);
int intra_clash_cg_ibond(double** co,double **p1co, double fact_co, double fact_ocb, int bond_number, int nr_atoms_loop,int atomshift, int *marker2, int reverse, int cg_mode);

/**
* Project a "co" minimal backbone representation (just N, CA, and C atoms) into a pre-existent map.
* (Mon made 27/8/2018)
* vol  = Input/Output volume
* res  = Resolution [Angstroms]
* nat  = Number of loop atoms (nr_atoms_loop)
* co   = N, CA, C, coordinates array
* mask = Mask of indices (array of integers) for fast cross-correlation. WARNING: A negative integer must indicate the array end!
*/
void co2map_real( vlVolume *vol, float res, int nat, double **co);

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
void co2map_real( vlVolume *vol, float res, int nat, double **co, double **pco=NULL, int cg=0, int *residuemarker=NULL);

// Required by "co2map_real"
int roundFloat_( float x );



void multiPDB(double **co,double **p1co,double **p2co,int number_loops,int number_residues,char *name, char **secu, int start);
void multiPDB_cg(double **co,double **p1co,int number_loops,int number_residues, char *name, char **secu, int start, int mode=0, char chain='A');
void multiPDB_write(double **co,int number_loops,int number_residues_loop, char *name);
// Compute RMSD (N,CA,C) from two "co" arrays
double total_rmsd_NCAC(double** co_i,double** co, int number_residues);
// Compute RMSD (N,CA,C,O) from two "co" and "p1co" arrays
double total_rmsd_NCACO(double** co_i,double** co,double** pco_i,double** pco, int number_residues);
//************************************************************
// Alex: Functions definitions
float *accunorm(float **pdf, int size); // Ramachandran neighbor dependent functions

// Now the following functions are in rama.cpp
// float *****read_dunbrack(char *file, int *p_size, float *p_step);
// float **MatrixProduct(float **M1, float **M2, int size);
// float **MatrixDivision(float **M1, float **M2, int size);
// float **map_gen(int C, int L, int R, float *****pdfs, int size);
// void norm_map(float **map, int size);
// void show_map(float **map,int size, char *file);
// Avoid quantization problems when "angle" is -180 or +180 (-M_PI or +M_PI)
// 	angle = ]-2*M_PI,+2*M_PI[
// 	returned index is bounded to [0,size-1]
// int angle2index(double angle, int size);

void rang_index2D_dun(float *pdf_accum, int size, int step, float *phi, float *psi);
float **map_bin(float **map, float thr, int size);
float **map_bin2(float **map, float thr, int size); // MON's normalization
float **map_bin3(float **map, float thr, int size); // PABLO's normalization
// Computes the cutoff value that accounts for the selected probability using a relative threshold ( thr=]0,1[ )
// Warning, map values are sorted first (not using quick-sort algorithm, thus not very efficient...)
// (it allocates output map memory)
float map_cutoff(float **map, float thr, int size);
// Set to zero those elements with values higher or lower than some threshold. (not allocates memory)
// less = true --> lower values than "thr" will be set to zero, otherwise higher values than "thr" will be set to zero
void map_setzero(float **map, int size, float thr, bool less);
// Copy a 2D map of "size" side elements. (it allocates output map memory)
float **map_copy(float **map, int size);
float **map_copyh3(float ***map, int k, int size);
float **map_mergeh3(float ***map, float **map2, int k, int size, float factor);

// Computes the dihedral angle around bond between at2 and at3
double get_dihedral(double *at1, double *at2, double *at3, double *at4);
// CHECK (and count) ALL DIHEDRAL ANGLES
// count_bad --> true= return the total number of aminoacids with bad dihedrals, false= return 1 as soon as any aminoacid has bad dihedrals (default)
int check_angles(float ***binmap, double *dihedral, int size, int nr_atoms_loop, int reverse, bool count_bad=false, int *marker=NULL, float *cutoffs=NULL);
// CHECK JUST A GIVEN DIHEDRAL ANGLES PAIR (returns 1 if bad angles or 0 if OK)
int check_angles_final(float **binmap, double ang_phi, double ang_psi, int size, int nr_atoms_loop);
// Quick screen dump of dihedral angles
void dump_dihedrals(char *name, double *dihedral_angle, int nr_atoms_loop);
// Quick screen dump of dihedral angles
void dump_IC(double *bondlength, double *valence,double *dihedral_angle, int nr_atoms_loop, char **secu,float FirstPhi,float LastPsi);
// Save into text file all dihedral angles (Phi, Psi and Omega, in degrees)
void save_dihedrals(double *dihedral_angle, int nr_atoms_loop, char *file, char *seq=NULL);


double index2angle(int index, int size);
double getoptphi(int indexphi, int indexpsi, int aux_index, float ***binmaps, double *dihedral_angle, int ind, int ibond, int size, int *opt_index, float cutoff=0.5);
double getoptpsi(int indexphi, int indexpsi, int aux_index, float ***binmaps, double *dihedral_angle, int ind, int ibond, int size, int *opt_index, float cutoff=0.5);
// Get the most likely index of Phi, if bestphi=true, (or Psi if bestphi=false) for a given Psi (or Phi) "index" from a PDF distribution "map" of side "size".
int best_index(float **map, bool bestphi, int index, int size);
// Ramachandran energy derived from neighbor-dependant PDFs from Dunbrack's paper.
//float rama_energy(double *dihedral_angle,int nr_atoms_loop,float ***maps,int size);
// Plots a 2D map and (optionally) phi and psi values.
void paint_map_point(float **map,int size, char *file, int indexphi = -1, int indexpsi = -1);
//float *accu_1D(float **pdf, int size, double angle, bool isphi);
void accu_1D(float **pdf, int size, double angle, bool isphi, float **p_acc);
// Input "index" instead of "angle" (this one does not allocate memory) (isphi = true, if the fixed angle is Phi)
void accu_1D(float **pdf, int size, int index, bool isphi, float *acc);
float rang_index1D_dun(float *pdf_accum, int size); // aleator --> rg
// Plot the Ramachandran maps for all aminoacids with all dihedral angles. G
// INPUT:  map --> Ramachandran map
//         size --> Map size
//         nres --> Number of residues including anchors (nr_loop+2)
//         base --> Output files base name
//         dihedral_angle --> Dihedral angles array
//         reverse --> 0= forward, 1= backward
// OUTPUT: One Ramachandran map per aminoacid with current dihedral angles represented therein.
void paint_maps_points(float ***maps, int size, int nres, char *base, double *dihedral_angle = NULL, int reverse = 0 );
// Plot the Ramachandran maps for all aminoacids with all dihedral angles. G
// INPUT:  map --> Ramachandran map
//         size --> Map size
//         dihedral_angle --> Dihedral angles array of the "nsol" solutions (nsol*nres*3)
//         nres --> Number of residues including anchors (nr_loop+2)
//         reverse --> 0= forward, 1= backward
//         base --> Output files base name
//         nsol --> Number of solutions
// OUTPUT: One Ramachandran map per aminoacid with current dihedral angles represented therein.
void paint_all_dihedrals(double *dihedral_angle, int size, int nres, int reverse, char *base, int nsol);
void show_points_angle(float anglephi, float anglepsi, int size, char *file);
// Given Phi (or Psi), it tells whether there is any valid Psi (or Phi) "in Ramachandran".
// angle --> Given Phi or Psi (fixed)
// isphi --> =true if "angle" is Phi, =false otherwise.
// count_bad --> =true, counts the number of valid Ramachandran map positions (default = false)
int check1angle(float **binmap, double angle, bool isphi, int size, int reverse, bool count_bad=false, float cutoff=0.5);
// Compute RMSD between two conformations of "nr_atoms_loop" atoms and "cg_mode" coarse graining
double rmsd_co_p1co(double **co1,double **p1co1, double **co2,double **p1co2, int nr_atoms_loop, int cg_mode=0);
// Quick-sort recursive routine to sort everything in between `low' <-> `high'
void quicksort(double *arr, int low, int high, int *index);
void quicksort(double **arr, int low, int high);


#endif /* INCLUDE_RCD_STUFF_H_ */
