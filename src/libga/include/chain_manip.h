#ifndef  GUARD_chain_manip_h
#define  GUARD_chain_manip_h

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "libvolume/include/vlvolume.h" // vlVolume is needed for the clash grids (independent of VOLUME_INCLUDED)
#include "libpdb/include/Macromolecule.h"
#include "libpdb/include/pdbIter.h"
#include "libpdb/include/Fragment.h"

#include "ga.h"
#include "statistics.h"
#include "own_io.h"
#include "chain_manip.h"

#define TINY 	 		1.0e-20

using namespace std;

// Import standard values for internal coordinate values but not for anchor extremities*/
void standic(double *bond_length,double *valence_angle,int nratoms,char *seq, char seq_ct1);

//void standic(double *bond_length,double *valence_angle,int nratoms,int* residuemarker);

//// Obsolete?
//void randic(double *bond_length,double *valence_angle,int nratoms,double lengthptb,double valenceptb,int *residuemarker,int randbond,int randval);
//void randic(double *bond_length,double *valence_angle,double *dihedral_angle_i,int nratoms,double lengthptb,double valenceptb,double *prearray1,double *prearray2,double *prearray3,double *prearray4);

int checkself(double** co,double **pco,int check, int index,int nr_atoms_loop,int atomshift,double clash_factor,int reverse);

//double rachprob(double dihedral, int phipsi, int marker);

void intersec(double *a,double *b,double c[3],double d[3]);
void orient(double v0[],double a[],double b[],double v[]);
double bestangle_ana(double **co, int bondnumber,double **co_i, int nr_atoms);

int	bondbreaksel(int start,int range,int *bond_previous,int *bond_preimage,int counter,int flag);

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
						);

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
						);


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
						);

// Set a constant CA-CB distance
void set_cacb_dist2(double **co, double **p1co, int natoms, double dist2, char *seq);

void isafrotatetreef_cg(double **co,				// loop coordinates
						double **p1co,				// Cb && O coordinates
						int number_atoms,		// atoms in loop
						double *bond,			// bond lengths loop
						double *dihedral,			// dihedral angles loop
						double *cosarrayd,		// preallocated -cosines
						double *sinarrayd,		// preallocated sines
						double *cosarrayd4,		// preallocated -cosines
						double *sinarrayd4,		// preallocated sines
						double *cosarray,		// preallocated -cosines
						double *sinarray,		// preallocated sines
						char *seq,				// sequence (one-letter code)
						double *full,
						int mode=0,		// mode of operation: 0=(N,CA,C), 1=(N,CA,C,O), 2=(N,CA,C,CB), 3=(N,CA,C,O,CB)
						double factor=-1			// CA-CB common scaling factor (set to <=0 for a default constant distance)
						);




void isafrotatetreef_cg_old(double **co,				// loop coordinates
						double **p1co,				// Cb && O coordinates
						int number_atoms,		// atoms in loop
						double *bond,			// bond lengths loop
						double *dihedral,			// dihedral angles loop
						double *cosarrayd,		// preallocated -cosines
						double *sinarrayd,		// preallocated sines
						double *cosarrayd4,		// preallocated -cosines
						double *sinarrayd4,		// preallocated sines
						double *cosarray,		// preallocated -cosines
						double *sinarray,		// preallocated sines
						char *seq,				// sequence (one-letter code)
						double *full,
						int mode=0,		// mode of operation: 0=(N,CA,C), 1=(N,CA,C,O), 2=(N,CA,C,CB), 3=(N,CA,C,O,CB)
						double factor=-1			// CA-CB common scaling factor (set to <=0 for a default constant distance)
						);


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
						);

void copyco(double **co1,double **co2, int size);
void copyco(double **co1,double **co2,int bond, int size);


void pose(double **co1,double **co2,double *quaternion);
//void    valencefind(double **,double *,int);

/*Internal coordinate related*/
double  dihedral(double **co);
void    findlength(double **, double **,double *,int);
void 	correctdihedral(double *);
void    correctdihedral(double);
void 	valencefind(double **,double *,int );
void 	findlength(double **,double *,int);
void 	finddihedral(double **,double *,int );
void 	dihedralfind(double **,double *,int);

/*rotation related functions*/
void 	rotatematrix(double*,double**);
void 	rotatematrix2(double**,double**);
void 	rotateframe(double**,int,double**);
void rotatetree(double** co,double **pco1,double **pco2, int bond_number, int number_atoms, double angle, double **part,double **part1,double **part2,double **matrix);
int rotatetree_check(vlVolume *vlm, double** co, int bond_number, double angle, int natoms);
// Mon: New rotate & check function (faster)
// WARNING: Coordinates CAN'T clash before "bond_number" rotable bond.
int rotatetree_check2(vlVolume *vlm, double **co, int bond_number, double angle, int natoms);

int rotatetree_check3(vlVolume *vlm, double **co, int bond_number, double angle, int natoms);
void rotatetree_co(double** co, int bond_number, int number_atoms, double angle);
void rotatetree_cg(double **co, double **p1co,int bond_number, int number_atoms, double angle, int reverse, int mode=0);
int rotatetree_cg_check(double **co, double **p1co,int bond_number, int number_atoms,
		double angle, int reverse, int mode, vlVolume *vlm,
		float posV[3], int stepV[3], float unitV[3], int limitV[3], int *marker2);


void 	rotate(double**, int, int, double,double **,double**);
void 	rotate_bw(double **,double *,double *, double *, int);
void    rotate_bw2(double **,double *,double *,double *, double*, double *, int);
//void 	rotate(double**,int,int,double);
void 	rotate_angles(double**, double*,int *, int,int, double *,double *,double *);
void 	rotate_angles2(double **, double *,int *, int,int, double *,double *,double *);
void 	rotate_angles3(double **, double *,int *, int,int, double *,double *,double *);
void 	rotate_angles4(double **, double *,int *, int,int, double *,double *,double *);
void 	rotate_anglesm(double **, double *,int *, int,int, double *,double *,double *);

void matrixrotate(double **co,
					const double *angles,
					double **posematrix,
					const int *bondnumbers,
					const int bondrange,
					const int number_atoms,
					const double *valence,
					const double *bond,
					double *dihedral,
					const double *cosarray,
					const double *sinarray,
					double **matrix1,
					double **matrix2,
					double **matrixtot);

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
						);

void ga_rotate_anglesm2tree(double **co,			/*loop coordinates*/
							double **p1co,			/*coordinates 1st parallel chain*/
							double **p2co,			/*coordinates 2nd parallel chain*/
						const double *angles,		/*angle variations from optimization*/
						double *posespinor,			/*pose spinor for anchor*/
						const int *bondnumbers,		/*physical bond numbers*/
						const int bondrange,		/*# of bonds rotating*/
						const int number_atoms,		/*# atoms in loop*/
						const double *valence,		/*valence angles loop*/
						const double *bond,			/*bond lengths loop*/
						double *dihedral,			/*dihedral angles loop*/
						double *p1dihedral,			/*dihedral angles parallel loop 1*/
						double *p2dihedral,			/*dihedral angles parallel loop 2*/
						int *marker,
						double *spinor1,			/*auxiliary spinor 1-preallocated*/
						double *spinor2,			/*auxiliary spinor 2-preallocated*/
						double *spinortot			/*total spinor-preallocated*/
						);

void spinorrotate(double **co,				/*loop coordinates*/
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
						);

void ga_rotate_anglesm4(double **co,				/*loop coordinates*/
						const double *angles,		/*angle variations from optimization*/
						const int *bondnumbers,		/*physical bond numbers*/
						const int bondrange,		/*# of bonds rotating*/
						const int number_atoms,		/*# atoms in loop*/
						double *spinor1,			/*auxiliary spinor 1-preallocated*/
						double *spinor2,			/*auxiliary spinor 2-preallocated*/
						double *spinortot			/*total spinor-preallocated*/
						);

void garotate(double **co,				/*loop coordinates*/
						double **bond,
						double **bondn,
						const double *length,
						const double *angles,		/*angle variations from optimization*/
						const int *marker,
						const int number_atoms,		/*# atoms in loop*/
						double *between,
						double *spinor1,			/*auxiliary spinor 1-preallocated*/
						double *spinor2,			/*auxiliary spinor 2-preallocated*/
						double *spinortot);			/*total spinor-preallocated*/


/*dihedral optimization related*/
void posespinor(	double **co,
					int number_atoms,
					double *spinortot,
					double *spinor1,			/*auxiliary spinor 1-preallocated*/
					double *spinor2);			/*auxiliary spinor 2-preallocated*/
void 	startspinor(double **co, double *posespinor);
void 	startmatrix(double **co,
					double **matrix,
					double **matrixrot);
double 	bestangle(double**, int, int, double*);
double 	bestangle_anchor(double**, int, int, double**);
double  bestangle_anchor(double** co,int bond_number,int number_atoms,double** anchor,
			double** matrix,double** polpoint,double** anchorpoint,double** result1,double** result2, double** result3);  //	Preallocated
double  bestangle_anchor_new(double** co,int bond_number,int number_atoms,double** anchor,
			double** matrix,double** polpoint,double** anchorpoint,double** result1,double** result2, double** result3);  //	Preallocated
void 	bestanglesm3(double** co,double *angles,int *bondsequence,int bondrange,int number_atoms,double** anchor,						//	Preallocated version
			double** matrix,double** polpoint,double** anchorpoint,double** result1,double** result2, double** result3) ;
void bestanglesm(double** co,		/*loop coordinates*/
				double *angles,		/*output dihedral angle variations*/
				double *angles2,		/*output dihedral angle variations*/
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
					);
void 	bestangles(double** co,double *angles,int *bondsequence,int bondrange,int number_atoms,double** anchor,
			double** matrix,double** polpoint,double** anchorpoint,double** result1,double** result2, double** result3);
//void 	bestanglesm(double**,double *,int *,int,int,double**,double**,double**,double**,double**,double**, double**);


void 	opt_rotation(double **, double **,double **,int);

void 	fix_anchor(double **,double **,double **,double **,int);
void    fix_cramer(double **,double **,double **,double **,int);
void 	gafixanchor(double **,double **,double **,int,double**,double**,double**);


int		bond_select3(int range);
int		bond_select4(int range,int *bond_previous);
int 	bond_select5(int range, int *bond_previous);
int		bond_select6(int,int *,int *);

void constructres(char *pdbname, int start, int end, int nratoms,int *marker);
void constructpsi(int *marker, int nratoms);
void    constructmarker(int *, int );
int	    bondselect(int ,int *,int *,int *);
int	   bondselect(int,int *,int *,int *, int*);
int		bond_select7(int ,int *);
int		bond_select_oddeven(int*,int,int*);

void	centroid(double **, double *,int);

void 	relativeco(double **,double **,int);

double  rmsd_anchor2(double**,double**);
double  rmsd_anchor3(double**,double**);
void 	movechain(double**,double**,int,int,int);
void 	selectpart(double**,int,int,double** part);
void 	opt_shift(double **, double**,double **,int);
void 	rand_shift(double **,int,double);
void    startspinor( double **co,		    /*loop coordinates*/
					 double **matrix,	    /*output matrix*/
					 double **matrixrot,	/*temporary matrix-preallocated*/
					 double *spinor);
//void    startspinor2( double **co,		    /*loop coordinates*/
//					  int number_atoms_loop,
//					  double *spinor);



#endif /* GUARD_chain_manip_h_ */

/*===============================================================================
 * CHAIN_MANIP.H
 *
 *  Created on: Dec 10, 2010
 *      Author: Pieter Chys
 *===============================================================================
 */
