/*
 * coutsias.h
 *
 *  Created on: Jun 17, 2016
 *      Author: mon
 *
 *  Here you will find all defines, function definitions, and globals from loop.h and strum.h files.
 *  Some day code must be cleaned and use only globals when strictly necessary (just for clarity and mental health...;-)
 *
 */

#ifndef INCLUDE_COUTSIAS_H_
#define INCLUDE_COUTSIAS_H_

//
// From loop.h
//
#define  pi 3.141592653589793238462643383279502884197e0
#define two_pi 2.0e0*pi
// MON: Not needed:
//#define deg2rad pi/180.0e0
//#define rad2deg 180.0e0/pi
#define Max(a,b) ((a) > (b))? (a) : (b)
#define Min(a,b) ((a) < (b))? (a) : (b)
int print_level = 0;
//  ! parameters for tripeptide loop (including bond lengths & angles)
double len0[6], b_ang0[7], t_ang0[2];
double aa13_min_sqr, aa13_max_sqr;
double delta[4], xi[3], eta[3], alpha[3], theta[3];
double cos_alpha[3], sin_alpha[3], cos_theta[3], sin_theta[3];
double cos_delta[4], sin_delta[4];
double cos_xi[3], cos_eta[3], sin_xi[3], sin_eta[3];
double r_a1a3[3], r_a1n1[3], r_a3c3[3];
double b_a1a3[3], b_a1n1[3], b_a3c3[3];
double len_na[3], len_ac[3], len_aa[3];
//  ! used for polynomial coefficients
double C0[3][3], C1[3][3], C2[3][3];
double Q[5][17], R[3][17];

#define MAX_SOLN 16
#define deg_pol 16
double dot_productL(double va[], double vb[]);
void matmulL(double ma[3][3], double mb[3], double mc[3]);
double signL(double a, double b);
void solve_3pep_poly(double r_n1[3], double r_a1[3], double r_a3[3], double r_c3[3], double r_soln_n[MAX_SOLN][3][3], double r_soln_a[MAX_SOLN][3][3], double r_soln_c[MAX_SOLN][3][3], int *n_soln);
void initialize_loop_closure(double b_len[], double b_ang[], double t_ang[]);
void get_input_angles(int *n_soln, double r_n1[], double r_a1[], double r_a3[], double r_c3[]);
void test_two_cone_existence_soln(double tt, double kx, double et, double ap, int *n_soln, char cone_type[]);
void get_poly_coeff(double poly_coeff[]);
void poly_mul_sub2(double u1[5][5], double u2[5][5], double u3[5][5], double u4[5][5], int p1[2], int p2[2], int p3[2], int p4[2], double u5[5][5], int p5[2]);
void poly_mul2(double u1[5][5], double u2[5][5], int p1[2], int p2[2], double u3[5][5], int p3[2]);
void poly_sub2(double u1[5][5], double u2[5][5], int p1[2], int p2[2], double u3[5][5], int p3[2]);
void poly_mul_sub1(double u1[], double u2[], double u3[], double u4[], int p1, int p2, int p3, int p4, double u5[], int *p5);
void poly_mul1(double u1[], double u2[], int p1, int p2, double u3[], int *p3);
void poly_sub1(double u1[], double u2[], int p1, int p2, double u3[], int *p3);
void coord_from_poly_roots(int *n_soln, double roots[MAX_SOLN], double r_n1[3], double r_a1[3], double r_a3[3], double r_c3[3], double r_soln_n[MAX_SOLN][3][3], double r_soln_a[MAX_SOLN][3][3], double r_soln_c[MAX_SOLN][3][3]);
double calc_t2(double t0);
double calc_t1(double t0, double t2);
void calc_dih_ang(double r1[], double r2[], double r3[], double *angle);
void calc_bnd_ang(double r1[], double r2[], double *angle);
void cross_tc(double p[], double q[], double s[]);
void quaternionL(double axis[], double quarter_ang, double p[]);
void rotation_matrixL(double q[4], double U[3][3]);

//
// From strum.h
//
#define PRINT_LEVEL 0
#define	MAX_ORDER 16
#define MAXPOW 32
#define	SMALL_ENOUGH 1.0e-18
typedef struct pol { // structure type for representing a polynomial
  int	ord;
  double coef[MAX_ORDER+1];
} poly;
double RELERROR;
int MAXIT, MAX_ITER_SECANT;
void initialize_sturm(double *tol_secant, int *max_iter_sturm, int *max_iter_secant);
void solve_sturm(int p_order, int *n_root, double *poly_coeffs, double *roots);
double hyper_tan(double a, double x);
static int modpL(poly *u, poly *v, poly *r);
int buildsturm(int ord, poly *sseq);
int numroots(int np, poly *sseq, int *atneg, int *atpos);
int numchanges(int np, poly *sseq, double a);
void sbisect(int np, poly *sseq, double min, double max, int atmin, int atmax, double *roots);
double evalpoly(int ord, double *coef, double x);
int modrf(int ord, double *coef, double	a, double b, double *val);


#endif /* INCLUDE_COUTSIAS_H_ */
