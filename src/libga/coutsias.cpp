/*
 * coutsias.cpp
 *
 *  Created on: Jun 17, 2016
 *      Author: mon
 *
 *  Here you will find all the meat of the Coutsias et al. analytic loop closure algorithm.
 *
 */

//!----------------------------------------------------------------------
//! Copyright (C) 2003
//!      Chaok Seok, Evangelos Coutsias, Matthew Jacobson, and Ken Dill
//!      UCSF and Univeristy of New Mexico
//! Witten by Chaok Seok 2003.
//!----------------------------------------------------------------------------
//!----------------------------------------------------------------------------
//!*************** Tripeptide Loop Closure Algorithm *****************
//! files to be compiled with:
//!    tripeptide_closure.f90
//!    sturm.c
//!
//!*******************************************************************
//! subroutine  initialize_loop_closure(b_len, b_ang, t_ang)
//!*******************************************************************
//! connectivity of atoms:
//!   N1-A1-C1-N2-A2-C2-N3-A3-C3
//!
//! input:
//!
//  * b_len(1:6): bond lengths (A1-C1, C1-N2, ..., N3-A3)
//!  * b_ang(1:7): bond angles (N1-A1-C1, A1-C1-N2, ..., N3-A3-C3)
//!  * t_ang(1:2): torsion angles (A1-C1-N2-A2, A2-C2-N3-A3)
//!*******************************************************************
//!
//!*******************************************************************
//! subroutine solv_3pep_poly(r_n1, r_a1, r_a3, r_c3, &
//!     r_soln_n, r_soln_a, r_soln_c, n_soln)
//!*******************************************************************
//! input:
//!  * r_n1(3), r_a1(3), r_a3(3), r_c3(3):
//!       Cartesian coordinates of N and CA atoms of the first residue and
//!        CA and C atoms of the last (third) residue.
//! output:
//!  * n_soln: number of alternative loop closure solutions.
//!  * r_soln_n(3,3,8), r_soln_a(3,3,8), r_soln_c(3,3,8):
//!       Cartesian coordinates of loop closure solutions.
//!       first dimension: x, y, z component
//!       second dim: residue number
//!       third dim: solution number
//!*******************************************************************

#include "coutsias.h"
#include "rcd_stuff.h"

//
// FUNCTIONS FROM loop.h
//

double dot_productL(double va[3], double vb[3])
{
 return va[0]*vb[0] + va[1]*vb[1] + va[2]*vb[2];
}

void matmulL(double ma[3][3], double mb[3], double mc[3])
 {
  int i, j;

  for(i=0;i<3;i++)
   {
    mc[i] = 0.;
    for(j=0;j<3;j++)
      mc[i] += ma[i][j]*mb[j];
   }
  return;
 }

double signL(double a, double b)
 {
  if (b>=0.)
    return fabs(a);
  else
    return -fabs(a);
 }

void solve_3pep_poly(double r_n1[3], double r_a1[3], double r_a3[3], double r_c3[3],
                     double r_soln_n[MAX_SOLN][3][3], double r_soln_a[MAX_SOLN][3][3],
                     double r_soln_c[MAX_SOLN][3][3], int *n_soln)
 {

  double poly_coeff[deg_pol+1], roots[MAX_SOLN];

  get_input_angles(n_soln, r_n1, r_a1, r_a3, r_c3);

  if (*n_soln == 0)
    return;

  get_poly_coeff(poly_coeff);

  solve_sturm(deg_pol, n_soln, poly_coeff, roots);

  if (*n_soln == 0) {
    //fprintf(stderr, "No hay ninguna solucion\n");
    return;
    //exit(1);
  }

  coord_from_poly_roots(n_soln, roots, r_n1, r_a1, r_a3, r_c3, r_soln_n, r_soln_a, r_soln_c);

  return;
 }

void initialize_loop_closure(double b_len[6], double b_ang[7], double t_ang[2])
 {

//!-----------------------------------------------------------------------
//! Input angles for the given bond lengths and angles
//!-----------------------------------------------------------------------
  double len1, len2, a_min, a_max;
  double axis[3], rr_a1[3], rr_c1[3], rr_n2[3], rr_a2[3], rr_n2a2_ref[3], rr_c1a1[3];
  double rr_a1a2[3], dr[3], bb_c1a1[3], bb_a1a2[3], bb_a2n2[3];
  double p[4], Us[3][3];
  double mulpro[3];
  double tmp_val[3];
  double tol_secant = 1.0e-15;
  int max_iter_sturm = 100;
  int max_iter_secant = 20;
  int i, j;

  initialize_sturm(&tol_secant, &max_iter_sturm, &max_iter_secant);

  for(i=0;i<6;i++)
    len0[i] = b_len[i];

  for(i=0;i<7;i++){
    b_ang0[i] = b_ang[i];
  }

  for(i=0;i<2;i++)
    t_ang0[i] = t_ang[i];

  for(i=0;i<3;i++)
    rr_c1[i] = 0.;

  axis[0] = 1.;
  axis[1] = 0.;
  axis[2] = 0.;

  for(i=0;i<2;i++)
   {
     rr_a1[0] = cos(b_ang0[3*i+1])*len0[3*i];
     rr_a1[1] = sin(b_ang0[3*i+1])*len0[3*i];
     rr_a1[2] = 0.0e0;

     rr_n2[0] = len0[3*i+1];
     rr_n2[1] = 0.0e0;
     rr_n2[2] = 0.0e0;

     for(j=0;j<3;j++)
       rr_c1a1[j] = rr_a1[j] - rr_c1[j];

     rr_n2a2_ref[0] = -cos(b_ang0[3*i+2])*len0[3*i+2];
     rr_n2a2_ref[1] = sin(b_ang0[3*i+2])*len0[3*i+2];
     rr_n2a2_ref[2] = 0.0e0;

     quaternionL(axis, t_ang0[i]*0.25e0, p);

     rotation_matrixL(p, Us);

     matmulL(Us, rr_n2a2_ref, mulpro);
     for(j=0;j<3;j++)
      {
       rr_a2[j] =  mulpro[j] + rr_n2[j];
       rr_a1a2[j] = rr_a2[j] - rr_a1[j];
       dr[j] = rr_a1a2[j];
      }

     len2 = dot_productL(dr, dr);

     len1 = sqrt(len2);
     len_aa[i+1] = len1;

     for(j=0;j<3;j++)
      {
       bb_c1a1[j] = rr_c1a1[j]/len0[3*i];
       bb_a1a2[j] = rr_a1a2[j]/len1;
       bb_a2n2[j] = (rr_n2[j] - rr_a2[j])/len0[3*i+2];
      }

//     ! xi
     for(j=0;j<3;j++)
       tmp_val[j] = -bb_a1a2[j];
     calc_bnd_ang(tmp_val, bb_a2n2, &xi[i+1]);

//     ! eta
     for(j=0;j<3;j++)
       tmp_val[j] = -bb_c1a1[j];
     calc_bnd_ang(bb_a1a2, tmp_val, &eta[i]);

//     ! delta: pi -  dih of N(1)CA(1)CA(3)C(3)
     calc_dih_ang(bb_c1a1, bb_a1a2, bb_a2n2, &delta[i+1]);
     delta[i+1] = pi - delta[i+1];
   }

  a_min = b_ang[3] - (xi[1] + eta[1]);
  a_max = Min((b_ang[3] + (xi[1] + eta[1])), pi);

//  ! min/max of base length
  aa13_min_sqr = pow(len_aa[1],2) + pow(len_aa[2],2) - 2.0e0*len_aa[1]*len_aa[2]*cos(a_min);
  aa13_max_sqr = pow(len_aa[1],2) + pow(len_aa[2],2) - 2.0e0*len_aa[1]*len_aa[2]*cos(a_max);

 }

void get_input_angles(int *n_soln, double r_n1[3], double r_a1[3], double r_a3[3], double r_c3[3])
 {
//!-----------------------------------------------------------------------
//! Input angles and vectors (later used in coordinates)
//!-----------------------------------------------------------------------
  double dr_sqr;
  double tmp_val[3];
  int i;
  char cone_type[2];
//!-----------------------------------------------------------------------

  *n_soln = MAX_SOLN;

//  ! vertual bond
  for(i=0;i<3;i++)
    r_a1a3[i] = r_a3[i] - r_a1[i];

  dr_sqr = dot_productL(r_a1a3,r_a1a3);

  len_aa[0] = sqrt(dr_sqr);

  if ((dr_sqr < aa13_min_sqr) || (dr_sqr > aa13_max_sqr))
   {
    *n_soln = 0;
    return;
   }

//  ! bond lengths
  for(i=0;i<3;i++)
    r_a1n1[i] = r_n1[i] - r_a1[i];

  len_na[0] = sqrt(dot_productL(r_a1n1,r_a1n1));
  len_na[1] = len0[2];
  len_na[2] = len0[5];

  for(i=0;i<3;i++)
    r_a3c3[i] = r_c3[i] - r_a3[i];

  len_ac[0] = len0[0];
  len_ac[1] = len0[3];
  len_ac[2] = sqrt(dot_productL(r_a3c3,r_a3c3));

//  ! unit vectors
  for(i=0;i<3;i++)
   {
    b_a1n1[i] = r_a1n1[i]/len_na[0];
    b_a3c3[i] = r_a3c3[i]/len_ac[2];
    b_a1a3[i] = r_a1a3[i]/len_aa[0];
   }

//  ! delta(3): dih of N(1)CA(1)CA(3)C(3)
  for(i=0;i<3;i++)
    tmp_val[i] = -b_a1n1[i];
  calc_dih_ang(tmp_val, b_a1a3, b_a3c3, &delta[3]);

  delta[0] = delta[3];

//  ! xi(1)
  for(i=0;i<3;i++)
    tmp_val[i] = -b_a1a3[i];
  calc_bnd_ang(tmp_val, b_a1n1, &xi[0]);

//  ! eta(3)
  calc_bnd_ang(b_a1a3, b_a3c3, &eta[2]);

  for(i=0;i<3;i++)
   {
     cos_delta[i+1] = cos(delta[i+1]);
     sin_delta[i+1] = sin(delta[i+1]);
     cos_xi[i] = cos(xi[i]);
     sin_xi[i] = sin(xi[i]);
     sin_xi[i] = sin(xi[i]);
     cos_eta[i] = cos(eta[i]);
     cos_eta[i] = cos(eta[i]);
     sin_eta[i] = sin(eta[i]);
     sin_eta[i] = sin(eta[i]);
   }

  cos_delta[0] = cos_delta[3];
  sin_delta[0] = sin_delta[3];

//  ! theta (N, CA, C) bond angle
  theta[0] = b_ang0[0];
  theta[1] = b_ang0[3];
  theta[2] = b_ang0[6];

  for(i=0;i<3;i++)
    cos_theta[i] = cos(theta[i]);

//  ! alpha
  cos_alpha[0] = -(pow(len_aa[0],2) + pow(len_aa[1],2) - pow(len_aa[2],2))/(2.0e0*len_aa[0]*len_aa[1]);
  alpha[0] = acos(cos_alpha[0]);
  sin_alpha[0] = sin(alpha[0]);
  cos_alpha[1] = (pow(len_aa[1],2) + pow(len_aa[2],2) - pow(len_aa[0],2))/(2.0e0*len_aa[1]*len_aa[2]);
  alpha[1] = acos(cos_alpha[1]);
  sin_alpha[1] = sin(alpha[1]);
  alpha[2] = pi - alpha[0] + alpha[1];
  cos_alpha[2] = cos(alpha[2]);
  sin_alpha[2] = sin(alpha[2]);


  if (print_level > 0)
   {
     printf("xi = %9.4f%9.4f%9.4f\n", xi[0]*rad2deg, xi[1]*rad2deg, xi[2]*rad2deg);
     printf("eta = %9.4f%9.4f%9.4f\n", eta[0]*rad2deg, eta[1]*rad2deg, eta[2]*rad2deg);
     printf("delta = %9.4f%9.4f%9.4f\n", delta[1]*rad2deg, delta[2]*rad2deg, delta[3]*rad2deg);
     printf("theta = %9.4f%9.4f%9.4f\n", theta[0]*rad2deg, theta[1]*rad2deg, theta[2]*rad2deg);
     printf("alpha = %9.4f%9.4f%9.4f\n", alpha[0]*rad2deg, alpha[1]*rad2deg, alpha[2]*rad2deg);
   }

//  ! check for existence of soln
  for(i=0;i<3;i++)
   {
//     call test_two_cone_existence_soln(theta(i), xi(i), eta(i), alpha(i), &
//          n_soln, cone_type)
    test_two_cone_existence_soln(theta[i], xi[i], eta[i], alpha[i], n_soln, cone_type);
    if (*n_soln == 0)
      return;
   }

  return;
 }


void test_two_cone_existence_soln(double tt, double kx, double et, double ap, int *n_soln, char cone_type[2])
 {

  double at, ex, abs_at, ap1, kx1, et1;
  double cos_tx1, cos_tx2, cos_te1, cos_te2, cos_ea1, cos_ea2, cos_xa1, cos_xa2;
  int s1, s2, t1, t2;
  int complicated = 0;

  *n_soln = MAX_SOLN;

  ap1 = ap;
  kx1 = kx;
  et1 = et;

  at = ap1 - tt;
  ex = kx1 + et1;
  abs_at = fabs(at);

  if (abs_at > ex)
   {
    *n_soln = 0;
    return;
   }

  if (complicated)
   {
    cos_tx1 = cos(tt+kx1);
    cos_tx2 = cos(tt-kx1);
    cos_te1 = cos(tt+et1);
    cos_te2 = cos(tt-et1);
    cos_ea1 = cos(et1+ap1);
    cos_ea2 = cos(et1-ap1);
    cos_xa1 = cos(kx1+ap1);
    cos_xa2 = cos(kx1-ap1);

    s1 = 0;
    s2 = 0;
    t1 = 0;
    t2 = 0;

    if ((cos_te1-cos_xa2)*(cos_te1-cos_xa1) <= 0.0e0)
      s1 = 0;

    if ((cos_te2-cos_xa2)*(cos_te2-cos_xa1) <= 0.0e0)
      s2 = 0;

    if ((cos_tx1-cos_ea2)*(cos_tx1-cos_ea1) <= 0.0e0)
      t1 = 0;

    if ((cos_tx2-cos_ea2)*(cos_tx2-cos_ea1) <= 0.0e0)
      t2 = 0;
   }

  return;
 }

void get_poly_coeff(double poly_coeff[deg_pol+1])
 {

  int i, j;
  double A0, A1, A2, A3, A4, A21, A22, A31, A32, A41, A42;
  double B0[3], B1[3], B2[3], B3[3], B4[3], B5[3], B6[3], B7[3], B8[3];
  double u11[5][5], u12[5][5], u13[5][5], u31[5][5], u32[5][5], u33[5][5];
  double um1[5][5], um2[5][5], um3[5][5], um4[5][5], um5[5][5], um6[5][5], q_tmp[5][5];
  int p1 [2], p3[2], p_um1[2], p_um2[2], p_um3[2], p_um4[2], p_um5[2], p_um6[2], p_Q[2];
  int p2, p4, p_f1, p_f2, p_f3, p_f4, p_f5, p_f6, p_f7, p_f8, p_f9;
  int p_f10, p_f11, p_f12, p_f13, p_f14, p_f15, p_f16, p_f17, p_f18;
  int p_f19, p_f20, p_f21, p_f22, p_f23, p_f24, p_f25, p_f26;
  int p_final;
  double f1[17], f2[17], f3[17], f4[17], f5[17], f6[17], f7[17], f8[17], f9[17];
  double f10[17], f11[17], f12[17], f13[17], f14[17], f15[17], f16[17], f17[17], f18[17];
  double f19[17], f20[17], f21[17], f22[17], f23[17], f24[17], f25[17], f26[17];

//  ! A0, B0
  for(i=0;i<3;i++)
   {
    A0 = cos_alpha[i]*cos_xi[i]*cos_eta[i] - cos_theta[i];
    A1 = -sin_alpha[i]*cos_xi[i]*sin_eta[i];
    A2 = sin_alpha[i]*sin_xi[i]*cos_eta[i];
    A3 = sin_xi[i]*sin_eta[i];
    A4 = A3*cos_alpha[i];

    j = i;

    A21 = A2*cos_delta[j];
    A22 = A2*sin_delta[j];
    A31 = A3*cos_delta[j];
    A32 = A3*sin_delta[j];
    A41 = A4*cos_delta[j];
    A42 = A4*sin_delta[j];

    B0[i] = A0 + A22 + A31;
    B1[i] = 2.0e0*(A1 + A42);
    B2[i] = 2.0e0*(A32 - A21);
    B3[i] = -4.0e0*A41;
    B4[i] = A0 + A22 - A31;
    B5[i] = A0 - A22 - A31;
    B6[i] = -2.0e0*(A21 + A32);
    B7[i] = 2.0e0*(A1 - A42);
    B8[i] = A0 - A22 + A31;
   }

//  ! C0i
  i = 0;

  C0[i][0] = B0[i];
  C0[i][1] = B2[i];
  C0[i][2] = B5[i];

  C1[i][0] = B1[i];
  C1[i][1] = B3[i];
  C1[i][2] = B7[i];

  C2[i][0] = B4[i];
  C2[i][1] = B6[i];
  C2[i][2] = B8[i];

  for(i=1;i<3;i++)
   {
    C0[i][0] = B0[i];
    C0[i][1] = B1[i];
    C0[i][2] = B4[i];

    C1[i][0] = B2[i];
    C1[i][1] = B3[i];
    C1[i][2] = B6[i];

    C2[i][0] = B5[i];
    C2[i][1] = B7[i];
    C2[i][2] = B8[i];
   }

//  ! first determinant
  for(i=0;i<3;i++)
   {
     u11[0][i] = C0[0][i];
     u12[0][i] = C1[0][i];
     u13[0][i] = C2[0][i];
     u31[i][0] = C0[1][i];
     u32[i][0] = C1[1][i];
     u33[i][0] = C2[1][i];
   }

  p1[0] = 2;
  p1[1] = 0;

  p3[0] = 0;
  p3[1] = 2;

  poly_mul_sub2(u32, u32, u31, u33, p3, p3, p3, p3, um1, p_um1);
  poly_mul_sub2(u12, u32, u11, u33, p1, p3, p1, p3, um2, p_um2);
  poly_mul_sub2(u12, u33, u13, u32, p1, p3, p1, p3, um3, p_um3);
  poly_mul_sub2(u11, u33, u31, u13, p1, p3, p3, p1, um4, p_um4);
  poly_mul_sub2(u13, um1, u33, um2, p1, p_um1, p3, p_um2, um5, p_um5);
  poly_mul_sub2(u13, um4, u12, um3, p1, p_um4, p1, p_um3, um6, p_um6);
  poly_mul_sub2(u11, um5, u31, um6, p1, p_um5, p3, p_um6, q_tmp, p_Q);

  for(i=0;i<5;i++)
    for(j=0;j<5;j++)
      Q[i][j] = q_tmp[i][j];

//  ! second determinant
  for(i=0;i<3;i++)
    for(j=0;j<17;j++)
      R[i][j] = 0.;

  for(i=0;i<3;i++)
   {
    R[0][i] = C0[2][i];
    R[1][i] = C1[2][i];
    R[2][i] = C2[2][i];
   }

  p2 = 2;
  p4 = 4;

  poly_mul_sub1(R[1], R[1], R[0], R[2], p2, p2, p2, p2, f1, &p_f1);
  poly_mul1(R[1], R[2], p2, p2, f2, &p_f2);
  poly_mul_sub1(R[1], f1, R[0], f2, p2, p_f1, p2, p_f2, f3, &p_f3);
  poly_mul1(R[2], f1, p2, p_f1, f4, &p_f4);
  poly_mul_sub1(R[1], f3, R[0], f4, p2, p_f3, p2, p_f4, f5, &p_f5);

  poly_mul_sub1(Q[1], R[1], Q[0], R[2], p4, p2, p4, p2, f6, &p_f6);
  poly_mul_sub1(Q[2], f1, R[2], f6, p4, p_f1, p2, p_f6, f7, &p_f7);
  poly_mul_sub1(Q[3], f3, R[2], f7, p4, p_f3, p2, p_f7, f8, &p_f8);
  poly_mul_sub1(Q[4], f5, R[2], f8, p4, p_f5, p2, p_f8, f9, &p_f9);

  poly_mul_sub1(Q[3], R[1], Q[4], R[0], p4, p2, p4, p2, f10, &p_f10);
  poly_mul_sub1(Q[2], f1, R[0], f10, p4, p_f1, p2, p_f10, f11, &p_f11);
  poly_mul_sub1(Q[1], f3, R[0], f11, p4, p_f3, p2, p_f11, f12, &p_f12);

  poly_mul_sub1(Q[2], R[1], Q[1], R[2], p4, p2, p4, p2, f13, &p_f13);
  poly_mul_sub1(Q[3], f1, R[2], f13, p4, p_f1, p2, p_f13, f14, &p_f14);
  poly_mul_sub1(Q[3], R[1], Q[2], R[2], p4, p2, p4, p2, f15, &p_f15);
  poly_mul_sub1(Q[4], f1, R[2], f15, p4, p_f1, p2, p_f15, f16, &p_f16);
  poly_mul_sub1(Q[1], f14, Q[0], f16, p4, p_f14, p4, p_f16, f17, &p_f17);

  poly_mul_sub1(Q[2], R[2], Q[3], R[1], p4, p2, p4, p2, f18, &p_f18);
  poly_mul_sub1(Q[1], R[2], Q[3], R[0], p4, p2, p4, p2, f19, &p_f19);
  poly_mul_sub1(Q[3], f19, Q[2], f18, p4, p_f19, p4, p_f18, f20, &p_f20);
  poly_mul_sub1(Q[1], R[1], Q[2], R[0], p4, p2, p4, p2, f21, &p_f21);
  poly_mul1(Q[4], f21, p4, p_f21, f22, &p_f22);
  poly_sub1(f20, f22, p_f20, p_f22, f23, &p_f23);
  poly_mul1(R[0], f23, p2, p_f23, f24, &p_f24);
  poly_sub1(f17, f24, p_f17, p_f24, f25, &p_f25);
  poly_mul_sub1(Q[4], f12, R[2], f25, p4, p_f12, p2, p_f25, f26, &p_f26);
  poly_mul_sub1(Q[0], f9, R[0], f26, p4, p_f9, p2, p_f26, poly_coeff, &p_final);

  if (p_final != 16)
   {
    printf("Error. Degree of polynomial is not 16!\n");
    exit(1);
   }

  if (poly_coeff[16] < 0.0e0)
    for(i=0;i<17;i++)
     poly_coeff[i] *= -1.0;

  if (print_level > 0)
   {
     printf("poly_coeff\n");
     for(i=0;i<17;i++)
        printf("%5d%15.6f\n", i, poly_coeff[i]);
   }

  return;
 }

void poly_mul_sub2(double u1[5][5], double u2[5][5], double u3[5][5], double u4[5][5], int p1[2], int p2[2], int p3[2], int p4[2], double u5[5][5], int p5[2])
 {

  double d1[5][5], d2[5][5];
  int pd1[2], pd2[2];

  poly_mul2(u1, u2, p1, p2, d1, pd1);
  poly_mul2(u3, u4, p3, p4, d2, pd2);
  poly_sub2(d1, d2, pd1, pd2, u5, p5);

 }

void poly_mul2(double u1[5][5], double u2[5][5], int p1[2], int p2[2], double u3[5][5], int p3[2])
 {
  int i1, j1, i2, j2, i3, j3, p11, p12, p21, p22;
  int i, j;
  double u1ij;

  for(i=0;i<2;i++)
    p3[i] = p1[i] + p2[i];
  for(i=0;i<5;i++)
    for(j=0;j<5;j++)
      u3[i][j] = 0.0e0;

  p11 = p1[0];
  p12 = p1[1];
  p21 = p2[0];
  p22 = p2[1];

  for(i1=0;i1<=p12;i1++)
   {
     for(j1=0;j1<=p11;j1++)
      {
        u1ij = u1[i1][j1];
        for(i2=0;i2<=p22;i2++)
         {
           i3 = i1 + i2;
           for (j2=0;j2<=p21;j2++)
            {
              j3 = j1 + j2;
              u3[i3][j3] = u3[i3][j3] + u1ij*u2[i2][j2];
            }
         }
      }
   }

 }

void poly_sub2(double u1[5][5], double u2[5][5], int p1[2], int p2[2], double u3[5][5], int p3[2])
 {

  int i, j, p11, p12, p21, p22;
  int i1_ok, i2_ok;

  p11 = p1[0];
  p12 = p1[1];
  p21 = p2[0];
  p22 = p2[1];
  p3[0] = Max(p11,p21);
  p3[1] = Max(p12,p22);

  for(i=0;i<=p3[1];i++)
   {
    i1_ok = (i > p12);
    i2_ok = (i > p22);
    for(j=0;j<=p3[0];j++)
     {
      if (i2_ok || (j > p21))
        u3[i][j] = u1[i][j];
      else if (i1_ok || (j > p11))
        u3[i][j] = -u2[i][j];
      else
        u3[i][j] = u1[i][j] - u2[i][j];
     }
   }

  return;
 }

void poly_mul_sub1(double u1[17], double u2[17], double u3[17], double u4[17], int p1, int p2, int p3, int p4, double u5[17], int *p5)
 {

  double d1[17], d2[17];
  int pd1, pd2;

  poly_mul1(u1, u2, p1, p2, d1, &pd1);
  poly_mul1(u3, u4, p3, p4, d2, &pd2);
  poly_sub1(d1, d2, pd1, pd2, u5, p5);

  return;
 }

void poly_mul1(double u1[17], double u2[17], int p1, int p2, double u3[17], int *p3)
 {

  int i, i1, i2, i3;
  double u1i;

  *p3 = p1 + p2;

  for(i=0;i<17;i++)
    u3[i] = 0.;

  for(i1=0;i1<=p1;i1++)
   {
    u1i = u1[i1];
    for(i2=0;i2<=p2;i2++)
     {
        i3 = i1 + i2;
        u3[i3] = u3[i3] + u1i*u2[i2];
     }
   }

  return;
 }

void poly_sub1(double u1[17], double u2[17], int p1, int p2, double u3[17], int *p3)
 {
  int i;

  *p3 = Max(p1, p2);

  for(i=0;i<=*p3;i++)
   {
    if (i > p2)
      u3[i] = u1[i];
    else if (i > p1)
      u3[i] = -u2[i];
    else
      u3[i] = u1[i] - u2[i];
   }

  return;
 }

void coord_from_poly_roots(int *n_soln, double roots[MAX_SOLN], double r_n1[3], double r_a1[3], double r_a3[3], double r_c3[3], double r_soln_n[MAX_SOLN][3][3], double r_soln_a[MAX_SOLN][3][3], double r_soln_c[MAX_SOLN][3][3])
 {

  double ex[3], ey[3], ez[3], b_a1a2[3], b_a3a2[3], r_tmp[3];
  double p_s[3][3], s1[3][3], s2[3][3], p_t[3][3], t1[3][3], t2[3][3];
  double p_s_c[3][3], s1_s[3][3], s2_s[3][3], p_t_c[3][3], t1_s[3][3], t2_s[3][3];
  double angle, sig1_init, half_tan[3];
  double cos_tau[4], sin_tau[4], cos_sig[3], sin_sig[3], ht, tmp, sig1;
  double r_s[3], r_t[3], r0[3], r_n[3][3], r_a[3][3], r_c[3][3], p[4], Us[3][3];
  int i_soln, i, j;
  double a1c1, c1n2, n2a2, a2c2, c2n3, n3a3, a1a2, a2a3;
  double rr_a1c1[3], rr_c1n2[3], rr_n2a2[3], rr_a2c2[3], rr_c2n3[3], rr_n3a3[3], rr_a1a2[3], rr_a2a3[3];
  double a3a1a2, a2a3a1, n1a1c1, n2a2c2, n3a3c3, a1c1n2a2, a2c2n3a3;
  double tmp_value, ex_tmp[3];
  double tmp_array[3], tmp_array1[3], tmp_array2[3], tmp_array3[3];
  double mat1[3], mat2[3], mat3[3], mat4[3], mat5[3];
  double mat11[3], mat22[3], mat33[3], mat44[3], mat55[3];

  if (*n_soln == 0)
    return;

//  ! Define body frame (ex, ey, ez)
  for(i=0;i<3;i++)
    ex[i] = b_a1a3[i];

  cross_tc(r_a1n1, ex, ez);

  tmp_value = sqrt(dot_productL(ez,ez));
  for(i=0;i<3;i++)
    ez[i] = ez[i]/tmp_value;

  cross_tc(ez, ex, ey);

//  ! vertual bond vectors in the reference plane
  for(i=0;i<3;i++)
   {
    b_a1a2[i] = -cos_alpha[0]*ex[i] + sin_alpha[0]*ey[i];
    b_a3a2[i] = cos_alpha[2]*ex[i] + sin_alpha[2]*ey[i];
   }

//  !! Define cone coordinates for each angle joint.
  for(i=0;i<3;i++)
   {
    p_s[0][i] = -ex[i];
    s1[0][i]  = ez[i];
    s2[0][i]  = ey[i];
    p_t[0][i] = b_a1a2[i];
    t1[0][i]  = ez[i];
    t2[0][i]  = sin_alpha[0]*ex[i] + cos_alpha[0]*ey[i];
   }

//  ! residue 2
  for(i=0;i<3;i++)
   {
    p_s[1][i] = -b_a1a2[i];
    s1[1][i]  = -ez[i];
    s2[1][i]  = t2[0][i];
    p_t[1][i] = -b_a3a2[i];
    t1[1][i]  = -ez[i];
    t2[1][i]  = sin_alpha[2]*ex[i] - cos_alpha[2]*ey[i];
   }

//  ! residue 3
  for(i=0;i<3;i++)
   {
    p_s[2][i] = b_a3a2[i];
    s2[2][i]  = t2[1][i];
    s1[2][i]  = ez[i];
    p_t[2][i] = ex[i];
    t1[2][i] =  ez[i];
    t2[2][i] = -ey[i];
   }

//  ! scale vectors
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
     {
      p_s_c[i][j] = p_s[i][j]*cos_xi[i];
      s1_s[i][j]  = s1[i][j]*sin_xi[i];
      s2_s[i][j]  = s2[i][j]*sin_xi[i];
      p_t_c[i][j] = p_t[i][j]*cos_eta[i];
      t1_s[i][j]  = t1[i][j]*sin_eta[i];
      t2_s[i][j]  = t2[i][j]*sin_eta[i];
     }

  for(i=0;i<3;i++)
    r_tmp[i] = (r_a1n1[i]/len_na[0] - p_s_c[0][i])/sin_xi[0];

  calc_bnd_ang(s1[0], r_tmp, &angle);
  sig1_init = signL(angle, dot_productL(r_tmp,s2[0]));

//  ! CA
  for(i=0;i<3;i++)
   {
    r_a[0][i] = r_a1[i];
    r_a[1][i] = r_a1[i] + len_aa[1]*b_a1a2[i];
    r_a[2][i] = r_a3[i];
    r0[i] = r_a1[i];
   }

  for(i_soln=0;i_soln<*n_soln;i_soln++)
   {
     half_tan[2] = roots[i_soln];
     half_tan[1] = calc_t2(half_tan[2]);
     half_tan[0] = calc_t1(half_tan[2], half_tan[1]);

     for(i=1;i<=3;i++)
      {
       ht = half_tan[i-1];
       tmp = 1.0e0 + ht*ht;
       cos_tau[i] = (1.0e0 - ht*ht)/tmp;
       sin_tau[i] = 2.0e0*ht/tmp;
      }

     cos_tau[0] = cos_tau[3];

     sin_tau[0] = sin_tau[3];

     for(i=0;i<3;i++)
      {
       cos_sig[i] = cos_delta[i]*cos_tau[i] + sin_delta[i]*sin_tau[i];
       sin_sig[i] = sin_delta[i]*cos_tau[i] - cos_delta[i]*sin_tau[i];
      }

    for(i=0;i<3;i++)
      for(j=0;j<3;j++)
       {
        r_s[j] = p_s_c[i][j] + cos_sig[i]*s1_s[i][j] + sin_sig[i]*s2_s[i][j];
        r_t[j] = p_t_c[i][j] + cos_tau[i+1]*t1_s[i][j] + sin_tau[i+1]*t2_s[i][j];
        r_n[i][j] = r_s[j]*len_na[i] + r_a[i][j];
        r_c[i][j] = r_t[j]*len_ac[i] + r_a[i][j];
       }

//     ! rotate back atoms by -(sig(1) - sig1_init) around -ex
     sig1 = atan2(sin_sig[0], cos_sig[0]);

     ex_tmp[0] = -ex[0];
     ex_tmp[1] = -ex[1];
     ex_tmp[2] = -ex[2];
     tmp_value = -(sig1-sig1_init)*0.25;
     quaternionL(ex_tmp, tmp_value, p);

     rotation_matrixL(p, Us);

     for(i=0;i<3;i++)
      {
       mat11[i] = r_c[0][i]-r0[i];
       mat22[i] = r_n[1][i]-r0[i];
       mat33[i] = r_a[1][i]-r0[i];
       mat44[i] = r_c[1][i]-r0[i];
       mat55[i] = r_n[2][i]-r0[i];
      }
     matmulL(Us,mat11,mat1);
     matmulL(Us,mat22,mat2);
     matmulL(Us,mat33,mat3);
     matmulL(Us,mat44,mat4);
     matmulL(Us,mat55,mat5);
     for(i=0;i<3;i++)
      {
       r_soln_n[i_soln][0][i] = r_n1[i];
       r_soln_a[i_soln][0][i] = r_a1[i];
       r_soln_c[i_soln][0][i] = mat1[i] + r0[i];
       r_soln_n[i_soln][1][i] = mat2[i] + r0[i];
       r_soln_a[i_soln][1][i] = mat3[i] + r0[i];
       r_soln_c[i_soln][1][i] = mat4[i] + r0[i];
       r_soln_n[i_soln][2][i] = mat5[i] + r0[i];
       r_soln_a[i_soln][2][i] = r_a3[i];
       r_soln_c[i_soln][2][i] = r_c3[i];
      }


     if (print_level > 0)
      {
        printf("roots: t0, t2, t1 %d\n", i_soln);
        printf("%15.6f %15.6f %15.6f\n", half_tan[2], half_tan[1], half_tan[0]);

        for(i=0;i<3;i++)
         {
          rr_a1c1[i] = r_soln_c[i_soln][0][i] - r_soln_a[i_soln][0][i];
          rr_c1n2[i] = r_soln_n[i_soln][1][i] - r_soln_c[i_soln][0][i];
          rr_n2a2[i] = r_soln_a[i_soln][1][i] - r_soln_n[i_soln][1][i];
          rr_a2c2[i] = r_soln_c[i_soln][1][i] - r_soln_a[i_soln][1][i];
          rr_c2n3[i] = r_soln_n[i_soln][2][i] - r_soln_c[i_soln][1][i];
          rr_n3a3[i] = r_soln_a[i_soln][2][i] - r_soln_n[i_soln][2][i];
          rr_a1a2[i] = r_soln_a[i_soln][1][i] - r_soln_a[i_soln][0][i];
          rr_a2a3[i] = r_soln_a[i_soln][2][i] - r_soln_a[i_soln][1][i];
         }

        a1c1 = sqrt(dot_productL(rr_a1c1, rr_a1c1));
        c1n2 = sqrt(dot_productL(rr_c1n2, rr_c1n2));
        n2a2 = sqrt(dot_productL(rr_n2a2, rr_n2a2));
        a2c2 = sqrt(dot_productL(rr_a2c2, rr_a2c2));
        c2n3 = sqrt(dot_productL(rr_c2n3, rr_c2n3));
        n3a3 = sqrt(dot_productL(rr_n3a3, rr_n3a3));
        a1a2 = sqrt(dot_productL(rr_a1a2, rr_a1a2));
        a2a3 = sqrt(dot_productL(rr_a2a3, rr_a2a3));

        printf("na: n2a2, n3a3 = %9.3f%9.3f%9.3f%9.3f\n", len0[2], n2a2, len0[5], n3a3);
        printf("ac: a1c1, a2c2 = %9.3f%9.3f%9.3f%9.3f\n", len0[0], a1c1, len0[3], a2c2);
        printf("cn: c1n2, c2n3 = %9.3f%9.3f%9.3f%9.3f\n", len0[1], c1n2, len0[4], c2n3);
        printf("aa: a1a2, a2a3 = %9.3f%9.3f%9.3f%9.3f\n", len_aa[1], a1a2, len_aa[2], a2a3);

        for(i=0;i<3;i++)
          tmp_array[i] = rr_a1a2[i]/a1a2;
        calc_bnd_ang(b_a1a3, tmp_array, &a3a1a2);

        for(i=0;i<3;i++)
          tmp_array[i] = rr_a2a3[i]/a2a3;
        calc_bnd_ang(tmp_array, b_a1a3, &a2a3a1);

        printf("alpha1, alpha3 = %9.3f%9.3f\n", (pi-a3a1a2)*rad2deg, (pi-a2a3a1)*rad2deg);

        for(i=0;i<3;i++)
          tmp_array[i] = rr_a1c1[i]/a1c1;
        calc_bnd_ang(b_a1n1, tmp_array, &n1a1c1);

        for(i=0;i<3;i++)
          tmp_array[i] = -rr_n3a3[i]/n3a3;
        calc_bnd_ang(b_a3c3, tmp_array, &n3a3c3);

        for(i=0;i<3;i++)
         {
          tmp_array1[i] = rr_a2c2[i]/a2c2;
          tmp_array2[i] = -rr_n2a2[i]/n2a2;
         }
        calc_bnd_ang(tmp_array1, tmp_array2, &n2a2c2);

        printf("ang_nac = %9.3f%9.3f\n", b_ang0[0]*rad2deg, n1a1c1*rad2deg);
        printf("ang_nac = %9.3f%9.3f\n", b_ang0[3]*rad2deg, n2a2c2*rad2deg);
        printf("ang_nac = %9.3f%9.3f\n", b_ang0[6]*rad2deg, n3a3c3*rad2deg);

        for(i=0;i<3;i++)
         {
          tmp_array1[i] = rr_a1c1[i]/a1c1;
          tmp_array2[i] = rr_c1n2[i]/c1n2;
          tmp_array3[i] = rr_n2a2[i]/n2a2;
         }
        calc_dih_ang(tmp_array1, tmp_array2, tmp_array3, &a1c1n2a2);

        for(i=0;i<3;i++)
         {
          tmp_array1[i] = rr_a2c2[i]/a2c2;
          tmp_array2[i] = rr_c2n3[i]/c2n3;
          tmp_array3[i] = rr_n3a3[i]/n3a3;
         }
        calc_dih_ang(tmp_array1, tmp_array2, tmp_array3, &a2c2n3a3);

        printf("t_ang1 = %9.3f%9.3f\n", t_ang0[0]*rad2deg, a1c1n2a2*rad2deg);
        printf("t_ang2 = %9.3f%9.3f\n", t_ang0[1]*rad2deg, a2c2n3a3*rad2deg);
      }

   }

  return;
 }

double calc_t2(double t0)
 {

  double tmp_value;
  double B0, B1, B2, A0, A1, A2, A3, A4, B2_2, B2_3;
  double K0, K1, K2, K3, t0_2, t0_3, t0_4;

  t0_2 = t0*t0;
  t0_3 = t0_2*t0;
  t0_4 = t0_3*t0;

  A0 = Q[0][0] + Q[0][1]*t0 + Q[0][2]*t0_2 + Q[0][3]*t0_3 + Q[0][4]*t0_4;
  A1 = Q[1][0] + Q[1][1]*t0 + Q[1][2]*t0_2 + Q[1][3]*t0_3 + Q[1][4]*t0_4;
  A2 = Q[2][0] + Q[2][1]*t0 + Q[2][2]*t0_2 + Q[2][3]*t0_3 + Q[2][4]*t0_4;
  A3 = Q[3][0] + Q[3][1]*t0 + Q[3][2]*t0_2 + Q[3][3]*t0_3 + Q[3][4]*t0_4;
  A4 = Q[4][0] + Q[4][1]*t0 + Q[4][2]*t0_2 + Q[4][3]*t0_3 + Q[4][4]*t0_4;

  B0 = R[0][0] + R[0][1]*t0 + R[0][2]*t0_2;
  B1 = R[1][0] + R[1][1]*t0 + R[1][2]*t0_2;
  B2 = R[2][0] + R[2][1]*t0 + R[2][2]*t0_2;

  B2_2 = B2*B2;
  B2_3 = B2_2*B2;

  K0 = A2*B2 - A4*B0;
  K1 = A3*B2 - A4*B1;
  K2 = A1*B2_2 - K1*B0;
  K3 = K0*B2 - K1*B1;

  tmp_value = (K3*B0 - A0*B2_3)/(K2*B2 - K3*B1);

  return tmp_value;
 }

double calc_t1(double t0, double t2)
 {

  double tmp_value;
  double U11, U12, U13, U31, U32, U33;
  double t0_2, t2_2;

  t0_2 = t0*t0;
  t2_2 = t2*t2;

  U11 = C0[0][0] + C0[0][1]*t0 + C0[0][2]*t0_2;
  U12 = C1[0][0] + C1[0][1]*t0 + C1[0][2]*t0_2;
  U13 = C2[0][0] + C2[0][1]*t0 + C2[0][2]*t0_2;
  U31 = C0[1][0] + C0[1][1]*t2 + C0[1][2]*t2_2;
  U32 = C1[1][0] + C1[1][1]*t2 + C1[1][2]*t2_2;
  U33 = C2[1][0] + C2[1][1]*t2 + C2[1][2]*t2_2;

  tmp_value = (U31*U13-U11*U33)/(U12*U33-U13*U32);

  return tmp_value;
 }

void calc_dih_ang(double r1[3], double r2[3], double r3[3], double *angle)
 {
//!-----------------------------------------------------------------------
//! r1=Rab, r2=Rbc, r3=Rcd : angle between planes abc and bcd
//!-----------------------------------------------------------------------
  double p[3], q[3], s[3];
  double arg;

  cross_tc(r1, r2, p);
  cross_tc(r2, r3, q);
  cross_tc(r3, r1, s);

  arg = dot_productL(p,q)/sqrt(dot_productL(p,p)*dot_productL(q,q));
  arg = signL(Min(fabs(arg),1.0e0),arg);

  *angle = signL(acos(arg), dot_productL(s,r2));

  return;
 }

void calc_bnd_ang(double r1[3], double r2[3], double *angle)
 {
//!-----------------------------------------------------------------------
//! assume that each vector is normalized
//! r1=Rba, r2=Rbc: angle between r1 and r2
//!-----------------------------------------------------------------------
  double arg;

  arg = dot_productL(r1, r2);
  arg = signL(Min(fabs(arg),1.0e0),arg);

  *angle = acos(arg);

  return;
 }

void cross_tc(double p[3], double q[3], double s[3])
 {
  s[0] = p[1]*q[2] - p[2]*q[1];
  s[1] = p[2]*q[0] - p[0]*q[2];
  s[2] = p[0]*q[1] - p[1]*q[0];

  return;
 }


void quaternionL(double axis[3], double quarter_ang, double p[4])
 {
//!-----------------------------------------------------------------------
//! calculate quaternion, given rotation axis and angle.
//!-----------------------------------------------------------------------
  double tan_w, tan_sqr, tan1, cosine, sine;

  tan_w = tan(quarter_ang);
  tan_sqr = tan_w * tan_w;
  tan1 = 1.0e0 + tan_sqr;
  cosine = (1.0e0 - tan_sqr)/tan1;
  sine = 2.0e0*tan_w/tan1;

  p[0] = cosine;
  p[1] = axis[0] * sine;
  p[2] = axis[1] * sine;
  p[3] = axis[2] * sine;

  return;

 }

void rotation_matrixL(double q[4], double U[3][3])
 {
//!-----------------------------------------------------------------------
//! constructs rotation matrix U from quaternion q.
//!-----------------------------------------------------------------------
  double q0,q1,q2,q3,b0,b1,b2,b3,q00,q01,q02,q03,q11,q12,q13,q22,q23,q33;

  q0 = q[0];
  q1 = q[1];
  q2 = q[2];
  q3 = q[3];

  b0 = 2.0e0*q0;
  b1 = 2.0e0*q1;

  q00 = b0*q0-1.0e0;
  q02 = b0*q2;
  q03 = b0*q3;

  q11 = b1*q1;
  q12 = b1*q2;
  q13 = b1*q3;

  b2 = 2.0e0*q2;
  b3 = 2.0e0*q3;

  q01 = b0*q1;
  q22 = b2*q2;
  q23 = b2*q3;
  q33 = b3*q3;

  U[0][0] = q00+q11;
  U[0][1] = q12-q03;
  U[0][2] = q13+q02;

  U[1][0] = q12+q03;
  U[1][1] = q00+q22;
  U[1][2] = q23-q01;

  U[2][0] = q13-q02;
  U[2][1] = q23+q01;
  U[2][2] = q00+q33;

  return;
 }
//
// END OF FUNCTIONS FROM: loop.h
//






//
// FUNCTIONS FROM: strum.h
//

// set termination criteria for polynomial solver
void initialize_sturm(double *tol_secant, int *max_iter_sturm, int *max_iter_secant)
{
  RELERROR = *tol_secant;
  MAXIT = *max_iter_sturm;
  MAX_ITER_SECANT = *max_iter_secant;
}

void solve_sturm(int p_order, int *n_root, double *poly_coeffs, double *roots)
{
  poly sseq[MAX_ORDER*2];
  double min, max;
  int order, i, j, nroots, nchanges, np, atmin, atmax;

  order = p_order;

  for (i = order; i >= 0; i--)
    {
      sseq[0].coef[i] = poly_coeffs[i];
    }

  if (PRINT_LEVEL > 0)
    {
      for (i = order; i >= 0; i--)
        {
          printf("coefficients in Sturm solver\n");
          printf("%d %lf\n", i, sseq[0].coef[i]);
        }
    }

  /*
   * build the Sturm sequence
   */
  np = buildsturm(order, sseq);

  if (PRINT_LEVEL > 0)
    {
      printf("Sturm sequence for:\n");
      for (i = order; i >= 0; i--)
        printf("%lf ", sseq[0].coef[i]);
      printf("\n\n");
      for (i = 0; i <= np; i++)
        {
          for (j = sseq[i].ord; j >= 0; j--)
            printf("%lf ", sseq[i].coef[j]);
          printf("\n");
        }
      printf("\n");
    }

  /*
   * get the number of real roots
   */

  nroots = numroots(np, sseq, &atmin, &atmax);


  if (nroots == 0)
    {
      //printf("solve: no real roots\n");
      *n_root = nroots;
      return;
    }

  if (PRINT_LEVEL > 0)
    printf("Number of real roots: %d\n", nroots);

  /*
   * calculate the bracket that the roots live in
   */
  min = -1.0;
  nchanges = numchanges(np, sseq, min);

  for (i = 0; nchanges != atmin && i != MAXPOW; i++) {
    min *= 10.0;
    nchanges = numchanges(np, sseq, min);
  }

  if (nchanges != atmin) {
    //printf("solve: unable to bracket all negative roots\n");
    atmin = nchanges;
  }

  max = 1.0;
  nchanges = numchanges(np, sseq, max);
  for (i = 0; nchanges != atmax && i != MAXPOW; i++) {
    max *= 10.0;
    nchanges = numchanges(np, sseq, max);
  }

  if (nchanges != atmax) {
    //printf("solve: unable to bracket all positive roots\n");
    atmax = nchanges;
  }

  nroots = atmin - atmax;

  /*
   * perform the bisection.
   */

  sbisect(np, sseq, min, max, atmin, atmax, roots);

  *n_root = nroots;

  /*
   * write out the roots...
   */
  if (PRINT_LEVEL > 0)
    {
      if (nroots == 1) {
        printf("\n1 distinct real root at x = %f\n", roots[0]);
      } else {
        printf("\n%d distinct real roots for x: \n", nroots);

        for (i = 0; i != nroots; i++)
          {
            printf("%f\n", roots[i]);
          }
      }
    }

  return;
}

/*
 * modp
 *
 *	calculates the modulus of u(x) / v(x) leaving it in r, it
 *  returns 0 if r(x) is a constant.
 *  note: this function assumes the leading coefficient of v
 *	is 1 or -1
 */

double hyper_tan(double a, double x)
{
  double exp_x1, exp_x2, ax;

  ax = a*x;
  if (ax > 100.0)
    {
      return(1.0);
    }
  else if (ax < -100.0)
    {
      return(-1.0);
    }
  else
    {
      exp_x1 = exp(ax);
      exp_x2 = exp(-ax);
      return (exp_x1 - exp_x2)/(exp_x1 + exp_x2);
    }
}

static int modpL(poly *u, poly *v, poly *r)
//	poly *u, *v, *z;
{
        int		k, j;
        double	*nr, *end, *uc;

        nr = r->coef;
        end = &u->coef[u->ord];

        uc = u->coef;
        while (uc <= end)
                        *nr++ = *uc++;

        if (v->coef[v->ord] < 0.0) {


                        for (k = u->ord - v->ord - 1; k >= 0; k -= 2)
                                r->coef[k] = -r->coef[k];

                        for (k = u->ord - v->ord; k >= 0; k--)
                                for (j = v->ord + k - 1; j >= k; j--)
                                        r->coef[j] = -r->coef[j] - r->coef[v->ord + k]
                                        * v->coef[j - k];
        } else {
                        for (k = u->ord - v->ord; k >= 0; k--)
                                for (j = v->ord + k - 1; j >= k; j--)
                                r->coef[j] -= r->coef[v->ord + k] * v->coef[j - k];
        }

        k = v->ord - 1;
        while (k >= 0 && fabs(r->coef[k]) < SMALL_ENOUGH) {
                r->coef[k] = 0.0;
                k--;
        }

        r->ord = (k < 0) ? 0 : k;

        return(r->ord);
}

/*
 * buildsturm
 *
 *	build up a sturm sequence for a polynomial in smat, returning
 * the number of polynomials in the sequence
 */
int buildsturm(int ord, poly *sseq)
//	int	ord;
//	poly	*sseq;
{
        int		i;
        double	f, *fp, *fc;
        poly	*sp;

        sseq[0].ord = ord;
        sseq[1].ord = ord - 1;

        /*
         * calculate the derivative and normalise the leading
         * coefficient.
         */
        f = fabs(sseq[0].coef[ord]*ord);

        fp = sseq[1].coef;
        fc = sseq[0].coef + 1;

        for (i = 1; i <= ord; i++)
          {
            *fp++ = *fc++ * i / f;
          }


        /*
         * construct the rest of the Sturm sequence
         */
        //	for (sp = sseq + 2; modp(sp - 2, sp - 1, sp); sp++) {
        for (sp = sseq + 2; modpL(sp - 2, sp - 1, sp); sp++) {

          /*
           * reverse the sign and normalise
           */

          f = -fabs(sp->coef[sp->ord]);
          for (fp = &sp->coef[sp->ord]; fp >= sp->coef; fp--)
                  *fp /= f;
        }

        sp->coef[0] = -sp->coef[0];	/* reverse the sign */

        return(sp - sseq);
}

/*
 * numroots
 *
 *	return the number of distinct real roots of the polynomial
 * described in sseq.
 */
int numroots(int np, poly *sseq, int *atneg, int *atpos)
//		int		np;
//		poly	*sseq;
//		int		*atneg, *atpos;
{
                int		atposinf, atneginf;
                poly	*s;
                double	f, lf;

                atposinf = atneginf = 0;


        /*
         * changes at positive infinity
         */
        lf = sseq[0].coef[sseq[0].ord];

        for (s = sseq + 1; s <= sseq + np; s++) {
                        f = s->coef[s->ord];
                        if (lf == 0.0 || lf * f < 0)
                                atposinf++;
                lf = f;
        }

        /*
         * changes at negative infinity
         */
        if (sseq[0].ord & 1)
                        lf = -sseq[0].coef[sseq[0].ord];
        else
                        lf = sseq[0].coef[sseq[0].ord];

        for (s = sseq + 1; s <= sseq + np; s++) {
                        if (s->ord & 1)
                                f = -s->coef[s->ord];
                        else
                                f = s->coef[s->ord];
                        if (lf == 0.0 || lf * f < 0)
                                atneginf++;
                        lf = f;
        }

        *atneg = atneginf;
        *atpos = atposinf;
        //	printf("atneginf, atposinf = %d %d\n", atneginf, atposinf);
        return(atneginf - atposinf);
}

/*
 * numchanges
 *
 *	return the number of sign changes in the Sturm sequence in
 * sseq at the value a.
 */
int numchanges(int np, poly *sseq, double a)
//	int		np;
//	poly	*sseq;
//	double	a;

{
        int		changes;
        double	f, lf;
        poly	*s;

        changes = 0;

        lf = evalpoly(sseq[0].ord, sseq[0].coef, a);

        for (s = sseq + 1; s <= sseq + np; s++) {
                        f = evalpoly(s->ord, s->coef, a);
                        if (lf == 0.0 || lf * f < 0)
                                changes++;
                        lf = f;
//			printf("lf %lf %d \n", f, changes);
        }

        //	printf("%d \n", changes);
        return(changes);
}

/*
 * sbisect
 *
 *	uses a bisection based on the sturm sequence for the polynomial
 * described in sseq to isolate intervals in which roots occur,
 * the roots are returned in the roots array in order of magnitude.
 */
void sbisect(int np, poly *sseq, double min, double max, int atmin, int atmax, double *roots)
//	int	np;
//	poly	*sseq;
//	double	min, max;
//	int		atmin, atmax;
//	double	*roots;
{
  double	mid=0.;
  int		n1 = 0, n2 = 0, its, atmid, nroot;

  if ((nroot = atmin - atmax) == 1)
    {

      /*
       * first try a less expensive technique.
       */
      //      printf("min max %lf, %lf \n", min, max);
      if (modrf(sseq->ord, sseq->coef, min, max, &roots[0]))
        return;

      //      printf("try hard way\n");
      /*
       * if we get here we have to evaluate the root the hard
       * way by using the Sturm sequence.
       */
      for (its = 0; its < MAXIT; its++)
        {
          mid = (min + max) / 2;

          atmid = numchanges(np, sseq, mid);

          if (fabs(mid) > RELERROR)
            {
              if (fabs((max - min) / mid) < RELERROR)
                {
                  roots[0] = mid;
                  return;
                }
            }
          else if (fabs(max - min) < RELERROR)
            {
              roots[0] = mid;
              return;
            }

          if ((atmin - atmid) == 0)
            min = mid;
          else
            max = mid;
        }

      if (its == MAXIT)
        {
          //fprintf(stderr, "sbisect: overflow min %f max %f\
                   //                     diff %e nroot %d n1 %d n2 %d\n",
                  //min, max, max - min, nroot, n1, n2);
          roots[0] = mid;
        }
      return;
    }

  /*
   * more than one root in the interval, we have to bisect...
   */

  for (its = 0; its < MAXIT; its++)
    {

      mid = (min + max) / 2;

      atmid = numchanges(np, sseq, mid);

      n1 = atmin - atmid;
      n2 = atmid - atmax;

      if (n1 != 0 && n2 != 0)
        {
          sbisect(np, sseq, min, mid, atmin, atmid, roots);
          sbisect(np, sseq, mid, max, atmid, atmax, &roots[n1]);
          break;
        }

      if (n1 == 0)
        min = mid;
      else
        max = mid;
    }

  if (its == MAXIT)
    {
      /*
      fprintf(stderr, "sbisect: roots too close together\n");
      fprintf(stderr, "sbisect: overflow min %f max %f diff %e\
                                nroot %d n1 %d n2 %d\n",
              min, max, max - min, nroot, n1, n2);
      */
      for (n1 = atmax; n1 < atmin; n1++)
        roots[n1 - atmax] = mid;
    }
}

/*
 * evalpoly
 *
 *	evaluate polynomial defined in coef returning its value.
 */
double evalpoly(int ord, double *coef, double x)
//	int		ord;
//	double	*coef, x;
{
        double	*fp, f;

        fp = &coef[ord];
        f = *fp;

        for (fp--; fp >= coef; fp--)
        f = x * f + *fp;

        return(f);
}


/*
 * modrf
 *
 *	uses the modified regula-falsi method to evaluate the root
 * in interval [a,b] of the polynomial described in coef. The
 * root is returned is returned in *val. The routine returns zero
 * if it can't converge.
 */
int modrf(int ord, double *coef, double	a, double b, double *val)
{
  int its;
  double fa, fb, x, fx, lfx;
  double *fp, *scoef, *ecoef;

  scoef = coef;
  ecoef = &coef[ord];

  fb = fa = *ecoef;
  for (fp = ecoef - 1; fp >= scoef; fp--) {
    fa = a * fa + *fp;
    fb = b * fb + *fp;
  }

  /*
   * if there is no sign difference the method won't work
   */
  if (fa * fb > 0.0)
    return(0);

  /*  commented out to avoid duplicate solutions when the bounds are close to the roots

  if (fabs(fa) < RELERROR)
    {
      *val = a;
      return(1);
    }

  if (fabs(fb) < RELERROR)
    {
      *val = b;
      return(1);
    }
  */

  lfx = fa;

  for (its = 0; its < MAX_ITER_SECANT; its++)
    {
      x = (fb * a - fa * b) / (fb - fa);

      // constrain that x stays in the bounds
      if (x < a || x > b)
        x = 0.5 * (a+b);

      fx = *ecoef;
      for (fp = ecoef - 1; fp >= scoef; fp--)
        fx = x * fx + *fp;

      if (fabs(x) > RELERROR)
        {
          if (fabs(fx / x) < RELERROR)
            {
              *val = x;
              //	      printf(" x, fx %lf %lf\n", x, fx);
              return(1);
            }
        }
      else if (fabs(fx) < RELERROR)
        {
          *val = x;
          //	  printf(" x, fx %lf %lf\n", x, fx);
          return(1);
        }

      if ((fa * fx) < 0)
        {
          b = x;
          fb = fx;
          if ((lfx * fx) > 0)
            fa /= 2;
        }
      else
        {
          a = x;
          fa = fx;
          if ((lfx * fx) > 0)
            fb /= 2;
        }

      lfx = fx;
    }

  //fprintf(stderr, "modrf overflow %f %f %f\n", a, b, fx);

  return(0);
}
//
// END OF FUNCTIONS FROM: strum.h
//

// Coutsia's LoopClosure driver routine for RCD.
float LoopClosure( double **co,  double **anchor_fixed,  double *bond_length, double *valence_angle, double *dihedral_angles, int nr_atoms_loop, int reverse)
{
	//!-----------------------------------------------------------------------
	//! This is a sample driver routine to reconstruct tripeptide loops
	//! from the coordinates of a pdb file using a canonical bond lengths and
	//! angles.
	//!-----------------------------------------------------------------------
	//  use in_out
	//  use tripep_closure
	//!-----------------------------------------------------------------------
	bool debug = false;
	char out_pdb[100];
	int n_soln, i, j, k;
	double r_soln_n[MAX_SOLN] [3] [3], r_soln_a[MAX_SOLN] [3] [3], r_soln_c[MAX_SOLN] [3] [3];
	double rmsd, sum, dr[3];
	double r0_n[3] [3], r0_a[3] [3], r0_c[3] [3];
	double r_sols_n[MAX_SOLN] [5] [3], r_sols_a[MAX_SOLN] [5] [3], r_sols_c[MAX_SOLN] [5] [3];

	int calc_rmsd, write_out_pdb;
	double b_len[6], b_ang[7], t_ang[2];

	// Para ordenar las soluciones segun su rmsd de mejor a peor
	int orden[16];
	double rmsds[16];
	int encontrado;

	double posespinor[4];
	double coordX[5][3],coordY[3][3];
	double *coord[5],*coord2[3];
	// "Chapa" to use some functions directly...
	for(int i=0; i<5; i++)
		coord[i] = &coordX[i][0];
	for(int i=0; i<3; i++)
		coord2[i] = &coordY[i][0];

	//TRANSLATE AND COPY COORDINATES TO REVERSE ARRAY
	for (int z = 0; z < 5; z++) // copy the coordinates of the last 5 atoms in reverse order.
		for (int x = 0; x < 3; x++)
			coord[z][x] = co[nr_atoms_loop-1-z][x] - co[nr_atoms_loop-1][x] + anchor_fixed[2][x];

	if(debug)
	{
		printf("Anchor:\n");
		for(int i=0; i<3; i++)
			printf (" %f %f %f\n", anchor_fixed[i][0], anchor_fixed[i][1], anchor_fixed[i][2]);

		printf("Coord:\n");
		for(int i=0; i<5; i++)
			printf (" %f %f %f\n", coord[i][0], coord[i][1], coord[i][2]);
	}

	//ROTATE REVERSE ANCHOR IN CORRECT PLACE
	for (int z = 0; z < 3; z++)
		for (int x = 0; x < 3; x++)
			coord2[z][x] = anchor_fixed[2 - z][x];

	pose(coord,coord2,posespinor);
	quatrotate_cg(coord,NULL,NULL,5,posespinor,0);

//	if(debug)
//	{
//		printf("Coord:\n");
//		for(int i=0; i<5; i++)
//			printf (" %f %f %f\n", coord[i][0], coord[i][1], coord[i][2]);
//	}
//	exit(0);

	//  !! initialize: set up input parameters
	//  ! bond lengths
	double r_n[4][3], r_a[4][3], r_c[4][3];
    // [0] anchor; [1][2] moving aa; [3] anchor; anchor [4] for compute dihedrals
	int nres=nr_atoms_loop/3;
	int indb;

	if (reverse==0)
	{
		// Forward
		if(debug)
			printf (" Forward\n");

		for (int j = 0; j < 3; j++)
			for (int i = 0; i < 3; i++)
			{
				// Get coordinates of terminal atoms (residues N-4 to N-2)
				r_n[j][i]=co[nr_atoms_loop-3*(4-j)-3][i];
				r_a[j][i]=co[nr_atoms_loop-3*(4-j)-2][i];
				r_c[j][i]=co[nr_atoms_loop-3*(4-j)-1][i];
			}
		// Get coordinates of terminal residue (N), i.e. the terminal anchor.
		// anchor terminal (fix, from PDB) (Nt or Ct, depending on "reverse")
		for (int i = 0; i < 3; i++)
		{
			r_c[3][i]=coord[3][i];
			r_a[3][i]=coord[4][i];
			r_n[3][i]=co[nr_atoms_loop-6][i]; // for visualization only
		}
//		for(int i=0; i<5; i++)
//			printf (" %f %f %f\n", coord[i][0], coord[i][1], coord[i][2]);

//		if(debug)
//		{
//			printf (" %f %f %f\n", r_n[3][0], r_n[3][1], r_n[3][2]);
//			printf (" %f %f %f\n", r_a[3][0], r_a[3][1], r_a[3][2]);
//			printf (" %f %f %f\n\n", r_c[3][0], r_c[3][1], r_c[3][2]);
//
//			printf (" %f %f %f\n", co[nr_atoms_loop-3][0], co[nr_atoms_loop-3][1], co[nr_atoms_loop-3][2]);
//			printf (" %f %f %f\n", co[nr_atoms_loop-2][0], co[nr_atoms_loop-2][1], co[nr_atoms_loop-2][2]);
//			printf (" %f %f %f\n\n", co[nr_atoms_loop-1][0], co[nr_atoms_loop-1][1], co[nr_atoms_loop-1][2]);
//		}

		indb=nr_atoms_loop-10;
//		// bond lenghts
//		b_len[0] = bond_length[indb];
//		b_len[1] = bond_length[indb+1];
//		b_len[2] = bond_length[indb+2];
//		b_len[3] = bond_length[indb+3];
//		b_len[4] = bond_length[indb+4];
//		b_len[5] = bond_length[indb+5];
//		// angles
//		b_ang[0] = valence_angle[indb];
//		b_ang[1] = valence_angle[indb+1];
//		b_ang[2] = valence_angle[indb+2];
//		b_ang[3] = valence_angle[indb+3];
//		b_ang[4] = valence_angle[indb+4];
//		b_ang[5] = valence_angle[indb+5];
//		b_ang[6] = valence_angle[indb+6];
		//  peptide torsion angles (omegas)
		t_ang[0] = dihedral_angles[nr_atoms_loop-8];
		t_ang[1] = dihedral_angles[nr_atoms_loop-5];
	}
	else // backward
	{
		for (int j = 2; j < 4; j++)
			for (int i = 0; i < 3; i++)
			{
				r_c[j][i]=co[nr_atoms_loop-3*j-3][i];
				r_a[j][i]=co[nr_atoms_loop-3*j-2][i];
				r_n[j][i]=co[nr_atoms_loop-3*j-1][i];
			}

		// anchor (residue N)
		for (int i = 0; i < 3; i++)
		{
			r_c[0][i]=anchor_fixed[0][i];
			r_a[0][i]=anchor_fixed[1][i];
			r_n[0][i]=anchor_fixed[2][i];
			r_n[1][i]=coord[3][i];
			r_a[1][i]=coord[4][i];
			r_c[1][i]=co[nr_atoms_loop-6][i];
		}

//		if(debug)
//		{
//			printf (" %f %f %f\n", r_n[0][0], r_n[0][1], r_n[0][2]);
//			printf (" %f %f %f\n", r_a[0][0], r_a[0][1], r_a[0][2]);
//			printf (" %f %f %f\n\n", r_c[0][0], r_c[0][1], r_c[0][2]);
//
//			printf (" %f %f %f\n", co[nr_atoms_loop-1][0], co[nr_atoms_loop-1][1], co[nr_atoms_loop-1][2]);
//			printf (" %f %f %f\n", co[nr_atoms_loop-2][0], co[nr_atoms_loop-2][1], co[nr_atoms_loop-2][2]);
//			printf (" %f %f %f\n", co[nr_atoms_loop-3][0], co[nr_atoms_loop-3][1], co[nr_atoms_loop-3][2]);
//		}

		indb=5;
//        // bond lengths
//		b_len[0] = bond_length[indb];
//		b_len[1] = bond_length[indb+1];
//	    b_len[2] = bond_length[indb+2];
//		b_len[3] = bond_length[indb+3];
//		b_len[4] = bond_length[indb+4];
//		b_len[5] = bond_length[indb+5];
//		// angles
//		b_ang[0] = valence_angle[indb];
//		b_ang[1] = valence_angle[indb+1];
//		b_ang[2] = valence_angle[indb+2];
//		b_ang[3] = valence_angle[indb+3];
//		b_ang[4] = valence_angle[indb+4];
//		b_ang[5] = valence_angle[indb+5];
//		b_ang[6] = valence_angle[indb+6];
		//  omega angles
//		t_ang[0] = dihedral_angles[4];
//		t_ang[1] = dihedral_angles[7];
		t_ang[0] = dihedral_angles[nr_atoms_loop-5];
		t_ang[1] = dihedral_angles[nr_atoms_loop-8];
	}

//	if(debug)
//		for (int j = 0; j < 6; j++)
//			printf("bond_lengh %8.6lf\n", b_len[j]   );
//
//	b_len[0] = bond_length2( r_a[0], r_c[0] );
//	b_len[1] = bond_length2( r_c[0], r_n[1] );
//	b_len[2] = bond_length2( r_n[1], r_a[1] );
//	b_len[3] = bond_length2( r_a[1], r_c[1] );
//	b_len[4] = bond_length2( r_c[1], r_n[2] );
//	b_len[5] = bond_length2( r_n[2], r_a[2] );

//	if(debug)
//	{
//		for (int j = 0; j < 6; j++)
//			printf("bond_lengh %8.6lf\n", b_len[j]   );
//		for (int j = 0; j < nr_atoms_loop/3; j++)
//			printf("bond_lengh j %d CON %f NCA %f CACO %f\n", j, bond_length[3*j], bond_length[3*j+1], bond_length[3*j+2] );
//	}

//    // ! bond angles
//	if(debug)
//		for (int j = 0; j < 7; j++)
//			printf("angles bk %8.6lf\n", b_ang[j] );

//	angle_bk( r_a[0], r_n[0], r_a[0], r_c[0], b_ang[0] );
//	angle_bk( r_c[0], r_a[0], r_c[0], r_n[1], b_ang[1] );
//	angle_bk( r_n[1], r_c[0], r_n[1], r_a[1], b_ang[2] );
//	angle_bk( r_a[1], r_n[1], r_a[1], r_c[1], b_ang[3] );
//	angle_bk( r_c[1], r_a[1], r_c[1], r_n[2], b_ang[4] );
//	angle_bk( r_n[2], r_c[1], r_n[2], r_a[2], b_ang[5] );
//	angle_bk( r_a[2], r_n[2], r_a[2], r_c[2], b_ang[6] );

//	if(debug)
//	{
//		for (int j = 0; j < 7; j++)
//			printf("angles bk %8.6lf\n", b_ang[j] );
//		for (int j = 0; j < nr_atoms_loop/3; j++)
//			printf("valence j %d CACON %f CONCA %f NCACO %f\n", j, valence_angle[3*j],valence_angle[3*j+1],valence_angle[3*j+2] );
//		printf("omega %8.6lf %8.6lf\n", t_ang[0],t_ang[1]);
//		for (int j = 0; j < nr_atoms_loop; j++)
//			printf("j= %d  Dihedral %8.6lf\n", j, dihedral_angles[j]);
//	}

	//for (int i = 0; i < 2; i++)
	//if (t_ang[i]<0)  t_ang[i]=2*M_PI+t_ang[i];


	//! connectivity of atoms:
	//!   N1-A1-C1-N2-A2-C2-N3-A3-C3
	//!
	//! input:
	//!
	//!  * b_len(1:6): bond lengths (A1-C1, C1-N2, ..., N3-A3)
	//!  * b_ang(1:7): bond angles (N1-A1-C1, A1-C1-N2, ..., N3-A3-C3)
	//!  * t_ang(1:2): torsion angles (A1-C1-N2-A2, A2-C2-N3-A3)

	//	initialize_loop_closure( b_len, b_ang, t_ang );
	initialize_loop_closure( bond_length+indb, valence_angle+indb, t_ang );

	solve_3pep_poly( r_n[1], r_a[1], r_a[3], r_c[3], r_soln_n, r_soln_a, r_soln_c, &n_soln ); //   ! call tripeptide loop closure routine

//	double dihedral_lc[6];
//	// Computes the dihedral angle around bond between at2 and at3
//	dihedral_lc[0] = get_dihedral(r_c[0], r_n[1], r_a[1], r_c[1]); // Phi_1
//	dihedral_lc[1] = get_dihedral(r_n[1], r_a[1], r_c[1], r_n[2]); // Psi_1
//	dihedral_lc[2] = get_dihedral(r_c[1], r_n[2], r_a[2], r_c[2]); // Phi_2
//	dihedral_lc[3] = get_dihedral(r_n[2], r_a[2], r_c[2], r_n[3]); // Psi_2
//	dihedral_lc[4] = get_dihedral(r_c[2], r_n[3], r_a[3], r_c[3]); // Phi_3
//	dihedral_lc[5] = get_dihedral(r_n[3], r_a[3], r_c[3], r_n[4]); // Psi_3

//	for(int i=0; i<6; i++)
//		printf("i= %d  dh= %lf\n", i, dihedral_lc[i]);

	if(debug)
		printf("No. of solutions: %d\n", n_soln);

	float min=9999;
	int k_min=0;
	for ( k = 0; k < n_soln; k++ )
	{
		// Copiar las coordenadas de los cinco atomos de las soluciones
		// para despues poder calcular los angulos dihedros
		for ( j = 0; j < 3; j++ )
		{
			r_sols_n[k] [0] [j] = r_n[0] [j];
			r_sols_a[k] [0] [j] = r_a[0] [j];
			r_sols_c[k] [0] [j] = r_c[0] [j];
		}

		for ( i = 0; i < 3; i++ )
			for ( j = 0; j < 3; j++ )
			{
				r_sols_n[k] [i + 1] [j] = r_soln_n[k] [i] [j];
				r_sols_a[k] [i + 1] [j] = r_soln_a[k] [i] [j];
				r_sols_c[k] [i + 1] [j] = r_soln_c[k] [i] [j];
			}

		if(reverse==0)
			for ( j = 0; j < 3; j++ )
			{
				r_sols_n[k][4][j] = anchor_fixed[0][j];
				r_sols_a[k][4][j] = anchor_fixed[1][j];
				r_sols_c[k][4][j] = anchor_fixed[2][j];
			}
		else
			for ( j = 0; j < 3; j++ )
			{
//				r_sols_n[k][4][j] = anchor_fixed[2][j];
//				r_sols_a[k][4][j] = anchor_fixed[1][j];
//				r_sols_c[k][4][j] = anchor_fixed[0][j];
				r_sols_n[k][4][j] = co[nr_atoms_loop-13][j];
				r_sols_a[k][4][j] = co[nr_atoms_loop-14][j];
				r_sols_c[k][4][j] = co[nr_atoms_loop-15][j];
			}

		// ! rmsd
		sum = 0.0;

//		for ( i = 1; i < 3; i++ )
//		{
//			for ( j = 0; j < 3; j++ )
//				dr[j] = r_sols_n[k] [i] [j] - r_n[i] [j];
//			sum += dot_productL( dr, dr );
//			for ( j = 0; j < 3; j++ )
//				dr[j] = r_sols_a[k] [i] [j] - r_a[i] [j];
//			sum += dot_productL( dr, dr );
//			for ( j = 0; j < 3; j++ )
//				dr[j] = r_sols_c[k] [i] [j] - r_c[i] [j];
//			sum += dot_productL( dr, dr );
//		}

		// Compute RMSD with the moving atoms (only 5!)
		if(reverse==0)
			for(j=0; j<3; j++)
			{
				sum += pow(r_sols_c[k][1][j] - co[nr_atoms_loop-3*(4-1)-1][j], 2);
				sum += pow(r_sols_n[k][2][j] - co[nr_atoms_loop-3*(4-2)-3][j], 2);
				sum += pow(r_sols_a[k][2][j] - co[nr_atoms_loop-3*(4-2)-2][j], 2);
				sum += pow(r_sols_c[k][2][j] - co[nr_atoms_loop-3*(4-2)-1][j], 2);
				sum += pow(r_sols_n[k][3][j] - co[nr_atoms_loop-3*(4-3)-3][j], 2);
			}
		else
			for(j=0; j<3; j++)
			{
				sum += pow(r_sols_c[k][1][j] - co[nr_atoms_loop-3*1-3][j], 2);
				sum += pow(r_sols_n[k][2][j] - co[nr_atoms_loop-3*2-1][j], 2);
				sum += pow(r_sols_a[k][2][j] - co[nr_atoms_loop-3*2-2][j], 2);
				sum += pow(r_sols_c[k][2][j] - co[nr_atoms_loop-3*2-3][j], 2);
				sum += pow(r_sols_n[k][3][j] - co[nr_atoms_loop-3*3-1][j], 2);
			}

		// rmsd = sqrt( sum / 5.0 );
		rmsd = sum / 5.0;
		encontrado = 0;
		if(debug)
			printf("Msd for solution #%d is %.3f\n", k+1, rmsd);
		if(rmsd < min)
		{
			min = rmsd;
			k_min = k;
		}
	}

//    float minimo;
//    if (n_soln>0)
//    	minimo=rmsds[orden[0]];
//    else
//    	minimo=9999;

	if(debug && n_soln > 0)
	{
		FILE *fpdb;
		fpdb=fopen("loopcl.pdb", "w");
		j=0;
		fprintf(fpdb, "MODEL %d\n",1);
		k=0;
		// RDC's loop
		if(reverse==0)
			for(int i = 0; i < 4; i++)
			{
				fprintf(fpdb, "ATOM   %4d  N   GLY  %4d    %8.3lf%8.3lf%8.3lf\n",
						k++, i, co[nr_atoms_loop-3*(4-i)-3][0], co[nr_atoms_loop-3*(4-i)-3][1], co[nr_atoms_loop-3*(4-i)-3][2]);

				fprintf(fpdb, "ATOM   %4d  CA  GLY  %4d    %8.3lf%8.3lf%8.3lf\n",
						k++,  i, co[nr_atoms_loop-3*(4-i)-2][0], co[nr_atoms_loop-3*(4-i)-2][1], co[nr_atoms_loop-3*(4-i)-2][2]);

				fprintf(fpdb, "ATOM   %4d  C   GLY  %4d    %8.3lf%8.3lf%8.3lf\n",
						k++,  i, co[nr_atoms_loop-3*(4-i)-1][0], co[nr_atoms_loop-3*(4-i)-1][1], co[nr_atoms_loop-3*(4-i)-1][2]);
			}
		else
			for(int i = 0; i < 4; i++)
			{
				fprintf(fpdb, "ATOM   %4d  N   GLY  %4d    %8.3lf%8.3lf%8.3lf\n",
						k++, i, co[nr_atoms_loop-3*i-1][0], co[nr_atoms_loop-3*i-1][1], co[nr_atoms_loop-3*i-1][2]);

				fprintf(fpdb, "ATOM   %4d  CA  GLY  %4d    %8.3lf%8.3lf%8.3lf\n",
						k++,  i, co[nr_atoms_loop-3*i-2][0], co[nr_atoms_loop-3*i-2][1], co[nr_atoms_loop-3*i-2][2]);

				fprintf(fpdb, "ATOM   %4d  C   GLY  %4d    %8.3lf%8.3lf%8.3lf\n",
						k++,  i, co[nr_atoms_loop-3*i-3][0], co[nr_atoms_loop-3*i-3][1], co[nr_atoms_loop-3*i-3][2]);
			}

		fprintf(fpdb, "ENDMDL\n");
		// Pre-Closure Loop
		fprintf(fpdb, "MODEL %d\n",2);
		for(int i = 0; i < 4; i++)
		{
			fprintf(fpdb, "ATOM   %4d  N   GLY  %4d    %8.3lf%8.3lf%8.3lf\n",
					k++, i, r_n[i][0], r_n[i][1], r_n[i][2]);

			fprintf(fpdb, "ATOM   %4d  CA  GLY  %4d    %8.3lf%8.3lf%8.3lf\n",
					k++,  i, r_a[i][0], r_a[i][1], r_a[i][2]);

			fprintf(fpdb, "ATOM   %4d  C   GLY  %4d    %8.3lf%8.3lf%8.3lf\n",
					k++,  i, r_c[i][0], r_c[i][1], r_c[i][2]);
		}
		fprintf(fpdb, "ENDMDL\n");
		// Loop-Closure Loops (Coutsias')
		for (int j = 0; j < n_soln; j++)
		{
			fprintf(fpdb, "MODEL %d\n",j+3);
			int k=0;
			for(int i = 0; i < 5; i++)
			{

				fprintf(fpdb, "ATOM   %4d  N   GLY  %4d    %8.3lf%8.3lf%8.3lf\n",
						k++, i, r_sols_n[j] [i] [0], r_sols_n[j] [i] [1], r_sols_n[j] [i] [2]);

				fprintf(fpdb, "ATOM   %4d  CA  GLY  %4d    %8.3lf%8.3lf%8.3lf\n",
						k++, i, r_sols_a[j] [i] [0], r_sols_a[j] [i] [1], r_sols_a[j] [i] [2]);

				fprintf(fpdb, "ATOM   %4d  C   GLY  %4d    %8.3lf%8.3lf%8.3lf\n",
						k++, i, r_sols_c[j] [i] [0], r_sols_c[j] [i] [1], r_sols_c[j] [i] [2]);

			}
			fprintf(fpdb, "ENDMDL\n");
		}
		fclose(fpdb);
	}


	// Warning: Loops must be larger than 4 residues!!!
	if(n_soln > 0)
		if(reverse==0) // Forward --> Ct end
		{
			dihedral_angles[nr_atoms_loop-10] = get_dihedral(r_sols_c[k_min][0], r_sols_n[k_min][1], r_sols_a[k_min][1], r_sols_c[k_min][1]); // Phi_1
			dihedral_angles[nr_atoms_loop-9] = get_dihedral(r_sols_n[k_min][1], r_sols_a[k_min][1], r_sols_c[k_min][1], r_sols_n[k_min][2]); // Psi_1
			correctdihedral(&dihedral_angles[nr_atoms_loop-10]);
			correctdihedral(&dihedral_angles[nr_atoms_loop-9]);

			dihedral_angles[nr_atoms_loop-7] = get_dihedral(r_sols_c[k_min][1], r_sols_n[k_min][2], r_sols_a[k_min][2], r_sols_c[k_min][2]); // Phi_2
			dihedral_angles[nr_atoms_loop-6] = get_dihedral(r_sols_n[k_min][2], r_sols_a[k_min][2], r_sols_c[k_min][2], r_sols_n[k_min][3]); // Psi_2
			correctdihedral(&dihedral_angles[nr_atoms_loop-7]);
			correctdihedral(&dihedral_angles[nr_atoms_loop-6]);

			dihedral_angles[nr_atoms_loop-4] = get_dihedral(r_sols_c[k_min][2], r_sols_n[k_min][3], r_sols_a[k_min][3], r_sols_c[k_min][3]); // Phi_3
			dihedral_angles[nr_atoms_loop-3] = get_dihedral(r_sols_n[k_min][3], r_sols_a[k_min][3], r_sols_c[k_min][3], r_sols_n[k_min][4]); // Psi_3
			correctdihedral(&dihedral_angles[nr_atoms_loop-4]);
			correctdihedral(&dihedral_angles[nr_atoms_loop-3]);
		}
		else // Backward --> Nt end
		{
			dihedral_angles[nr_atoms_loop-3] = get_dihedral(r_sols_c[k_min][0], r_sols_n[k_min][1], r_sols_a[k_min][1], r_sols_c[k_min][1]); // Phi_1
			dihedral_angles[nr_atoms_loop-4] = get_dihedral(r_sols_n[k_min][1], r_sols_a[k_min][1], r_sols_c[k_min][1], r_sols_n[k_min][2]); // Psi_1
			correctdihedral(&dihedral_angles[nr_atoms_loop-3]);
			correctdihedral(&dihedral_angles[nr_atoms_loop-4]);

			dihedral_angles[nr_atoms_loop-6] = get_dihedral(r_sols_c[k_min][1], r_sols_n[k_min][2], r_sols_a[k_min][2], r_sols_c[k_min][2]); // Phi_2
			dihedral_angles[nr_atoms_loop-7] = get_dihedral(r_sols_n[k_min][2], r_sols_a[k_min][2], r_sols_c[k_min][2], r_sols_n[k_min][3]); // Psi_2
			correctdihedral(&dihedral_angles[nr_atoms_loop-6]);
			correctdihedral(&dihedral_angles[nr_atoms_loop-7]);

			dihedral_angles[nr_atoms_loop-9] = get_dihedral(r_sols_c[k_min][2], r_sols_n[k_min][3], r_sols_a[k_min][3], r_sols_c[k_min][3]); // Phi_3
			dihedral_angles[nr_atoms_loop-10] = get_dihedral(r_sols_n[k_min][3], r_sols_a[k_min][3], r_sols_c[k_min][3], r_sols_n[k_min][4]); // Psi_3
			correctdihedral(&dihedral_angles[nr_atoms_loop-9]);
			correctdihedral(&dihedral_angles[nr_atoms_loop-10]);
		}

	// exit(0);
	return min;
}

// Mon: this code was cut-pasted from the main code of RCD (it should be a guide for re-implementation ;-)
////
//// ANALYTIC LOOP-CLOSURE (Coutsias et al.)
////
//// Buffering anchor
//double *anchor_buf[3];
//anchor_buf[0] = (double *)malloc(sizeof(double)*3);
//anchor_buf[1] = (double *)malloc(sizeof(double)*3);
//anchor_buf[2] = (double *)malloc(sizeof(double)*3);
//// Tri-peptide loop closure
//if( enable_lc && rmsd < rmsd_crit )
//{
//	n_almost_closed++;
//
//	// fprintf(stderr,"BEFORE k= %d  rmsd= %f\n",k,rmsd);
//
//	// Checking intra-chain clashes here...
//	if( intra_clash_cg(co->el,p1co->el,factor_co,factor_ocb,nr_atoms_loop,atomshift,marker2,reverse,cg_mode) == 0 )
//	{
//		for(int i=0; i<nr_atoms_loop; i++)
//			dihedral_angle_aux[i] = dihedral_angle[i]; // Backup dihedral angles
//
////						dump_dihedrals("Before LoopClosure:",dihedral_angle,nr_atoms_loop);
////						multiPDB_cg(co->el,p1co->el,1,nr_loop+2,"test_lc0.pdb",secu,1,0);
////						if(check_angles(binmaps, dihedral_angle, size, nr_atoms_loop, reverse, false, marker) != 0)
////						{
////							fprintf(stderr,"Check_angles before found!\n");
////						}
//
//		// ***************************
//		// Mon's shake:
//		// ***************************
//		check2 = 1;
//		float shake = 0.02;
//		int n_current_lc = 0;
//		int n_max = 20;
//		int lc_bad_rama = 0;
//		int n_bad_lc = 0;
//		while(check2 != 0 && n_current_lc < n_max)
//		{
//			n_LoopClosure++;
//			n_current_lc++;
//
//			for(int k=0; k<3; k++)
//			{
//				anchor_buf[0][k] = anchor_end->el[0][k] + shake * (rg->Random()-0.5);
//				anchor_buf[1][k] = anchor_end->el[1][k] + shake * (rg->Random()-0.5);
//				anchor_buf[2][k] = anchor_end->el[2][k] + shake * (rg->Random()-0.5);
//			}
//
//			float min = LoopClosure( co->el, anchor_buf , bond_length, valence_angle, dihedral_angle, nr_atoms_loop, reverse);
//			printf("LC: k=%d n_lc= %d  rmsd=%lf  total_counter=%ld  min=%f  rama_fail= %d  rand_repairs=%d  rmsd_count= %d\n", k, n_current_lc, rmsd, total_counter, min, n_rama_fail, n_rand_repairs, rmsd_count);
//
//			if(min < 1000) // If valid solutions
//			{
//				check2 = check_angles(binmaps, dihedral_angle, size, nr_atoms_loop, reverse, false, marker);
//				if(check2 != 0) // if bad angles
//				{
//					for(int i=0; i<nr_atoms_loop; i++)
//						dihedral_angle[i] = dihedral_angle_aux[i]; // Restore dihedral angles (not neccessary...)
//					lc_bad_rama++;
//				}
//				else // Good angles, then update coordinates and leave...
//				{
//					if(reverse==0) // Forward
//					{
//						isafrotatetreef_cg(co->el,p1co->el,nr_atoms_loop,bond_length,dihedral_angle,prearray3,prearray4,prearray5,prearray6,prearray1,prearray2,full,cg_mode);
//					}
//					else // Backward
//					{
//						//FORWARD COORDINATES ARE NOW CREATED
//						for (int z = 1; z < nr_atoms_loop-2; z++) // Initialize reverse dihedrals
//							dihedral_angle2[z+2]	=	dihedral_angle[nr_atoms_loop-z];
//						dihedral_angle2[2] = FirstPhi; // forward
//						dihedral_angle2[0] = LastPsi; // forward
//
//						co_dummy->copy(0,0,&co_i->el[0],nr_atoms_loop,3);	// Needed the first 3 coordinates in co matrix (forward style)
//						isafrotatetreef_cg(co_dummy->el,p1co_dummy->el,nr_atoms_loop,bond_length,dihedral_angle2,prearray3,prearray4,prearray5,prearray6,prearray1,prearray2,full,cg_mode);
//
//						// Convert "FORWARD" coordinates into "REVERSE" in "co" and "p1co"
//						// Note: Now re-alignment it is not needed because all ICs are OK perfect !!!
//						for(int i = 0; i < nr_atoms_loop ; i++)																	/*Store randomized loop conformation*/
//							for(int j = 0; j < 3; j++)
//							{
//								co->el[i][j] =  co_dummy->el[nr_atoms_loop-1-i][j];
//								if(cg_mode>0)
//									p1co->el[i][j] =  p1co_dummy->el[nr_atoms_loop-1-i][j];						/*store parallel side chain 1*/
//							}
//					}
//
//					rmsd = rmsd_anchor3(co->el+last,anchor_end->el);
//					//fprintf(stderr,"After LoopClosure RMSD: %f\n",rmsd);
//					rmsd_count=0;
//					//fprintf(stderr,"I successfully closed the loop! k=%d\n",k);
//				}
//			}
//			else
//				n_bad_lc++;
//			// Additional check for Analytical Loop-Closure...
//		}
//		fprintf(stderr,"k= %d  n_current_lc= %d  rmsd= %f  rmsd2= %f  rmsd3= %f\n",k,n_current_lc,rmsd,rmsd_anchor3(co->el+last,anchor_end->el),rmsd_anchor3(co->el+last,anchor_buf));
//		if(n_current_lc == n_max)
//		{
//			total_counter = max_path; // Otherwise, it does save as valid solution!
//			rmsd = rmsd_crit + 1e-6; // to continue in RCD loop <-- SUPER SLOW !!!
//			fprintf(stderr,"BAD: (k=%d) n_current_lc= %d  lc_bad_rama= %d  n_bad_lc= %d\n",k,n_current_lc,lc_bad_rama,n_bad_lc);
//		}
//		else
//			fprintf(stderr,"GOOD: (k=%d) n_current_lc= %d  lc_bad_rama= %d  n_bad_lc= %d\n",k,n_current_lc,lc_bad_rama,n_bad_lc);
//		// End Mon's shake
//		// ***************************
//
//
//
////							float min = LoopClosure( co->el, anchor_end->el , bond_length, valence_angle, dihedral_angle, nr_atoms_loop, reverse);
////							printf("LC: k=%d rmsd=%lf total_counter=%ld min=%f  rama_fail= %d  rand_repairs=%d  rmsd_count= %d\n", k, rmsd, total_counter, min, n_rama_fail, n_rand_repairs, rmsd_count);
////							// if(min < 0.005)
////							if(check_angles(binmaps, dihedral_angle, size, nr_atoms_loop, reverse, false, marker) != 0)
////							{
////								isafrotatetreef_cg(co->el,p1co->el,nr_atoms_loop,bond_length,dihedral_angle,prearray3,prearray4,prearray5,prearray6,prearray1,prearray2,full,cg_mode);
////								multiPDB_cg(co->el,p1co->el,1,nr_loop+2,"test_lc.pdb",secu,1,0);
////								dump_dihedrals("New LoopClosure:",dihedral_angle,nr_atoms_loop);
////								fprintf(stderr,"Check_angles after= %d\n",check_angles(binmaps, dihedral_angle, size, nr_atoms_loop, reverse, false, marker));
////								paint_maps_points(binmaps, size, nr_loop+2, "lc_maps", dihedral_angle, reverse);
////								exit(0);
////							}
////							if(min < 1000 && check_angles(binmaps, dihedral_angle, size, nr_atoms_loop, reverse, false, marker) != 0)
////							{
////								for(int i=0; i<nr_atoms_loop; i++)
////									dihedral_angle[i] = dihedral_angle_aux[i]; // Restore dihedral angles
////								rmsd_count++;
////								//							n_rama_fail++;
////								rmsd = rmsd_crit + 1e-6; // to continue in RCD loop
////							}
////							else // If good dihedrals...
////							{
////								if(min > lc_msd) // Minimum MSD is above threshold (Bad), then continue in RCD loop.
////								{
////									rmsd = rmsd_crit + 1e-6; // to continue in RCD loop
////									rmsd_count++;
////									for(int i=0; i<nr_atoms_loop; i++)
////										dihedral_angle[i] = dihedral_angle_aux[i]; // Restore dihedral angles
////									n_bad_closures++;
////								}
////								else // Good analytic loop-closure
////								{
////									n_good_closures++;
////									if(reverse==0) // Forward
////									{
////										isafrotatetreef_cg(co->el,p1co->el,nr_atoms_loop,bond_length,dihedral_angle,prearray3,prearray4,prearray5,prearray6,prearray1,prearray2,full,cg_mode);
////									}
////									else // Backward
////									{
////										//FORWARD COORDINATES ARE NOW CREATED
////										for (int z = 1; z < nr_atoms_loop-2; z++) // Initialize reverse dihedrals
////											dihedral_angle2[z+2]	=	dihedral_angle[nr_atoms_loop-z];
////										dihedral_angle2[2] = FirstPhi; // forward
////										dihedral_angle2[0] = LastPsi; // forward
////
////										co_dummy->copy(0,0,&co_i->el[0],nr_atoms_loop,3);	// Needed the first 3 coordinates in co matrix (forward style)
////										isafrotatetreef_cg(co_dummy->el,p1co_dummy->el,nr_atoms_loop,bond_length,dihedral_angle2,prearray3,prearray4,prearray5,prearray6,prearray1,prearray2,full,cg_mode);
////
////										// Convert "FORWARD" coordinates into "REVERSE" in "co" and "p1co"
////										// Note: Now re-alignment it is not needed because all ICs are OK perfect !!!
////										for(int i = 0; i < nr_atoms_loop ; i++)																	/*Store randomized loop conformation*/
////											for(int j = 0; j < 3; j++)
////											{
////												co->el[i][j] =  co_dummy->el[nr_atoms_loop-1-i][j];
////												if(cg_mode>0)
////													p1co->el[i][j] =  p1co_dummy->el[nr_atoms_loop-1-i][j];						/*store parallel side chain 1*/
////											}
////									}
////
////									rmsd = rmsd_anchor3(co->el+last,anchor_end->el);
////									//fprintf(stderr,"After LoopClosure RMSD: %f\n",rmsd);
////									rmsd_count=0;
////									//fprintf(stderr,"I successfully closed the loop! k=%d\n",k);
////								}
////							}
////							// Additional check for Analytical Loop-Closure...
////							if(check2==0 && enable_lc && rmsd < rmsd_crit )
////							{
////								if(check2 == 0)
////									check2 = clash_cg(vol,co->el,p1co->el,nr_atoms_loop,reverse,cg_mode);
////							}
//
//	}
//	else
//		self_clash_count++;
//}


//			//
//			// "Little vibrator"
//			//
//			int max_shakes = 1000;
//			float shake_rmsd = 9999;
//			float shake_ang = 5 * M_PI/180;
//			if(enable_vib && rmsd < rmsd_crit )
//			{
//				if (reverse == 1)
//					for (int z = 0; z < nr_atoms_loop; z++) // Mon
//						marker[z] =	marker_r[z]; // Marker re-initialization reversely
//				else
//					for (int z = 0; z < nr_atoms_loop; z++) // Mon
//						marker[z] = marker_i[z];
//
//				for(int i=0; i<nr_atoms_loop; i++)
//					dihedral_angle_aux[i] = dihedral_angle[i]; // Backup dihedral angles
//
//				int n_shakes=0;
//				while (rmsd > 0.1 && n_shakes < max_shakes)
//				{
//					// Randomize all mobile torsions
//					for(int i = 2; i < nr_atoms_loop-1; i++ )
//					{
//						if(marker[i] == 1)
//						{
//							dihedral_angle_aux[i + 1] = shake_ang * ((float)i/nr_atoms_loop) * (rg->Random()-0.5) + dihedral_angle[i+1]; // Mon: Check this!
//							correctdihedral(&dihedral_angle_aux[i + 1]);
//						}
//					}
//					isafrotatetreef_cg(co->el,p1co->el,nr_atoms_loop,bond_length,dihedral_angle_aux,prearray3,prearray4,prearray5,prearray6,prearray1,prearray2,full,cg_mode);
//					rmsd = rmsd_anchor3(co->el+last,anchor_end->el);
//					if(rmsd < shake_rmsd)
//					{
//						shake_rmsd = rmsd;
//						for(int i=0; i<nr_atoms_loop; i++)
//							dihedral_angle[i] = dihedral_angle_aux[i]; // Backup dihedral angles
//					}
//					// fprintf(stderr,"Shake rmsd= %f\n", rmsd);
//					n_shakes++;
//				}
//
//				// fprintf(stderr,"Minimum rmsd= %f  n_shakes= %d\n", shake_rmsd, n_shakes);
//			}

