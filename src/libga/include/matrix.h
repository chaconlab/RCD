/*
 * matrix.h
 *
 *  Created on: November 18, 2011
 *      Author: Pieter Chys
 */

#ifndef MATRIX_H_
#define MATRIX_H_



class Matrix {
	public:
		Matrix(int i,int j);
		~Matrix();												//	virtual normally in front

		void copy(int x,int y,double **input, int diff_x, int diff_y);
		void dihedral(int i,int j);
		//void print(int k,char *text);

		double **el;
		int 	size_x;
		int     size_y;
		//double **bond;
		//double **bondnorm;


		//double sin;
		//double cos;
		//double 	phi;
		//double 	teta;
		//int 	size;

		//void cl_length();
		//void cl_valence();
		//void cl_dihedral();

		//void initialize();

	private:


};



#endif /* MATRIX_H_ */
