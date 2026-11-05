/*
 * Stat.h
 *
 *  Created on: June 20, 2011
 *      Author: Pieter Chys
 */

#ifndef STAT_H_
#define STAT_H_


#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <vector>
#include <cmath>

#include "statistics.h"



class Stat {
public:
	Stat(int i);												//	Constructor
	~Stat();													//	Destructor

	double *array;											// 	Data array
	double *x_array;										//	x-values Data array

	double minimum;												//	Minimum value in array	--- altered from private
	double average;												//	Average value in array  --- altered from private
	double deviation;											//	Standard deviation values in array

	void read(double *);
	void max();
	void min();
	void avg();
	void std();
	void skew();
	void kurt();
	void imax();
	void imin();
	void compute();
	void calculate();
	void print2(char *text,int num_format, FILE *output);
	void io(double *);
	void histogram(int size_histogram,char* xlabel,char* ylabel);
	void plot2d(char* xlabel,char* ylabel,int cycle_limit);
	void polar(int size_histogram,char* xlabel,char* ylabel);
	void datafile(FILE *);
	void printsize(FILE *);
	// Mon's methods
	int nless(double thr);

private:


    double maximum;												//	Maximum value in array
	double skewness;											//	Skewness array values
	double kurtosis;											//	Kurtosis array values
	int imaximum;												//	Integer marker maximum value
	int iminimum;												//	Integer marker minimum value
	int size;													//	Size array
	int index;													//	index active element


};



#endif /* STAT_H_ */
