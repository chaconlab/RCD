#ifndef GUARD_own_io_h
#define GUARD_own_io_h

//own_io.h
#include <iostream>
#include <cstdlib>
#include <cstring>
#include "libpdb/include/Macromolecule.h"
#include "libpdb/include/pdbIter.h"
#include "libpdb/include/Fragment.h"

using namespace std;



void pointeralloc(double ***,int, int);
//void pointeralloc(double**,int);
void freepointer(double **,int);

/*message output*/
void show3(char *,
		   char *,
		   double **,
		   int);

void show3(char *,
		   char *,
		   double *,
		   double *,
		   double *,
		   int);

void show1(char *,
		   char *,
		   double *,
		   int);

void show1(char *,
		   char *,
		   int *,
		   int);

void show1(char *,
		   int);

void show1(char *,
		   double);

void message(char *,
		     int);

//void indicator(char *text,
//			   int index,
//			   int size_loop);


#endif /* OWN_IO_H_ */

/*===============================================================================
 * own_io.h
 *
 *  Created on: Dec 7, 2010
 *      Author: Pieter Chys
 *===============================================================================
 */
