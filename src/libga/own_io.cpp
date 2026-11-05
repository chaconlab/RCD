/***********************************************************************
 * 									   								   *
 * COPYRIGHT: 2011,2012 Pieter Jozef Chys							   *
 *   	   															   *
 ***********************************************************************
 * 																	   *
 * FILE: 				OWN_IO.CPP								   	   *
 *																	   *
 * DEPENDS ON:		    - linux										   *
 *  					- libpdb/Macromolecule.h					   *
 *  					- libpdb/pdbIter.h					           *
 *                      - libpdb/Fragment.h                            *
 * 																	   *
 * DATE OF CREATION: 	December 10th, 2011							   *
 * LAST MODIFIED:	 	April 18th, 2012							   *
 * AUTHOR (S):  		Pieter Jozef Chys							   *
 *      															   *
 ***********************************************************************
 *																	   *
 * DESCRIPTION:	cpp-file which contains the input/output routines	   *
 * 				relavant in the use of the library MY_LIBRARY          *
 *																	   *
 **********************************************************************/



#include <iostream>
#include <cstdlib>
#include <cstring>

#include "libpdb/include/Macromolecule.h"
#include "libpdb/include/pdbIter.h"
#include "libpdb/include/Fragment.h"


using namespace std;


//:

/*======================================================================
 *MESSAGE prints a text message to the standard output and has an input
 *argument which allows to parameterize the execution the function.
 *This argument has Boolean functioning:if button is equal to 1 we
 *the message is printed, if equal to 0 the output is not printed.
 *----------------------------------------------------------------------
 *INPUT:- char *TITLE which is the text to be displayed
 * 	    - int BUTTON which determines whether output is produced (= 1)
 *        or not (= 0)
 *----------------------------------------------------------------------
 *OUTPUT: - explicitly VOID
 *        - output TITLE to standard output (normally screen) depending
 *          on the integer value of BUTTON (= 0 or  = 1)
 *======================================================================*/
void message(char *title,int button)
{
	if (button == 1)
	{
		fprintf(stdout,"%-s\n",title);
		fflush(stdout);

	}	else if (button == 0) {

		/*do not output to standard output*/
	}
}


/*======================================================================
 *SHOW3 prints the XYZ coordinates present in a 2D array[0][3] to the
 *standard output.
 *----------------------------------------------------------------------
 *INPUT:- char *TITLE which is the title to the numerical columns to be
 *		  displayed
 * 	    - char *FORMATTING which the text that must go together with
 *        the numerical data (most left column)
 *      - double **DATA which is of format [][3] and contains the
 *        numerical data (normally XYZ coordinates)
 *     -  int SIZE which is the size of the array
 *----------------------------------------------------------------------
 *OUTPUT: - explicitly VOID
 *        - output TITLE and the numerical data to the standard output
 *          (normally screen)
 *----------------------------------------------------------------------
 *REMARKS: - function closely related sister variants
 *======================================================================*/
void show3(char *title,char *formatting,double** data,int size)
{
	fprintf(stdout,"\n");
	fprintf(stdout,"%s\n",title);
	fprintf(stdout,"\n");
	fflush(stdout);

	for(int i = 0; i < size; i++ )
	{
		fprintf(stdout,"%-s %8.3f %8.3f %8.3f\n",formatting,data[i][0],data[i][1],data[i][2]);
		fflush(stdout);
	}

	fprintf(stdout,"\n");
	fflush(stdout);

}


/*======================================================================
 *This variant of SHOW3 prints the XYZ coordinates present in three distinct
 *1D arrays to the standard output (normally screen). This variant is
 *identical to SHOW3 with a single 2D array but uses three single 1D arrays
 *instead.
 *----------------------------------------------------------------------
 *INPUT:- char *TITLE which is the title to the numerical columns to be
 *		  displayed
 * 	    - char *FORMATTING which is the format specification as should
 *        be used in PRINTF and FPRINTF functions for the three numerical
 *        columns
 *      - double *DATA1 which is of format [] and contains the
 *        numerical data to be printed to the 1st column
 *      - double *DATA2 which is of format [] and contains the
 *        numerical data to be printed to the 2nd column
 *      - double *DATA3 which is of format [] and contains the
 *        numerical data to be printed to the 3rd column
 *      - int SIZE which is the size of the array
 *----------------------------------------------------------------------
 *OUTPUT: - explicitly VOID
 *        - output TITLE and the numerical data to the standard output
 *          (normally screen) from the three data arrays
 *----------------------------------------------------------------------
 *REMARKS: - function closely related sister variants
 *======================================================================*/
void show3(char *title,char *formatting,double* data1,double* data2,double* data3,int size)
{
	fprintf(stdout,"\n");
	fprintf(stdout,"%s\n",title);
	fprintf(stdout,"\n");
	fflush(stdout);

	for(int i = 0; i < size; i++ )
	{
		fprintf(stdout,"%-s %8.3f %8.3f %8.3f\n",formatting,data1[i],data2[i],data3[i]);
		fflush(stdout);
	}
	fprintf(stdout,"\n");
	fflush(stdout);

}


/*======================================================================cpp
 *SHOW1 prints numerical data present in a 1D array[] to the standard
 *output. It the variant from SHOW3 which only contains 1 column output
 *and outputs data of type DOUBLE.
 *----------------------------------------------------------------------
 *INPUT:- char *TITLE which is the title to the numerical columns to be
 *		  displayed
 * 	    - char *FORMATTING which is the text that goes together with
 *        the numerical data (most left column)
 *      - double *DATA which is of format [] and contains the
 *        numerical data
 *      - int SIZE which refers to the size of the array
 *----------------------------------------------------------------------
 *OUTPUT: - explicitly VOID
 *        - output TITLE and the numerical data to the standard output
 *          (normally screen)
 *----------------------------------------------------------------------
 *REMARKS: - function closely related sister variants
 *         - output for DOUBLE numerical data
 *======================================================================*/
void show1(char *title,char *formatting,double* data,int size)
{
	fprintf(stdout,"\n");
	fprintf(stdout,"%s\n",title);
	fprintf(stdout,"\n");
	fflush(stdout);

	for(int i = 0; i < size; i++ )
	{
		fprintf(stdout,"%-s %8.3f\n",formatting,data[i]);
		fflush(stdout);
	}

	fprintf(stdout,"\n");
	fflush(stdout);

}


/*======================================================================
 *SHOW1 prints numerical data present in a 1D array[] to the
 *standard output. It the variant from SHOW3 which only contains 1
 *column output and outputs data from type INT.
 *----------------------------------------------------------------------
 *INPUT:- char *TITLE which is the title to the numerical columns to be
 *		  displayed
 * 	    - char *FORMATTING which is the text that goes together with
 *        the numerical data (most left column)
 *      - double *DATA which is of format [] and contains the
 *        numerical data
 *      - int SIZE which refers to the size of the array
 *----------------------------------------------------------------------
 *OUTPUT: - explicitly VOID
 *        - output TITLE and the numerical data to the standard output
 *          (normally screen)
 *----------------------------------------------------------------------
 *REMARKS: - function closely related sister variants
 *         - output for DOUBLE numerical data
 *======================================================================*/
void show1(char *title,char *formatting,int* data,int size)
{
	fprintf(stdout,"\n");
	fprintf(stdout,"%s\n",title);
	fprintf(stdout,"\n");
	fflush(stdout);

	for(int i = 0; i < size; i++ )
	{
		fprintf(stdout,"%-s %i\n",formatting,data[i]);
		fflush(stdout);
	}

	fprintf(stdout,"\n");
	fflush(stdout);

}


/*======================================================================
 *SHOW1 prints a numerical INT number to the standard output. It is
 *the variant from SHOW3 which only contains 1 element and outputs
 *data from type INT.
 *----------------------------------------------------------------------
 *INPUT:- char *TITLE which is the title to the numerical columns to be
 *		  displayed
 * 	    - char *FORMATTING which is the text that goes together with
 *        the numerical data (most left column)
 *      - int DATA which is an integer number
 *----------------------------------------------------------------------
 *OUTPUT: - explicitly VOID
 *        - output TITLE and integer number to the standard output
 *          (normally screen)
 *----------------------------------------------------------------------
 *REMARKS: - function closely related sister variants
 *         - output for INT
 *======================================================================*/
void show1(char *title,int data)
{
	fprintf(stdout,"%-s %i\n",title,data);
	fflush(stdout);
}


/*======================================================================
 *SHOW1 prints a numerical DOUBLE number to the standard output. It the variant
 *from SHOW3 which only contains 1 element and outputs data from type INT.
 *----------------------------------------------------------------------
 *INPUT:- char *TITLE which is the title to the numerical columns to be
 *		  displayed
 * 	    - char *FORMATTING which is the text that goes together with
 *        the numerical data (most left column)
 *      - double DATA which is a double number
 *----------------------------------------------------------------------
 *OUTPUT: - explicitly VOID
 *        - output TITLE and double number to the standard output
 *          (normally screen)
 *----------------------------------------------------------------------
 *REMARKS: - function closely related sister variants
 *         - output for DOUBLE
 *======================================================================*/
void show1(char *title,double data)
{
	fprintf(stdout,"%-s %f\n",title,data);
	fflush(stdout);
}


/*============================================================================
 *READINPDB reads in a macromolecule from a PDB file and uses to that purpose
 *some functionality of the libary "Macromolecule.h" by means of its classes.
 *A macromolecule object is created and the pdb-file read in into that
 *macromolecule.  Member functions are used from this class to retrieve
 *the backbone atom coordinates (N, CA and C). Hereafter the coordinates are
 *returned to a 2D array returned as a pointer. The original macromolecules
 *destroyed.
 *------------------------------------------------------------------------------
 *INPUT: - char *FILENAME which refers to the pdb-file from which the atom
 *         coordinates are going to be selected and read in
 *       - macromolecule *MOL which is going to be used as auxiliary object 1
 *       - macromolecule *MOL2 which is going to be used as auxiliary object 2
 *       - int *STARTRES will contain the number of rows that will have to be
 *         skipped in the 2D array CO to get the atom coordinate of the start
 *         residue
 *       - int START is input residue position of the start residue
 *       - float **CO which will contain the coordinates
 *-------------------------------------------------------------------------------
 *OUTPUT: - explicitly VOID
 *        - double** CO which will contain the atom coordinates of the loop
 *          conformations stacked on top of each other: it refers to the
 *          backbone coordinates
 *-------------------------------------------------------------------------------
 *REMARKS: - uses functionality of "Macromolecule.h"
 *============================================================================= */



/*======================================================================
 *POINTERALLOC creates an allocated 2-dimensional array of doubles with
 *specified X dimension X_SIZE and Y dimension Y_SIZE. In addition, the
 *2D array is initialized to zero values.Thus (X_SIZE.Y_SIZE)*sizeof(double)
 *bytes are allocated covering*(X_SIZE X Y_SIZE) elements which are set
 *to zero. The deallocation has to be done with the sister function
 *FREEPOINTER.  The purpose of the function is mainly in allowing
 *compactification of the code if a lot of allocation needs to be done.
 *----------------------------------------------------------------------
 *INPUT:- explicitly VOID
 * 		- double ***MATRIX
 * 	    - int X_size = number elements in X dimension (number rows)
 *      - int Y_size = number elements in Y dimension (number columns)
 *----------------------------------------------------------------------
 *OUTPUT: - double ***MATRIX which now points to allocated memory
 *----------------------------------------------------------------------
 *REMARK: - the triple pointer must be passed from a double **POINTER
 *          MATRIX as &MATRIX; as such no triple pointer is ever created
 * 			in the work code
 *        - sister function FREEPOINTER deallocates the allocated memory
 *======================================================================*/
void pointeralloc(double*** matrix,int x_size, int y_size)
{
	*matrix = (double**) malloc(sizeof(double *)* x_size);

	for(int i = 0; i < x_size; i++)
	{
		(*matrix)[i] = (double *) malloc(sizeof(double) * y_size);
	}

	for(int i = 0; i < x_size; i++)
	{	for(int j = 0; j < y_size; j++)
		{
			(*matrix)[i][j] =  0.0;
		}
	}
}


 /*======================================================================
 *FREEPOINTER deallocates a 2-dimensional memory array of doubles with
 *specified X dimension X_SIZE. The array must have been created with
 *the brother function POINTERALLOC. The purpose of the function is
 *similar to POINTERALLOC: compactification of the code if a lot of
 *allocation needs to be done.
 *----------------------------------------------------------------------
 *INPUT:- double **MATRIX
 * 	    - int X_size = number elements in X dimension (number rows)
 *----------------------------------------------------------------------
 *OUTPUT: - explictly VOID
 *----------------------------------------------------------------------
 *REMARK: - the memory array must be passed with the original double
 *          **MATRIX pointer
 *        - sister function POINTERALLOC allocates the original memory
 *          array
 *======================================================================*/
void freepointer(double** matrix, int x_size)
{
	for(int i = 0; i < x_size; i++)
	{
		free(matrix[i]);
	}

	free(matrix);
}


/*======================================================================
 *INDICATOR shows at the terminal screen the progression of a number
 *of cycles
 *----------------------------------------------------------------------
 *	INPUT:  - int SIZE_LOOP indicating the integer value that corresponds
 *            with a 100% bar
 * 			- int INDEX which is the current integer index and will be
 * 		      converted to a percentage from the 100% bar
 *          - char *pointer TEXT with text which will be in front of bar
 *----------------------------------------------------------------------
 *  OUTPUT: - explicitly VOID
 * 			- output to terminal screen upon execution
 *----------------------------------------------------------------------
 *  REMARK: - in the computational cycles nothing should be calculated
 *            since otherwise the bar representation is interrupted
 *          - uses C++ strings
 *          - SIZE_BAR could be made variable and hence is possible source
 *            of BUGs
 *=====================================================================*/
//void indicator(char *text,int index,int size_loop)
//{
//	int bar;
//	int size_bar = 40;
//	float ratio;
//	ratio = (float) index/ (float) size_loop;
//	bar = (int) floor( ratio * (float) size_bar );
//
//	fprintf(stdout,"\r%s [",text);
//
//	for(int i=0; i<bar; i++)
//		fprintf(stdout,"=");
//	fprintf(stdout,">");
//
//	for(int i=bar+1; i<size_bar; i++)
//		fprintf(stdout,"-");
//
//	if(bar == size_bar)
//		fprintf(stdout,"] %3.0f%% Finished!\n",ratio*100);
//	else
//		fprintf(stdout,"] %3.0f%%",ratio*100);
//	fflush(stdout);
//}






