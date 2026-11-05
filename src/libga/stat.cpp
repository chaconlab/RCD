/*
 * stat.cpp
 *
 *  Created on: Jun 20, 2011
 *      Author: Pieter Chys
 */

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <vector>
#include <cmath>

#include "stat.h"
#include "statistics.h"




Stat::Stat(int i)
{
	// TODO Auto-generated constructor stub
	size  		= 	i;
	//array 		= 	new double[size];
	//x_array		=	new double[size];
	array = (double *) malloc( sizeof(double) * size);
	x_array = (double *) malloc( sizeof(double) * size);
	//	Initialize array values to zero
	for(int j = 0; j < size; j++)
	{
		array[j] 	= 0.0;
		x_array[j]	=	j;										//	initialize as natural numbers
	}

	maximum   = 	0.0;
	minimum   = 	0.0;
	average   = 	0.0;
	deviation = 	0.0;
	skewness  = 	0.0;
	kurtosis  =	0.0;
	imaximum  =	0;
	iminimum  = 	0;
	index	  =	0;

}

// Destructor of class 		Stat
Stat::~Stat()
{
	// TODO Auto-generated destructor stub
	//delete[] array;
	//delete[] x_array;
	free(array);
	free(x_array);
}

//	Stat::read method to read in data for array
void Stat::read(double *data)
{
	for(int j = 0; j < size; j++ )
	{
		array[j] = data[j];
	}
}



// Stat::mat method to obtain maximum value in array
void Stat::max()
{
	maximum = array[0];

	for (int i = 0; i < size; i++)
	{
		if (array[i] > maximum)
		{
			maximum 	= 	array[i];
			imaximum 	= 	i;
		}
	}
}



// Stat::min method to obtain minimum value in array
void Stat::min()
{
	max();
	minimum = maximum;
	for (int i = 0; i < size; i++)
	{
		if (array[i] < minimum)
		{
			minimum 	= 	array[i];
			iminimum 	= 	i;
		}
	}
}



//	Stat::avg method to calculate average value in array
void Stat::avg()
{
	double sum = 0.0;
	for (int i = 0; i < size; i++)
	{
		sum += array[i];
	}

	average = sum/size;
}



//	Stat::std method to calculate standard deviation in array
void Stat::std()
{
	double sum = 0.0;
	avg();
	for (int i = 0; i < size; i++)
	{
		sum  +=  pow((double) ( array[i] - average),2);
	}

	if (size != 1)
	{
		deviation = sqrt( double (sum/(size - 1)) );

	} else if (size == 1)
	{
		//printf("\n");
		//printf("Number of cycles for valuable statistics should be > 30\n");
		//printf("\n");

		deviation = 0.0;

	}
}



//	Stat::skew method to calculate skewness of array values
void Stat::skew()
{

	double m3 = 0.0;
	double m2 = 0.0;
	avg();

	for (int i = 0; i < size; i++)
	{
		m3 	+=  powl((array[i]  - average),3);			//	power function for 	double		type
		m2  +=  powl((array[i]  - average),2);
	}

	m3  /=  size;
	m2  /=  size;
	m2  =   powl(sqrtl(m2),3);

	//IN DEBUG MODE THE NUMBERS GO OUT OF SCALE ---- ALGORITHM MUST BE CHANGED TO NOT SUM COMPLETELY BUT DO IT STEPWISE ADDITIVELY
	//printf("Size is %15i\n",size);
	//printf("m2 is equal to %15.3lf\n",m2);
	//printf("m3 is equal to %15.3lf\n",m3);
	//fflush(stdout);
	//getchar();

	if(size != 1)
	{
		skewness = m3/m2;

	} else if (size == 1)
	{
		skewness = 0.0;			//:REMARK: in principel skewness is undefined: 0/0 is undefined
	}
}



//	Stat::kurt method to calculate kurtosis of array values
void Stat::kurt()
{
	double m4 = 0.0;
	double m2 = 0.0;
	avg();

	for (int i = 0; i < size; i++)
	{
		m4 	+=  powl((double) ( array[i] - average),4);
		m2  +=  powl((double) ( array[i] - average),2);
	}
	m4 /= size;
	m2 /= size;
	m2 = powl(m2,2);

	if(size != 1)
	{
		kurtosis =  m4/m2;

	} else if (size == 1)
	{
		kurtosis = 0.0;			//:REMARK: in principel curtosis is undefined: 0/0 is undefined
	}

}



//	Stat::imax method to give index of maximum value of array values
void Stat::imax()
{

	maximum = array[0];

	for (int i = 0; i < size; i++)
	{
		if (array[i] > maximum)
		{
			maximum = array[i];
			imaximum =	i;

		} else {	}
	}

}



//	Stat::imin method to give index of minimum value of array values
void Stat::imin()
{
	max();

	minimum = maximum;

	for (int i = 0; i < size; i++)
	{
		if (array[i] < minimum)
		{
			minimum = array[i];
			iminimum = i;

		} else { 	}
	}
}



//	Stat::compute  calculates all statistical parameters for an 1-dimensional index
void Stat::compute()
{
	max();
	min();
	avg();
	std();
	skew();
	kurt();
	imax();
	imin();

}



//	Stat::compute  calculates all statistical parameters for an 1-dimensional index
void Stat::io(double *data)
{
	read(data);
	max();
	min();
	avg();
	std();
	skew();
	kurt();
	imax();
	imin();

}



//	Stat::compute  calculates all statistical parameters for an 1-dimensional index
void Stat::print2(char *text, int num_format,FILE *stream)
{


	std::string line(80,'-');

	if (num_format == 1)
	{
		fprintf(stream,"%s\n"," ");
		fprintf(stream,"%-30s\n",text);
		fprintf(stream,"%s\n"," ");
		fprintf(stream,"%s\n",line.c_str());
		fprintf(stream,"%-30s %15s %-5s %15s %-5s\n","Parameter","Exact", " ","Rounded"," ");
		fprintf(stream,"%s\n",line.c_str());
		fprintf(stream,"%s\n"," ");
		fprintf(stream,"%-30s %15.6E %-5s %15.3E %-5s\n","Average",average," ",average," ");
		fprintf(stream,"%-30s %15.6E %-5s %15.3E %-5s\n","Standard deviation",deviation," ",deviation," ");
		fprintf(stream,"%-30s %15.6E %-5s %15.3E %-5s\n","Maximum",maximum," ",maximum," ");
		fprintf(stream,"%-30s %15.6E %-5s %15.3E %-5s\n","Minimum",minimum," ",minimum," ");
		fprintf(stream,"%-30s %15.6E %-5s %15.3E %-5s\n","Skewness",skewness," ",skewness," ");
		fprintf(stream,"%-30s %15.6E %-5s %15.3E %-5s\n","Kurtosis",kurtosis," ",kurtosis," ");

	} else if (num_format == 2)					/*Output for doubles*/
	{
		fprintf(stream,"%s\n"," ");
		fprintf(stream,"%-30s\n",text);
		fprintf(stream,"%s\n"," ");
		fprintf(stream,"%s\n",line.c_str());
		fprintf(stream,"%-30s %15s %-5s %15s %-5s\n","Parameter","Exact (micros)", " ","Rounded (micros)"," ");
		fprintf(stream,"%s\n",line.c_str());
		fprintf(stream,"%s\n"," ");
		fprintf(stream,"%-30s %15.6lf %-5s %15.0lf %-5s\n","Average",average," ",average," ");
		fprintf(stream,"%-30s %15.6lf %-5s %15.0lf %-5s\n","Standard deviation",deviation," ",deviation," ");
		fprintf(stream,"%-30s %15.6lf %-5s %15.0lf %-5s\n","Maximum",maximum," ",maximum," ");
		fprintf(stream,"%-30s %15.6lf %-5s %15.0lf %-5s\n","Minimum",minimum," ",minimum," ");
		fprintf(stream,"%-30s %15.6lf %-5s %15.3lf %-5s\n","Skewness",skewness," ",skewness," ");
		fprintf(stream,"%-30s %15.6lf %-5s %15.3lf %-5s\n","Kurtosis",kurtosis," ",kurtosis," ");

	} else if (num_format == 3)					/*Output for doubles but with ceil function for some values*/
	{											/*LIKELY REDUNDANT WITH FORMAT 2 DESPITE EXPLICIT CEILING ->    *                 	 	 	 	 	 	 	 	 	 	 *implementation of printf specifies this               	 	*/
		fprintf(stream,"%s\n"," ");
		fprintf(stream,"%-30s\n",text);
		fprintf(stream,"%s\n"," ");
		fprintf(stream,"%s\n",line.c_str());
		fprintf(stream,"%-30s %15s %-5s %15s %-5s\n","Parameter","Exact", " ","Rounded"," ");
		fprintf(stream,"%s\n",line.c_str());
		fprintf(stream,"%s\n"," ");
		fprintf(stream,"%-30s %15.6lf %-5s %15.0lf %-5s\n","Average",average," ",  ceil(average)," ");
		fprintf(stream,"%-30s %15.6lf %-5s %15.0lf %-5s\n","Standard deviation",deviation," ", ceil(deviation)," ");
		fprintf(stream,"%-30s %15.6lf %-5s %15.0lf %-5s\n","Maximum",maximum," ", ceil(maximum)," ");
		fprintf(stream,"%-30s %15.6lf %-5s %15.0lf %-5s\n","Minimum",minimum," ", ceil(minimum)," ");
		fprintf(stream,"%-30s %15.6lf %-5s %15.3lf %-5s\n","Skewness",skewness," ",skewness," ");
		fprintf(stream,"%-30s %15.6lf %-5s %15.3lf %-5s\n","Kurtosis",kurtosis," ",kurtosis," ");

	}  else if (num_format == 4)
	{
		fprintf(stream,"%s\n"," ");
		fprintf(stream,"%-30s\n",text);
		fprintf(stream,"%s\n"," ");
		fprintf(stream,"%s\n",line.c_str());
		fprintf(stream,"%-30s %15s %-5s %15s %-5s\n","Parameter","Exact (Angstrom)", " ","Rounded (Angstrom)"," ");
		fprintf(stream,"%s\n",line.c_str());
		fprintf(stream,"%s\n"," ");
		fprintf(stream,"%-30s %15.6lf %-5s %15.3lf %-5s\n","Average",average," ",average," ");
		fprintf(stream,"%-30s %15.6lf %-5s %15.3lf %-5s\n","Standard deviation",deviation," ",deviation," ");
		fprintf(stream,"%-30s %15.6lf %-5s %15.3lf %-5s\n","Maximum",maximum," ",maximum," ");
		fprintf(stream,"%-30s %15.6lf %-5s %15.3lf %-5s\n","Minimum",minimum," ",minimum," ");
		fprintf(stream,"%-30s %15.6lf %-5s %15.3lf %-5s\n","Skewness",skewness," ",skewness," ");
		fprintf(stream,"%-30s %15.6lf %-5s %15.3lf %-5s\n","Kurtosis",kurtosis," ",kurtosis," ");

	} else if (num_format == 5)
	{
          fprintf(stream,"%s avg %.3lf+/-%.0lf \tMax %.0lf \tMin %.0lf\n",text, average, deviation, maximum, minimum);
	} else if (num_format == 6) // rmsd
	{
          fprintf(stream,"%s avg %.3lf+/-%.3lf \tMax %.3lf \tMin %.3lf\n",text, average, deviation, maximum, minimum);
	}




}

void Stat::printsize(FILE * stream)
{
	//fprintf(stream,"%s\n"," ");
	fprintf(stream,"%-40s %25d\n","Number samples",size);
	//fprintf(stream,"%s\n"," ");
}

void Stat::histogram(int size_histogram,char* xlabel,char* ylabel)
{

	double *xhist;
	double *yhist;

	xhist = (double *) malloc(sizeof(double) * size_histogram);
	yhist = (double *) malloc(sizeof(double) * size_histogram);
	for (int i = 0; i < size_histogram; i++) {
		xhist[i] = 0.0;
		yhist[i] = 0.0;
	}

	max();
	min();
	makehistogram(array,size,(double) minimum, (double) maximum,xhist,yhist,size_histogram);
 /*
	plot2d_write(xhist,yhist,xlabel,ylabel,size_histogram);	 	//	some problems before but seems to be settled out
	plot2d_screen(xhist,yhist,xlabel,ylabel,size_histogram);
*/
	//	Deallocation of arrays
	free(xhist);			//	lines added but not yet completely debugged
	free(yhist);			//	lines added but not yet debugged

}

void Stat::plot2d(char* xlabel,char* ylabel,int cycle_limit)
{

	//double *xhist;
	//double *yhist;

	//xhist = (double *) malloc(sizeof(double) * size_histogram);
	//yhist = (double *) malloc(sizeof(double) * size_histogram);
	//for (int i = 0; i < size_histogram; i++) {
		//xhist[i] = 0.0;
		//yhist[i] = 0.0;
	//}


	//max();
	//min();
	//makehistogram(array,size,minimum,maximum,xhist,yhist,size_histogram);

	//plot2d_write(xhist,yhist,xlabel,ylabel,size_histogram);	 	//	some problems before but seems to be settled out
	//plot2d_screen(x_array,array,xlabel,ylabel,cycle_limit);

}


void Stat::polar(int size_histogram,char* xlabel,char* ylabel){	//	PROBLEM FOR PLOTTIN SMALL VALUES

	double *xhist;
	double *yhist;

	xhist = (double *) malloc(sizeof(double) * size_histogram);
	yhist = (double *) malloc(sizeof(double) * size_histogram);

	for (int i = 0; i < size_histogram; i++) {
		xhist[i] = 0.0;
		yhist[i] = 0.0;
	}

	max();
	min();
	makehistogram(array,size, (double) minimum, (double) maximum,xhist,yhist,size_histogram);

	for (int i = 0; i > size; i++){
		//	Correct negative angles to PI - 2*Pi range
		if(xhist[i] < 0.0){
			xhist[i] += 2.0*M_PI;
		} else { };

		//	scale to 360 degrees
		xhist[i] *= 180.0/M_PI;
		//yhist[i] /= (double) -log((double) (yhist[i]));		//	extra function to make small values more visible // does not work that well
		//yhist[i] = (double) pow((double) yhist[i],0.0001);		//	extra function to make small values more visible // does not work that well
	}



//	polar_write(xhist,yhist,xlabel,ylabel,size_histogram);	 	//	some problems before but seems to be settled out
	//plot2d_write(xhist,yhist,xlabel,ylabel,size_histogram);	 //	some problems before but seems to be settled out
	//plot2d_screen(xhist,yhist,xlabel,ylabel,size_histogram);

	//	Deallocation of arrays
	free(xhist);			//	lines added but not yet completely debugged
	free(yhist);			//	lines added but not yet debugged

}


void Stat::datafile(FILE* p_file)
{
	//FILE* data_file;													//	pointer to file
	//p_file	= fopen(filename,"w");		 							//	read the benchmark.txt

	for(int i = 0; i < size; i++)
	{
		fprintf(p_file,"%15.5lf\n", array[i]);
	}

	//fclose(p_file);
}

// Mon added methods:

// Returns the number of values with less or equal values than "thr"
int Stat::nless(double thr)
{
	int n = 0;
	for(int i = 0; i < size; i++)
		if(array[i] <= thr)
			n++;
	return n;
}



