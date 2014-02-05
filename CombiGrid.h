#include<stdio.h>

#ifndef COMBIGRID_H_
	#define COMBIGRID_H_

	class CombiGrid {
		double * val; // array to store coefficients/function values
		int d; // dimension
		int * l; // level
		int * n; // number of grid points in each dimension
		int n0Aligned; // number of values in first (0) dimension with padding
		long unsigned int size; // number of grid points
		long unsigned int arraySize; // number of values with padding
	public:
		CombiGrid(int dim, int* levels);
		CombiGrid(int dim, int* levels, int alignment );

		~CombiGrid();

		void append (int alignment); // alignment in bytes

		int getSize();
		int getArraySize();
		double getValue(int pos);
		void printSize();
		void print();
		void printValues2DArr(int size, int offset, int n0);
		void printValues();
		void printValues(int alignment);

		void compare(CombiGrid * cc); // compare to combination grids point by point

		void setValues(double (* func) (double *) );
		void setValues(double (* func) (double *) , int alignment);

		void hierarchize1DUnoptimized(int start, int stride, int size, int dim);
		double hierarchizeUnoptimized(); // returns time needed for hierarchization in seconds

		void hierarchize1DOptimized(int start, int stride, int size, int dim, int blockSize);
		double hierarchizeOptimized(int blockSize); // returns time needed for hierarchization in seconds
	};

#endif /* COMBIGRID_H_ */
