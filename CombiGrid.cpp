#define mypow2(b) (1 << b)
#define mydiv2(b) (b >> 1)

#include "CombiGrid.h"
#include <sstream>
#include <iostream>
//#include <malloc.h>
#include "func.cpp"
#include "immintrin.h"
#include "omp.h"
#include <chrono>

// wieder raus?
#include <math.h>


using namespace std;

CombiGrid::CombiGrid(int dd, int * ll) {
	d = dd;
	l = new int[d];
	n = new int[d];
	int i;
	size = 1;
	for (i = 0; i < d; i++){
		l[i]= ll[i];
		n[i] = mypow2(ll[i])-1; // number of grid points in dimension i
		size *= n[i]; // total number of grid points
	}
	arraySize = size;
	val = new double[size];
	if (val == 0) exit(35);
	n0Aligned = n[0];
}

CombiGrid::CombiGrid(int dd, int * ll, int alignment) {
	// alignment in bytes
	// alignment multiple of 32 bytes for AVX
	//grid needs to be aligned for the blocked (optimized) version of the code
	d = dd;
	l = new int[d];
	n = new int[d];
	int i;
	// lengthen 1st dimensions
	n[0] = mypow2(ll[0])-1;
	l[0]= ll[0];
	n0Aligned = ceil((double) (mypow2(ll[0])-1)/alignment*sizeof(double)) * alignment/sizeof(double);
	size = n[0];
	arraySize = n0Aligned;
	for (i = 1; i < d; i++){
		l[i]= ll[i];
		n[i] = mypow2(ll[i])-1;
		size *= n[i];
		arraySize *= n[i];
	}
	val = (double * )_mm_malloc(arraySize*sizeof(double),alignment);
	if (val == 0) {cout << "EXIT: malloc Konstruktor"<< endl; exit(2);}
}


void CombiGrid::append(int alignment){ // if grid has been initalized without alignment this function pads/aligns it a posteriori
	// alignment in bytes
	// alignment multiple of 32 bytes for AVX
	// grid needs to be aligned for the blocked (optimized) version of the code
	int alignmentDoubles = alignment/sizeof(double);
	n0Aligned = ceil((double) n[0]/alignmentDoubles) * alignmentDoubles;
	arraySize = size/n[0]*n0Aligned;
	double * newVal = (double * )_mm_malloc(arraySize*sizeof(double),alignment);
	if (newVal == 0) {
		cout << "EXIT: malloc align"<< endl;
		exit(2);
	}
	 int i;
	 int ctr = 0;
	 int chunkSize = n0Aligned;
	 int nbrOfChunks = arraySize / chunkSize;
	 for (int chunkCtr = 0; chunkCtr < nbrOfChunks; chunkCtr++) {
		 for (i = 0; i < n[0]; i++) {
			 newVal[chunkCtr*chunkSize+ i] = val[ctr];
			 ctr++;
		 }
		 newVal[chunkCtr*chunkSize+i] = INFINITY;
	 }
	 delete [] val;
	 val = newVal;
}


CombiGrid::~CombiGrid(){
	if (size == arraySize){
		delete [] l;
		delete [] n;
		delete  [] val;
	} else {
		delete [] l;
		delete [] n;
		_mm_free(val);
	}
}

int CombiGrid::getSize(){
	return size;
}

int CombiGrid::getArraySize(){
	return arraySize;
}
double CombiGrid::getValue(int pos){
	return val[pos];
}

void CombiGrid::printSize(){
	cout << "Approx. GridSize: "<< 8*((double) arraySize)/(1024*1024*1024) << " GB\n";
}

void CombiGrid::print() {
	cout << "Dimension: " << d << "\nLevel Vektor: ";
	for (int i=0; i < d; i++)
		cout << l[i] << " ";
	cout << "\nGrid points per Dim: ";
	for (int i=0; i< d; i++)
		cout << n[i] << " " ;
	cout << "\nTotal number of grid points: "<< size << "\n\n";
	return;
}



 void CombiGrid::printValues2DArr(int size, int offset, int n0) {
	for (long unsigned int  i=1; i <= size; i++){
		cout << val[offset+i-1]<< "\t";
		if (i % n0 == 0)
			cout << endl;
		}
	cout << endl;
	return;
}


void CombiGrid::printValues() {
	printf("\n");
	if (d <= 2) printValues2DArr(size, 0, n[0]);
	else {
		int * currentLevels = new int [d];
		int chunkSize = n[0]*n[1];
		int nbrOfChunks = size / chunkSize;
		for (int ctr = 0; ctr < nbrOfChunks; ctr++){
			cout << endl;
			printValues2DArr(chunkSize, ctr*chunkSize, n[0]);
			// rest is only for formatting - make an additioanl new line when increasing the level for d >= 2
			currentLevels[2]++;
			for (int dd = 2; dd < d; dd++){
				if (currentLevels[dd] == n[dd]) {
					currentLevels[dd] = 0;
					currentLevels[dd+1]++;
					cout << endl;
				}
			}
		}
	}
	return;
}

void CombiGrid::printValues(int alignment) {
	printf("\n");
	if (d <= 2) printValues2DArr(arraySize, 0, n0Aligned);
	else {
		int * currentLevels = new int [d];
		int chunks = n0Aligned*n[1];
		int nbrOfChunks = arraySize / chunks;
		for (int ctr = 0; ctr < nbrOfChunks; ctr++){
			cout << endl;
			printValues2DArr(chunks, ctr*chunks, n0Aligned);
			// rest is only for formatting - make an additioanl new line when increasing the level for d >= 2
			currentLevels[2]++;
			for (int dd = 2; dd < d; dd++){
				if (currentLevels[dd] == n[dd]) {
					currentLevels[dd] = 0;
					currentLevels[dd+1]++;
					cout << endl;
				}
			}
		}
	}
	return;
}

void CombiGrid::compare(CombiGrid * cc){
	  bool err = false;
	  if (this->getArraySize() != cc->getArraySize()){
		  cerr << "Grid sizes do not agree!"<< endl;
		  err = true;
		  return;
	  }
	  int i;
	  for ( i = 0; i < this->getArraySize(); i++){
		  if (fabs(this->getValue(i) - cc->getValue(i)) >1e-10){
			  cerr << "Grid point i = " << i << " does not agree"<< endl;
			  cerr << "the values are " << this->getValue(i) << " and " << cc->getValue(i) << endl;
			  err = true;
		  }
	  }
	  if (!err){
		  cout << "The combigrids agree."<< endl;
	  }
}

void CombiGrid::setValues(double (* func) (double *)){
	int ctr;
	double * stepsize = new double[d];
	double * x = new double [d];
	long unsigned int * levelSets = new long unsigned int [d];
	levelSets[0] = 1;
	int * dimCtr = new int [d];

	for (int dd = 1; dd < d; dd++){
		levelSets[dd] = levelSets[dd-1]*n[dd-1];
	}

	for (int dd = 0; dd < d; dd++){
		stepsize[dd] = pow(2,-l[dd]);
		dimCtr[dd] = 1;
	}

	for (ctr = 0; ctr< size; ctr++){
		for (int dd = 0; dd < d; dd++){
			if (dimCtr[dd] > n[dd] && dd == d-1) {
				cout << "at this point we should be done setting function values: CTR "<< ctr <<endl;
				exit (43);
			}
			if (dimCtr[dd] > n[dd] && dd < d-1) {
				dimCtr[dd] = 1;
				dimCtr[dd+1]++;
			}
		}
		for (int dd = 0; dd < d; dd ++){
			x[dd] = dimCtr[dd]*stepsize[dd];
		}
		val[ctr] = func(x);
		dimCtr[0]++;
	}
	return ;
}

void CombiGrid::setValues(double (* func) (double *), int alignment){
	alignment = alignment / 8; // from bytes to doubles
	int ctr, pos;
	double * stepsize = new double[d];
	double * x = new double [d];
	long unsigned int * levelSets = new long unsigned int [d];
	levelSets[0] = 1;
	int * dimCtr = new int [d];

	for (int dd = 1; dd < d; dd++){
		levelSets[dd] = levelSets[dd-1]*n[dd-1];
	}

	for (int dd = 0; dd < d; dd++){
		stepsize[dd] = pow(2,-l[dd]);
		dimCtr[dd] = 1;
	}
	pos = 0;
	for (ctr = 0; ctr< size; ctr++){
		for (int dd = 0; dd < d; dd++){
			if (dimCtr[dd] > n[dd] && dd == d-1) {
				cout << "at this point we should be done setting function values: CTR "<< ctr <<endl;
				exit (43);
			}
			if (dd ==0 && dimCtr[dd] > n[dd]) { // we are in first dim (padded!) and need to increment pos
				val[pos] = INFINITY; // pos points to padded point
				pos++;
				dimCtr[dd] = 1;
				dimCtr[dd+1]++;
			}
			if (dimCtr[dd] > n[dd] && dd < d-1) {
				dimCtr[dd] = 1;
				dimCtr[dd+1]++;
			}
		}
		for (int dd = 0; dd < d; dd ++){
			x[dd] = dimCtr[dd]*stepsize[dd];
		}
		val[pos] = func(x);
		dimCtr[0]++;
		pos++;
	}
	// set very last padded value
	val[arraySize-1] = INFINITY;
	return ;
}



inline void CombiGrid::hierarchize1DUnoptimized(int start, int stride, int size, int dim){
			int ll;
			int steps;
			int ctr;
			int offset, parentOffset;
			int stepsize;
			int parOffsetStrided;

			// ssa variables
			double val1, val2, val3, parL, parR;

			ll = l[dim];
			steps = mypow2(ll-1);
			offset = 0;
			stepsize = 2;
			parentOffset =  1;

			for (ll--; ll > 1; ll--){
				parOffsetStrided = parentOffset*stride;
				val[start+offset*stride] -= 0.5*val[start+offset*stride+parOffsetStrided];
				offset += stepsize;
				parL= 0.5*val[start+offset*stride -parOffsetStrided];;
				for (ctr = 1; ctr < steps-1; ctr++){
						val1 = val[start+offset*stride];
						parR = 0.5*val[start+offset*stride+parOffsetStrided];
						val2 = val1 - parL;
						val3 = val2 - parR;
						val[start+offset*stride] = val3;
						parL = parR;
						offset += stepsize;
				}
				val[start+offset*stride] -= parR;
				steps = steps >> 1;
				offset = mypow2(l[dim]- ll )-1;
				parentOffset =  stepsize;
				stepsize = stepsize << 1;
			}
			//	 level = 2 seperate
			// parRight is now used to store the root for both vertices of level 2
			parR =0.5*val[start+(offset+parentOffset)*stride];
			val[start+offset*stride] -= parR;
			offset += stepsize;
			val[start+offset*stride] -= parR;
			return ;
}

double CombiGrid::hierarchizeUnoptimized(){ // unoptimized but parallelized version of hierarchization
	int dim;
	int start;
	int stride =1 ;
	int ndim;
	int nbrOfPoles;
	int jump;
	div_t divresult;

	std::chrono::time_point<std::chrono::system_clock> startChrono, endChrono;
	std::chrono::duration<double> elapsed_seconds;
	startChrono = std::chrono::system_clock::now();

	//	 dimension 1 separate as start of each pole is easier to calculate
	ndim = n[0];
	nbrOfPoles = size /ndim;
	#pragma omp parallel for schedule(static) firstprivate (ndim, nbrOfPoles, start)
		for (int kk = 0; kk < nbrOfPoles; kk++){
			start = kk *ndim;
			hierarchize1DUnoptimized(start, 1, ndim,0);
		}
	// end dimension 1

	for (dim = 1; dim < d; dim ++){ // hierarchize for all dims
		stride *=ndim;
		ndim = n[dim];
		jump = stride*ndim;
		nbrOfPoles = size/ndim;
		#pragma omp parallel for schedule(static) firstprivate (divresult, stride, dim, ndim, nbrOfPoles, start, jump)
		for (int nn = 0; nn< nbrOfPoles; nn++){ // integer operations form bottleneck here -- nested loops are twice as slow
			divresult = div (nn,stride);
			start = divresult.quot*jump +divresult.rem;
			hierarchize1DUnoptimized(start, stride, ndim,dim);
		}
	} // end loop over dimension 2 to d

	endChrono = std::chrono::system_clock::now();
	elapsed_seconds = endChrono-startChrono;
	return 	elapsed_seconds.count();
}



inline void CombiGrid::hierarchize1DOptimized(int start, int stride, int size, int dim, int unroll){
			// only applicable for d >= 1  -- assert dim >= 1 ?
			int ll;
			int steps;
			int ctr;
			int offset, parentOffset;
			int stepsize;
			int parOffsetStrided;

			__m256d _valOrg, _valTemp, _valRes, _parR, _parL, _parR05, _parL05, _const05;
			_const05 = _mm256_set_pd(-0.5,-0.5,-0.5,-0.5);

			ll = l[dim];

			steps = mypow2(ll-1);
			offset = 0;
			stepsize = 2;
			parentOffset =  1;
			for (ll--; ll > 1; ll--){
				parOffsetStrided = parentOffset*stride;
				for (int poleLoop = 0; poleLoop < unroll; poleLoop+=4){
					_valOrg = _mm256_load_pd(val+start +offset*stride+poleLoop);
					_parR = _mm256_load_pd(val+start +offset*stride+parOffsetStrided+poleLoop);
					_parR05 = _mm256_mul_pd(_const05, _parR);
					_valRes = _mm256_add_pd(_valOrg, _parR05);
					_mm256_store_pd(val+start +offset*stride+poleLoop, _valRes);

				} // end poleLoop
				offset += stepsize;

				for (ctr = 1; ctr < steps-1; ctr++){
					for (int poleLoop = 0; poleLoop < unroll; poleLoop+=4){
						_parL = _mm256_load_pd(val+start +offset*stride-parOffsetStrided+poleLoop);
						_parL05 = _mm256_mul_pd(_const05, _parL);
						_valOrg = _mm256_load_pd(val+start +offset*stride+poleLoop);
						_parR = _mm256_load_pd(val+start +offset*stride+parOffsetStrided+poleLoop);
						_parR05 = _mm256_mul_pd(_const05, _parR);
						_valTemp = _mm256_add_pd(_valOrg, _parL05);
						_valRes = _mm256_add_pd(_valTemp, _parR05);
						_mm256_store_pd(val+start +offset*stride+poleLoop, _valRes);
					}
				offset += stepsize;
				}
				for (int poleLoop = 0; poleLoop < unroll; poleLoop+=4){
					_parL = _mm256_load_pd(val+start +offset*stride-parOffsetStrided+poleLoop);
					_parL05 = _mm256_mul_pd(_const05, _parL);
					_valOrg = _mm256_load_pd(val+start +offset*stride+poleLoop);
					_valRes = _mm256_add_pd(_valOrg, _parL05);
					_mm256_store_pd(val+start +offset*stride+poleLoop, _valRes);
				}
				steps = steps >> 1;
				offset = mypow2(l[dim]- ll )-1;
				parentOffset =  stepsize;
				stepsize = stepsize << 1;
		}	// end loop over levels
		// level = 2 seperate
		for (int poleLoop = 0; poleLoop < unroll; poleLoop += 4){
			// parRight is now used to store the root for both vertices of level 2
			_parR = _mm256_load_pd(val+start+(offset+parentOffset)*stride+poleLoop);
			_parR05 = _mm256_mul_pd(_const05, _parR);
			_valOrg = _mm256_load_pd(val+start +offset*stride+poleLoop); //checkmark2
			_valRes = _mm256_add_pd(_valOrg, _parR05);
			_mm256_store_pd(val+start +offset*stride+poleLoop, _valRes);
		}
		offset += stepsize;
		for (int poleLoop = 0; poleLoop < unroll; poleLoop += 4){
			_parL = _mm256_load_pd(val+start+(offset-parentOffset)*stride+poleLoop);
			_parL05 = _mm256_mul_pd(_const05, _parL);
			_valTemp = _mm256_load_pd(val+start +offset*stride+poleLoop);
			_valRes = _mm256_add_pd(_valTemp, _parL05);
			_mm256_store_pd(val+start +offset*stride+poleLoop, _valRes);
		} // end PoleLoop for level 2
}


double CombiGrid::hierarchizeOptimized(int blockSize){
	int dim;
	int start;
	int stride =1 ;
	int ndim;
	int nbrOfPoles;
	int jump;
	div_t divresult;

	std::chrono::time_point<std::chrono::system_clock> startChrono, endChrono;
	std::chrono::duration<double> elapsed_seconds;
	startChrono = std::chrono::system_clock::now();

	//	 dim 0
	ndim = n0Aligned;
	nbrOfPoles = arraySize /ndim;
#pragma omp parallel for schedule(static) firstprivate (ndim, nbrOfPoles, start)
		for (int kk = 0; kk < nbrOfPoles; kk++){
			start = kk *ndim;
			hierarchize1DUnoptimized(start, 1, ndim,0);
		}

	for (dim = 1; dim < d; dim ++){ // hierarchize d >=1
		stride *=ndim;
		ndim = n[dim];
		jump = stride*ndim;
		nbrOfPoles = arraySize/ndim;// do loop over first dim in 1d Parts
		#pragma omp parallel for schedule(static) firstprivate (jump, stride, divresult, ndim, nbrOfPoles, start, dim)// start1,start2,start3,start4)
		for (int nn = 0; nn< nbrOfPoles; nn+=blockSize){ // integer operations form bottleneck here -- nested loops are twice as slow
			divresult = div (nn,stride);
			start = divresult.quot*jump +divresult.rem;
			hierarchize1DOptimized(start, stride, ndim,dim, blockSize);
		}
	}

	endChrono = std::chrono::system_clock::now();
	elapsed_seconds = endChrono-startChrono;
	return 	elapsed_seconds.count();
}

