#include "CombiGrid.h"
#include "CombiGrid.cpp"
#include <sstream>

long static flops(int d, int * l){ // get number of flops performed (calculated by formula given in paper)
	  long flop = 0;
	  int i, j;
	  long temp;;
	  for (i = 0 ; i < d; i++){
		  temp = 1;
		  for (j = 0; j < d; j++){
			  if (j != i){
				  temp *= (1 << l[j])-1;
			  }
		  }
		  flop+=temp * ( ( 1 << (l[i]+1))-2*l[i]-2 );
	  }
	  return 2*flop;//fuer mults
}

string static composeFileName(string prefix, string method, int d, int lMin, int lMax, int procs, int alignment, int unroll, string schedule){
 		std::stringstream ss;
 	 	ss <<prefix<<method<<"D" << d << "L";
 	 	ss << lMin << "to" << lMax <<"P"<<procs<<"Alignment"<<alignment<<"Bytes-Unrolled"<<unroll<<"doublesSchedule"<<schedule<<".txt";
 	 	cout << ss.str() << endl;
 	 	return ss.str();
}



void static timeUnoptimized(int d, int * ll, int repetitions, FILE * output){
 	 		 double time, timeMin, timeMax, timeAcc;
 	 		 timeAcc = 0;
 	 		 timeMax = 0;
 	 		 timeMin = INFINITY;

 	 		 for (int reps = 0; reps < repetitions; reps++) {
				 CombiGrid * cc  = new CombiGrid (d,ll);
				 if (cc  == 0 ) exit(37);
				 cc->setValues(allOnes);

				 // warm up cache
				 cc->hierarchizeUnoptimized();
				 // time
				 time = cc->hierarchizeUnoptimized();
				 timeAcc = timeAcc+time;
				 if (time < timeMin){
					 timeMin = time;
				 }
				 if(time > timeMax){
					 timeMax = time;
				 }
				 delete cc;
 	 		 }
 	 		 timeAcc = timeAcc/repetitions;
 	 		fprintf(output, "%f %f %f\n", timeAcc,timeMin,timeMax);
 }

 	 void static timeOptimized(int d, int * ll, int alignment, int blockSize, int repetitions, FILE * output){
 	 	 		 double time, timeMin, timeMax, timeAcc;
 	 	 		 timeAcc = 0;
 	 	 		 timeMax = 0;
 	 	 		 timeMin = INFINITY;
 	 	 		 for (int reps = 0; reps < repetitions; reps++) {
					 CombiGrid * cc  = new CombiGrid (d,ll, alignment);
					 if (cc  == 0 ) exit(37);
					 cc->setValues(allOnes);
					 // warm up cache
					 cc->hierarchizeOptimized(blockSize);
					 // time
					 time = cc->hierarchizeOptimized(blockSize); // requires alignment
					 timeAcc = timeAcc+time;
					 if (time < timeMin){
						 timeMin = time;
					 }
					 if(time > timeMax){
						 timeMax = time;
					 }
					 delete cc;
 	 	 		 }
 	 	 		 timeAcc = timeAcc/repetitions;
 	 	 		fprintf(output, "%f %f %f\n", timeAcc,timeMin,timeMax);
 	 }


 	 void static varyDimAndLevelMain(int minDim, int maxDim,   int * minLevels, int * maxLevels, int * blockSize, int alignment, int procs, int reps, string prefix) {
 		 omp_set_num_threads(procs);
 		 omp_set_schedule(omp_sched_static,0);
 		 string scheduleStr = "ScheduleStatic0";

 		int blockSizeTmp;
 		FILE * outputFile;
 		string str;

 		 for (int dim = minDim; dim <= maxDim; dim++ ){
 			 str = composeFileName(prefix+"Increaselevel", "Optimized", dim, minLevels[dim-1], maxLevels[dim-1], procs, alignment, blockSize[dim-1], scheduleStr);
 			 outputFile = fopen((char*)str.c_str(), "w");
 			 if (outputFile == NULL) {
 			 cerr <<  "I couldn't open file for writing.\n";
 			    exit(0);
 			 }
 	 		 int * levels = new int[dim];
 			 for (int lCtr = minLevels[dim-1]; lCtr <= maxLevels[dim-1]; lCtr++){
				 // set level vector for combiGrid
 				 for (int i =0; i < dim; i++){
					 levels[i] = lCtr;
				 }
 				blockSizeTmp = blockSize[dim-1];
 	 			 if (blockSizeTmp > (1<<lCtr)){ // blockSize larger than number of grid points in first dimension
 	 				 blockSizeTmp = (1 << lCtr);
 	 			 }
 				 timeOptimized(dim,levels,alignment, blockSizeTmp, reps, outputFile);
 				 fflush(outputFile);
 			 }
 			 fclose(outputFile);
 		 }
 	 }

 	 void static scalingMain(int minDim, int maxDim, int * levelsDim, int minProcs, int maxProcs, int * blockSize, int alignment, int reps,  string prefix) { // hierarchize only for d>=2 parallelized
		 omp_set_schedule(omp_sched_static,0);
		 string scheduleStr = "Static0";

		 int blockSizeTmp;
		 FILE * outputFile;
		 string str;


 		// scaling for optimized version
 		 for (int dim = minDim; dim <= maxDim; dim++ ){
 			 int lCtr = levelsDim[dim-1];
 			 str = composeFileName(prefix+"ScalingP"+std::to_string(minProcs)+"-"+std::to_string(maxProcs), "Optimized", dim, lCtr, lCtr, 0, alignment, blockSize[dim-1], scheduleStr);
 			 outputFile = fopen((char*)str.c_str(), "w");
 			 if (outputFile == NULL) {
 			 cerr <<  "I couldn't open file for writing.\n";
 			    exit(0);
 			 }
			 // set level vector for combiGrid
 	 		 int * levels = new int[dim];
 			 for (int i =0; i < dim; i++){
				 levels[i] = lCtr;
			 }
 			blockSizeTmp = blockSize[dim-1];
 	 		 if (blockSize[dim-1]> (1<<lCtr)){ // blockSize larger than number of grid points in first dimension
 	 			blockSizeTmp = (1 << lCtr);
 	 		 }
 	 		 for (int p= minProcs; p <= maxProcs; p++){
				 omp_set_num_threads(p);
				 timeOptimized(dim,levels,alignment, blockSizeTmp, reps, outputFile);
				 fflush(outputFile);
 	 		 	 }
 	 		 fclose(outputFile);
 		 }
 		 // scaling for unoptimized version
 		 for (int dim = minDim; dim <= maxDim; dim++ ){
 			 int lCtr = levelsDim[dim-1];
 			 str = composeFileName(prefix+"ScalingP"+std::to_string(minProcs)+"-"+std::to_string(maxProcs), "Unoptimized", dim, lCtr, lCtr, 0, alignment, 0, scheduleStr);
 			 outputFile = fopen((char*)str.c_str(), "w");
 			 if (outputFile == NULL) {
 			 cerr <<  "I couldn't open file for writing.\n";
 			    exit(0);
 			 }
		// set level vector for combiGrid
 	 		 int * levels = new int[dim];
 			 for (int i =0; i < dim; i++){
				 levels[i] = lCtr;
			 }
 	 		 for (int p= minProcs; p <= maxProcs; p++){
				 omp_set_num_threads(p);
				 timeUnoptimized(dim,levels, reps, outputFile);
				 fflush(outputFile);
 	 		 	 }
 	 		 fclose(outputFile);
 		 }
 	 }

 	void static anisotropicMain(int dd, int avgL, int stdL, int maxL, int procs, int maxUnroll, int alignment, int reps, string prefix) {
 	 		 omp_set_schedule(omp_sched_static,0);
 	 		 string scheduleStr = "static0";
 	 		 omp_set_num_threads(procs);

 	 		 int blockSizeTmp;
 	 		 int * levels = new int[dd];
 	 		 FILE * outputFile;
 	 		 string str;
 			 str = composeFileName(prefix+"AnisotropicD"+std::to_string(dd), "Optimized", 0, stdL, maxL, procs, alignment, maxUnroll, scheduleStr);
 			 outputFile = fopen((char*)str.c_str(), "w");
 			 if (outputFile == NULL) {
 			 cerr <<  "I couldn't open file for writing.\n";
 			    exit(0);
 			 }

 			 // base case of isotropic, regular refined grid.
  	 			 for (int i =0; i < dd; i++){
  					 levels[i] = avgL;
  				 }
  	 	 		 blockSizeTmp = maxUnroll;
  	 	 		 if (blockSizeTmp > (1<<levels[0])){ // blockSize larger than number of grid points in first dimension
  	 	 			blockSizeTmp = (1 << levels[0]);
  	 	 		 }
  	 	 		 timeOptimized(dd,levels,alignment, blockSizeTmp, reps, outputFile);
  				 fflush(outputFile);

 	 		// ansiotropic grids
 	 		for (int dim = 0; dim < dd; dim ++){
 		  // set level vector for anisotropic grid
 	 			 for (int i =0; i < dd; i++){
 					 levels[i] = stdL;
 				 }
 	 			 levels[dim] = maxL;
  	 	 		 blockSizeTmp = maxUnroll;
  	 	 		 if (blockSizeTmp > (1<<levels[0])){ // blockSize larger than number of grid points in first dimension
  	 	 			blockSizeTmp = (1 << levels[0]);
  	 	 		 }
 	 	 		 timeOptimized(dd,levels,alignment, blockSizeTmp, reps, outputFile);
 				 fflush(outputFile);
 	 		}
 	 	 	fclose(outputFile);
 	 	 }


int main (int argc, char *argv[] ) {
	string prefix;
    if (argc == 1){
		prefix = "";
 	}
	else {
		prefix = argv[1];
	}

	int alignment = 32; // in bytes
	int repititions = 10;
	int procs = 3; // 3 processor fixed as these achieved best result

	int minDim =1; // if minDim is set >1 then the blockSize, minLevels and maxLevels Arrays need to have dummy values in the first (minDim-1) positions
	int maxDim = 5;
	int blockSize[5] = {1,16384,512,128,64};
	int minLevels[5] = {3,3,3,3,3};
	// make sure main memory is large enough to hold those grids
//	int maxLevels[5] = {30,15,10,7,6};
	int maxLevels[5] = {20,10,7,5,4};
	/* ********************************************************************************************************************
	 * Test Case 1: Incresing level for different dimensions.
	********************************************************************************************************************* */
	varyDimAndLevelMain(minDim, maxDim, minLevels, maxLevels, blockSize, alignment, procs, repititions, prefix);


	if (minDim <2) { // code is only parallelized for d >= 2
		minDim =2;
	}
	int minProcs =1;
	int maxProcs =4;
	/* ********************************************************************************************************************
	 * Test Case 2: Scaling for the optimized and unoptimized code.
	********************************************************************************************************************* */
	scalingMain(minDim, maxDim, maxLevels, minProcs, maxProcs, blockSize, alignment, repititions, prefix);



	int dd = 6;
	int stdL = 3;
	// make sure main memory is large enough to hold those grids
//	int avgL = 5;
//	int maxL = 15;
	int avgL = 4;
	int maxL = 8;
    int maxUnroll = 16384;
	/* ********************************************************************************************************************
	 * Test Case 3: Comparing the isotropic baseline with anisotropic grids.
	********************************************************************************************************************* */
	anisotropicMain(dd, avgL, stdL, maxL, procs, maxUnroll, alignment, repititions, prefix);


	return 0;
}
