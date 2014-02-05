#include <math.h>

static double allOnes (double * x){
	  return 1;
 }

static double testFunc2D (double * x){
	  return 1+sin(x[0]*x[1])+cos(x[0]*x[0])+x[1]*x[1]*x[1];
 }

static double testFunc1D (double * x){
	  return 1+x[0]-fabs(x[0] -0.5)+3*x[0]*(x[0]-0.3);
 }

static double testFunc5D (double * x){
	  return 1+sin(x[0]*x[1])*cos(x[2]+x[3]*x[4])*sin(4*x[3])+cos(x[0]*x[0]+x[2])*sin(x[4])*sin(x[4])*cos(2*x[2]+x[0])*sin(x[3])*cos(x[4]);
 }

static double testFunc4D (double * x){
	  return 1+8*(x[0]*x[1])*cos(x[2]+x[3])*(4*x[3])+cos(x[0]*x[0]+x[2])*sin(x[3])*sin(x[3])*cos(2*x[2]+x[0])*sin(x[3]);
 }


 static double parabel1D (double * x){
	  return 4* x[0]*(1-x[0]);
 }
 static double parabel2D (double * x){
	  return 4*4* x[0]*(1-x[0])*x[1]*(1-x[1]);
 }

 static double parabel3D (double * x){
	  return 64* x[0]*(1-x[0])*x[1]*(1-x[1])*x[2]*(1-x[2]);
 }

 static double parabel4D (double * x){
	  return 64*4* x[0]*(1-x[0])*x[1]*(1-x[1])*x[2]*(1-x[2])*x[3]*(1-x[3]);
 }




