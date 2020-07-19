#include <stdio.h>
#include <math.h>

#define _USE_MATH_DEFINES

#include "./header/function.h"
#include "./header/rootfind.h"

double g(double x,...){
	return 10.0*x/(1.0+x*x);}
double h(double x,...){
	return 10.0*x*x*exp(-1.0*x*x)/(1.0+x*x*x*x);}
double f(double x,...){
	return sin(g(x)*sin(h(x)));}

int main(int argc, char** argv){
	
	double ansatz[4] = {0.35,0.6,0.8,1.3};
	double sol[5] = {0,};
	func_1d ftn = {f};

	printf("\nI developed simple proof of 'there is no solution larger then 2' on supplementary.\nThus I gave 4 ansatz from observation.(please see obs.png(command: make plot) for justfication.)\n\n");
	
	sol[0] = 0.0; // there is trivial solution zero.

	for(int i=0;i<4;i++){
		sol[i+1] = NR_1d(ftn,ansatz[i],0.00001,1e-11);
	}
	for(int i=0;i<5;i++){
		printf("root%d: %012.11e\n",i+1,sol[i]);
	}

	printf("\nSince there is No solution larger then 2,\nthere is only 5 solution on 0 to inf\n\n");
	return 0;
}
