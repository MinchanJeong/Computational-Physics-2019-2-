#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "function.h"
#include "random.h"
#include "linalgb.h"
#include "monte.h"

int monte_debug_flag = 0;

int monte_debug(int flag){ monte_debug_flag = flag; return 0;}

int init_monte_random(){

	unsigned long long init[4]  = {(unsigned long long)time(NULL), 0x23456ULL, 0x34567ULL, 0x45678ULL};
	int length = 4;
	init_by_array64(init, length);
	return 0;
}

int set_rand_vector(matrix* M, double* domain){
	
	int N = M->cdim;
	int i;
	for(i=0;i<N;i++){
		M->mat[i]=genrand64_real1()*(domain[2*i+1]-domain[2*i])+domain[2*i];
	}
	return 0;
}

double montecarlo(func_multi ftn, double *domain, int sampleno){
	
	int N = ftn.dim;
	double (*f)(matrix,...) = ftn.fpart;
	
	matrix *dummy = new_0matrix(N,1);
	double *values = (double*)malloc(sizeof(double)*sampleno);
	
	//initial t0
	set_rand_vector(dummy,domain);
	values[0] = f(*dummy);

	//make sequence
	int piv = 0;
	double d1,d2,p;
	while(piv<(sampleno-1)){

		set_rand_vector(dummy,domain);
			values[piv+1] = f(*dummy);
		piv++;
	}

	// summation all f(ti) and makes integral
	double I=0;
	for(int i=0;i<sampleno;i++){
		I += values[i];
	}
	I = I/sampleno;
	
	//now multiply volume
	double V=1;
	for(int i=0;i<N;i++){
		V *= (-domain[2*i]+domain[2*i+1]);
	}

	I *= V;
	
	//freeing allocated memory
	freematrix(dummy);
	free(values);
	return I;
}
double metropolis(func_multi ftn,func_multi prob_w, double *domain, int sampleno){
	
	int N = ftn.dim;
	double (*f)(matrix,...) = ftn.fpart;
	double (*w)(matrix,...) = prob_w.fpart;
	
	matrix *dummy = new_0matrix(N,1);
	double *weights = (double*)malloc(sizeof(double)*sampleno);
	double *values = (double*)malloc(sizeof(double)*sampleno);
	
	//initial t0
	set_rand_vector(dummy,domain);
	weights[0] = w(*dummy);
	values[0] = f(*dummy);

	//make sequence
	int piv = 0;
	bool sign = true;
	double d1,d2,p;
	while(piv<(sampleno-1)){

		sign = false; d1 = weights[piv];
		set_rand_vector(dummy,domain);
		d2 = w(*dummy);

		if(d1<0){
			fprintf(stderr,"\nMTP: WEIGHTS SIGN ERROR!!\n");
			break;
		}
		else if(d1 < 1e-300){sign = true;}
		else if((d2/d1 - 1.0) >= -1e-30){sign = true;}
		else if(d2/d1 < 1.0){
			p = genrand64_real1();
			if(p < d2/d1){sign = true;}
		}

		if(sign){
			weights[piv+1] = d2;
			values[piv+1] = f(*dummy);
		}
		else{
			weights[piv+1] = weights[piv]; // = d1
			values[piv+1] = values[piv];
		}
		piv++;

		if(monte_debug_flag==1){
			matrix_print(*dummy,13,11,true);
			printf("%11.13f\n",weights[piv]);
			getchar();
		}

	}

	// summation all f(ti) and makes integral
	double I=0;
	for(int i=0;i<sampleno;i++){
		I += values[i]/weights[i];

		if(monte_debug_flag==2){
				printf("%11.13f\n",I);
				getchar();
		}

	}
	I = I/sampleno;
	
	I = I * montecarlo(prob_w,domain,sampleno); // normalize weight
	
	//freeing allocated memory
	freematrix(dummy);
	free(values);
	free(weights);
	return I;
}

