#include <stdio.h>
#include <math.h>
#include <string.h>

#include "function.h"
#include "calculus.h"
#include "linalgb.h"

/* Integration */

int coef_trep[3] = {1,2,1};
int coef_simp[3] = {1,4,1};
int coef_booles[5] = {14,64,24,64,14};
numrcoefs trep = {coef_trep,2,3};
numrcoefs simp = {coef_simp,3,3};
numrcoefs booles = {coef_booles,45,5};

numrcoefs qdr_str_to_coefs(const char* option){

	int p = 10;
	if(strcmp(option,"trep") == 0){p = 0;}
	else if(strcmp(option, "simp") == 0){p = 1;}
	else if(strcmp(option, "booles") == 0){p = 2;}

	numrcoefs co;
	switch(p){
		case 0: co = trep; break;
		case 1: co = simp; break;
		case 2: co = booles; break;
		default: fprintf(stderr,"\nOption not found!!\n"); exit(3); break;
	}

	return co;
}

double qdr_calc_base(double* f_arr, int N, double h, numrcoefs co){

	int *c = co.coefs;
	int d = co.divider;
	int gap = co.size-1;

	int i,j;
	double I = 0;
	double dum = 0;
	for(i = 0 ; i+gap < N ; i += gap){

		dum = 0;
		for(j = 0 ; j <= gap ; j++){
			dum += f_arr[i+j] * c[j];
		}
		I += dum / d;
	}

	I = I * h;
	return I;

}

double qdr_calc_1d(func_1d ftn, numrcoefs co, double a, double b, double h, int N){
	
	double (*f)(double,...) = ftn.fpart;

	if ( a == b ){ return 0;}
	else if ( fabs(b-a) < h ){ fprintf(stderr,"\n interval error!!\n"); exit(1); }
	else if ( h <= 0 ){ fprintf(stderr,"\n h value error!! \n"); exit(2); }

	double f_arr[N] = {0,};
	int i = 0;
	for(i=0;i<N;i++){
		f_arr[i] = f(a + i*h);
	}
	
	return qdr_calc_base(f_arr, N, h, co);
}

double quadr(func_1d ftn, double a, double b, double h, const char *option){
	
	int N = abs(b - a) / h + 1;
	numrcoefs co = qdr_str_to_coefs(option);

	return qdr_calc_1d(ftn,co,a,b,h,N);
}

/*
double quadr_N(func_1d ftn, double a, double b, int N, const char *option){
	
	if(N <= 1){ fprintf(stderr, "\nN value (>1) error!!\n");}
	double h = abs(b - a) / (N-1);
	numrcoefs co = qdr_str_to_coefs(option);

	return qdr_calc_1d(ftn,co,a,b,N);
}
*/


double quadr_2d(func_multi ftn, double* dom, double* h, const char* option){
	//dom[4] = {a,b,c,d}. we consider [a,b]x[c,d].
	//h[2] = {hx,hy}
	numrcoefs co = qdr_str_to_coefs(option);
	double (*f)(matrix,...) = ftn.fpart;

	double a1,b1,a2,b2,h1,h2;
	a1= dom[0]; b1 = dom[1]; a2 = dom[2]; b2 = dom[3];
	h1 = h[0]; h2 = h[1];

	int N1 = abs(b1-a1) / h1 + 1;
	int N2 = abs(b2-a2) / h2 + 1;
	
	matrix v; v.cdim = 2; v.rdim =1;
	v.mat = (double*)malloc(sizeof(double)*2);
	double *fx_y0_arr = (double*)malloc(sizeof(double)*N1);
	double *fy_arr = (double*)malloc(sizeof(double)*N2);
	
	int i,j;
	for(j=0;j<N2;j++){

		v.mat[1] = a2 + j*h2;

		for(i=0;i<N1;i++){
			v.mat[0] = a1 + i*h1;
			fx_y0_arr[i] = f(v);
		}
		fy_arr[j] = qdr_calc_base(fx_y0_arr,N1,h1,co);
	}

	double I = qdr_calc_base(fy_arr,N2,h2,co);

	free(fx_y0_arr);
	free(fy_arr);
	free(v.mat);
	return I;
}

/* Diffrentiation*/

int coef_sq[1] = {1}; //{-1,1}
int coef_qd[2] = {8,-1}; // {1,-8,8,-1}
int coef_hx[3] = {45,-9,1};
numrcoefs hsq = {coef_sq,2,1};
numrcoefs hqd = {coef_qd,12,2};
numrcoefs hhx = {coef_hx,60,3};

numrcoefs diff_str_to_coefs(const char* option){

	int p = 10;
	if(strcmp(option,"h^2") == 0){p = 0;}
	else if(strcmp(option, "h^4") == 0){p = 1;}
	else if(strcmp(option, "h^6") == 0){p = 2;}

	numrcoefs co;
	switch(p){
		case 0: co = hsq; break;
		case 1: co = hqd; break;
		case 2: co = hhx; break;
		default: fprintf(stderr,"\nOption not found!!\n"); exit(3); break;
	}

	return co;
}

double diff_calc_base(double* f_arr,double h, numrcoefs co){

	int *c = co.coefs;
	int d = co.divider;
	int N = co.size;

	double Df=0;
	double _h=h/N;

	int i=1;
	for(i=1;i<=N;i++){
		Df += ((f_arr[2*N+1 - i] - f_arr[i-1]))*c[N-i];
	}

	Df = Df / (_h * d);
	return Df;
}

double diff_calc_1d(func_1d ftn,double x0, double h, numrcoefs co){
	
	double (*f)(double,...) = ftn.fpart;
	int N = co.size;
	double _h = h / N;
	
	double f_arr[2*N+1];

	int i;
	for(i=0;i<=2*N;i++){
		f_arr[i] = f(x0 + (i - N ) *_h);
	}

	return diff_calc_base(f_arr,h, co);

}

double diff_1d(func_1d ftn,double x0,double h, const char* option){

	numrcoefs co = diff_str_to_coefs(option);

	return diff_calc_1d(ftn, x0, h, co);
}

double diff_calc_multi(func_multi ftn, matrix v0, int k, double h,  numrcoefs co){

	double (*f)(matrix,...) = ftn.fpart;

	matrix v;
	v.cdim = v0.cdim;
	v.rdim = 1; 
	v.mat = (double*)malloc(sizeof(double)*(v.cdim));
	matrix_copy(&v,v0);

	int l = co.size;
	double _h = h / l;
	
	double f_arr[2*l+1];
	
	int i;
	for(i=0;i<=2*l;i++){
		v.mat[k-1] = v0.mat[k-1] + (i - l)*_h;
		f_arr[i] = f(v);
	}
	
	free(v.mat);
	return diff_calc_base(f_arr,h, co);
}

double diff_multi(func_multi ftn, matrix v0, int k, double h, const char* option){

	numrcoefs co = diff_str_to_coefs(option);

	return diff_calc_multi(ftn, v0, k, h, co);

}

matrix* evaluate(func_multi* vectorftn, matrix* rtn_vec){
	
	int N = rtn_vec->cdim;

	matrix* dummy = new_0matrix(N,1);

	for(int i=0;i<N;i++){
	dummy->mat[i] = vectorftn[i].fpart(*rtn_vec);
	}
	matrix_copy(rtn_vec,*dummy);
	
	freematrix(dummy);
	return rtn_vec;
}

matrix* Jacobian(matrix* rtn,func_multi* vectorftn ,int N,matrix v0,double *h,const char* option){
	
	int i,j;
	
	for(i=0;i<N;i++){
		
		if(vectorftn[i].dim != N){
			fprintf(stderr,"\nJACOBIAN: FUNCTION DIMMENSION ERROR!!\n");
			return NULL;
		}
	}

	
	if(rtn != NULL){
	if(rtn->cdim != N || rtn->rdim != N){
		fprintf(stderr,"\njacobian: return matrix size error!!\n");
		return NULL;
	}}

	matrix *dummy = new_0matrix(N,N);

	for(i=0;i<N;i++){
	for(j=0;j<N;j++){

		dummy->mat[i+N*j] = diff_multi(vectorftn[i],v0,j+1,h[j],option);
		
	}}
	if(rtn != NULL){
	matrix_copy(rtn,*dummy);
	freematrix(dummy);
	return rtn;
	}
	else{
	write_mat_allocdlist(dummy);
	return dummy;
	}

}
