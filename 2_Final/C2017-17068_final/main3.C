#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#define _USE_MATH_DEFINES

#include "./header/function.h"
#include "./header/calculus.h"
#include "./header/linalgb.h"
#include "./header/monte.h"
#include "./header/bessel.h"

#define eps 1e-30
#define R M_PI
#define a_to_b true
#define c_to_d true
#define trial_a 10

double B(matrix k_mat){
	double *k = k_mat.mat;
	if(fabs(k[0])<eps&&fabs(k[1])<eps&&fabs(k[2]<eps)&&fabs(k[3]<eps)){return 0;}
	return  0.25/(sin(k[0]/2)*sin(k[0]/2)+sin(k[1]/2)*sin(k[1]/2)
				+sin(k[2]/2)*sin(k[2]/2)+sin(k[3]/2)*sin(k[3]/2));
}

double real_intgrand(matrix k_mat,int *q_no){
	double *k = k_mat.mat;
	double y = pow(2*M_PI,-4)*cos(k[0]*q_no[0]+k[1]*q_no[1]+k[2]*q_no[2]+k[3]*q_no[3])*B(k_mat);
	return y;
}
double imag_intgrand(matrix k_mat,int *q_no){
	double *k = k_mat.mat;
	double y = pow(2*M_PI,-4)*sin(k[0]*q_no[0]+k[1]*q_no[1]+k[2]*q_no[2]+k[3]*q_no[3])*B(k_mat);
	return y;
}

double weightftn(matrix k_mat,...){
	double *k = k_mat.mat;
	double r2 =  k[0]*k[0]+k[1]*k[1]+k[2]*k[2]+k[3]*k[3];
	return 1/(1e-6+r2);
}

double reducedftn(double x,...){
	
	if(fabs(x)<eps){return 0;}
	va_list ap;
	va_start(ap,x);
	int pivot;
	
	double y=0;
	int n[4] = {0,};
	for(int i=0;i<4;i++){
		 pivot = va_arg(ap,int);
		 n[i]=pivot;
	}
	va_end(ap);

	y = exp(-8.0*x)*bessel_I(n[0],2*x)*bessel_I(n[1],2*x)*bessel_I(n[2],2*x)*bessel_I(n[3],2*x);
	return y;
}

int mu[4] = {0,0,0,0};
int set_mu_null(int *arr){
	for(int i=0;i<4;i++){arr[i]=0;}
	return 0;
}

int main(int argc, char** argv){
	
	init_monte_random();
	int i,j;
	double dum1, dum2; dum1=0; dum2=0;
	
if(a_to_b){
	
	// 3 - (a)

	printf("\n3-(a) integral using montecarlo method.(trial# = %d)\n\n",trial_a);

	double domain[8] = {-R,R,-R,R,-R,R,-R,R};
	
	func_multi real; real.dim = 4;
	func_multi imag; imag.dim = 4;

	double r_value[5][trial_a] = {0,};
	double i_value[5][trial_a] = {0,};

	double r_avg[5] = {0,};
	double i_avg[5] = {0,};

	int n = 2e6;
	set_mu_null(mu);
	for(i=0;i<5;i++){
	
		dum1=0; dum2=0;

		if(i>0){mu[i-1] = 1;}
		double Z_real(matrix,...); real.fpart = Z_real;
		double Z_imag(matrix,...); imag.fpart = Z_imag;

		for(j=0;j<trial_a;j++){
		r_value[i][j] = montecarlo(real,domain,n);
		i_value[i][j] = montecarlo(imag,domain,n);
	
		printf("Z%d%d%d%d %dth trial; %13.11e + %013.11e j\n",
				mu[0],mu[1],mu[2],mu[3],j+1,r_value[i][j],i_value[i][j]);
		dum1 += r_value[i][j]; dum2 += i_value[i][j];
		}
		r_avg[i] = dum1/trial_a; i_avg[i] = dum2/trial_a;
		printf("Z%d%d%d%d MC method average; %13.11e + %013.11e j\n",
				mu[0],mu[1],mu[2],mu[3],r_avg[i],i_avg[i]);
		printf("\n");
	}

}

if(c_to_d){
	// 3 - (c)
	printf("\n3-(c) Integral using quadratrue metnod\n\n");

	func_1d reduced;
	
	double value[5] = {0,};
	double trepvalue[5] = {0,};

	double h1 = 1.1e-6; double h2 = 1.6e-5; double h3 = 4e-5;
	set_mu_null(mu);
	for(i=0;i<5;i++){

		if(i>0){mu[i-1] = 1;}
		double Z_reduced(double,...); reduced.fpart = Z_reduced;

		trepvalue[i] = quadr(reduced,0,1,h1,"trep") +quadr(reduced,1,15,h2,"trep") + quadr(reduced,15,44.1,h3,"trep");
		printf("Z%d%d%d%d in trep method   ; %13.11e.\n",
				mu[0],mu[1],mu[2],mu[3],trepvalue[i]);
		value[i] = quadr(reduced,0,1,h1,"simp") + quadr(reduced,1,15,h2,"simp") +quadr(reduced,15,44.1,h3,"simp");
		printf("Z%d%d%d%d in simpson method; %13.11e.\n",
				mu[0],mu[1],mu[2],mu[3],value[i]);
		printf("\n");
	}

	printf("\n3-(d): Please see '2e3d_readme'\n\n");
}

	return 0;
}
double Z_real(matrix M,...){ return real_intgrand(M,mu);}
double Z_imag(matrix M,...){ return imag_intgrand(M,mu);}
double Z_reduced(double x,...){ return reducedftn(x,mu[0],mu[1],mu[2],mu[3]);}

