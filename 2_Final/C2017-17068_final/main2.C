#include <stdio.h>
#include <math.h>
#define _USE_MATH_DEFINES

#include "./header/function.h"
#include "./header/calculus.h"
#include "./header/linalgb.h"
#include "./header/monte.h"

#define eps 1.0e-30
#define a 1.26 // for quadrature method
#define R 3.5
#define _a 0.001 // for Montecarlo method
#define a_to_c true
#define d_to_e true
#define trial 3

double p(matrix v0){
	double *v = v0.mat;
	return v[0]*v[0]+v[0]*v[1]+v[1]*v[1];}
double j(matrix v0){
	double *v = v0.mat;
	return (v[0]*v[0]*v[0]*v[0])+(v[0]*v[0]*v[1]*v[1])+(v[1]*v[1]*v[1]*v[1]);}
double k(matrix v0){
	double *v = v0.mat;
	return pow(v[0]*v[0]+v[1]*v[1],0.5);}
double g(matrix v0){
	double y;

	if(fabs(k(v0))<1e-14){y = 0;}//since singularity 'on integral' is removable
	else{
	y= cos(j(v0)*cos(j(v0)))*sin(k(v0)*sin(k(v0)))/((j(v0)*j(v0)+1)*pow(k(v0),3.5));}
	
	return y;
}

double f_ftn(matrix v0,...){
	double y = g(v0)*exp(-1.0*p(v0));
	//double y = k(v0);
	return y;
}
double chi(matrix v0,...){ // characteristic function of 1/4 disk
	double y;
	if(k(v0)>=a){y = 1;}
	else{ y = 0; }
	return y;
}
double fp_ftn(matrix v0,...){
	if(k(v0)<eps){return 0;}
	return (1-chi(v0))*(f_ftn(v0) - 1/pow(k(v0),1.5));
}
double fq_ftn(matrix v0,...){
	return (chi(v0))*f_ftn(v0);
}
double f_weight_ftn(matrix v0,...){
	return exp(-1*p(v0))/(1+j(v0)*j(v0));
}

//Code testing function
double A_ftn(matrix v0,...){
	double *v = v0.mat;
	return exp(v[0])*sin(v[1]);
}
double A_weight(matrix v0,...){
	//double *v = v0.mat;
	return 1;
}


int main(int argc, char** argv){
	
	func_multi f = {f_ftn,2};
	func_multi f_p = {fp_ftn,2};
	func_multi f_q = {fq_ftn,2};

	double I,Ip,Iq;
	double dom_p[4] = {0,a,0,a};
	double dom_q[4] = {0,R,0,R};
	
	double h1,h2; h1 = 2e-4; h2 = 1.0e-3;
	double h_p[2] = {h1,h1};
	double h_q[2] = {h2,h2};


	// 2-(a),(b),(c)
if(a_to_c){
	printf("\nCalculuation strategy is on Supplementary.pdf.\n");
	printf("our h on domain 1 is %5.2e, on domain 2 is %5.2e\n",h1,h2);
	printf("our a = %5.2e, R = %5.2e.\n",a,R);
	
	printf("\n2-(a): Calculating....\n");
	Ip = quadr_2d(f_p,dom_p,h_p,"trep");
	Iq = quadr_2d(f_q,dom_q,h_q,"trep");
	I = Ip + Iq + M_PI*pow(a,.5);

	printf("Integral using trep rule is %013.11e + %013.11e = [%013.11e]\n",Ip,Iq,I);

	printf("\n2-(b): Calculating....\n");
	Ip = quadr_2d(f_p,dom_p,h_p,"simp");
	Iq = quadr_2d(f_q,dom_q,h_q,"simp");
	I = Ip+ Iq + M_PI*pow(a,.5);

	printf("Integral using simpson rule is %013.11e + %013.11e = [%013.11e]\n",Ip,Iq,I);
	
	printf("\n2-(c): they are quitely same, but not much for scale of 1e-10.(1e-5~1e-6)\n");
}

if(d_to_e){
	double dom_1[4] = {0,_a,0,_a};
	double dom_2[4] = {0,_a,_a,R};
	double dom_3[4] = {_a,R,0,_a};
	double dom_4[4] = {_a,R,_a,R};

	double I1,I2,I3,I4; 

	func_multi f_weight = {f_weight_ftn,2};

	printf("\n\n2-(d):integral in montecarlo method.(trial# = %d)\n\n",trial);

	for(int i=0;i<trial;i++){
	printf("Calculating....(%d / %d)\n",i+1,trial);
	monte_debug(0); I1 = metropolis(f,f_weight,dom_1,2e7); monte_debug(0);
	I2 = metropolis(f,f_weight,dom_2,2e7);
	I3 = metropolis(f,f_weight,dom_3,2e7);
	I4 = metropolis(f,f_weight,dom_4,7e7);
	
	I = I1+I2+I3+I4;
	printf("Integral is %13.11e + %13.11e + \n%13.11e + %13.11e = [%13.11e] (%d/%d)\n\n",I1,I2,I3,I4,I,i+1,trial);
	}

	printf("\n2-(e): Please see '2e3d_readme'\n\n");

}
	return 0;
}
