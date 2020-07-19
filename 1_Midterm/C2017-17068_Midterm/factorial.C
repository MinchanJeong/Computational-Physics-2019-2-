#include <stdio.h>
#include <stdlib.h>

#include "common.h"

// exfac: reconstruct 'deg' to represent number multiplied with n!.
int exfac(deg *no, int n)
{
	if (n<=0){
		return 0;}
	int i=1;
	for(i=1;i<=n;i++){
		deg_x_uint(no,i,0);
	}
	return 0;
}

int exfacfac(deg *no, int n)
{
	if (n<=0){
		return 0;}
	int i=1;
	for(i=n;i>0;i -= 2){
		deg_x_uint(no,i,0);
	}
	return 0;
}

//like exfac(), bi()  convert array 'no' to  represent 'no * (nCr)'
int bi(deg *no,int n,int r){

	if((n <= r) || (r == 0))
	{return 1;}

	int i=0;

	exfac(no,n);

	for(i=1;i<=(n-r);i++){
	deg_d_uint(no,i);
	}

	for(i=1;i<=r;i++){
	deg_d_uint(no,i);
	}

	return 0;
}

int coef(deg *no,int N,int n){
	
	if(N <= n){
		return 1;}

	int i=0;

	for(i=(n+1);i<N;i++){

	deg_x_uint(no,2*i+1,0);
	}
	return 0;
}

//reduction of fraction
int isdivd(deg *a,deg *b,int m){

	///dummy declaration
	deg *dum_a;
	deg *dum_b;
	dum_a = (deg*)malloc(100*sizeof(deg));
	dum_b = (deg*)malloc(100*sizeof(deg));
	*dum_a = {1,1,NULL};
	*dum_b = {1,1,NULL};

	int rtn = 0;

	copy(dum_a,a);
	copy(dum_b,b);

	int ra = deg_d_uint(dum_a,m);
	int rb = deg_d_uint(dum_b,m);

	if((ra==0)&&(rb==0)){
		copy(a,dum_a);
		copy(b,dum_b);
		rtn = 0;
	}
	else{
		rtn = 1;
	}
	
	free(dum_a);
	free(dum_b);
	return rtn;

}

// fac: print n! by double with inner dynamic memory allocation using exfac.
double fac(int n){

	if (n<=0){
		return 0;}

	deg *largenum;
	largenum = (deg*)malloc(1000*sizeof(deg));
	*largenum = {1,1,NULL};

	exfac(largenum,n);
	
	double output = deg_to_double(largenum);
		
	free(largenum);
	return output;
}

//binomial coeff using fac()
double dbbi(int m,int n){

	if(m<n)
	{return 0;}
	
	double output = 1;

	output = fac(m) / (fac(n) * fac(m-n)) ;

	return output;
	
}
