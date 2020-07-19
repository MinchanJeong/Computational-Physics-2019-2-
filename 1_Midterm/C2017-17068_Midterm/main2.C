#include <stdio.h>
#include <stdlib.h>

#include "common.h"

int main(int argc,char** argv){

	int n=1; // output = n!
	double result =0;

    printf("\n<Want to Calculate n!>\n\nWhat is 'n'?:");
	if(scanf("%d",&n)==1){}
	getchar();
	printf("\n");

	if(n<0){
		printf("\nI have no such option...\n\n");
		exit(2);
	}
	else if(n==0){
		printf("\n0! = 0 conventionally.\n\n");
		exit(1);
	}

	///print by full precision

	// initial setting for calculating factorial
	deg *bignum;
	bignum = (deg*)malloc(1000*sizeof(deg));
	*bignum = {1,1,NULL};

	/// now bignum represent n!.
	exfac(bignum,n); 

	// printing
	printf("%d! by %d-nary notation(exact precision):\n\n%d! = ",n,onemax,n);
	printdeg(bignum);

	printf("\n\n");

	/// print by double precision

	result = deg_to_double(bignum);

	printf("%d! by double representation:\n",n);
	printf("\n%d! = %lf\n\n",n,result);

	free(bignum);
	return 0;
}
