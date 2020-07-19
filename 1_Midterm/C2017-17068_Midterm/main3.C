#include <stdio.h>
#include <stdlib.h>

#include "common.h"

int main(int argc, char** argv){

	int i = 1;
	int j = 1;

	//# 3-(a)
	double B2n_1[11]={1,0,};

	for(i=1;i<11;i++){

		B2n_1[i] = i - 0.5;
		for(j=1;j<i;j++){

			B2n_1[i] -= dbbi(2*i+1,2*j)*B2n_1[j];
		}
		B2n_1[i] = B2n_1[i] / (double)(2*i+1);
	}

	printf("\n# 3 - (a)\n\n");
	for(i=0;i<11;i++){
		printf("B_%02d = %lf\n",2*i,B2n_1[i]);
	}
	printf("\n---------------\n");

	//# 3-(b) (I used N->N+1 version of (5))
	double B2n_2[11]={1,0,};

	for(i=1;i<11;i++){

		B2n_2[i] = i;
		for(j=1;j<i;j++){

			B2n_2[i] -= dbbi(2*i+2,2*j)*B2n_2[j];
		}
		B2n_2[i] = B2n_2[i] / dbbi(2*i+2,2*i);
	}

	printf("\n# 3 - (b)\n\n");
	for(i=0;i<11;i++){
		printf("B_%02d = %lf\n",2*i,B2n_2[i]);
	}
	printf("\n---------------\n");

	// #3-(c)
	//Since I will use equation '2*(2N-1)!!*(eq(4))' for convention, I will calculate 2*(B_2n) first. Thus, below B2n array also represent twice of B-number for a while.

	quot B2n[11];
	
	//memory allocating
	for(i=0;i<11;i++){
		B2n[i].n = (deg*)malloc(50*sizeof(deg));
		B2n[i].d = (deg*)malloc(50*sizeof(deg));

		*(B2n[i].n) = {1,1,NULL};
		*(B2n[i].d) = {1,1,NULL};
	}

	//I will initially define denominator of B(2n) by (2n+1)!!.
	//Then Below calculation is based on (eq(4))*2*(2N-1)!!, which is integer equality by simple induction.

	//Define denominator(n>=1)
	for(i=1;i<11;i++){
		exfacfac(B2n[i].d,2*i+1);
	}
	
	//Calculate B2n[i].n, using (eq(4))*2*(2N-1)!!.
	for(i=1;i<11;i++){

		exfacfac(B2n[i].n,2*i-1);
		deg_x_uint(B2n[i].n,2*i-1,0);
		for(j=1;j<i;j++){

			deg *dumy1;
			dumy1 = (deg*)malloc(50*sizeof(deg));
			*dumy1 = {1,0,NULL}; // dummy1 = 0
			copy(dumy1,B2n[j].n);

			bi(dumy1,2*i+1,2*j);
			coef(dumy1,i,j);

			deg_m_deg(B2n[i].n,dumy1);

			free(dumy1);
		}
		// now make results half by twices denominator
		deg_x_uint(B2n[i].d,2,0);
	}

	printf("\n# 3 - (c)\n\n");

	for(i=0;i<11;i++){
		printf("B_%02d = ",2*i);
		printdeg(B2n[i].n);
		printf("/ ");
		printdeg(B2n[i].d);
		printf("\n");
	}
	printf("\n---------------\n");

	// #3-(d)
	// as I noted on 3-(c), denominator is (2N+1)!! (N: 1 ~ 10). Thus we only have to check 2,3,5,7,9,11,13,17,19.(prime factorization).
	int prime[8]={2,3,5,7,11,13,17,19};
	
	for(i=0;i<11;i++){
	for(j=0;j<8;j++){
		while(isdivd(B2n[i].n,B2n[i].d,prime[j])==0){
		}
	}
	}

	printf("\n# 3 - (d)\n\n");
	for(i=0;i<11;i++){
		printf("B_%02d = ",2*i);
		printdeg(B2n[i].n);
		printf("/ ");
		printdeg(B2n[i].d);
		printf("\n");
	}

	//freeing memory
	for(i=0;i<11;i++){
		free(B2n[i].n); 	
		free(B2n[i].d);	
	}

	printf("\n------------------\n\n");
}

