#include <stdio.h>
#include <math.h>

#include "common.h"

int indexmax = 1000;

float avg(int *xi,float *f){
	int j=0;
	int size = 0;
	float sum = 0;

	for(j=0;j<indexmax;j++){
		if(xi[j]!=-1){
			sum += f[xi[j]];
			size++;
		}
		else{
			break;
		}
	}
	return (sum/size);
}

float err(int size,int *xi,float *f){
	int j=0;
	float sqsum = 0;
	float m = avg(xi,f);

	for(j=0;j<size;j++){
		sqsum += ((f[xi[j]] - m) * (f[xi[j]] - m));
	}

	return (sqsum / ((size * (size-1))));
}

float mean(int size,float *a){
	float sum = 0;

	int i=0;
	for(i=0;i<size;i++){
		sum += a[i];
	}
	return (sum/size);
}


