#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "common.h"

int onemax = 10000; // maximum number+1 in one array

int lastlistno(deg *no){
	int i=0;
	while(((no+i)->upper) != NULL){
		i++;
	}

	return i;
}

//copy
int copy(deg *obj , deg *stuf){
	
	int a = lastlistno(stuf);
	int b = lastlistno(obj); 
	// obj will be new duplicated.

	// make length equal
	int i=0;
	if(a>b){

		for(i=b;i<a;i++){
			((obj+i) -> upper) = obj+(i+1);
			*((obj+i)->upper) = {0,0,NULL};
		}
	}
	else if(a<b){

		for(i=b;i>=a;i--){
			*(obj+i) = {0,0,NULL};
		}
	}
	
	//duplicate
	i=0;
	for(i=0;i<=a;i++){
		(obj+i)->numpart = (stuf+i)->numpart;
		(obj+i)->sign = (stuf+i)->sign;
	}
	obj->sign = stuf->sign;

	return 0;
}

int printdeg(deg *no){

	int i = lastlistno(no);

	if((no->sign) == -1){
		printf("-");}

	printf("%d ",(no+i)->numpart);
	i--;

	while(i>=0){
	printf("%04d ",(no+i)->numpart);
	i--;
	}
return 0;

}

//I used pow()
double deg_to_double(deg *no){
	int n = lastlistno(no);
	int i = 0;

	double output = 0;

	for(i=0;i<=n;i++){
		output += ((no+i)->numpart)*pow(onemax,i);
	}
	if((no->sign) == -1){
		output *= -1;
	}

	return output;
}


//on user input(not referencial), third integer is 0.
int deg_x_uint(deg *no, int mtplier, int updig){

int q,n,r;

n = (no->numpart)*mtplier + updig;
q = n/onemax;
r = n%onemax;

(no->numpart) = r;

if(((no->upper)==NULL) && (q==0)){
	return 1;
}
else if(((no->upper)==NULL) && (q>0)){
	(no->upper) = no+1;
	*(no->upper) = {0,0,NULL};
}

deg_x_uint((no->upper),mtplier,q);

return 0;

}

int deg_p_deg(deg *no1, deg *no2){
	//make same number of digit
	int a = lastlistno(no1);
	int b = lastlistno(no2);
	
	int sign1 = no1->sign;
	int sign2 = no2->sign;

	int c = 0;
	if (a>b){
		c = a;}
	else{
		c = b;}
	int i = 0;
	if(a>b){

		for(i=b;i<a;i++){
			((no2+i) -> upper) = no2+(i+1);
			*((no2+i)->upper) = {0,0,NULL};
		}
	}
	else if(a<b){

		for(i=a;i<b;i++){
			((no1+i) -> upper) = no1+(i+1);
			*((no1+i)->upper) = {0,0,NULL};
		}
	}
	no1->sign = sign1;
	no2->sign = sign2;

	// we consider case 1 : same sign / case 2: different sign
	//case 1
	if((no1->sign)==(no2->sign)){
		// sign will be conserved

		// add compotentwise
		for(i=0;i<=c;i++){
			(no1+i)->numpart = ((no1+i)->numpart)+ ((no2+i)->numpart);
		}
		// deg-ization
		for(i=0;i<c;i++){
			if(((no1+i)->numpart) > onemax){
				(no1+(i+1))->numpart += 1;
				(no1+i)->numpart -= onemax;
			}
		}
		if(((no1+c)->numpart) > onemax){
			((no1+c) -> upper) = no1+(c+1);
			*((no1+c) -> upper) = {0,1,NULL};

			(no1+c)->numpart -= onemax;
		}
	}
	//case 2(with declare new sign, define new leftest)
	else{
		//signing;
		for(i=0;i<=c;i++){
			(no1+i)->numpart = ((no1+i)->numpart) * (no1->sign);
			(no2+i)->numpart = ((no2+i)->numpart) * (no2->sign);
		}
		// add compotentwise. each component has -onemax+1 ~ onemax-1
		for(i=0;i<=c;i++){
			(no1+i)->numpart = ((no1+i)->numpart)+ ((no2+i)->numpart);
		}
		// deg-ization
		for(i=0;i<c;i++){
			if(((no1+i)->numpart) < 0){
				(no1+(i+1))->numpart -= 1;
				(no1+i)->numpart += onemax;
			}
		}
		if(((no1+c)->numpart) < 0){
			no1->sign = -1;
			(no1+c)->numpart += onemax;
			
			int m;
			for(i=0;i<=c;i++){

			m = (no1+i)->numpart; 

			if(m != 0){
				(no1+i)->numpart = onemax - ((no1+i)->numpart);
				if(i<c){
				(no1+(i+1))->numpart += 1;}
			}
			}
		}
		else{
			no1->sign = 1;
		}
		// define new 'leftest of no1'
		for(i=c;i>0;i--){
			if (((no1+i)->numpart) != 0){
				break;
			}
			else{
				// cut link 'no+i-1->upper = no+i'
				((no1+(i-1))->upper) = NULL;
			}
		}
		// bring numpart of no2 back to positive
		for(i=0;i<=b;i++){
			(no2+i)->numpart = ((no2+i)->numpart) * (no2->sign);		
		}
	}

	// bring 'leftest of no2' back
	if(a>b){
		for(i=b;i<=a;i++){
			((no2+i) -> upper) = NULL;
		}
	}

	return 0;

}

int deg_m_deg(deg *no1, deg *no2){

	no2->sign *= -1;

	deg_p_deg(no1, no2);

	no2->sign *= -1;

	return 0;
}

int deg_d_uint(deg *n, int m){
	
	if(m==1){
		return 1;
	}

	int r=0; //remainder
	
	deg *q;
	q = (deg*)malloc(100*sizeof(deg));
	*q = {(n->sign),1,NULL};

	///from 'int copy(.,.)'
	int a = lastlistno(n);
	int b = lastlistno(q); //( = 0 )
	// obj is initialized state. this stuf's listno is no less then obj's.

	// make length equal
	int i=0;
	if(a>b){

		for(i=b;i<a;i++){
			((q+i) -> upper) = q+(i+1);
			*((q+i) -> upper) = {0,0,NULL};
		}
	}
	///
	
	///dummy declaration
	deg *dum;
	dum = (deg*)malloc(100*sizeof(deg));
	*dum = {1,1,NULL};
	for(i=1;i<=a;i++){
		((dum+i) -> upper) = dum+(i+1);
		*((dum+i) -> upper) = {0,0,NULL};
	}

	//compotentwise deviding
	for(i=0;i<=a;i++){

		(q+i)->numpart = ((n+i)->numpart)/m;
		(dum+i)->numpart = ((n+i)->numpart)%m;
	}

	// making q
	int x=0;
	for(i=a;i>=1;i--){

		x=((dum+i)->numpart)*onemax + ((dum+(i-1))->numpart);
		
		(q+(i-1))->numpart += x/m;

		(dum+i)->numpart = 0;
		(dum+(i-1))->numpart = x%m;
	}

	r = (dum->numpart);

	// deg-ization
	for(i=0;i<a;i++){
		if(((q+i)->numpart) > onemax){
			(q+(i+1))->numpart += 1;
			(q+i)->numpart -= onemax;
		}
	}

	if((a!=0) && (((q+a)->numpart) == 0)){
		((q+(a-1))->upper) = NULL;
	}

	if(((q+a)->numpart) > onemax){
		((q+a) -> upper) = (q+(a+1));
		*((q+a) -> upper) = {0,1,NULL};

		(q+a)->numpart -= onemax;
	}
	

	copy(n,q);

	free(q);
	free(dum);
	return r;
}
