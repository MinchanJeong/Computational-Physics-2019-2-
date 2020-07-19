#include <stdio.h>
#include <math.h>

#include "common.h"

int main(int argc, char** argv){

	/// # 1-(a)
	printf("\n# 1 - (a)\n");
	FILE *fp;
	fp = fopen("data.20191105","r");

	int length=0;
	
	// Convention: index starting with 0.
	
	int x_[7][indexmax]={0,};
	// array initialized with 0.
	// x_[x0][i] = 1 if and onlyif x_i=x0
	// x_[x0][i] = -1 if i = length of file
	float f[indexmax]={0,}; //f[i]: value of ftn at index_i

	float u[indexmax]={0,}; //x_i for calculating regression b
		
	int n=0;
	float g=0;
	int i=0;
	while(feof(fp)==0){
		if(fscanf(fp,"%d	%f\n",&n,&g)==1){}

		x_[n][i]=1;

		f[i]=g;
		u[i]=n;

		i++;
	}
	printf("\nDataset Determined..");

	length = i;
	printf("\nLength of data is %d.",length);

	for(i=0;i<6;i++){
		x_[i][length]=-1;
	}//marking

	int state = fclose(fp);
	if(state!=0){
		printf("closing error!");
	}
	printf("\n\n---------------\n");

	/// # 1-(b) : functions are defined at func1.C. I defined new array x for calculation.
	int x[6][length]={0,}; // x[x0] = { .. index having x0 ..,-1(end),0,0...}

	int xsize[6]={0,}; // how many 0,1,2,3,4,5 in data?

	int j=0;
	int k=0;

	for(i=0;i<6;i++){
		k=0;
		for(j=0;j<=length;j++){
			if(x_[i][j] == 1){
				x[i][k]=j;
				k++;
			}
		}
		x[i][k]=-1;
		xsize[i]=k;
	}

	/// # 1-(c)
	printf("\n# 1 - (c)\n");
	printf("\n<average and stderr of f[x]>\n\n");
	for(i=0;i<6;i++){
	printf("x = %d : avg(%d) = %e  err(%d) = %e\n",i,i,avg(x[i],f),i,err(xsize[i],x[i],f));
	}
	printf("\n---------------\n");

	/// # 1-(d): I will use linear regression method
	
	float a,b;
	float dum=0;

	float xmean = mean(length,u);
	float fmean = mean(length,f);

	//b?
	b=0;
	for(i=0;i<length;i++){
		b+=((u[i] - xmean)*(f[i] - fmean));
	}
	for(i=0;i<length;i++){
		dum+=((u[i] - xmean)*(u[i] - xmean));
	}
	b = b / dum ;
	
	//a?
	a = fmean - b *xmean;

	//print
	printf("\n# 1 - (d)\n");
	printf("\n<Interpretation>\n");
	printf("\nLet us our functional is f(x) = ax + b.\n\nThen ");
	printf("a= %f, b= %f.\n",a,b);
	printf("\n---------------\n");

	///1-(e). ref:  stats.stackexchange.com/question/138938/

	float omega[7] ={0,}; // w_i = sigma_i ^ (-.5)
	for(i=0;i<6;i++){
	omega[i] = pow(err(xsize[i],x[i],f),-0.5);
	}

	float sqerr_b=0; // sqaured error of b

	dum = 0;
	float dum2=0;
	float dum3=0;

	for(i=0;i<length;i++){
		dum += u[i]*u[i]*omega[(int)u[i]];
	}
	for(i=0;i<length;i++){
		dum2 += u[i]*omega[(int)u[i]];
	}
	for(i=0;i<length;i++){
		dum3 += omega[(int)u[i]];
	}

	dum = dum - dum2*dum2/dum3;
	dum = 1/dum;

	dum2=0;
	for(i=0;i<length;i++){
		dum2 += omega[(int)u[i]]*pow((f[i]-a-b*u[i]),2);
	}
	dum = (dum * dum2 / (length -2));

	sqerr_b = dum;

	printf("\n# 1 - (e)\n");
	printf("\n<statistical error>\n");
	printf("\nsquared error of b = %f.\n",sqerr_b);
	printf("\n---------------\n");
	return 0;

}
