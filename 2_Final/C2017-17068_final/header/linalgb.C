#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "common.h"
#include "linalgb.h"

#define epsln 1e-9
#define errmessage false

int matrix_print(matrix M,int dgt,int dcml,bool inrow){

	int m = M.cdim;
	int n = M.rdim;
	
	if(inrow == true){
		m = m+n; n = m - n; m = m - n;
	}

	int i,j;
	for(i=0;i<m;i++){
		
		printf("| ");
		for(j=0;j<n;j++){
			printf("%0*.*e ",dgt,dcml,M.mat[i+m*j]);
		}
		printf("|");
		if(i!=m-1){printf("\n");}
	}

	return 0;

}

void write_mat_allocdlist(matrix* M){
	
	write_allocdlist(M->mat);
	write_allocdlist(M);

}

int matrix_copy(matrix *obj, matrix origin){
	
	int m = origin.cdim;
	int n = origin.rdim;
	if((obj->cdim != m || obj->rdim != n) && errmessage){
		fprintf(stderr,"\nmatrix copy size error!!\n");
		return 1;
	}
	
	int i;
	for(i=0;i<m*n;i++){
			obj->mat[i] = origin.mat[i];
	}

	return 0;
}

matrix *new_0matrix(int m,int n){
	matrix *M = (matrix*)malloc(sizeof(matrix));
	M->cdim = m;
	M->rdim = n;
	M->mat = (double*)malloc(m *n * sizeof(double));
	int i;
	for (i=0;i<m*n;i++){
		M->mat[i]=0.0;
	}
	return M;
}

int freematrix(matrix* M){
	free(M->mat);
	free(M);
	return 0;
}

matrix *mat_x_mat(matrix* rtn, matrix mat_A, matrix mat_B){
	
	if((mat_A.rdim != mat_B.cdim) && errmessage){
		fprintf(stderr,"\nmatrix product size error!!\n");
		return NULL;
	}

	int m = mat_A.cdim;
	int n = mat_B.rdim;

	if(rtn != NULL){
	if((rtn->cdim != m || rtn->rdim != n) && errmessage){
		fprintf(stderr,"\nreturn matrix size error!!\n");
		return NULL;
	}}

	matrix *dummy = new_0matrix(m,n);

	int i,j,k;
	int l = mat_A.rdim;
	for(i=0;i<m;i++){
	for(j=0;j<n;j++){
		dummy->mat[i+m*j]=0;
		for(k=0;k<l;k++){
			(dummy->mat)[i+m*j] += mat_A.mat[i+m*k] * mat_B.mat[k+m*j];
		}
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

matrix *mat_p_mat(matrix* rtn, matrix mat_A, matrix mat_B){
	
	if(mat_A.rdim != mat_B.rdim || mat_A.cdim != mat_B.cdim){
		fprintf(stderr,"\nmatrix addition size error!!\n");
		return NULL;
	}

	int m = mat_A.cdim;
	int n = mat_A.rdim;

	if(rtn != NULL){
	if(rtn->cdim != m || rtn->rdim != n){
		fprintf(stderr,"\nreturn matrix size error!!\n");
		return NULL;
	}}

	matrix *dummy = new_0matrix(m,n);

	int i,j;
	for(i=0;i<m;i++){
	for(j=0;j<n;j++){
		dummy->mat[i+m*j]= mat_A.mat[i+m*j] + mat_B.mat[i+m*j];
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

matrix *const_x_mat(matrix* rtn, double c){
	
	int m = rtn->cdim;
	int n = rtn->rdim;

	int i,j;
	for(i=0;i<m;i++){
	for(j=0;j<n;j++){
		rtn->mat[i+m*j] *= c;
	}}
	
	return rtn;
}

//THIS FUNCTION ALLOCS 1 matrix* & 2 matrix.
matrix** matrix_to_LU(matrix M)//Dolittle Algorithm,used in det and inv
{
	if(M.rdim != M.cdim){
		fprintf(stderr,"\nUnsupported size for LU decomposition.\n");
		return NULL;
	}
	int N = M.rdim;


	matrix **LU = (matrix**)malloc(2*sizeof(matrix*));//LU[0] = L, LU[1] = U
	matrix *L = new_0matrix(N,N);
	matrix *U = new_0matrix(N,N);
	
	LU[0] = L, LU[1] = U;
	
	int i,j,k;
	double sum;
	for(i=0;i<N;i++){
		for(k=i;k<N;k++){
			sum = 0;
			for(j=0;j<i;j++){sum += ((L->mat[i+N*j]) * (U->mat[j+N*k])); }
			(U->mat)[i+N*k] = (M.mat)[i+N*k] - sum ; 
		}

		for(k=i;k<N;k++){
			if(i==k){L->mat[i+N*i] = 1;}
			else{
			sum = 0;
			for(j=0;j<i;j++){sum += (L->mat[k+N*j]) * (U->mat[j+N*i]); }
			
			if(!is_num_proper(U->mat[i+N*i])){
				
				if(errmessage){
				fprintf(stderr,"\nLU: MATRIX ELEMENT ERROR!\n");}
				U->cdim = -1;
				return LU;
			}
			else if(is_num_proper(U->mat[i+N*i])&&fabs(U->mat[i+N*i])<epsln){
				
				if(errmessage)
				{fprintf(stderr,"\nLU: MATRIX RANK WARNING!\n");}
				U->mat[i+N*i] = epsln;
			}
				
		
			L->mat[k+N*i] = (M.mat[k+N*i] - sum) / U->mat[i+N*i];
			}
		}
	}
	
	return LU;
}

// THIS FUNCTION ALLOCS 1 matrix "if return is not NULL"
matrix* matrix_U_inv(matrix* mat_U){
	
	int N = mat_U->cdim;
	
	double *U = mat_U->mat;
	int i;

	if(N == -1 && errmessage){
			fprintf(stderr,"\nUpper Triangular is not invertible!!!!\n");
	}

	matrix* Uinv = new_0matrix(N,N);
	for(i=0;i<N;i++){
		Uinv->mat[i+N*i]=1;
	}
	
	int j,k;
	for(i=N-1;i>=0;i--){
	for(j=i-1;j>=0;j--){
		for(k=0;k<N;k++){
			Uinv->mat[j+N*k] -= (U[j+N*i]/U[i+N*i])*Uinv->mat[i+N*k];
		}
	
	}}

	for(i=0;i<N;i++){
	for(k=0;k<N;k++){
		Uinv->mat[i+N*k] /= U[i+N*i];
	}}
	
	return Uinv;
}

// THIS FUNCTION ALLOCS 1 matrix "if return is not NULL"
matrix* matrix_L_inv(matrix* mat_L){
	
	int N = mat_L->cdim;
	
	double *L = mat_L->mat;
	
	int i;

	matrix* Linv = new_0matrix(N,N);
	for(i=0;i<N;i++){
		Linv->mat[i+N*i]=1;
	}
	
	int j,k;
	for(i=0;i<N;i++){
	for(j=i+1;j<N;j++){
		for(k=0;k<N;k++){
			Linv->mat[j+N*k] -= (L[j+N*i]/L[i+N*i])*Linv->mat[i+N*k];
		}
	
	}}
	
	return Linv;
}

matrix* matrix_inv(matrix* rtn,matrix M){
	
	if(M.rdim != M.cdim){
		fprintf(stderr,"\nINVERSE: MATRIX NOT SQUARE!!\n");
		return NULL;
	}

	matrix** LU = matrix_to_LU(M);
	if(LU[1]->cdim != -1){
	
		matrix* Linv = matrix_L_inv(LU[0]);
		matrix* Uinv = matrix_U_inv(LU[1]);
	
		int N = M.cdim;
		matrix* dummy = new_0matrix(N,N);
		mat_x_mat(dummy,*Uinv,*Linv);

		freematrix(Linv);
		freematrix(Uinv); 
		freematrix(LU[0]); 
		freematrix(LU[1]);
		free(LU);

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

	else{
		if(errmessage){
		fprintf(stderr,"\nINVERSE DOES NOT EXIST!!\n");}
		freematrix(LU[0]);
		freematrix(LU[1]);
		free(LU);
		return NULL;
		}
	}


double absmean(matrix M, double *w){
	
	double L2 = 0;

	int m = M.cdim;
	int n = M.rdim;

	int i;
	for(i=0;i<(m*n);i++){
		L2 += fabs(M.mat[i])*w[i];
	}
	L2 /= (m*n);	

	return L2;
}

bool is_num_proper(double num){

	bool flag = false;
	if(isinf(num)==0 && isnan(num)==0)
	{ flag = true;} 
	//else{printf("%f\n",num);}
	return flag;
}

bool is_mat_proper(matrix M){

	int m = M.cdim;
	int n = M.rdim;
	double *arr = M.mat;
	
	bool flag = true;

	for(int i=0;i<m*n;i++){	
		if(!is_num_proper(arr[i]))
		{//printf("\n%f",arr[i]);
		flag = false; break;}
	}
	return flag;
}
