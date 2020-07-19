#include <stdio.h>
#include <math.h>

#include "common.h"
#include "function.h"
#include "calculus.h"
#include "rootfind.h"
#include "linalgb.h"

#define stableNo 4
#define breakNo 50

bool debug = false;

int setdebug(bool flag){

	debug = flag;
	return 0;
}

double NR_1d(func_1d ftn, double x0, double h, double e){

	double x = x0;
	double delta,Df;
	double (*f)(double,...) = ftn.fpart;
	
	bool sign = false;
	int count = 0;
	int loopNo = 0;
	do{
		if(loopNo>breakNo){sign = true; break;}
		loopNo ++;
		Df = diff_1d(ftn,x,h,"h^6");
		delta = -1 * f(x)/Df;
			
		if(!is_num_proper(Df) || !is_num_proper(delta)){sign = true; break;}

		x += delta;

		if(fabs(delta)<e){count++;}

	}while(count<stableNo);
	
	if(sign == false)
	{return x;}
	else
	{return 0.0/0.0;}
}

matrix* NR_multi(matrix* rtn,func_multi* vectorftn, matrix v0, matrix y, double *filter_arr, double *h, double *w, double e){
	
	int N = v0.cdim;

	matrix *root = new_0matrix(N,1);

	matrix *dummy_root = new_0matrix(N,1);
	matrix *dummy_jac = new_0matrix(N,N);
	matrix *filter = new_0matrix(N,N);
	for(int i=0;i<N;i++){filter->mat[i+N*i] = filter_arr[i];}

	matrix_copy(root,v0);
	matrix_copy(dummy_root,v0);

	int sign = false; // true: jacobian or result is improper
	matrix *ptr;
	char ch;

	int count = 0;
	int loopNo = 0;
	do{
		loopNo++;
		if(loopNo>breakNo){sign = true; break;}

		Jacobian(dummy_jac,vectorftn,N,*root,h,"h^4");
		
		if(debug==true){
		getchar(); matrix_print(*dummy_jac,7,3,false);printf("\n");}

		ptr = matrix_inv(dummy_jac,*dummy_jac);
		if(ptr == NULL){
			sign = true;
			break;
		}

		matrix_copy(dummy_root,*root);
		evaluate(vectorftn, dummy_root);
		const_x_mat(dummy_root,-1);
		mat_p_mat(dummy_root,*dummy_root,y);
		
		mat_x_mat(dummy_root,*dummy_jac,*dummy_root);
		mat_x_mat(dummy_root,*filter,*dummy_root);
		
		if(is_mat_proper(*dummy_root)&&(absmean(*dummy_root,w) < e)){
			count++;
		}

		mat_p_mat(root,*root,*dummy_root);
		
		if(debug == true){
			getchar(); 
			matrix_print(*root,7,3,true); printf("\n\n");
		}
		if(is_mat_proper(*root)==false){
			sign = true; break;
		}
	}while(count<stableNo);
	freematrix(dummy_root);
	freematrix(dummy_jac);
	freematrix(filter);

	if(sign == true){freematrix(root); return NULL;}
	
	if(rtn == NULL){
		write_mat_allocdlist(root);
		return root;
	}
	else{
		matrix_copy(rtn,*root);
		freematrix(root);
		return rtn;
	}
}

root_list *init_root_list(){
	
	root_list *list = (root_list*)malloc(sizeof(root_list));
	list->rootptr = NULL;
	list->next = NULL;

	write_allocdlist(list);
	
	return list;
}

bool is_root_overlapped(matrix root, root_list *list, double *w, double err){
	
	root_list* ptr = list;
	
	matrix *_root = new_0matrix(root.cdim,1);
	matrix_copy(_root,root);

	double dum[root.cdim] = {0,};
	matrix dummy = {dum,root.cdim,1};

	const_x_mat(_root,-1.0);

	if(ptr->rootptr != NULL){
	int count=0;
	do{	
		if(count != 0){ptr = ptr->next;}
		count ++;
		//printf("%d\n",count);
		if(absmean(*(mat_p_mat(&dummy,*_root,*(ptr->rootptr))),w)<2*err){
			freematrix(_root);
			return true;
		}

	}while(ptr->next != NULL);

	ptr->next = (root_list*)malloc(sizeof(root_list));
	write_allocdlist(ptr->next);
	ptr = ptr->next;
	}

	const_x_mat(_root,-1.0);
	ptr->rootptr = _root;
	write_mat_allocdlist(_root);
	ptr->next = NULL;
	return false;
}
