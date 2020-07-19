#ifndef LINALG_H
#define LINALG_H

struct matrix_struct{
	double *mat;// mat_ij = mat[i+j*cdim]
	int cdim; // size of column
	int rdim; // size of row
};
typedef matrix_struct matrix;

extern int matrix_print(matrix M,int dgt,int dcml,bool inrow);

extern int matrix_copy(matrix*, matrix);

extern void write_mat_allocdlist(matrix*);
extern matrix* new_0matrix(int cdim, int rdim);
extern int freematrix(matrix*);

extern matrix* mat_x_mat(matrix* rtn,matrix,matrix);
extern matrix* mat_p_mat(matrix* rtn, matrix, matrix);
extern matrix* const_x_mat(matrix* rtn, double); // call by reference only
//
extern matrix** matrix_to_LU(matrix); // return *p s.t p[0] = L p[1] = U;
extern matrix* matrix_L_inv(matrix* L);
extern matrix* matrix_U_inv(matrix* U);
extern matrix* matrix_inv(matrix*,matrix);
//extern double matrix_det(matrix);

extern double absmean(matrix,double*);
extern bool is_num_proper(double);
extern bool is_mat_proper(matrix);;

#endif
