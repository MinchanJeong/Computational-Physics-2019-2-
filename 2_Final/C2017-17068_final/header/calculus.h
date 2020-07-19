#ifndef QUAD_H
#define QUAD_H

#include "function.h"

/* Integration */
extern numrcoefs qdr_str_to_coefs(const char* option);//
extern double qdr_calc_base(double* f_arr, int N, double h, numrcoefs co);

extern double qdr_calc_1d(func_1d,numrcoefs,double a,double b,double h,int N);

extern double quadr(func_1d, double a,double b, double h, const char* option);
extern double quadr_N(func_1d, double a, double b, int N, const char* option);

extern double quadr_2d(func_multi,double *dom, double *h, const char* option);
#endif


#ifndef DIFFERENTIAL_H
#define DIFFERENTIAL_H

#include "function.h"
#include "linalgb.h"

/* Differentiation */
extern numrcoefs diff_str_to_coefs(const char* option);
extern double diff_calc_base(double* f_arr, double h, numrcoefs co);

extern double diff_calc_1d(func_1d, double x0, double h, numrcoefs co);
extern double diff_1d(func_1d, double x0, double h,  const char* option);

extern double diff_calc_multi(func_multi, matrix v0, int k, double h, numrcoefs co);
extern double diff_multi(func_multi, matrix v0, int k, double h, const char* option);

#endif


#ifndef MULTIDIFF_H
#define MULTIDIFF_H

#include "function.h"
#include "linalgb.h"

extern matrix* evaluate(func_multi* vectorftn, matrix* rtn); // call for ref 
extern matrix* Jacobian(matrix* rtn, func_multi*, int dim, matrix v0, double *h, const char* option);

#endif

