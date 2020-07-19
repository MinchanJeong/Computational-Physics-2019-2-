#ifndef NEWTONRAPHSON_H
#define NEWTONRAPHSON_H

#include "function.h"
#include "linalgb.h"

extern bool debug;
extern int setdebug(bool);

extern double NR_1d(func_1d, double guess,double h, double epsilon);

extern matrix* NR_multi(matrix* rtn,func_multi*, matrix v0, matrix y, double *filter, double *h, double *w, double err);

struct root_list_struct{
	matrix *rootptr;
	struct root_list_struct *next;
};
typedef root_list_struct root_list;
//last 
extern root_list *init_root_list();
extern bool is_root_overlapped(matrix root, root_list *list, double *w ,double err);

#endif
