#ifndef FUNC_H

#define FUNC_H

struct func_1d_struct{
//  const char* name;
	double (*fpart)(double,...);
};
typedef func_1d_struct func_1d;

struct numrcoefs_struct{
	int* coefs;
	int divider; // divide coeff with this value
	int size;
};
typedef numrcoefs_struct numrcoefs;

#include "linalgb.h"

struct func_multi_struct{
	double (*fpart)(matrix,...);
	int dim;
};
typedef func_multi_struct func_multi;

#endif
