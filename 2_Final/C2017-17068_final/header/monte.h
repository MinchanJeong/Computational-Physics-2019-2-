#ifndef MONTE_H
#define MONTE_H
#include "random.h"
#include "linalgb.h"

extern int monte_debug(int);
extern int init_monte_random();
extern int set_rand_vector(matrix* M, double *domain);

extern double montecarlo(func_multi,double *domain,int sampleno);
extern double metropolis(func_multi,func_multi weight,double *domain,int sampleno);

#endif


