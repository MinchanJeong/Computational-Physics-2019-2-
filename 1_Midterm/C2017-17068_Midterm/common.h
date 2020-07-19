//#1
#ifndef COMMON_1
#define COMMON_1

/*function in #1 (func1.C) */

extern int indexmax;

extern float avg (int*,float*);
extern float err (int,int*,float*);

extern float mean(int,float*);

#endif

//#2,#3
#ifndef COMMON_2to3
#define COMMON_2to3

struct degree {
	int sign;
	// only -1 or 1 on first list(digit). otherwise '*.sign = 0'
	// if zero, sign is 1.
	int numpart; //numpart has value of 0~onemax-1(defined at func2.C)
	struct degree *upper;
};
	// (leftest list) -> upper = NULL : "must" be satisfied.

typedef struct degree deg;

extern int onemax;//'maximal # + 1' of number in one list.

/*Basic function with deg (degcore.C)*/

extern int lastlistno(deg*); //size of deg
extern int copy(deg*,deg*); //copy second to first
extern int printdeg(deg*); //print deg type exact number
extern double deg_to_double(deg*); // return dobule type representation of deg by using pow()

/*Operation with deg (degcore.C)*/ 
	//operation with deg act like 'first <- <first,second>'
	//i.e operations change 'value at given pointer'.

extern int deg_x_uint(deg*, int, int); //deg = deg*(small poitive integer)
extern int deg_p_deg(deg*, deg*); // first = first + second
extern int deg_m_deg(deg*, deg*); // first = first - second
extern int deg_d_uint(deg*,int); 
	// first = second/int, return is remainder. 
	// devide with small integer ( 1 ~ 100)
	// first deg is initialized, i.e 0 or 1.

/*Related to Factorial and Bernoulli no (factorial.C)*/

struct quotient {

	deg *n; // numerator
	deg *d; // denominator
};

typedef struct quotient quot;

extern int exfac(deg*, int); // deg -> deg*n!
extern int exfacfac(deg*, int); // n!! = n(n-2)(n-4)....
extern int bi(deg*, int, int);//binomial coeff by exfac()
extern int coef(deg*, int, int);//coefficient in equality deg operation

extern int isdivd(deg*,deg*,int); 
	// if two deg types are both devided by the int, return 0 and two deg type will be devided by the int. Otherwise, return 1 and there is no change for two deg types.

extern double fac(int); //factorial in double using exfac(), deg_to_double()
extern double dbbi(int , int );// binomial coeff by double using fac()

#endif

