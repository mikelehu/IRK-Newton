/*------------------------------------------------------------------------------*/
/*										*/
/*                                GaussCommon.h					*/
/*										*/
/* ----------------------------------------------------------------------------*/

#include <stdio.h> 
#include <stdlib.h> 
#include <string.h>
#include <math.h>
#include <def.h>
#include <sys/stat.h>
#include <GaussUserProblem.h> 
#include <GaussCoefficients.h>
#include <omp.h>
#include <time.h>
#include <sys/types.h>
#include <mathlink.h>
#include <cblas.h>
#include <lapacke.h>


val_type NormalizedDistance0
(int neq,int ns,toptions *options,val_type0 *z,val_type0 *zold
 );


int Yi_update0
(solution *u, val_type0 *z0,ode_sys *system,
 gauss_method *method,solver_stat *thestatptr
);

void RemoveDigitsFcn0
(ode_sys *system,gauss_method *gsmethod,val_type0 *z0, int m
);

void UpdateDMin0
(ode_sys *system,gauss_method *method,
 int *D0,bool *cont,val_type0 *DMin,val_type0 *Y, val_type0 *Yold
);


void Newton_it_Mix
(ode_sys *system, solution *u, val_type tn,val_type h, 
 toptions *options,gauss_method *method,solver_stat *thestatptr
);


void Newton_it_Mix2
(ode_sys *system, solution *u, val_type tn,val_type h, 
 toptions *options,gauss_method *method,solver_stat *thestatptr
);



void MyKroneckerProd0
(int alpha,lowfloat0 *A, lowfloat0 *B, 
 int beta, lowfloat0 *C,int i0,
 int n1, int m1, int n2, int m2
);


void MMfun0 
(int neq, gauss_method *method,
 lowfloat0 *Jac0,lowfloat0 *MM0
);


void MMfunOsoa0 
(int neq, gauss_method *method,val_type tn,
 ode_sys *system,solver_stat *thestatptr,lowfloat0 *MM
);

int Newton_Step0
(ode_sys *system, solution *u, val_type tn,val_type h, 
 gauss_method *method,solver_stat *thestatptr,
 lowfloat0 *MM, int *ipiv
);

