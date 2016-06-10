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
#include <GaussCommon0.h>


void print_u
(int neq,val_type *u
);

void InitStat
(ode_sys *system,gauss_method *gsmethod, solver_stat *thestatptr
 );

val_type NormalizedDistance
(int neq,int ns,toptions *options,val_type *z,val_type *zold
 );

int Li_init
(solution *u, val_type *z,ode_sys *system,
 gauss_method *method,solver_stat *thestatptr, toptions *options
);

int Yi_update
(solution *u, val_type *z,ode_sys *system,
 gauss_method *method,solver_stat *thestatptr
);

int statlinit
(ode_sys *system,gauss_method *method,
 solver_stat *thestatptr
);

void RemoveDigitsFcn
(ode_sys *system,gauss_method *gsmethod,val_type *z, int m
);

void UpdateDMin
(ode_sys *system,gauss_method *method,
 int *D0,bool *cont,val_type *DMin,val_type *Y, val_type *Yold
);


void Newton_it
(ode_sys *system, solution *u, val_type tn,val_type h, 
 toptions *options,gauss_method *method,solver_stat *thestatptr
);


void MyKroneckerProd
(int alpha,lowfloat *A, lowfloat *B, 
 int beta, lowfloat *C,int i0,
 int n1, int m1, int n2, int m2
);


void MMfun 
(int neq, gauss_method *method,
 lowfloat *Jac,lowfloat *MM
);

void MMfunOsoa 
(int neq, gauss_method *method,val_type tn,
 ode_sys *system,solver_stat *thestatptr,lowfloat *MM
);


int Newton_Step
(ode_sys *system, solution *u, val_type tn,val_type h, 
 gauss_method *method,solver_stat *thestatptr,
 lowfloat *MM, int *ipiv
);

void TheOutput
(ode_sys *system,gauss_method *method,val_type t,solution *u,
 solver_stat *thestatptr, parameters *params,toptions *options,FILE *loga
);

void TheOutput2
(ode_sys *system,gauss_method *method,val_type t,solution *u,solution *u2,
 solver_stat *thestatptr,parameters *params,toptions *options,FILE *loga
);


void RKG 
(gauss_method *gsmethod,
 solution *u,
 ode_sys *system, toptions *options,
 void RKG_Step (), solver_stat *thestatptr
);


void RKG2 
(gauss_method *gsmethod, gauss_method *gsmethod2,
 solution *u, solution *u2,
 ode_sys *system,toptions *options,toptions *options2,
 void RKG_Step (), solver_stat *thestatptr, solver_stat *thestatptr2
);


void select_gauss
(gauss_method *gsmethodptr, gauss_method *gsmethod2ptr,
 solution *uptr,solution *u2ptr,ode_sys *systemptr,
 toptions *optionsptr,toptions *options2ptr,
 solver_stat *thestatptr,solver_stat *thestat2ptr
);




