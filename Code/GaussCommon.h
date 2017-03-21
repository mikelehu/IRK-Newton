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


void print_u
(int neq,val_type *u
);

void InitStat
(ode_sys *system,gauss_method *gsmethod, solver_stat *thestatptr
 );

val_type NormalizedDistance
(int neq,int ns,toptions *options,val_type *z,val_type *zold
 );

void RemoveDigitsFcn
(val_type *x,int m
);

int Li_init
(solution *u, val_type *z,ode_sys *system,
 gauss_method *method,solver_stat *thestatptr, toptions *options
);

int Yi_update
(solution *u, val_type *z,ode_sys *system,
 gauss_method *method,solver_stat *thestatptr
);

int Yi_update_Classic
(solution *u, val_type *z,ode_sys *system,
 gauss_method *method,solver_stat *thestatptr
);

int statlinit
(ode_sys *system,gauss_method *method,
 solver_stat *thestatptr
);


void StopCriterion
(ode_sys *system,gauss_method *method,
 int *D0,bool *cont,val_type *DMin,val_type *L, val_type *Lold
);

void StopCriterionFloat         
(ode_sys *system,gauss_method *method,
 int *D0,bool *cont,val_type *DMin,val_type *L, val_type *Lold
);


void NSS_Step
(ode_sys *system, solution *u, val_type tn,val_type h, 
 toptions *options,gauss_method *method,solver_stat *thestatptr
);


void NSS_Step_plus		
(ode_sys *system, solution *u, val_type tn,val_type h, 
 toptions *options,gauss_method *method,solver_stat *thestatptr
);


void NSS_MIX_Step
(ode_sys *system, solution *u, val_type tn,val_type h, 
  toptions *options,gauss_method *method,solver_stat *thestatptr
);

void MatAdd
(val_type alpha,val_type *A, int lda, 
 val_type beta,val_type *B, int ldb, 
 val_type gamma, val_type *C, int ldc,
 int m, int n
);

void MatAddID
(val_type aa,val_type bb,val_type *B, int ldb, 
 val_type cc, val_type *C, int ldc, int n
);


void MMNewtonSS
(int neq, gauss_method *method, val_type tn, val_type h, 
 ode_sys *system, val_type *u, val_type *MM
);


int NSS_Solve
(ode_sys *system, solution *u, val_type tn,val_type h, 
 gauss_method *method,solver_stat *thestatptr,
 val_type *MM, int *ipiv
);


int NSS_Solve_plus 		
(ode_sys *system, solution *u, val_type tn,val_type h, 
 gauss_method *method,solver_stat *thestatptr,
 val_type *IDM, int *ipiv
);


void Compute_MM 
(ode_sys *system,val_type h,val_type *MM, 
 gauss_method *method,solver_stat *thestatptr
);


void Compute_R            
(ode_sys *system, val_type h,val_type *g, val_type *RR, 
 gauss_method *method,solver_stat *thestatptr
);


void Compute_Z 	  	  
(ode_sys *system, val_type h, val_type *IDM,val_type *RR, 
 val_type *ZZ,gauss_method *method,solver_stat *thestatptr,
 int *ipiv
);


void Compute_W1         
(ode_sys *system, val_type h, val_type *RR,val_type *ZZ,
 val_type *W1, gauss_method *method,solver_stat *thestatpt
);


void Compute_W2       
(ode_sys *system, val_type h, val_type *W1,val_type *ZZ,
 val_type *r, val_type *W2, gauss_method *method,
 solver_stat *thestatptr
);


void Compute_DL		
(ode_sys *system, val_type *W1,val_type *W2,
 val_type *DL, gauss_method *method,
 solver_stat *thestatptr
);

void Compute_GG        
(ode_sys *system, val_type h, 
 val_type *g, val_type *DL, val_type *GG,
 gauss_method *method, solver_stat *thestatptr
);


void TheOutput
(ode_sys *system,gauss_method *method,val_type t,solution *u,
 solver_stat *thestatptr, parameters *params,toptions *options,FILE *loga
);


void CompensatedSummation 
(gauss_method *gsmethod,
 solution *u,
 ode_sys *system, toptions *options,
 solver_stat *thestatptr
);


void RKG 
(gauss_method *gsmethod,
 solution *u,
 ode_sys *system, toptions *options,
 void RKG_Step (), solver_stat *thestatptr
);


void select_gauss
(gauss_method *gsmethodptr, 
 solution *uptr,ode_sys *systemptr,
 toptions *optionsptr, solver_stat *thestatptr
);

void select_odefun
(int codfun, ode_sys *system
);



