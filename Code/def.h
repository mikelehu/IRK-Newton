/*------------------------------------------------------------------------------*/
/*										*/
/*                                def.h    					*/
/*										*/
/* -----------------------------------------------------------------------------*/
#include "prec.h"
#include <stdbool.h>

/* -----------------------------------------------------------------------------*/
/*										*/
/*	Parameters								*/
/*										*/
/* -----------------------------------------------------------------------------*/

#ifndef DEFH_
#define DEFH_

#define PI atan(1.) * 4.   
#define INF 99999.
#define MAXNEQ	200				// Maximum neq.
#define MAXPARAM 30				// Number of maximum parameters.
#define MAXIT 50				// Maximum number of fixed point iterations 50
#define MAXKSW 10                               // Max number of steps to change second 
                                                // integration initialization mode.
#define SUCCESS 0
#define FAIL -1
#define STRMAX 40				// Filename maximum string.
#define RTOL pow(10.,-12);			// pow(2.,-40);
#define ATOL pow(10.,-12);
#define ESTERRTHRS pow(10.,16)			// EstimatedErrorThreshold. pow(10.,12)
#define IOUT
#define PARALLEL				// Active opem-mpi parallel execution     

/* -----------------------------------------------------------------------------*/
/*										*/
/*	PREC: FLOAT,DOUBLEPRECISION,QUADRUPLEPRECISION				*/
/*										*/
/* -----------------------------------------------------------------------------*/

//#define PREC 2        to be specified in the Makefile file ( -D PREC=value)
//1=DOUBLEPRECISION    				// Specify by: gcc -D DOUBLE ...
//2=QUADRUPLEPRECISION			  
//3=FLOAT

#if PREC ==2 
#define QUADRUPLEPRECISION 
typedef __float128 val_type;
typedef double lowfloat;
#define GETRF LAPACKE_dgetrf
#define GETRS LAPACKE_dgetrs
#define GEMM cblas_dgemm
typedef double val_type0;
typedef double lowfloat0;			//float or double
#define GETRF0 LAPACKE_dgetrf			//LAPACKE_sgetrf,LAPACKE_dgetrf
#define GETRS0 LAPACKE_dgetrs			//LAPACKE_sgetrs,LAPACKE_dgetrs
#define GEMM0 cblas_dgemm			//cblas_sgemm,cblas_dgemm

#elif PREC ==3
#define FLOAT
typedef float val_type;
typedef float lowfloat;
typedef float val_type0;
typedef float lowfloat0;
#define GETRF LAPACKE_sgetrf
#define GETRS LAPACKE_sgetrs
#define GEMM cblas_sgemm
#define RTOL pow(2.,-16);
#define ATOL pow(2.,-16);
#define GETRF0 LAPACKE_sgetrf
#define GETRS0 LAPACKE_sgetrs
#define GEMM0 cblas_sgemm

#else
#define DOUBLEPRECISION
typedef double val_type;
typedef double lowfloat;			 //float or double
typedef float val_type0;
typedef float lowfloat0;
#define GETRF LAPACKE_dgetrf			//LAPACKE_sgetrf,LAPACKE_dgetrf
#define GETRS LAPACKE_dgetrs			//LAPACKE_sgetrs,LAPACKE_dgetrs
#define GEMM cblas_dgemm			//cblas_sgemm,cblas_dgemm
#define GETRF0 LAPACKE_sgetrf
#define GETRS0 LAPACKE_sgetrs
#define GEMM0 cblas_sgemm
#endif


#ifdef QUADRUPLEPRECISION

/* val_type functions */

#define POW(x, y)     powq(x, y)
#define SQRT(x)       sqrtq(x)

#define EXP(x)        expq(x)
#define LOG(x)        logq(x)

#define SIN(x)        sinq(x)
#define COS(x)        cosq(x)
#define TAN(x)        tanq(x)

#define FABS(x)	      fabsq(x)
#define FMAX(x,y)     fmaxq(x,y)

/* val_type0 functions */

#define POW0(x, y)     pow(x, y)
#define SQRT0(x)       sqrt(x)

#define EXP0(x)        exp(x)
#define LOG0(x)        log(x)

#define SIN0(x)        sin(x)
#define COS0(x)        cos(x)
#define TAN0(x)        tan(x)

#define FABS0(x)       fabs(x)
#define FMAX0(x,y)     fmax(x,y)


#else //DOUBLEPRECISION or FLOAT

/* val_type functions */

#define POW(x, y)      pow(x, y)
#define SQRT(x)        sqrt(x)

#define EXP(x)        exp(x)
#define LOG(x)        log(x)

#define SIN(x)        sin(x)
#define COS(x)        cos(x)
#define TAN(x)        tan(x)

#define FABS(x)	      fabs(x)
#define FMAX(x,y)     fmax(x,y)

/* val_type0 functions */

#define POW0(x, y)     pow(x, y)
#define SQRT0(x)       sqrt(x)

#define EXP0(x)        exp(x)
#define LOG0(x)        log(x)

#define SIN0(x)        sin(x)
#define COS0(x)        cos(x)
#define TAN0(x)        tan(x)

#define FABS0(x)       fabs(x)
#define FMAX0(x,y)     fmax(x,y)

#endif


/* -----------------------------------------------------------------------------*/
/*										*/
/*	General definitions							*/
/*										*/
/* -----------------------------------------------------------------------------*/

typedef struct gauss_method
   {
     int ns;
     val_type *m;	 			 // mij=aij/bj and mij+mji-1=0.
     val_type *c,*b;	 			 // c,b coefficients.
     val_type *hc;       			 // hc=h*c.
     val_type *hb;       			 // hb=h*b.
     val_type *nu;                               // interpolate coeficcients (a*/bj).
     val_type *munu;				 // interpolate coeficcients (Li formulation).
     int *orderedindices;			 // ascending order indices for bi coeficients

     val_type0 *m0;	 			 // mij=aij/bj and mij+mji-1=0.
     val_type0 *c0,*b0;	 			 // c,b coefficients.
     val_type0 *hc0;       			 // hc=h*c.
     val_type0 *hb0;       			 // hb=h*b.
     val_type0 *nu0;                             // interpolate coeficcients (a*/bj).
     val_type0 *munu0;				 // interpolate coeficcients (Li formulation).

   } gauss_method;


typedef struct solution
  { 
     val_type *uu,*ee;

  } solution;


typedef struct toptions
  {
     val_type h; 				  // stepsize
     val_type t0;
     val_type t1;
     val_type *rtol,*atol;
     int algorithm;			  
     int sampling;
     int approximation;
     int rdigits,mrdigits;
     int (*iteration[2])();			 // Iteration : Jacobi1, Jacobi2.
   } toptions; 


typedef struct parameters
   { 
     int numrpar;				 // Number of real parameters.
     val_type *rpar;				 // Variables for specifying odefun real parameters.
     int numipar;				 // Number of int parametes. 
     int *ipar;     			         // Variables for specifying odefun integer parameters.
					         //  ipar[0]: which part of differential equation must be evaluated
     val_type0 *rpar0;
   } parameters;


typedef struct ode_sys       
   {
     int problem;
     int neq;					// number of equations.
     int n;					// n-body.
     void (*f)();				// odefun.
     void (*jac)();				// jacobian.
     val_type (*ham)();				// hamiltonian
     int cod[2];				// tO specify which part of ODE must be evaluated.
     parameters params;

     void (*f0)();				// odefun.
     void (*jac0)();				// jacobian.
     val_type0 (*ham0)();			// hamiltonian


    } ode_sys;


typedef struct solver_stat
    {

    int convergence;				  // SUCCESS or FAIL. 			        
    bool laststep;                                

    /* auxiliar*/
    val_type *z,*li,*liold,*lit0;
    val_type0 *z0,*li0;

    val_type E0;     		      		  // Initial Energy
    val_type MaxDE;		       		  // MaxDE=Abs(Ei-E0/E0)    
		 

    /* stadistics*/
    int stepcount;
    int itcount;
    int totitcount;
    int maxitcount;			
    int itzero;				
    int fcn;
    int *initqlty;				// ns*neq matrix. Quality of initialization of Li stages    

    /* output filename */
    char filename[STRMAX];			// Integration filename.
    int nout;                                   // number of output values.

    } solver_stat; 



	
#endif /*DEFH_*/
