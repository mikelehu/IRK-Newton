/*----------------------------------------------------------------------------*/
/*									      */
/*                                def.h    				      */
/*									      */
/* ---------------------------------------------------------------------------*/

#include <stdbool.h>

/* ---------------------------------------------------------------------------*/
/*									      */
/*	Parameters							      */
/*									      */
/* ---------------------------------------------------------------------------*/

#ifndef DEFH_
#define DEFH_

#define INF 99999.
#define MAXIT 1000		   // Maximum number of fixed point iterations 

#define SUCCESS 0
#define FAIL -1
#define STRMAX 256	           // Filename maximum string.
#define RTOL pow(10.,-12);	   // pow(2.,-40);
#define ATOL pow(10.,-12);
#define IOUT
//#define PARALLEL 		   // Active opem-mpi parallel execution     
#define RDIGITS
//#define TESTJAC                  // Test Jacobian is correct.

#define DIR_TERM "../CoefficientsData/"    // Path Coefficients (terminal)
#define DIR_MATH "../../CoefficientsData/" // Path Coefficients (mathematica)

#define DOUBLEPRECISION
typedef double val_type;


#define GETRF LAPACKE_dgetrf			
#define GETRS LAPACKE_dgetrs			
#define GBTRF LAPACKE_dgbtrf
#define GBTRS LAPACKE_dgbtrs
#define GETRI LAPACKE_dgetri
#define GEMM  cblas_dgemm			
#define GEMV cblas_dgemv			

//DOUBLEPRECISION or FLOAT

#define POW(x, y)     pow(x, y)
#define SQRT(x)       sqrt(x)

#define EXP(x)        exp(x)
#define LOG(x)        log(x)

#define SIN(x)        sin(x)
#define COS(x)        cos(x)
#define TAN(x)        tan(x)

#define FABS(x)	      fabs(x)
#define FMAX(x,y)     fmax(x,y)


/* ---------------------------------------------------------------------------*/
/*									      */
/*	General definitions						      */
/*									      */
/* ---------------------------------------------------------------------------*/

typedef struct gauss_method
   {
     int ns;
     int mm, smm;                       // 01-02-2017
     val_type *m;	 		// mij=aij/bj and mij+mji-1=0.
     val_type *c,*b,*a;	 		// c,b,a coefficients.
     val_type *hc;       		// hc=h*c.
     val_type *hb;       		// hb=h*b.
     val_type *nu;                      // interpolate coeficcients (a*/bj).
     val_type *munu;		        // interpolate coeficcients (Li formulation).
  
     val_type *Q1,*Q2;                  // 01-02-2017
     val_type *sigma2;                  // 01-02-2017 sigma2=sigma**2 !!!!
     val_type *alpha,*alpha2;           // 01-02-2017 alpha2=alpha**2 !!!!      
     val_type *hDQ2T,*hDT;              // 01-02-2017
     val_type *BQ1,*BQ2;                // 01-02-2017
     val_type *hBAB;			// 10-09-2016 Newton
     int *orderedindices;	        // ascending order indices for bi coeficients

   } gauss_method;


typedef struct solution
  { 
     val_type *uu,*ee;

  } solution;


typedef struct toptions
  {
     val_type *rtol,*atol;		  
     int sampling;
     int rdigits,mrdigits;
     char filename[STRMAX];    			// Output filename.
     void (*TheOutput)();       		// Output function.
     void (*StageInitFn)();
     void (*IRKNEWTON_Step_Fn) ();
   } toptions; 



typedef struct parameters
   { 
     int numrpar;				 // Number of real parameters.
     val_type *rpar;				 // Variables for specifying odefun real parameters.
     int numipar;				 // Number of int parametes. 
     int *ipar;     			         // Variables for specifying odefun integer parameters.
     int eval;	                                 // specify which part of differential equation 
                                                 //  must be evaluated
   } parameters;



typedef struct ode_sys       
   {
     int neq;					// number of equations.
     void (*f)();				// odefun.
     void (*jac)();				// jacobian.
     val_type (*ham)();				// hamiltonian
     int cod[2];				// tO specify which part of ODE must be evaluated.
     parameters params;

    } ode_sys;



typedef struct solver_stat
    {

    int convergence;				  // SUCCESS or FAIL. 			        
    bool laststep;                                

    /* auxiliar*/
    val_type *z,*li,*lit0,*lik;
    val_type *jac;				  //01-09-2016  jac=df/dy(y_n)
    val_type *jac_i;				  //01-11-2016  jac_i=df/dy(Y_i), i=1,...,s.
    val_type *INVIDJ2;                            //01-09-2016
    val_type *DL;				  //21-07-2016

    val_type E0;     		      		  // Initial Energy
    val_type MaxDE;		       		  // MaxDE=Abs(Ei-E0/E0)    
		 

    /* stadistics*/
    int stepcount;
    int itcount;
    long int totitcount;
    long int totitcountzero;	     		// 28-11-2016.
    int maxitcount;			
    int itzero;	

    int totitDL[3];                             // 08-03-2017
    int itzeroDL[3];                            // 08-03-2017		
    int fcn;
    int *initqlty;				// ns*neq matrix. Quality of initialization of Li stages   
    int nout;                                   // number of output values.

    } solver_stat; 

  
	
#endif /*DEFH_*/
