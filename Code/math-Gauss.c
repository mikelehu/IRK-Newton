/*----------------------------------------------------------------------------*/
/*									      */
/*                                math-Gauss.c				      */
/*									      */
/* ---------------------------------------------------------------------------*/

#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <mathlink.h>
#include <def.h>
#include <stdio.h>
#include <GaussUserProblem.h> 
#include <GaussCommon.h>
#include <GaussCoefficients.h>
#include <time.h>
#include <stdbool.h>

/* Global variables */

int thread_count;


void mathGauss (int neq, int n, int ns, double t0, double tend,
                double *u0,long ulen, double *e0,long elen, 
                double h, double *rpar,long rlen, int *ipar,long ilen,
                int approximation, int threads,int algorithm,               
                const char *myfilename,int sampling,int codfun)

{

/*------ declarations --------------------------------------------------------*/

    int i,is,initmean,totmean;
    solution u; 
    solver_stat thestat;
    gauss_method gsmethod;
    ode_sys system;
    parameters params;
    toptions options;
    int aux[2];
 
    u.uu = (val_type *)malloc(ulen*sizeof(val_type));
    u.ee = (val_type *)malloc(ulen*sizeof(val_type));


    options.rtol=malloc(ulen*sizeof(val_type));
    options.atol=malloc(ulen*sizeof(val_type));


    clock_t clock0, clock1; 
    time_t  wtime0,wtime1;


/* ----------- implementation  -----------------------------------------------*/    

    params.rpar =(val_type *)malloc(rlen*sizeof(val_type));
    params.ipar =(int *)malloc(ilen*sizeof(int));
    params.numrpar=rlen;
  
    thread_count=threads;
    options.t0=t0;
    options.t1=tend;
    options.h=h;
    options.algorithm=algorithm;
    options.sampling=sampling;
    options.approximation=approximation;
    gsmethod.ns=ns;
    gsmethod.mm=(ns+1)/2;
    gsmethod.smm=ns-gsmethod.mm;
    system.neq=neq;
    system.n=n;

    for (i=0; i<rlen; i++) params.rpar[i]=rpar[i];
    for (i=0; i<ilen; i++) params.ipar[i]=ipar[i];

    select_odefun (codfun,&system);

    system.params.rpar=&params.rpar[0];
    system.params.ipar=&params.ipar[0];

    for (i=0; i<system.neq; i++)        
    {
          options.rtol[i]=RTOL;
          options.atol[i]=ATOL;
    }
  

    GaussCoefficients (DIR_MATH,&gsmethod,&options);
    GaussCoefficientsNewton (DIR_MATH,&system,&gsmethod,&options);

    strncpy(thestat.filename, myfilename,STRMAX); 
    InitStat(&system,&gsmethod,&thestat);

    for (i=0; i<neq;i++)
    { 
          u.uu[i]=u0[i];
          u.ee[i]=e0[i];					
    }


    wtime0= time(NULL);
    clock0= clock();

    system.cod[0]=0;
    select_gauss(&gsmethod, &u, &system, &options,&thestat);

    clock1=clock();
    wtime1= time(NULL);


/* --------- Results ---------------------------------------------------------*/

    totmean=100;
    if (thestat.stepcount>1)
       for (is=0; is<gsmethod.ns; is++)
       {
            for (i=0; i<system.neq; i++) 
            {
                  initmean=(int)thestat.initqlty[is*system.neq+i]/(thestat.stepcount-1);
                  if (initmean < totmean) totmean=initmean;
            }               
       }
    else totmean=0;

    switch (options.algorithm)
    {  
    case 1: case 2: case 3: 

          if (!MLPutFunction(stdlink, "List",16)) MLErrorMessage(stdlink);
          if (!MLPutReal64( stdlink,thestat.MaxDE)) MLErrorMessage(stdlink); 
          if (!MLPutReal64List( stdlink, u.uu, neq))  MLErrorMessage(stdlink); 
          if (!MLPutReal64List( stdlink, u.ee, neq))  MLErrorMessage(stdlink); 
          if (!!MLPutInteger32( stdlink,thestat.stepcount)) 
             MLErrorMessage(stdlink);
          if (!!MLPutInteger32( stdlink,thestat.totitcount)) 
             MLErrorMessage(stdlink); 
          if (!!MLPutInteger32( stdlink,thestat.maxitcount)) 
             MLErrorMessage(stdlink); 
          if (!!MLPutInteger32( stdlink,thestat.itzero)) MLErrorMessage(stdlink); 
          if (!!MLPutInteger32( stdlink,thestat.fcn)) MLErrorMessage(stdlink); 
          if (!MLPutReal64( stdlink,(float) (clock1 - clock0)/CLOCKS_PER_SEC)) 
             MLErrorMessage(stdlink);
          if (!MLPutReal64( stdlink,(float) (wtime1 - wtime0))) 
             MLErrorMessage(stdlink); 
          if (!MLPutInteger32( stdlink,thestat.convergence)) 
             MLErrorMessage(stdlink);  
          if (!MLPutInteger32( stdlink,thestat.nout)) 
             MLErrorMessage(stdlink); 
          if (!MLPutInteger32( stdlink,totmean)) 
             MLErrorMessage(stdlink);    
          if (!MLPutInteger32( stdlink,thestat.totitcountzero)) 
             MLErrorMessage(stdlink);
          if (!!MLPutInteger32List (stdlink,thestat.totitDL,3)) MLErrorMessage(stdlink); 
          if (!!MLPutInteger32List (stdlink,thestat.itzeroDL,3)) MLErrorMessage(stdlink); 
          if(! MLEndPacket(stdlink)) MLErrorMessage(stdlink);
          if(! MLFlush(stdlink)) MLErrorMessage(stdlink); 

    break;


    default:
           printf("Error: incorrect algorithm\n");
           
    break;

     }

    free(params.rpar); free(params.ipar);
    free(u.uu); free(u.ee);

    free(options.rtol); free(options.atol);

    free(gsmethod.m);
    free(gsmethod.c); 
    free(gsmethod.b); 
    free(gsmethod.a);
    free(gsmethod.hc); 
    free(gsmethod.hb);      
    free(gsmethod.nu); 
    free(gsmethod.munu);     

    free(gsmethod.Q1);
    free(gsmethod.Q2);
    free(gsmethod.sigma2);
    free(gsmethod.alpha);
    free(gsmethod.alpha2);
    free(gsmethod.hDQ2T);
    free(gsmethod.hDT);
    free(gsmethod.BQ1);
    free(gsmethod.BQ2);
    free(gsmethod.hBAB);

    free(gsmethod.orderedindices);

    free(thestat.z); 
    free(thestat.li);
    free(thestat.lit0);
    free(thestat.lik);
    free(thestat.jac);
    free(thestat.jac_i);
    free(thestat.INVIDJ2);
    free(thestat.DL);
    free(thestat.initqlty);      
    


}


int main(int argc, char *argv[])
{
     return MLMain(argc, argv);
}


