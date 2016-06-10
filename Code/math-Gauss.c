/*------------------------------------------------------------------------------*/
/*										*/
/*                                math-Gauss.c					*/
/*										*/
/* -----------------------------------------------------------------------------*/

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
#include <GaussUserProblem0.h> 
#include <GaussCommon0.h>

/* Global variables */

int thread_count;


void mathGauss (int neq, int n, int ns, double t0, double tend,
                double *u0,long ulen, double *e0,long elen, 
                double h, double *rpar,long rlen, int *ipar,long ilen,
                int approximation, int threads,int algorithm,
                int rdigits1,int rdigits2,
                const char *myfilename,int sampling,int codfun)

{

//------ declarations --------------------------------------------------*

    int i,is,initmean,totmean;
    solution u,u2; 
    solver_stat thestat,thestat2;
    gauss_method gsmethod,gsmethod2;
    ode_sys system;
    parameters params;
    toptions options,options2;
    int aux[2];
 
    u.uu = (val_type *)malloc(ulen*sizeof(val_type));
    u.ee = (val_type *)malloc(ulen*sizeof(val_type));
    u2.uu = (val_type *)malloc(ulen*sizeof(val_type));
    u2.ee = (val_type *)malloc(ulen*sizeof(val_type));

    options.rtol=malloc(ulen*sizeof(val_type));
    options.atol=malloc(ulen*sizeof(val_type));
    options2.rtol=malloc(ulen*sizeof(val_type));
    options2.atol=malloc(ulen*sizeof(val_type));

    clock_t clock0, clock1; 
    time_t  wtime0,wtime1;

    options.rdigits=rdigits1;    
    options2.rdigits=rdigits2;   


// ----------- implementation  --------------------------------------*     

    params.rpar =(val_type *)malloc(rlen*sizeof(val_type));
    params.rpar0 =(val_type0 *)malloc(rlen*sizeof(val_type0));
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
    system.neq=neq;
    system.n=n;

    for (i=0; i<rlen; i++)
    {
          params.rpar[i]=rpar[i];
          params.rpar0[i]=rpar[i];
    }
    for (i=0; i<ilen; i++) params.ipar[i]=ipar[i];

    switch (codfun)
    { 
    case 1: 
          system.f = Ode1;
          system.ham= Ham1;
          system.jac= Jac1;     
          system.f0 = Ode1low;
          system.ham0= Ham1low;
          system.jac0= Jac1low;                        
    break;

    case 2: 
          system.f = Ode2;
          system.ham= Ham2;
          system.jac= Jac2;
          system.f0 = Ode2low;
          system.ham0= Ham2low;
          system.jac0= Jac2low;  
    break;

    case 3: 
          system.f = Ode3;
          system.ham= Ham3;
          system.jac= Jac3;       
          system.f0 = Ode3low;
          system.ham0= Ham3low;
          system.jac0= Jac3low;                 
    break;

    case 4: 
          system.f = Ode4;
          system.ham= Ham4;
          system.jac= Jac4;
          system.f0 = Ode4low;
          system.ham0= Ham4low;
          system.jac0= Jac4low;                        
    break;

    case 5: 
          system.f = Ode5;
          system.ham= Ham5;
          system.jac= Jac5;         
          system.f0 = Ode5low;
          system.ham0= Ham5low;
          system.jac0= Jac5low;               
    break;

    case 6: 
          system.f = Ode6;
          system.ham= Ham6;
          system.jac= Jac6;          
          system.f0 = Ode6low;
          system.ham0= Ham6low;
          system.jac0= Jac6low;              
    break;

    case 7: 
          system.f = Ode7;
          system.ham= Ham7;
          system.jac= Jac7;  
          system.f0 = Ode7low;
          system.ham0= Ham7low;
          system.jac0= Jac7low;                      
    break;

    case 8: 
          system.f = Ode8;
          system.ham= Ham8;
          system.jac= Jac8;    
          system.f0 = Ode8low;
          system.ham0= Ham8low;
          system.jac0= Jac8low;                    
    break;

    case 9: 
          system.f = Ode9;
          system.ham= Ham9;
          system.jac= Jac9;
          system.f0 = Ode9low;
          system.ham0= Ham9low;
          system.jac0= Jac9low;                       
    break;

    case 10: 
          system.f = Ode10;
          system.ham= Ham10;
          system.jac= Jac10; 
          system.f0 = Ode10low;
          system.ham0= Ham10low;
          system.jac0= Jac10low;                       
    break;

    case 11: 
          system.f = Ode11;
          system.ham= Ham1;
    break;

    case 12: 
          system.f = Ode12;
          system.ham= Ham2;
    break;

    default:
          printf("error. codfun\n");
    break;
    }
 
    system.params.rpar=&params.rpar[0];
    system.params.ipar=&params.ipar[0];

    for (i=0; i<params.numrpar; i++) params.rpar0[i]=params.rpar[i];
    system.params.rpar0=&params.rpar0[0];

    for (i=0; i<system.neq; i++)        
    {
          options.rtol[i]=RTOL;
          options.atol[i]=ATOL;
    }
  

    if (options.rdigits>0) options.mrdigits=pow(2,options.rdigits);
    if (options2.rdigits>0) options2.mrdigits=pow(2,options2.rdigits);

    GaussCoefficients(&gsmethod,&options);
    strncpy(thestat.filename, myfilename,STRMAX); 
    InitStat(&system,&gsmethod,&thestat);

    for (i=0; i<neq;i++)
    { 
          u.uu[i]=u0[i];
          u.ee[i]=e0[i];					
    }


    wtime0= time(NULL);
    clock0= clock();


    select_gauss(&gsmethod, &gsmethod2, &u, &u2,
           &system, &options, &options2,
           &thestat, &thestat2);

    clock1=clock();
    wtime1= time(NULL);

    totmean=100;
    for (is=0; is<gsmethod.ns; is++)
    {
            for (i=0; i<system.neq; i++) 
            {
                  initmean=(int)thestat.initqlty[is*system.neq+i]/(thestat.stepcount-1);
                  if (initmean < totmean) totmean=initmean;
            }
                
    }

    switch (options.algorithm)
    {  
    case 1: case 2: case 11 : case 12: case 13:

          if (!MLPutFunction(stdlink, "List",13)) MLErrorMessage(stdlink);
          if (!MLPutReal64( stdlink,thestat.MaxDE)) MLErrorMessage(stdlink); 
          if (!MLPutReal64List( stdlink, u.uu, neq))  MLErrorMessage(stdlink); 
          if (!MLPutReal64List( stdlink, u.ee, neq))  MLErrorMessage(stdlink); 
          if (!!MLPutInteger32( stdlink,thestat.stepcount)) MLErrorMessage(stdlink);
          if (!!MLPutInteger32( stdlink,thestat.totitcount)) MLErrorMessage(stdlink); 
          if (!!MLPutInteger32( stdlink,thestat.maxitcount)) MLErrorMessage(stdlink); 
          if (!!MLPutInteger32( stdlink,thestat.itzero)) MLErrorMessage(stdlink); 
          if (!!MLPutInteger32( stdlink,thestat.fcn)) MLErrorMessage(stdlink); 
          if (!MLPutReal64( stdlink,(float) (clock1 - clock0)/CLOCKS_PER_SEC)) MLErrorMessage(stdlink);
          if (!MLPutReal64( stdlink,(float) (wtime1 - wtime0))) MLErrorMessage(stdlink); 
          if (!MLPutInteger32( stdlink,thestat.convergence)) MLErrorMessage(stdlink);  
          if (!MLPutInteger32( stdlink,thestat.nout)) MLErrorMessage(stdlink); 
          if (!MLPutInteger32( stdlink,totmean)) MLErrorMessage(stdlink);    
          if(! MLEndPacket(stdlink)) MLErrorMessage(stdlink);
          if(! MLFlush(stdlink))  MLErrorMessage(stdlink); 

    break;

    case 21: case 22:

          aux[0]=thestat.totitcount;
          aux[1]=thestat2.totitcount;

          if (!MLPutFunction(stdlink, "List",13)) MLErrorMessage(stdlink);
          if (!MLPutReal64( stdlink,thestat.MaxDE)) MLErrorMessage(stdlink); 
          if (!MLPutReal64List( stdlink, u.uu, neq)) MLErrorMessage(stdlink); 
          if (!MLPutReal64List( stdlink, u.ee, neq))  MLErrorMessage(stdlink);
          if (!!MLPutInteger32( stdlink,thestat.stepcount)) MLErrorMessage(stdlink);
          if (!!MLPutInteger32List (stdlink,aux,2)) MLErrorMessage(stdlink); 
          if (!!MLPutInteger32( stdlink,thestat.maxitcount)) MLErrorMessage(stdlink); 
          if (!!MLPutInteger32( stdlink,thestat.itzero)) MLErrorMessage(stdlink);
          if (!!MLPutInteger32( stdlink,thestat.fcn)) MLErrorMessage(stdlink); 
          if (!MLPutReal64( stdlink,(float) (clock1 - clock0)/CLOCKS_PER_SEC)) MLErrorMessage(stdlink); 
          if (!MLPutReal64( stdlink,(float) (wtime1 - wtime0))) MLErrorMessage(stdlink); 
          if (!MLPutInteger32( stdlink,thestat.convergence)) MLErrorMessage(stdlink);   
          if (!MLPutInteger32( stdlink,thestat.nout)) MLErrorMessage(stdlink);  
          if (!MLPutInteger32( stdlink,totmean)) MLErrorMessage(stdlink); 
          if(! MLEndPacket(stdlink)) MLErrorMessage(stdlink);
          if(! MLFlush(stdlink)) MLErrorMessage(stdlink);
    break;

    default:
           printf("Error: incorrect algorithm\n");
           
    break;

     }

    free(params.rpar); free(params.ipar);
    free(params.rpar0);
    free(u.uu); free(u.ee);
    free(u2.uu); free(u2.ee);

    free(options.rtol); free(options.atol);
    free(options2.rtol); free(options2.atol);

    free(gsmethod.m); 
    free(gsmethod.b); 
    free(gsmethod.hb);      
    free(gsmethod.c); 
    free(gsmethod.hc); 
    free(gsmethod.nu); 
    free(gsmethod.munu); 
    free(gsmethod.orderedindices);

    free(thestat.z); 
    free(thestat.li);
    free(thestat.liold);
    free(thestat.lit0);
    free(thestat.initqlty);

}


int main(int argc, char *argv[])
{
     return MLMain(argc, argv);
}


