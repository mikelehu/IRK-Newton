/*------------------------------------------------------------------------------*/
/*										*/
/*                                GaussTerminal.c				*/
/*										*/
/* -----------------------------------------------------------------------------*/

#include <GaussTerminal.h>

/* Global variables */

int thread_count=1;

int main()
{

/*------ declarations --------------------------------------------------*/    

     int i,is,initmean,totmean,codfun;   
	
     gauss_method gsmethod;
     solution u;
     toptions options;
     parameters params;
     ode_sys system;   
     solver_stat thestat;

     clock_t clock0, clock1; 
     time_t  wtime0,wtime1;

     params.rpar =(val_type *)malloc(MAXPARAM*sizeof(val_type));
     params.ipar =(int *)malloc(MAXPARAM*sizeof(int));

     u.uu = (val_type *)malloc(MAXNEQ*sizeof(val_type));
     u.ee = (val_type *)malloc(MAXNEQ*sizeof(val_type));

     options.rtol=malloc(MAXNEQ*sizeof(val_type));
     options.atol=malloc(MAXNEQ*sizeof(val_type));
  
#    if PREC ==2  //QUADRUPLEPRECISION
     int n;
     int width = 46;
     char buf[128];
#    endif

/* ----------- implementation  --------------------------------------------*/     


/* ----------- integration parameters-------------------------------------*/  

     gsmethod.ns = 6;     				 //	Stages.
     gsmethod.mm=(gsmethod.ns+1)/2;
     gsmethod.smm=gsmethod.ns-gsmethod.mm;

     strncpy(thestat.filename, "Output.bin",STRMAX);     //     Output filename.

     options.approximation=0;   			 //     Approximation: Y^[0] (GaussCommon.c/Li_init()).

     options.algorithm=3;    				 //	   
     options.rdigits=0;       

/* ----------- problem parameters-----------------------------------------*/

/* Double pendulum stiff */ 

     codfun=1;						 //      Ode system of Double Pendulum Stiff.
     system.problem =3;  				 //	 Initial values : DP Stiff C=0.  (GaussInitData.c).
     system.problem =4;  				 //	 Initial values : DP Stiff C=32. (GaussInitData.c).
     system.problem =5;  				 //	 Initial values : DP Stiff C=64. (GaussInitData.c).

     options.h = POW(2,-7);	       			 //	 Stepsize.
     options.sampling=POW(2,10);   			 	 


/* ----------- execution  ----------------------------------------------------*/

    printf("Begin execution \n");


    InitialData (&options,&u,&params,&system); 

    select_odefun (codfun,&system);

    system.params.rpar=&params.rpar[0];
    system.params.ipar=&params.ipar[0];

    printf("method=%i, problem=%i, algorithm=%i\n",
            gsmethod.ns,system.problem,options.algorithm);
    printf("approximation=%i,sampling=%i\n",
            options.approximation,options.sampling);

    printf("options.h=%lg\n", options.h);        
    printf("----------------\n");         
                     
     for (i=0; i<system.neq; i++)
     {
          options.rtol[i]=RTOL;
          options.atol[i]=ATOL;
     }

     if (options.rdigits>0) options.mrdigits=pow(2,options.rdigits);

     GaussCoefficients(DIR_TERM,&gsmethod,&options);  
     GaussCoefficientsNewton(DIR_TERM,&system,&gsmethod,&options);

     InitStat(&system,&gsmethod,&thestat);	


     wtime0= time(NULL);
     clock0= clock();

     system.cod[0]=0;
     select_gauss(&gsmethod, &u, &system, &options,&thestat);
  
     clock1=clock();
     wtime1= time(NULL);

/* --------- Results ---------------------------------------------------------*/
     printf("End execution \n");   


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
     printf("\ntotmean=%i\n",totmean);

        
     if (thestat.convergence==SUCCESS)
     {
           printf("Execution Correct\n");
           printf("convergence=%i.  (=0 correct;=-1 incorrect)\n",thestat.convergence);        

           print_u (system.neq,u.uu);
           printf("Energy MaxDE=%.20lg\n",thestat.MaxDE);

           printf ("\nCPU time:%lg\n", (double) (clock1 - clock0)/CLOCKS_PER_SEC);
	   printf ("Elapsed wall clock time: %ld\n", (wtime1 - wtime0));
           printf ("Elapsed wall clock time: %lg\n", difftime(wtime1, wtime0));

           printf("\n");
           
           printf("stepcount=%i\n",thestat.stepcount);
           printf("nout=%i\n",thestat.nout);
           
           printf("fixed-point iterations=%i\n",thestat.totitcount);
           printf("max fixed-point iterations=%i\n",thestat.maxitcount);        
           printf("deltaZero iterations=%i\n",thestat.itzero); 
 
           printf("\n");

           printf ("function ebaluations\n");
           printf("fcn=%i\n",thestat.fcn);
           printf ("\n");

           printf ("Quality init\n");
           


     }
     else
     {
           printf("Execution InCorrect\n");
           printf("convergence=%i.  (=0 correct;=-1 incorrect)\n",thestat.convergence);
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

     free(thestat.z); 
     free(thestat.li);
     free(thestat.lit0);
     free(thestat.lik);
     free(thestat.jac);
     free(thestat.jac_i);
     free(thestat.INVIDJ2);
     free(thestat.DL);
     free(thestat.initqlty);      
    
    
     exit(0);

}





