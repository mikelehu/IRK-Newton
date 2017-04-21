/*----------------------------------------------------------------------------*/
/*									      */
/*    GaussCommon.c							      */
/*									      */
/*	Functions: 							      */
/*	 print_u():							      */
/*	 InitStat():							      */
/*	 NormalizedDistance():						      */
/*	 RemoveDigitsFcn():						      */
/*	 Default_Stage_init();						      */
/*	 Interpolated_Stage_init();					      */
/*	 Yi_update():							      */
/*	 Yi_update_Classic():						      */
/*	 statlinit():							      */
/*	 StopCriterion();						      */
/*       StopCriterionFloat();						      */
/*	 NSS_Step():							      */
/*	 NSS_Step_plus():						      */
/*	 NSS_MIX_Step():						      */
/*       MatAdd();                                                            */
/*       MatAddID();                                                          */
/*	 MMNewtonSS():							      */
/*	 NSS_Solve():							      */
/*	 NSS_Solve_plus():						      */
/*       Compute_MM();                                                        */
/*       Compute_R();                                                         */
/*       Compute_Z();                                                         */
/*       Compute_W1();                                                        */
/*       Compute_W2();                                                        */
/*       Compute_DL();                                                        */
/*       Compute_GG();                                                        */
/*	 MyOutput():							      */
/*	 CompensatedSummation():					      */
/*	 IRKNEWTON():							      */
/*                                                                            */
/* ---------------------------------------------------------------------------*/

#include <Common-IRKNEWTON.h>

/******************************************************************************/
/* 					   				      */
/* print_u 		         	  				      */
/*                                        				      */
/*									      */
/******************************************************************************/

void print_u (int neq,val_type *u)
{
     int i;


     for (i = 0;i<neq;i++) 
     {
           if (i<neq-1)
                 printf("%.20lg,", u[i]);
           else
                 printf("%.20lg\n", u[i]);
     }     


}


/******************************************************************************/
/* 					   				      */
/* InitStat 		         	  				      */
/*                                        				      */
/*									      */
/******************************************************************************/
      
void InitStat (ode_sys *system,gauss_method *gsmethod, solver_stat *thestatptr)
{

     int i,is,neq,ns,mm;

     ns=gsmethod->ns;
     neq=system->neq;
     mm=gsmethod->mm;     

     thestatptr->laststep = false;    
     thestatptr->stepcount = 0;
     thestatptr->itcount=0;
     thestatptr->totitcount=0;
     thestatptr->totitcountzero=0;
     thestatptr->maxitcount=0;
     thestatptr->itzero=0;
     thestatptr->fcn=0;
     thestatptr->convergence=SUCCESS;
     thestatptr->nout=0;
     thestatptr->MaxDE=0.;

     for (i=0; i<3; i++)
     {
          thestatptr->totitDL[i]=0;
          thestatptr->itzeroDL[i]=0;

     }
      
     thestatptr->z = (val_type *)malloc(neq*ns*sizeof(val_type));
     thestatptr->li = (val_type *)malloc(neq*ns*sizeof(val_type));
     thestatptr->lit0 = (val_type *)malloc(neq*ns*sizeof(val_type));
     thestatptr->initqlty = (int *)malloc(neq*ns*sizeof(int));
     thestatptr->lik = (val_type *)malloc(neq*ns*sizeof(val_type));
     thestatptr->DL = (val_type *)malloc(neq*ns*sizeof(val_type));
     thestatptr->jac = (val_type *)malloc(neq*neq*sizeof(val_type));
     thestatptr->INVIDJ2= (val_type *) malloc(mm*(neq*neq)*sizeof(val_type));
     thestatptr->jac_i = (val_type *)malloc(ns*(neq*neq)*sizeof(val_type));

    
     
     for (is=0; is<ns; is++)
          for (i=0; i<neq; i++)
          {
               thestatptr->li[is*neq+i]=0.;
               thestatptr->initqlty[is*neq+i]=0;
          }


     return;

}


/******************************************************************************/
/* 					   				      */
/* NormalizedDistance: 	         	  				      */
/*                                        				      */
/*									      */
/******************************************************************************/

val_type NormalizedDistance (int neq,int ns,
                             toptions *options,val_type *z,val_type *zold)
{


/*---------------- Declarations ----------------------------------------------*/

     int i,is,ix;
     val_type maxi,mayi,relerrors;
     val_type maxz,maxzold;

/* --------------- Implementation --------------------------------------------*/
   	
     maxi=0.;

     for (i=0; i<neq; i++)
     {  
          maxz=0.;
          maxzold=0.;
            
          for (is=0; is<ns;is++)
          {
               ix=neq*is+i;
               maxz=FMAX(maxz,FABS(z[ix]));
               maxzold=FMAX(maxzold,FABS(zold[ix]));
          }

          mayi=(maxz+maxzold)/2;
          relerrors=0.;
      
          for (is=0; is<ns; is++)
          {
               ix=neq*is+i;
               relerrors=FMAX(relerrors,
                         FABS(z[ix]-zold[ix])/(mayi*options->rtol[i]+options->atol[i]));
          }

          maxi=FMAX(maxi,relerrors);

     }

     return maxi;

}

/******************************************************************************/
/* 					   				      */
/* RemoveDigitsFcn	         	  				      */
/*                                        				      */
/*									      */
/******************************************************************************/

void RemoveDigitsFcn (val_type *x, int m)
{


/*---------------- Declarations ----------------------------------------------*/

     val_type aux,mx,mxx;

/* --------------- Implementation --------------------------------------------*/

     aux=*x;
     mx=m*aux;
     mxx=mx+aux;
     *x=mxx-mx;

     return;

}

/******************************************************************************/
/* 					   				      */
/* Default_Stage_init: 		         	  				      */
/*                                        				      */
/*									      */
/******************************************************************************/

void Default_Stage_init(solution *u, val_type *z,ode_sys *system,
            gauss_method *method,solver_stat *thestatptr, toptions *options)
{


/* ---------- First initializations ------------------------------------------*/

     int neq,ns;

     neq=system->neq;
     ns=method->ns;

/*------------- declarations -------------------------------------------------*/ 

     int i,is;
     val_type *li,*lit0;


/* ----------- implementation  -----------------------------------------------*/

     li=thestatptr->li;
     lit0=thestatptr->lit0;

     for (is = 0; is<ns; is++)
          for (i = 0; i<neq; i++)
          {
               li[is*neq+i]=0.;
               lit0[is*neq+i]=0.;       
           }
 
     return;

}


/******************************************************************************/
/* 					   				      */
/* Interpolated_Stage_init: 			         	  	      */
/*                                        				      */
/*									      */
/******************************************************************************/

void Interpolated_Stage_init(solution *u, val_type *z,ode_sys *system,
            gauss_method *method,solver_stat *thestatptr, toptions *options)
{


/* ---------- First initializations ------------------------------------------*/

     int neq,ns;

     neq=system->neq;
     ns=method->ns;

/*------------- declarations -------------------------------------------------*/ 

     int i,is,in,js,jsn;
     val_type *coef,*li,*lit0;
     val_type sum;

     val_type liold[neq*ns];

/* ----------- implementation  -----------------------------------------------*/
 

     for (is=0; is<ns;is++)
            for (i=0; i<neq; i++) liold[neq*is+i]=thestatptr->li[neq*is+i];

     li=thestatptr->li;
     lit0=thestatptr->lit0;
     coef=method->munu;     
  
     for (is = 0; is<ns; is++)
     {
          for (i = 0; i<neq; i++)
          {  
                in=neq*is+i;
                sum=0.;

                for (js = 0; js<ns; js++)
                {
                    jsn=ns*is+js;
                    sum+=liold[neq*js+i]*(coef[jsn]);
                }

                li[in]=sum;
                lit0[in]=li[in];
           }
     }  

  
     return;

}



/******************************************************************************/
/* 					   				      */
/* Yi_update: 		         	  				      */
/*                                        				      */
/*									      */
/******************************************************************************/

int Yi_update(solution *u, val_type *z,ode_sys *system,
            gauss_method *method,solver_stat *thestatptr)
{


/* ---------- First initializations ------------------------------------------*/

     int neq,ns;

     neq=system->neq;
     ns=method->ns;

/*------------- declarations -------------------------------------------------*/ 

     int i,is,in,js,jsn;
     val_type *coef,*li;
     val_type sum;

/* ----------- implementation  -----------------------------------------------*/


    li=thestatptr->li;
    coef=method->m;


    for (is = 0; is<ns; is++)
    {
                for (i = 0; i<neq; i++)
                {  
                     in=neq*is+i;
                     sum=coef[ns*is]*li[i]+u->ee[i];

                     for (js = 1; js<ns; js++)
                     {
                          jsn=ns*is+js;
                          sum+=(coef[jsn])*li[neq*js+i];
                     }

                     z[in]=u->uu[i]+sum;
                }
     }  
 
   
     return(0);

}


/******************************************************************************/
/* 					   				      */
/* Yi_update_Classic: 		    	  				      */
/*                                        				      */
/*									      */
/******************************************************************************/

int Yi_update_Classic (solution *u, val_type *z,ode_sys *system,
                       gauss_method *method,solver_stat *thestatptr)
{


/* ---------- First initializations ------------------------------------------*/

     int neq,ns;

     neq=system->neq;
     ns=method->ns;

/*------------- declarations -------------------------------------------------*/ 

     int i,is,in;
     int lda, ldb, ldc;
     val_type *coef,*li;
     val_type alpha, beta;

/* ----------- implementation  -----------------------------------------------*/


    li=thestatptr->li;
    coef=method->m;


    for (is = 0; is<ns; is++)
    {
                for (i = 0; i<neq; i++)
                {  
                     in=neq*is+i;                    
                     z[in]=u->uu[i];
                }
     }  


     alpha = 1.; beta = 1.;
     lda=ns; ldb=neq; ldc=neq;   

     GEMM(CblasRowMajor,CblasNoTrans,CblasNoTrans, ns,neq,ns,
          alpha, coef, lda, li, ldb, beta, z, ldc );

   
     return(0);

}




/******************************************************************************/
/* 					   				      */
/* statlinit: 		         	  				      */
/*                                        				      */
/*									      */
/******************************************************************************/

int statlinit (ode_sys *system,gauss_method *method,solver_stat *thestatptr)
{


/* ---------- First initializations ------------------------------------------*/

     int neq,ns;

     neq=system->neq;
     ns=method->ns;

/*------------- declarations -------------------------------------------------*/ 

     int i,is,in;
     val_type *li,*lit0;
     int *initqlty;
     int prec;

/* ----------- implementation  -----------------------------------------------*/


     li=thestatptr->li;
     lit0=thestatptr->lit0;
     initqlty=thestatptr->initqlty;
     prec=16;

     for (is = 0; is<ns; is++)
           for (i = 0; i<neq; i++)
           {                      
                 in=neq*is+i;
                 if (li[in]==lit0[in]) {initqlty[in]+=prec;}
                       else 
                       {
                        initqlty[in]+=(int)round(fabs(log10(fabs(li[in]-lit0[in])/fabs(lit0[in]))));
                       }
           }
                    
     return(0);

}



/******************************************************************************/
/* 									      */
/*   StopCriterion  						              */
/* 									      */
/* 									      */
/******************************************************************************/

void StopCriterion (ode_sys *system, gauss_method *method,
                    int *D0,bool *cont,val_type *DMin, 
                    val_type *Y, val_type *Yold)
{
/*----------------  First: initialization ------------------------------------*/
     int neq,ns;

     neq=system->neq;
     ns=method->ns;

/*------ declarations --------------------------------------------------------*/
     int i,is,isn;
     val_type dY;

/* ----------- implementation  -----------------------------------------------*/ 


     bool plusIt;

     if (*D0<0) plusIt=false;
         else plusIt=true;

     *D0=0;
     *cont=false;


     for (is=0; is<ns; is++)
     {
          for (i=0; i<neq; i++) 
          { 
               isn=is*neq+i;
               dY=FABS(Y[isn]-Yold[isn]);

               if (dY>0.)
               {                     
                   if (dY<DMin[isn])
                    {
                         DMin[isn]=dY;
                         *cont=true;
                    }

               }
               else
               {
                    *D0=*D0+1;
               } 
          }
     
     }

  
    if (*cont==false && *D0<(ns*neq) && plusIt)
    {
          *D0=-1;
          *cont=true;
    }


     return;

}



/******************************************************************************/
/* 									      */
/*   StopCriterionFloat						              */
/* 									      */
/* 									      */
/******************************************************************************/

void StopCriterionFloat (ode_sys *system, gauss_method *method,
                         int *D0,bool *cont,val_type *DMin, 
                         val_type *Y, val_type *Yold)
{
/*----------------  First: initialization ------------------------------------*/
     int neq,ns;

     neq=system->neq;
     ns=method->ns;

/*------ declarations --------------------------------------------------------*/
     int i,is,isn;
     float dY,YY,YYold;

/* ----------- implementation  -----------------------------------------------*/ 


     bool plusIt;

     if (*D0<0) plusIt=false;
         else plusIt=true;

     *D0=0;
     *cont=false;


     for (is=0; is<ns; is++)
     {
          for (i=0; i<neq; i++) 
          { 
               isn=is*neq+i;
               YY=Y[isn];
               YYold=Yold[isn];
               dY=FABS(YY-YYold);

               if (dY>0.)
               {                     
                   if (dY<DMin[isn])
                    {
                         DMin[isn]=dY;
                         *cont=true;
                    }

               }
               else
               {
                    *D0=*D0+1;
               } 
          }
     
     }

  
    if (*cont==false && *D0<(ns*neq) && plusIt)
    {
          *D0=-1;
          *cont=true;
    }


     return;

}



/******************************************************************************/
/* 							          	      */
/*   NSS_Step    							      */
/* 									      */
/* 									      */
/******************************************************************************/


void NSS_Step   ( ode_sys *system, solution *u, val_type tn,val_type h, 
                  toptions *options,gauss_method *method,solver_stat *thestatptr)

{ 

/*----------------  First: initialization ------------------------------------*/
     int neq,ns;

     neq=system->neq;
     ns=method->ns;

/*------ declarations --------------------------------------------------------*/

     int i,is;
     bool iter0;
     int D0;			
     val_type difftest,DMin[neq*ns];
     val_type *z;
     val_type zold[neq*ns];
     val_type *MM;

     int n,m,lda,info;
     int ipiv[neq*ns];

     
/* ----------- implementation  -----------------------------------------------*/ 

     MM=(val_type *)malloc((neq*(ns))*(neq*(ns))*sizeof(val_type));
     z=thestatptr->z;

     iter0=true;
     D0=0;
     for (is=0; is<ns; is++)
          for (i=0; i<neq; i++) DMin[is*neq+i]=INF;

     thestatptr->itcount=0;


/* ------LU factorization of a general n-by-m matrix A------------------------*/

     n=neq*(ns);
     m=n;
     lda=n;

     MMNewtonSS (neq,method,tn,h,system,u->uu,MM);

     info = GETRF(LAPACK_ROW_MAJOR, n, m, MM, lda, ipiv);
     if( info > 0 )
     {
	printf( "The diagonal element of the triangular factor of A,\n" );
	printf( "U(%i,%i) is zero, so that A is singular;\n", info, info );
	printf( "the solution could not be computed.\n" );
	exit( 1 );
     }


     while (iter0 && thestatptr->itcount<MAXIT)
     {            
           Yi_update (u,z,system,method,thestatptr);
              
           NSS_Solve (system,u,tn,h,method,thestatptr,MM,ipiv); 
           StopCriterion(system,method,&D0,&iter0,DMin,thestatptr->li,thestatptr->lik);
           thestatptr->itcount++;    

     }


     if (thestatptr->itcount==MAXIT)
     { 
           printf("Break: step(MAXIT)=%i\n",thestatptr->itcount);
           thestatptr->convergence=FAIL; 
     }

     else
           if (D0<(ns*neq))
           {
                 for (is=0; is<ns;is++) for (i=0; i<neq; i++)    
                             zold[neq*is+i]=thestatptr->z[neq*is+i];
                 Yi_update (u,z,system,method,thestatptr);
                 difftest=NormalizedDistance(neq,ns,options,z,zold);
                 if (difftest>1.)
                 {
                       thestatptr->convergence=FAIL;
                       printf("Lack of convegence of Newton iteration:\
                               step=%i,iteration=%i,",
                               thestatptr->stepcount,thestatptr->itcount);  

                       printf("difftest=%lg\n",difftest);                     
                 }
           }
           else
           {
                (thestatptr->totitcountzero)+=(thestatptr->itcount);
                thestatptr->itzero++;
           }


     free(MM);
     return;

}



/******************************************************************************/
/* 							          	      */
/*   NSS_Step_plus							      */
/* 									      */
/* 									      */
/******************************************************************************/


void NSS_Step_plus  ( ode_sys *system, solution *u, val_type tn,val_type h, 
                      toptions *options,gauss_method *method,solver_stat *thestatptr)

{ 

/*----------------  First: initialization ------------------------------------*/
     int neq,ns;

     neq=system->neq;
     ns=method->ns;

/*------ declarations --------------------------------------------------------*/


     int i,is;
     bool iter0;
     int D0;			
     val_type difftest,DMin[neq*ns];
     val_type *z;
     val_type *jac;
     parameters *params;
     val_type zold[neq*ns];

     int info;

     val_type MM[neq*neq];
     int ipiv[neq];

    
/* ----------- implementation  -----------------------------------------------*/ 

     z=thestatptr->z;
     jac=thestatptr->jac;
     params=&system->params;

     iter0=true;
     D0=0;
     for (is=0; is<ns; is++)
          for (i=0; i<neq; i++) DMin[is*neq+i]=INF;

     thestatptr->itcount=0;

     system->jac(neq,tn+h/2,u->uu,jac,params);

     /* -------- Compute M matrix ------------- */

     Compute_MM (system,h,MM,method,thestatptr);
 
     /* -------- LU factorization of (MM) matrix ----------- */

     info = GETRF(LAPACK_ROW_MAJOR, neq, neq, MM, neq, ipiv);
     if( info > 0 )
     {
 	   printf( "The diagonal element of the triangular factor of A,\n" );
	   printf( "U(%i,%i) is zero, so that A is singular;\n", info, info );
	   printf( "the solution could not be computed.\n" );
	   exit( 1 );
     }     


     while (iter0 && thestatptr->itcount<MAXIT)
     { 

           Yi_update (u,z,system,method,thestatptr);
            
           NSS_Solve_plus (system,u,tn,h,method,thestatptr,MM,ipiv); 
           StopCriterion(system,method,&D0,&iter0,DMin,thestatptr->li,thestatptr->lik);
           thestatptr->itcount++;
           
     }


     if (thestatptr->itcount==MAXIT)
     { 
           printf("Break: step(MAXIT)=%i\n",thestatptr->itcount);
           thestatptr->convergence=FAIL; 
     }

     else
           if (D0<(ns*neq))
           {
                 for (is=0; is<ns;is++) for (i=0; i<neq; i++)     
                             zold[neq*is+i]=thestatptr->z[neq*is+i];
                 Yi_update (u,z,system,method,thestatptr);
                 difftest=NormalizedDistance(neq,ns,options,z,zold);
                 if (difftest>1.)
                 {
                       thestatptr->convergence=FAIL;
                       printf("Lack of convegence of Newton iteration:\
                               step=%i,iteration=%i,",
                               thestatptr->stepcount,thestatptr->itcount);  
                       printf("difftest=%lg\n",difftest);                    
                 }
           }
           else
           {
                (thestatptr->totitcountzero)+=(thestatptr->itcount);
                thestatptr->itzero++;
           }


     return;

}


/******************************************************************************/
/* 							          	      */
/*   NSS_MIX_Step							      */
/* 									      */
/* 									      */
/******************************************************************************/


void NSS_MIX_Step   ( ode_sys *system, solution *u, val_type tn,val_type h, 
                      toptions *options,gauss_method *method,solver_stat *thestatptr)

{ 

/*----------------  First: initialization ------------------------------------*/
     int neq,ns;

     neq=system->neq;
     ns=method->ns;
     int mm=method->mm;
     int smm=method->smm;

/*------ declarations --------------------------------------------------------*/


     int i,is,j,k,l,isn;
     bool iter0;
     int D0;			
     val_type difftest,DMin[neq*ns];
     val_type *z,*li,*lik;
     val_type *DL,*jac;
     parameters *params;   

     int info;

     val_type MM[neq*neq];
     int ipiv[neq];

     val_type fz[ns*neq];
     val_type g[ns*neq];
     val_type RR[mm*neq];
     val_type ZZ[neq],W1[mm*neq],W2[smm*neq];
     val_type GG[ns*neq],DLl[ns*neq];  

     val_type DLOLD[ns*neq];
     val_type zold[neq*ns];

     val_type aa,bb;
     int lda, incx,incy;
    
/* ----------- implementation  ----------------------------------------------*/ 

     z=thestatptr->z;
     li=thestatptr->li;
     lik=thestatptr->lik;
     jac=thestatptr->jac;
     params=&system->params;
     DL=thestatptr->DL;


     thestatptr->itcount=0;

     system->jac(neq,tn+h/2,u->uu,jac,params);

     /* -------- Compute M matrix ------------- */

     Compute_MM (system,h,MM,method,thestatptr);
  
     /* -------- LU factorization of (MM) matrix ----------- */

     info = GETRF(LAPACK_ROW_MAJOR, neq, neq, MM, neq, ipiv);
     if( info > 0 )
     {
 	   printf( "The diagonal element of the triangular factor of A,\n" );
	   printf( "U(%i,%i) is zero, so that A is singular;\n", info, info );
	   printf( "the solution could not be computed.\n" );
	   exit( 1 );
     }     

/*------------ 1-Step: Newton-simplified iterations (Float) -----------------*/

     k=0;
     iter0=true; 
     D0=0; 
     for (is=0; is<ns; is++)
          for (i=0; i<neq; i++) DMin[is*neq+i]=INF;
     while (iter0  && thestatptr->itcount<MAXIT)
     { 

           thestatptr->itcount++; 
           k++;
           Yi_update_Classic (u,z,system,method,thestatptr);

           /* ---------- g = hB F(L)- L -------------------*/
#ifdef PARALLEL
#      pragma omp parallel for num_threads(thread_count) private(isn)
#endif       
           for (is = 0; is<ns; is++)
           {
                 isn=neq*is;
                 params->eval=system->cod[0];
                 thestatptr->fcn++;
                 system->f(neq,tn+method->hc[is],&z[isn],&fz[isn],params); 
                 for (i=0; i<neq; i++)  g[isn+i]=fz[isn+i]*method->hb[is]-li[isn+i];
           }


           Compute_R (system, h, g, RR, method, thestatptr);
           Compute_Z (system, h, MM, RR, ZZ, method, thestatptr, ipiv);
           Compute_W1 (system, h, RR, ZZ, W1, method, thestatptr);
           Compute_W2 (system, h, W1, ZZ, g, W2, method,thestatptr);
           Compute_DL (system,W1,W2,DL,method,thestatptr);

           /* ---- Lb=L+Dl ------------*/

           for (is=0; is<ns; is++)
           {
                 isn=neq*is;
                 for (i=0; i<neq; i++)
                 {
                       lik[isn+i]=li[isn+i];
                       li[isn+i]=li[isn+i]+DL[isn+i];
                 }
           }

           StopCriterionFloat (system,method,&D0,&iter0,DMin,li,lik);
         
     }

      
     if (thestatptr->itcount==MAXIT)
     { 
           printf("Break: step(MAXIT)=%i\n",thestatptr->itcount);
           thestatptr->convergence=FAIL; 
     }
     else
        if (D0==(ns*neq))
        {
                thestatptr->itzeroDL[0]++;
        }

    thestatptr->totitDL[0]+=k;
  

/*------------ 2-Step: Compute Jacobian matrices ----------------------------*/

     Yi_update_Classic (u,z,system,method,thestatptr);

     for (is=0; is<ns; is++)
          system->jac(neq,tn+method->hc[is],
                      &z[is*neq],&thestatptr->jac_i[is*(neq*neq)],params);

/*------------ 3-Step: Inner iterations -------------------------------------*/

     l=0;
     iter0=true;
     D0=0;
     for (is=0; is<ns; is++)
          for (i=0; i<neq; i++) DMin[is*neq+i]=INF;

     while (iter0 && l<MAXIT)
     {

          l++;

          Compute_GG (system, h, g, DL, GG, method, thestatptr);
          Compute_R (system, h, GG, RR, method, thestatptr);
          Compute_Z (system, h, MM, RR, ZZ, method, thestatptr, ipiv);
          Compute_W1 (system, h, RR, ZZ, W1, method, thestatptr);
          Compute_W2 (system, h, W1, ZZ, GG, W2, method,thestatptr);
          Compute_DL (system,W1,W2,DLl,method,thestatptr);
          
          for (is=0; is<ns; is++)
              for (j=0; j<neq; j++)
              { DLOLD[neq*is+j]=DL[neq*is+j];
                DL[is*neq+j]+=DLl[is*neq+j];
              }

          StopCriterionFloat (system,method,&D0,&iter0,DMin,DL,DLOLD);  

     }

    if (D0==(ns*neq))  thestatptr->itzeroDL[1]++;
    thestatptr->totitDL[1]+=l;

     /* ---- L^[k]=L^[k-1]+DL^[k,l]------------*/

     for (is=0; is<ns; is++)
     {
           isn=neq*is;
           for (i=0; i<neq; i++)
           {
                thestatptr->li[isn+i]=thestatptr->lik[isn+i]+DL[isn+i];   
           }
     }


/*------------ 4-Step: Inexact Newton iterations ----------------------------*/


     thestatptr->itcount++;
     Yi_update_Classic (u,z,system,method,thestatptr);

     /* ---------- g = hB F(L)- L -------------------*/
     #ifdef PARALLEL
     #      pragma omp parallel for num_threads(thread_count) private(isn)
     #endif    
     bb=1.;   
     lda=neq; incx=1; incy=1;
     for (is = 0; is<ns; is++)
     {
           isn=neq*is;
           params->eval=system->cod[0];
           thestatptr->fcn++;
           system->f(neq,tn+method->hc[is],&z[isn],&fz[isn],params); 
           for (i=0; i<neq; i++)  g[isn+i]=fz[isn+i]*method->hb[is]-li[isn+i];
           aa=method->hb[is];            
           GEMV (CblasRowMajor,CblasNoTrans,neq,neq,aa,	
                 &thestatptr->jac_i[is*(neq*neq)], lda, &u->ee[0], incx, bb, &g[isn], incy); 
     }
              

     /* ------DL^{[k,0]}=(I-hBAB x J)^{-1} g^{[k]} ---*/

     Compute_R (system, h, g, RR, method, thestatptr);
     Compute_Z (system, h, MM, RR, ZZ, method, thestatptr, ipiv);
     Compute_W1 (system, h, RR, ZZ, W1, method, thestatptr);
     Compute_W2 (system, h, W1, ZZ, g, W2, method,thestatptr);
     Compute_DL (system,W1,W2,DL,method,thestatptr);

     /*------ Inner iterations -----------------------*/

     iter0=true;
     D0=0;
     for (is=0; is<ns; is++)
          for (i=0; i<neq; i++) DMin[is*neq+i]=INF;
     l=0;
     while (iter0 && l<MAXIT)
     {
           l++;

           Compute_GG (system, h, g, DL, GG, method, thestatptr);
           Compute_R (system, h, GG, RR, method, thestatptr);
           Compute_Z (system, h, MM, RR, ZZ, method, thestatptr, ipiv);
           Compute_W1 (system, h, RR, ZZ, W1, method, thestatptr);
           Compute_W2 (system, h, W1, ZZ, GG, W2, method,thestatptr);
           Compute_DL (system,W1,W2,DLl,method,thestatptr);
           
           for (is=0; is<ns; is++)
                 for (j=0; j<neq; j++)
                 {     DLOLD[neq*is+j]=DL[neq*is+j];
                       DL[is*neq+j]+=DLl[is*neq+j];
                 }

           StopCriterionFloat (system,method,&D0,&iter0,DMin,DL,DLOLD);  

     }
   
     if (D0==(ns*neq))  thestatptr->itzeroDL[2]++;
     thestatptr->totitDL[2]+=l;


     /* ---- Lb=L+Dl ------------*/

     for (is=0; is<ns; is++)
     {
           isn=neq*is;
           for (i=0; i<neq; i++)
           {
                 lik[isn+i]=li[isn+i];
                 li[isn+i]=li[isn+i]+DL[isn+i];
           }
      }


/*------------ End: Check convergence ---------------------------------------*/

     for (is=0; is<ns;is++) for (i=0; i<neq; i++) zold[neq*is+i]=thestatptr->z[neq*is+i];
     Yi_update_Classic (u,z,system,method,thestatptr);
     difftest=NormalizedDistance(neq,ns,options,z,zold);
     if (difftest>1.)
     {
         thestatptr->convergence=FAIL;
         printf("Lack of convegence of Newton iteration:\
                 step=%i,iteration=%i,",
                 thestatptr->stepcount,thestatptr->itcount);  
                 printf("difftest=%lg\n",difftest);                 
     }

    return;

}



/******************************************************************************/
/* 								              */
/*   MMNewtonSS: Newton Super Simplified				      */
/* 									      */
/* 									      */
/******************************************************************************/


void MMNewtonSS (int neq, gauss_method *method, val_type tn, val_type h, 
                 ode_sys *system, val_type *u, val_type *MM)
{ 

/*----------------  First: initialization ------------------------------------*/

     int ns;
     ns=method->ns;


/*-------------declarations --------------------------------------------------*/

     int i,j,k,l;
     int n1,n2,m1,m2;
     int ci,cj;

     val_type Aij,Bkl;
     val_type Jac[neq*neq];
     parameters *params;

/* ----------- implementation  -----------------------------------------------*/ 

     params=&system->params;

     n1=ns;
     m1=ns;
     n2=neq;
     m2=neq;

     system->jac(neq,tn+h/2,u,Jac,params);

     for (i=0; i<n1; i++)
     {
           for (k=0; k<n2; k++)           
           {
                 ci=i*n2+k;

                 for (j=0; j<m1; j++)
                 {
                       Aij=method->hBAB[i*m1+j];
                       
                       for (l=0; l<m2; l++)
                       {
                             cj=m2*j+l;
                             Bkl=Jac[k*m2+l];
                             
                             if (ci==cj)
                                 MM[ci*(m1*m2)+cj]=1.-Aij*Bkl;
                             else
                                 MM[ci*(m1*m2)+cj]=-Aij*Bkl;

                       }
                 }
            }
     }

     
     return;

}




/******************************************************************************/
/* 									      */
/*   MatAdd     							      */
/* 									      */
/* 									      */
/******************************************************************************/

void MatAdd (val_type aa,val_type *A,int lda, val_type bb,val_type *B, int ldb, 
             val_type cc, val_type *C, int ldc,int m, int n)
{
/*----------------------------------------------------------------------------*/
/* C=aa*A+bb*B+gamma*C                                                        */
/* where A,B,C is an m-by-n matrix                                            */
/*----------------------------------------------------------------------------*/
/*-------------declarations --------------------------------------------------*/

     int i,j;
     
/* ----------- implementation  -----------------------------------------------*/ 

     for (i=0; i<m; i++)
           for (j=0; j<n; j++) 
               if (cc==0.)
                   C[i*ldc+j]=aa*A[i*lda+j]+bb*B[i*ldb+j];
               else
                   C[i*ldc+j]=aa*A[i*lda+j]+bb*B[i*ldb+j]+cc*C[i*ldc+j];

     return;

}

/******************************************************************************/
/* 									      */
/*   MatAddID     							      */
/* 									      */
/* 									      */
/******************************************************************************/

void MatAddID (val_type aa,val_type bb,val_type *B, int ldb, 
               val_type cc, val_type *C, int ldc, int n)
{
/*----------------------------------------------------------------------------*/
/* C=aa*Id+bb*B+gamma*C                                                       */
/* where Id,B,C is an n-by-n matrix                                           */
/* and Id is identity matrix                                                  */
/*----------------------------------------------------------------------------*/
/*-------------declarations --------------------------------------------------*/

     int i,j;
     
/* ----------- implementation  -----------------------------------------------*/ 

     for (i=0; i<n; i++)
     {
           /* diagonal */
           j=i;
           if (cc==0.)
                   C[i*ldc+j]=aa+bb*B[i*ldb+j];
               else
                   C[i*ldc+j]=aa+bb*B[i*ldb+j]+cc*C[i*ldc+j];

           /* no diagonal */
           
           for (j=i+1; j<n; j++) 
           {
               if (cc==0.)
               {
                   C[i*ldc+j]=bb*B[i*ldb+j];
                   C[j*ldc+i]=bb*B[j*ldb+i];
               }
               else
               {
                   C[i*ldc+j]=bb*B[i*ldb+j]+cc*C[i*ldc+j];
                   C[j*ldc+i]=bb*B[j*ldb+i]+cc*C[j*ldc+i];
               }
           }
      }
      return;

}




/******************************************************************************/
/*									      */
/*      NSS_Solve							      */
/*									      */
/******************************************************************************/


int NSS_Solve     (ode_sys *system, solution *u, val_type tn,val_type h, 
                   gauss_method *method,solver_stat *thestatptr,
                   val_type *MM, int *ipiv)
{


# ifdef PARALLEL
     int extern thread_count;
# endif

/* ---------- First initializations ------------------------------------------*/

     int neq,ns;
     parameters *params;
     val_type *z,*li,*lik;
     val_type *DL;
    
     neq=system->neq;
     ns=method->ns;
     params=&system->params;

     val_type fz[ns*neq];
     z=thestatptr->z;
     li=thestatptr->li;
     DL=thestatptr->DL;
     lik=thestatptr->lik;

/*------ declarations --------------------------------------------------------*/

     int i,is,isn;
     int info;

     int n,lda,ldb,nrhs;
     char trans;
  
/* ----------- implementation  -----------------------------------------------*/


#ifdef PARALLEL
#      pragma omp parallel for num_threads(thread_count) private(isn)
#endif  
     for (is = 0; is<ns; is++)
     {
           isn=neq*is;
           params->eval=system->cod[0];
           thestatptr->fcn++;
           system->f(neq,tn+method->hc[is],&z[isn],&fz[isn],params); 
           for (i=0; i<neq; i++) DL[isn+i]=fz[isn+i]*method->hb[is]-li[isn+i];	
     }


     trans = 'N';
     n = ns*neq;         
     nrhs=1;             
     lda = (ns)*neq;
     ldb = 1;

     info= GETRS (LAPACK_ROW_MAJOR,trans,n,nrhs,MM,lda,ipiv,DL,ldb); 

     if( info < 0 )
     {
	printf( "The execution _getrs is unsuccesfull. info=%i,\n",info );
	printf( "the solution could not be computed.\n" );
	exit( 1 );
     }

     /* ---- Lb=L+Dl ------------*/

     for (is=0; is<ns; is++)
     {
           isn=neq*is;
           for (i=0; i<neq; i++)
           {
                 lik[isn+i]=li[isn+i];
                 li[isn+i]=li[isn+i]+DL[isn+i];
           }
     }

     return(0);

}	




/******************************************************************************/
/*									      */
/*      NSS_Solve_plus   			   			      */
/*									      */
/******************************************************************************/


int NSS_Solve_plus   (ode_sys *system, solution *u, val_type tn,val_type h, 
                      gauss_method *method,solver_stat *thestatptr,
                      val_type *MM, int *ipiv)
{

# ifdef PARALLEL
     int extern thread_count;
# endif

/* ---------- First initializations ------------------------------------------*/

     int neq,ns,mm,smm;
     parameters *params;
     val_type *z,*li,*lik;
     val_type *DL;
      
     neq=system->neq;
     ns=method->ns;
     mm=method->mm;
     smm=method->smm;

     params=&system->params;

     val_type fz[ns*neq];
     val_type g[ns*neq];
     val_type RR[mm*neq];
     val_type ZZ[neq],W1[mm*neq],W2[smm*neq];
     z=thestatptr->z;
     li=thestatptr->li;
     DL=thestatptr->DL;
     lik=thestatptr->lik;


/*------ declarations --------------------------------------------------------*/

     int i,is,isn;
  
/* ----------- implementation  -----------------------------------------------*/


/* ---------- g = hB F(L)- L -------------------*/
#ifdef PARALLEL
#      pragma omp parallel for num_threads(thread_count) private(isn)
#endif       
     for (is = 0; is<ns; is++)
     {
           isn=neq*is;
           params->eval=system->cod[0];
           thestatptr->fcn++;
           system->f(neq,tn+method->hc[is],&z[isn],&fz[isn],params); 
           for (i=0; i<neq; i++)  g[isn+i]=fz[isn+i]*method->hb[is]-li[isn+i];
     }


/* ------- R , Z, W1, W2  ---------------------------------------------------- */

     Compute_R (system, h, g, RR, method, thestatptr);
     Compute_Z (system, h, MM, RR, ZZ, method, thestatptr, ipiv);
     Compute_W1 (system, h, RR, ZZ, W1, method, thestatptr);
     Compute_W2 (system, h, W1, ZZ, g, W2, method,thestatptr);

/* ---------- Compute DL --------------------------*/

     Compute_DL (system,W1,W2,DL,method,thestatptr);

     /* ---- Lb=L+Dl ------------*/

     for (is=0; is<ns; is++)
     {
           isn=neq*is;
           for (i=0; i<neq; i++)
           {
                 lik[isn+i]=li[isn+i];
                 li[isn+i]=li[isn+i]+DL[isn+i];
           }
      }

     return(0);

}	




/******************************************************************************/
/*									      */
/*      Compute_MM 							      */
/*									      */
/******************************************************************************/

void Compute_MM   (ode_sys *system,val_type h,val_type *MM, 
                  gauss_method *method,solver_stat *thestatptr)             

{

/*----------------  First: initialization ------------------------------------*/
     int neq;

     neq=system->neq;
     int mm=method->mm;

/*------ declarations --------------------------------------------------------*/

     int i,j,ii,jj;
     int ind0;
     int info;

     val_type C1[neq*neq],SUM[neq*neq];     

     val_type aa,bb,cc;
     int lda,ldb,ldc;
     int ipiv[neq];
     
/* ----------- implementation  -----------------------------------------------*/ 


     /* ---- Init: MM--------------------------------*/
 
     for (i=0; i<neq; i++)
           for (j=0; j<neq; j++)
           {
                if (i==j) MM[i*neq+j]=1.;
                else MM[i*neq+j]=0.;

                SUM[i*neq+j]=0.;
           }

     /* --- Id + h^2 sigma_i^2 J^2 ------------ */

     aa = h*h; bb = 0.;
     lda=neq; ldb=neq; ldc=neq;
     GEMM (CblasRowMajor,CblasNoTrans,CblasNoTrans,neq,neq,neq,aa,	
           thestatptr->jac, lda, thestatptr->jac, ldb, bb, C1, ldc);

     aa=1.;
     cc=0.;
     for (i=0; i<mm; i++)
     {

           ind0=i*(neq*neq);
           
           bb=method->sigma2[i];           
           MatAddID (aa,bb,C1,neq,cc,&thestatptr->INVIDJ2[ind0],neq,neq);

           info = GETRF(LAPACK_ROW_MAJOR, neq, neq,
                        &thestatptr->INVIDJ2[ind0], neq, ipiv);

           if( info > 0 )
           {
 	         printf( "The diagonal element of the triangular factor of A,\n" );
	         printf( "U(%i,%i) is zero, so that A is singular;\n", info, info );
	         printf( "the solution could not be computed.\n" );
	         exit( 1 );
           }

           info= GETRI (LAPACK_ROW_MAJOR,neq,
                        &thestatptr->INVIDJ2[ind0],neq,ipiv); 

           if( info < 0 )
           {
	         printf( "The execution (Compute IDM) _getri is unsuccesfull. info=%i,\n",info );
	         printf( "the solution could not be computed.\n" );
	         exit( 1 );
           }


           for (ii=0; ii<neq; ii++)
                 for (jj=0; jj<neq; jj++)
                      SUM[ii*neq+jj]+=method->alpha2[i]*thestatptr->INVIDJ2[ind0+ii*neq+jj];
 
      }    


     aa=-h/2; 
     bb=1.;
     lda=neq; ldb=neq; ldc=neq;
     GEMM (CblasRowMajor,CblasNoTrans,CblasNoTrans,neq,neq,neq,aa,	
           thestatptr->jac, lda, SUM, ldb, bb, MM, ldc);   


/*
     printf("\n(M). i=%i,\n", i);
     for (i=0; i<neq; i++)
     {      for (j=0; j<neq; j++) printf("%lg,",MM[i*neq+j]);
            printf("\n");
     }          

*/

     return;


}




/******************************************************************************/
/*									      */
/*      Compute_R  28-02-2017						      */
/*									      */
/******************************************************************************/

void Compute_R (ode_sys *system, val_type h,val_type *g, val_type *RR, 
                gauss_method *method,solver_stat *thestatptr)

{

/*----------------  First: initialization ------------------------------------*/
     int neq,ns;

     neq=system->neq;
     ns=method->ns;
     int mm=method->mm;

/*------ declarations --------------------------------------------------------*/

  
     val_type aa,bb;
     int lda,ldb,ldc;
     val_type C1[mm*neq];
     
/* ----------- implementation  -----------------------------------------------*/ 


/*   RR=(Q1T g+ (h DQ2T) g J^T) */
     aa=1.;
     bb=0.;
     lda=mm; ldb=neq; ldc=neq;
     GEMM (CblasRowMajor,CblasTrans,CblasNoTrans,mm,neq,ns,aa,	
           method->Q1, lda, g, ldb, bb, RR, ldc);

     aa=1.;
     bb=0.;
     lda=ns; ldb=neq; ldc=neq;
     GEMM (CblasRowMajor,CblasNoTrans,CblasNoTrans,mm,neq,ns,aa,	
           method->hDQ2T, lda, g, ldb, bb, C1, ldc);

     aa=1.;
     bb=1.;
     lda=neq; ldb=neq; ldc=neq;
     GEMM (CblasRowMajor,CblasNoTrans,CblasTrans,mm,neq,neq,aa,	
           C1, lda,thestatptr->jac, ldb, bb, RR, ldc);
     
/*

     int i,j;

     printf("\nRR\n");
     for (i=0; i<mm; i++) 
     {
          for (j=0; j<neq; j++) printf("%lg,",RR[i*neq+j]);
          printf("\n");
     }
     printf("\n");
*/

     return;


}





/******************************************************************************/
/*									      */
/*      Compute_Z   		         				      */
/*									      */
/******************************************************************************/

void Compute_Z (ode_sys *system, val_type h, val_type *MM,val_type *RR, 
               val_type *ZZ,gauss_method *method,solver_stat *thestatptr,
               int *ipiv)

{

/*----------------  First: initialization ------------------------------------*/
     int neq;

     neq=system->neq;
     int mm=method->mm;

/*------ declarations --------------------------------------------------------*/

     int i;
     int ind0;
     int info;
     val_type SUM[neq];  
  
     val_type aa,bb;
     int n,lda,ldb,nrhs;
     int incx,incy;
     char trans;
     
/* ----------- implementation  -----------------------------------------------*/ 
   

     /* ---- compute d= h J sum alpha_i ( Id+h^2*sigma^2*J^2 )^-1 --*/

     i=0;
     aa=method->alpha[i];
     bb=0.;
     lda=neq; incx=1; incy=1;

     GEMV (CblasRowMajor,CblasNoTrans,neq,neq,aa,	
           &thestatptr->INVIDJ2[0], lda, &RR[i*neq], incx, bb, SUM, incy);

     bb=1.;
     for (i=1; i<mm; i++)
     {

           ind0=i*(neq*neq);
           aa=method->alpha[i];
           
           GEMV (CblasRowMajor,CblasNoTrans,neq,neq,aa,	
           &thestatptr->INVIDJ2[ind0], lda, &RR[i*neq], incx, bb, SUM, incy);

     }

     aa=h ; 
     bb=0.;

     GEMV (CblasRowMajor,CblasNoTrans,neq,neq,aa,	
           thestatptr->jac, lda, SUM, incx, bb, ZZ, incy); 

     /* -----  Solve M Dz= d --------*/

     trans = 'N';
     n = neq;        
     nrhs=1;         
     lda = neq;
     ldb = 1;

     info= GETRS (LAPACK_ROW_MAJOR,trans,n,nrhs,MM,lda,ipiv,ZZ,ldb);

     if( info < 0 )
     {
	printf( "The execution (Compute z) _getrs is unsuccesfull. info=%i,\n",info );
	printf( "the solution could not be computed.\n" );
	exit( 1 );
     }


/*
     printf("\nz\n");
     for (i=0; i<neq; i++) printf("%lg,",ZZ[i]);
     printf("\n");
*/

     return;

}




/******************************************************************************/
/*									      */
/*      Compute_W1    							      */
/*									      */
/******************************************************************************/

void Compute_W1 (ode_sys *system, val_type h, val_type *RR,val_type *ZZ,
                val_type *W1, gauss_method *method,solver_stat *thestatptr)

{

/*----------------  First: initialization ------------------------------------*/
     int neq;

     neq=system->neq;
     int mm=method->mm;

/*------ declarations --------------------------------------------------------*/

     int i,j,ind;
     val_type C1[neq]; 
     val_type aux;
  
     val_type aa,bb;
     int lda;
     int incx,incy;
     
/* ----------- implementation  -----------------------------------------------*/ 


     lda=neq; incx=1; incy=1;
     aa=1.;
     bb=0.;   

     for (i=0; i<mm; i++)
     {
           ind=i*(neq*neq);
    
           aux=method->alpha[i]/2;
           for (j=0; j<neq; j++) C1[j]=(RR[i*neq+j]+aux*ZZ[j]);

           GEMV (CblasRowMajor,CblasNoTrans,neq,neq,aa,	
                 &thestatptr->INVIDJ2[ind], lda, C1, incx, bb, &W1[i*neq], incy);

     }

/*
     printf("\nW1\n");
     for (i=0; i<mm; i++) 
     {
          for (j=0; j<neq; j++) printf("%lg,",W1[i*neq+j]);
          printf("\n");
     }
     printf("\n");
*/


     return;
}




/******************************************************************************/
/*									      */
/*      Compute_W2    							      */
/*									      */
/******************************************************************************/

void Compute_W2 (ode_sys *system, val_type h,val_type *W1,val_type *ZZ,
                 val_type *g, val_type *W2, gauss_method *method,
                 solver_stat *thestatptr)

{

/*----------------  First: initialization ------------------------------------*/
     int neq,ns;

     neq=system->neq;
     ns=method->ns;
     int mm=method->mm;
     int smm=method->smm;

/*------ declarations --------------------------------------------------------*/

     val_type aa,bb;
     int lda,ldb,ldc;
     val_type C1[smm*neq];
     
/* ----------- implementation  -----------------------------------------------*/ 
    

     aa=1.;
     bb=0.;
     lda=smm; ldb=neq; ldc=neq;
     GEMM (CblasRowMajor,CblasTrans,CblasNoTrans,smm,neq,ns,aa,	
           method->Q2, lda, g, ldb, bb, W2, ldc);                

     aa=-1.;
     bb=0.;
     lda=mm; ldb=neq; ldc=neq;
     GEMM (CblasRowMajor,CblasNoTrans,CblasNoTrans,smm,neq,mm,aa,	
           method->hDT, lda, W1, ldb, bb, C1, ldc);

     aa=1.;
     bb=1.;
     lda=neq; ldb=neq; ldc=neq;
     GEMM (CblasRowMajor,CblasNoTrans,CblasTrans,smm,neq,neq,aa,	
           C1, lda,thestatptr->jac, ldb, bb, W2, ldc); 


/*

     int i,j;

     printf("\nW2\n");
     for (i=0; i<smm; i++) 
     {
          for (j=0; j<neq; j++) printf("%lg,",W2[i*neq+j]);
          printf("\n");
     }
     printf("\n");
*/

     return;

}




/******************************************************************************/
/*									      */
/*      Compute_DL 		         				      */
/*									      */
/******************************************************************************/

void Compute_DL (ode_sys *system, val_type *W1,val_type *W2,
                 val_type *DL, gauss_method *method,
                 solver_stat *thestatptr)

{

/*----------------  First: initialization ------------------------------------*/
     int neq,ns;

     neq=system->neq;
     ns=method->ns;

     int mm=method->mm;
     int smm=method->smm;


/*------ declarations --------------------------------------------------------*/

     val_type aa,bb;
     int lda,ldb,ldc;
     
/* ----------- implementation  -----------------------------------------------*/ 

     aa=1.;
     bb=0.;
     lda=mm; ldb=neq; ldc=neq;
     GEMM (CblasRowMajor,CblasNoTrans,CblasNoTrans,ns,neq,mm,aa,	
           method->BQ1, lda, W1, ldb, bb, DL, ldc); 

     aa=1.;
     bb=1.;
     lda=smm; ldb=neq; ldc=neq;
     GEMM (CblasRowMajor,CblasNoTrans,CblasNoTrans,ns,neq,smm,aa,	
           method->BQ2, lda, W2, ldb, bb, DL, ldc); 


/*

     int i,j;

     printf("\nDL\n");
     for (i=0; i<ns; i++) 
     {
          for (j=0; j<neq; j++) printf("%lg,",DL[i*neq+j]);
          printf("\n");
     }
     printf("\n");
*/


     return;


}


/******************************************************************************/
/*									      */
/*      Compute_GG 							      */
/*									      */
/******************************************************************************/

void Compute_GG (ode_sys *system, val_type h, 
                 val_type *g, val_type *DL, val_type *GG,
                 gauss_method *method, solver_stat *thestatptr)
{

/*----------------  First: initialization ------------------------------------*/
     int neq,ns;

     neq=system->neq;
     ns=method->ns;

/*------ declarations --------------------------------------------------------*/

     int i,is,in;
     int ind0,indi;

     val_type alpha,beta;
     int lda,ldb,ldc;
     int incx,incy;

     val_type SUM[ns*neq];   
     val_type aux1,aux2;

     
/* ----------- implementation  -----------------------------------------------*/ 


    alpha = 1.; beta = 0.;
    lda=ns; ldb=neq; ldc=neq;  
    GEMM(CblasRowMajor,CblasNoTrans,CblasNoTrans, ns,neq,ns,
          alpha, method->hBAB, lda, DL, ldb, beta, SUM, ldc );

    alpha = 1.; beta = 1.;
    lda=neq;
    incx=1; incy=1;

    for (is = 0; is<ns; is++)
    {
           ind0=is*(neq*neq);
           indi=is*neq;

           for (i = 0; i<neq; i++)
           {   
                 in=indi+i;
                 aux1=g[in];
                 aux2=DL[in];
                 GG[in]=(aux1-aux2);
           }

           GEMV (CblasRowMajor,CblasNoTrans,neq,neq,alpha,	
                 &thestatptr->jac_i[ind0], lda, &SUM[indi], incx, beta, &GG[indi], incy); 

     }  


/*
     printf("\nGG\n");
     for (i=0; i<ns; i++) 
     {
          for (j=0; j<neq; j++) printf("%lg,",GG[i*neq+j]);
          printf("\n");
     }
     printf("\n");
*/

     return;


}



/******************************************************************************/
/*									      */
/*      MyOutput (RKG)							      */
/*									      */
/******************************************************************************/
 
void MyOutput(ode_sys *system,gauss_method *method,val_type t,solution *u,
               solver_stat *thestatptr,
               parameters *params,toptions *options,FILE *myfile)
{

/* ---------- First initializations ------------------------------------------*/

     int neq;
     neq=system->neq;


/*------ declarations --------------------------------------------------------*/
     val_type DH;
     int i,node;
  
     struct rec
     {
     val_type t;
     val_type uu[neq];
     val_type ee[neq];
     };
  
     struct rec my_record;

/* ----------- implementation  -----------------------------------------------*/ 

     node=neq;
     DH=(system->ham(node,u,params)-thestatptr->E0)/thestatptr->E0;
     if (fabs(DH)>thestatptr->MaxDE) thestatptr->MaxDE=FABS(DH);

     if (((thestatptr->stepcount % options->sampling) == 0) || (thestatptr->laststep))
     {
           thestatptr->nout++;
           printf("%i,%i, %lg\n",thestatptr->stepcount,thestatptr->itcount,DH);

           my_record.t=t;
           for (i=0; i<neq;i++) 
           {
                 my_record.uu[i]=u->uu[i];
                 my_record.ee[i]=u->ee[i];
           }
           fwrite(&my_record, sizeof(struct rec), 1, myfile);
  
     } 

    if (thestatptr->stepcount>1) statlinit(system,method,thestatptr);

  return;

}



/******************************************************************************/
/* 									      */
/*   Compensated Summation:						      */
/* 									      */
/* 									      */
/******************************************************************************/

void CompensatedSummation (gauss_method *gsmethod,
                           solution *u,ode_sys *system,
                           toptions *options,solver_stat *thestatptr)

{
/* ---------- First initializations ------------------------------------------*/

     int neq,ns;
     
     neq=system->neq;
     ns=gsmethod->ns;

/*------ declarations --------------------------------------------------------*/

     int i,ix,is,isn;
     val_type aux;
     val_type s0,s1,delta,eli,ee;
     val_type *lik;
     lik=thestatptr->lik;

/* ----------- implementation  -----------------------------------------------*/ 


          for (i = 0; i<neq; i++)
          {
              eli= thestatptr->DL[i]; 

              for (is=1; is<ns; is++)
              {
                      isn=neq*is+i;
                      eli+=thestatptr->DL[isn];                      
              }

              u->ee[i]+=eli;
       
          }

     /* ----------------yn+1=yn+(Sum Li^{[k-1]}+e_n) 21-07-2016---------------------------*/

          for (i = 0; i<neq; i++)
          { 

               s0=u->uu[i];
               ee=u->ee[i];
		
               for (ix =0 ; ix<ns; ix++)
               {
                    is=gsmethod->orderedindices[ix];
                    isn=neq*is+i;                
                    delta=lik[isn]+ee;
                    if (options->rdigits>0) 
                        RemoveDigitsFcn(&delta,options->mrdigits);
                    s1=s0;
                    s0=s1+delta;
                    aux=s1-s0;
                    ee=aux+delta;
               } 

               u->uu[i]=s0;
               u->ee[i]=ee;
      	       	       
          }
     


     return;
}



/******************************************************************************/
/* 								              */
/*   RGK:  								      */
/* 									      */
/* 									      */
/*******************************************************************************/

void IRKNEWTON 
(val_type t0, val_type t1, val_type h,
 gauss_method *gsmethod, solution *u,
 ode_sys *system,toptions *options,
 solver_stat *thestatptr)

{

/* ---------- First initializations ------------------------------------------*/

     int neq;
     parameters *params;

     neq=system->neq;
     params=&system->params;

/*------ declarations --------------------------------------------------------*/

     FILE *myfile;
  
     int istep,nstep;
     val_type tn;             
     val_type *z;

     z=thestatptr->z;


/* ----------- implementation  -----------------------------------------------*/ 

     if (options->TheOutput != 0)
         myfile = fopen(options->filename,"wb"); 

/* ----------- initial energy (E0)   -----------------------------------------*/

     thestatptr->E0=system->ham(neq,u,params);
     printf("Initial energy=%lg\n", thestatptr->E0);

     tn=t0;
     nstep=(t1-t0)/h;

     if (options->TheOutput != 0)
          options->TheOutput(system,gsmethod,tn,u,thestatptr,params,options,myfile);

     for(istep=0; istep<nstep; istep++) 
     {     
          
          options->StageInitFn (u,z,system,gsmethod,thestatptr,options);              

          options->IRKNEWTON_Step_Fn (system,u,tn,h,options,gsmethod,thestatptr);

          if (thestatptr->convergence==FAIL)
          { 
               printf("Stop Fail. step=%i\n", istep);
               nstep=istep;
          } 	      

          else
          {
               CompensatedSummation (gsmethod,u,system,options,thestatptr);

               tn=(istep+1)*h;         

               thestatptr->stepcount++;		
               (thestatptr->totitcount)+=(thestatptr->itcount); 
               if ((thestatptr->itcount)>(thestatptr->maxitcount))
                   (thestatptr->maxitcount)=(thestatptr->itcount);

               if (options->TheOutput != 0)
                   options->TheOutput(system,gsmethod,tn,u,thestatptr,params,options,myfile);


               if ((tn+h)>=t1)
               { 
                     h=t1-tn;
                     thestatptr->laststep=true;
               } 

          }

     }


     fclose(myfile);

     return;

}


