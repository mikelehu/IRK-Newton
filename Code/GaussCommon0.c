/*------------------------------------------------------------------------------*/
/*										*/
/*    GaussCommon0.c								*/
/*										*/
/*	Functions: 								*/
/*	 Newton_it_Mix():							*/
/*	 Newton_it_Mix2():							*/
/*	 MyKroneckerProd0():							*/
/*	 MMfun0():								*/
/*	 MMfunOsoa0():								*/
/*	 Newton_Step0():							*/
/* -----------------------------------------------------------------------------*/

#include <GaussCommon0.h>
#include <quadmath.h>

/************************************************************************************/
/* 					   					    */
/* NormalizedDistance: 	         	  					    */
/*                                        					    */
/*										    */
/************************************************************************************/

val_type NormalizedDistance0 (int neq,int ns,
                             toptions *options,val_type0 *z,val_type0 *zold)
{


/*---------------- Declarations -----------------------------------------*/

     int i,is,ix;
     val_type0 maxi,mayi,relerrors;
     val_type0 maxz,maxzold;

/* --------------- Implementation --------------------------------------*/
   	
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
               relerrors=FMAX(relerrors,FABS(z[ix]-zold[ix])/(mayi*options->rtol[i]+options->atol[i]));
          }

          maxi=FMAX(maxi,relerrors);

     }

     return maxi;

}



/************************************************************************************/
/* 					   					    */
/* Yi_update: 		         	  					    */
/*                                        					    */
/*										    */
/************************************************************************************/

int Yi_update0 (solution *u, val_type0 *z0,ode_sys *system,
                gauss_method *method,solver_stat *thestatptr)
{


/* ---------- First initializations -------------------------------------*/

     int neq,ns;

     neq=system->neq;
     ns=method->ns;

/*------------- declarations -----------------------------------------*/ 

     int i,is,in,js,jsn;
     val_type0 *coef,*li0;
     val_type0 sum;

/* ----------- implementation  --------------------------------------*/


    li0=thestatptr->li0;
    coef=method->m0;

    for (is = 0; is<ns; is++)
    {
                for (i = 0; i<neq; i++)
                {  
                     in=neq*is+i;
                     sum=0.;

                     for (js = 0; js<ns; js++)
                     {
                          jsn=ns*is+js;
                          sum+=li0[neq*js+i]*(coef[jsn]);
                     }

                     z0[in]=u->uu[i]+sum;
                }
     }  
      
     return(0);


}


/************************************************************************************/
/* 					   					    */
/* RemoveDigitsFcn	         	  					    */
/*                                        					    */
/*										    */
/************************************************************************************/

void RemoveDigitsFcn0 (ode_sys *system,gauss_method *method,val_type0 *z0, int m)
{

/* ---------- First initializations ------------------------------------*/

     int neq,ns;
    
     neq=system->neq;
     ns=method->ns;

/*---------------- Declarations -----------------------------------------*/

     int i,is,ix;
     val_type0 mY,mYY;

/* --------------- Implementation --------------------------------------*/

     for (i=0; i<neq; i++)
          for (is=0; is<ns; is++)
          { 
               ix=neq*is+i;
               mY=m*z0[ix];
               mYY=mY+z0[ix];
               z0[ix]=mYY-mY;
           }

}


/****************************************************************************************/
/* 											*/
/*   UpdateDmin  								        */
/* 											*/
/* 											*/
/****************************************************************************************/
void UpdateDMin0 (ode_sys *system,gauss_method *method,
                  int *D0,bool *cont,val_type0 *DMin,val_type0 *Y,val_type0 *Yold)
{
/*----------------  First: initialization -------------------------------*/
     int neq,ns;

     neq=system->neq;
     ns=method->ns;

/*------ declarations --------------------------------------------------*/
     int i,is;
     val_type0 dY;

/* ----------- implementation  --------------------------------------*/ 

    *D0=0;
     *cont=false;

     for (is=0; is<ns; is++)
          for (i=0; i<neq; i++) 
          { 
               dY=FABS(Y[is*neq+i]-Yold[is*neq+i]);
               if (dY>0.)
               { 
                    
                   if (dY<DMin[is*neq+i])
                    {
                         DMin[is*neq+i]=dY;
                         *cont=true;
                    }
                   else
                    {
                         DMin[is*neq+i]=-1;
                    } 
               }
               else
               {
                    DMin[is*neq+i]=0;
                    *D0=*D0+1;
               } 
          }
  
     return;

}


/****************************************************************************************/
/* 											*/
/*   Newton_it_Mix  								        */
/* 											*/
/*        Azkeneko iterazioa doitasun altuan						*/
/* 											*/
/****************************************************************************************/


void Newton_it_Mix  ( ode_sys *system, solution *u, val_type tn,val_type h, 
                      toptions *options,gauss_method *method,solver_stat *thestatptr)

{ 

/*----------------  First: initialization -------------------------------*/
     int neq,ns;

     neq=system->neq;
     ns=method->ns;

/*------ declarations --------------------------------------------------*/

     int i,is,j;
     bool iter0;
//     bool D0;			
     int D0;
     val_type0 difftest,DMin[neq*ns];
     val_type0 *z0,*li0;
     parameters *params;
     val_type0 zold0[neq*ns];
     lowfloat0 *MM0;

     int n,m,lda,info;
     int ipiv[neq*(ns+1)];
     lowfloat0 Jac0[neq*neq];

     val_type *z,*li;
     val_type mynorm;
     
/* ----------- implementation  --------------------------------------*/ 


     MM0=(lowfloat0 *)malloc((neq*(ns+1))*(neq*(ns+1))*sizeof(lowfloat0));

     z0=thestatptr->z0;
     li0=thestatptr->li0;
     params=&system->params;
     z=thestatptr->z;
     li=thestatptr->li;

     for (is=0; is<ns; is++)
          for (i=0; i<neq; i++)
          {
                 DMin[is*neq+i]=INF;
                 z0[is*neq+i]=z[is*neq+i];
                 li0[is*neq+i]=li[is*neq+i];
          }  
 

     iter0=true;
     thestatptr->itcount=0;

/* ------LU factorization of a general n-by-m matrix A--------------*/

     n=neq*(ns+1);
     m=n;
     lda=n;


     switch (options->algorithm)
     { 
     case 11: 
           system->jac0(neq,tn,u->uu,Jac0,params);
           MMfun0(neq,method,Jac0,MM0);
     break;

     case 12:
	   MMfunOsoa0(neq,method,tn,system,thestatptr,MM0);
     break;

     default:
          printf("error. algorithm\n");
     break;
    }


    info = GETRF0(LAPACK_ROW_MAJOR, n, m, MM0, lda, ipiv);

     if( info > 0 )
     {
	printf( "The diagonal element of the triangular factor of A,\n" );
	printf( "U(%i,%i) is zero, so that A is singular;\n", info, info );
	printf( "the solution could not be computed.\n" );
	exit( 1 );
     }

     while (iter0 && thestatptr->itcount<MAXIT)
     { 
          
           for (is=0; is<ns;is++)
                for (i=0; i<neq; i++) zold0[neq*is+i]=thestatptr->z0[neq*is+i];
              
           options->iteration[0](system,u,tn,h,method,thestatptr,MM0,ipiv); 
           if (options->rdigits>0) RemoveDigitsFcn0(system,method,z0,options->mrdigits);
           UpdateDMin0(system,method,&D0,&iter0,DMin,z0,zold0);
           thestatptr->itcount++;
     }


     if (thestatptr->itcount==MAXIT)
     { 
           printf("Break: step(MAXIT)=%i\n",thestatptr->itcount);
           thestatptr->convergence=FAIL; 
     }
 
     if (D0<(ns*neq))
     {
           difftest=NormalizedDistance0(neq,ns,options,z0,zold0);
           if (difftest>1.)
           {
                 printf("Lack of convegence of fixed-point iteration:\
                         step=%i,difftest=%lg\n",thestatptr->itcount,difftest);
                 thestatptr->convergence=FAIL;
           }
     }
     else
     {
           thestatptr->itzero++;
     }


/* Last iteration in HighPrecision */


     for (is=0; is<ns; is++)
          for (i=0; i<neq; i++) {li[is*neq+i]=li0[is*neq+i];}

     Yi_update (u,z,system,method,thestatptr);

     if (thestatptr->stepcount==0)			// Only at the first step recalculate MM
     {
           switch (options->algorithm)
           { 
            case 11: 
                 system->jac0(neq,tn,u->uu,Jac0,params);
                 MMfun0(neq,method,Jac0,MM0);
           break;

           case 12:case 13:
	        MMfunOsoa0(neq,method,tn,system,thestatptr,MM0);
           break;

           default:
               printf("error. algorithm\n");
           break;
           }

           info = GETRF0(LAPACK_ROW_MAJOR, n, m, MM0, lda, ipiv);
           
     }

 /*  
     printf("Jakobiarra\n");
     mynorm=0.;     
     for (i=0; i<n; i++)
     {
            for (j=0; j<n; j++) 
             {
              mynorm+=MM0[i*n+j]*MM0[i*n+j];
              }
     }
     printf("\nNorma=%.20lg ***********",sqrt(mynorm));
*/

     if (thestatptr->convergence==SUCCESS)
      { options->iteration[1](system,u,tn,h,method,thestatptr,MM0,ipiv);
        thestatptr->itcount++;
      }  


     free(MM0);
     return;

}

/****************************************************************************************/
/* 											*/
/*   Newton_it_Mix2  								        */
/* 											*/
/*        Azkeneko iterazioa doitasun altuan						*/
/* 											*/
/****************************************************************************************/


void Newton_it_Mix2  ( ode_sys *system, solution *u, val_type tn,val_type h, 
                       toptions *options,gauss_method *method,solver_stat *thestatptr)

{ 

/*----------------  First: initialization -------------------------------*/
     int neq,ns;

     neq=system->neq;
     ns=method->ns;

/*------ declarations --------------------------------------------------*/

     int i,is,j;
     bool D0,iter0;			
     val_type0 difftest,DMin[neq*ns];
     val_type0 *z0,*li0;
     parameters *params;
     val_type0 zold0[neq*ns];
     lowfloat0 *MM0;

     int n,m,lda,info;
     int ipiv[neq*(ns+1)];
     lowfloat0 Jac0[neq*neq];

     val_type *z,*li;
     val_type mynorm;
     
/* ----------- implementation  --------------------------------------*/ 


     MM0=(lowfloat0 *)malloc((neq*(ns+1))*(neq*(ns+1))*sizeof(lowfloat0));

     z0=thestatptr->z0;
     li0=thestatptr->li0;
     params=&system->params;
     z=thestatptr->z;
     li=thestatptr->li;

     for (is=0; is<ns; is++)
          for (i=0; i<neq; i++)
          {
                 DMin[is*neq+i]=INF;
                 z0[is*neq+i]=z[is*neq+i];
                 li0[is*neq+i]=li[is*neq+i];
          }      

     iter0=true;
     thestatptr->itcount=0;

/* ------LU factorization of a general n-by-m matrix A--------------*/

     n=neq*(ns+1);
     m=n;
     lda=n;


     switch (options->algorithm)
     { 
     case 11: 
           system->jac0(neq,tn,u->uu,Jac0,params);
           MMfun0(neq,method,Jac0,MM0);
     break;

     case 12: case 13:
	   MMfunOsoa0(neq,method,tn,system,thestatptr,MM0);
     break;

     default:
          printf("error. algorithm\n");
     break;
    }


     info = GETRF0(LAPACK_ROW_MAJOR, n, m, MM0, lda, ipiv);

     if( info > 0 )
     {
	printf( "The diagonal element of the triangular factor of A,\n" );
	printf( "U(%i,%i) is zero, so that A is singular;\n", info, info );
	printf( "the solution could not be computed.\n" );
	exit( 1 );
     }

     while (iter0 && thestatptr->itcount<MAXIT)
     { 
          
           for (is=0; is<ns;is++)
                for (i=0; i<neq; i++) zold0[neq*is+i]=thestatptr->z0[neq*is+i];
              
           options->iteration[0](system,u,tn,h,method,thestatptr,MM0,ipiv); 
           if (options->rdigits>0) RemoveDigitsFcn0(system,method,z0,options->mrdigits);
           UpdateDMin0(system,method,&D0,&iter0,DMin,z0,zold0);
           thestatptr->itcount++;
     }


     if (thestatptr->itcount==MAXIT)
     { 
           printf("Break: step(MAXIT)=%i\n",thestatptr->itcount);
           thestatptr->convergence=FAIL; 
     }
 
     if (D0==false)
     {
           difftest=NormalizedDistance0(neq,ns,options,z0,zold0);
           if (difftest>1.)
           {
                 printf("Lack of convegence of fixed-point iteration:\
                         step=%i,difftest=%lg\n",thestatptr->itcount,difftest);
                 thestatptr->convergence=FAIL;
           }
     }
     else
     {
           thestatptr->itzero++;
     }


/* Last iteration in HighPrecision */


     for (is=0; is<ns; is++)
          for (i=0; i<neq; i++) {li[is*neq+i]=li0[is*neq+i];}

     Yi_update (u,z,system,method,thestatptr);


     switch (options->algorithm)
     { 
     case 11: 
           system->jac0(neq,tn,u->uu,Jac0,params);
           MMfun0(neq,method,Jac0,MM0);
     break;

     case 12:case 13:
	   MMfunOsoa0(neq,method,tn,system,thestatptr,MM0);
     break;

     default:
          printf("error. algorithm\n");
     break;
    }

     info = GETRF0(LAPACK_ROW_MAJOR, n, m, MM0, lda, ipiv);
  
     if( info > 0 )
     {
	printf( "The diagonal element of the triangular factor of A,\n" );
	printf( "U(%i,%i) is zero, so that A is singular;\n", info, info );
	printf( "the solution could not be computed.\n" );
	exit( 1 );
     }

/*
     printf("Jakobiarra\n");
     mynorm=0.;     
     for (i=0; i<n; i++)
     {
            for (j=0; j<n; j++) 
             {
              mynorm+=MM0[i*n+j]*MM0[i*n+j];
              }
     }
     printf("\nNorma=%.20lg ***********",sqrt(mynorm));
*/

     if (thestatptr->convergence==SUCCESS)
      { options->iteration[1](system,u,tn,h,method,thestatptr,MM0,ipiv);
        thestatptr->itcount++;
      }  


     free(MM0);
     return;

}


/****************************************************************************************/
/* 											*/
/*   MMfun	  								        */
/* 											*/
/* 											*/
/****************************************************************************************/


void MMfun0 (int neq, gauss_method *method,lowfloat0 *Jac0,lowfloat0 *MM0)

{ 

/*----------------  First: initialization -------------------------------*/
     int ns;
     int row,col;

     ns=method->ns;

     row=(ns+1)*neq;
     col=(ns+1)*neq;

/*-------------declarations --------------------------------------------*/

     int i,j;
     lowfloat0 hB[ns*ns],C[ns*(ns+1)];
     lowfloat0 Id[neq*neq],Ones[ns+1];
     lowfloat0 mm[ns*ns];		// behinbehinekoa

     lowfloat0 alpha,beta;
     int lda,ldb,ldc,i0;

/* ----------- implementation  ----------------------------------------*/ 

     /* ---- Init: MM --------------------------------*/
     for (i=0; i<row; i++)
           for (j=0; j<col; j++) MM0[i*row+j]=0.;
     for (i=0; i<(row-neq); i++) MM0[i*row+i]=1.;

     /* ---- Init: Id --------------------------------*/
     for (i=0; i<neq; i++)
           for (j=0; j<neq; j++)
                if (i==j) Id[i*neq+j]=1.;
                else Id[i*neq+j]=0.;

     /* ---- Init: Ones -----------------------------*/
     for (i=0; i<ns; i++) Ones[i]=1./2;
     Ones[ns]=-1.;

     /* ---- Init: hB, C ----------------------------*/
     for (i=0; i<ns; i++)
           for (j=0; j<ns; j++)
           {      
                 if (i==j) hB[i*ns+j]=method->hb0[i];
                     else  hB[i*ns+j]=0.;
                 C[i*(ns+1)+j]=-method->hb0[i]/2;
           }

     /* ---- j=ns column-----------------------------*/
     for (i=0; i<ns; i++)
     {
           j=ns;
           C[i*(ns+1)+j]=method->hb0[i];
     }   


     /* ---- calculate ------------------------------*/
     alpha = 1.; beta = 1.;
     lda=ns; ldb=ns; ldc=ns+1;

     for (i=0; i<(ns*ns); i++) mm[i]=method->m0[i];   //behinbehinekoa. 


     GEMM0(CblasRowMajor,CblasNoTrans,CblasNoTrans,ns,ns,ns,alpha,	
            hB, lda, mm, ldb, beta, C, ldc);

     alpha= -1.; beta= 1.;
     i0=0;
     MyKroneckerProd0(alpha,C,Jac0,beta,MM0,i0,ns,ns+1,neq,neq);
     
     alpha= 1.; beta= 1.;
     i0=ns;
     MyKroneckerProd0(alpha,Ones,Id,beta,MM0,i0,1,ns+1,neq,neq);
     

     return;

}


/****************************************************************************************/
/* 											*/
/*   MMfunOsoa	  								        */
/* 											*/
/* 											*/
/****************************************************************************************/


void MMfunOsoa0 (int neq, gauss_method *method,val_type tn,
                 ode_sys *system,solver_stat *thestatptr,lowfloat0 *MM0)
{ 

/*----------------  First: initialization -------------------------------*/
     int ns;
     int row,col;

     ns=method->ns;

     row=(ns+1)*neq;
     col=(ns+1)*neq;

/*-------------declarations --------------------------------------------*/

     int i,is,j;
     lowfloat0 hB[ns*ns],C[ns*(ns+1)];
     lowfloat0 Id[neq*neq],Ones[ns+1];
     lowfloat0 mm[ns*ns];		// behinbehinekoa
     lowfloat0 Jac[neq*neq];
     val_type0 *z0;
     parameters *params;

     lowfloat0 alpha,beta;
     int lda,ldb,ldc,i0;


/* ----------- implementation  ----------------------------------------*/ 

     z0=thestatptr->z0;
     params=&system->params;

     /* ---- Init: MM --------------------------------*/
     for (i=0; i<row; i++)
           for (j=0; j<col; j++) MM0[i*row+j]=0.;
     for (i=0; i<(row-neq); i++) MM0[i*row+i]=1.;

     /* ---- Init: Id --------------------------------*/
     for (i=0; i<neq; i++)
           for (j=0; j<neq; j++)
                if (i==j) Id[i*neq+j]=1.;
                else Id[i*neq+j]=0.;

     /* ---- Init: Ones -----------------------------*/
     for (i=0; i<ns; i++) Ones[i]=1./2;
     Ones[ns]=-1.;

     /* ---- Init: hB, C ----------------------------*/
     for (i=0; i<ns; i++)
           for (j=0; j<ns; j++)
           {      
                 if (i==j) hB[i*ns+j]=method->hb0[i];
                     else  hB[i*ns+j]=0.;
                 C[i*(ns+1)+j]=-method->hb0[i]/2;
           }

     /* ---- j=ns column-----------------------------*/
     for (i=0; i<ns; i++)
     {
           j=ns;
           C[i*(ns+1)+j]=method->hb0[i];
     }   


     /* ---- calculate ------------------------------*/
     alpha = 1.; beta = 1.;
     lda=ns; ldb=ns; ldc=ns+1;

     for (i=0; i<(ns*ns); i++) mm[i]=method->m0[i];   //behinbehinekoa. 


     GEMM0 (CblasRowMajor,CblasNoTrans,CblasNoTrans,ns,ns,ns,alpha,	
            hB, lda, mm, ldb, beta, C, ldc);

     alpha= -1.; beta= 1.;
     for (is=0; is<ns; is++)
     {
           system->jac0(neq,tn,&z0[is*neq],Jac,params);
           i0=is;
           MyKroneckerProd0(alpha,&C[is*(ns+1)],Jac,beta,MM0,i0,1,ns+1,neq,neq);

     }
     
     alpha= 1.; beta= 1.;
     i0=ns;
     MyKroneckerProd0(alpha,Ones,Id,beta,MM0,i0,1,ns+1,neq,neq);
 
     return;

}



/****************************************************************************************/
/* 											*/
/*   MyKroneckerProd  								        */
/* 											*/
/* 											*/
/****************************************************************************************/

void MyKroneckerProd0 (int alpha,lowfloat0 *A, lowfloat0 *B, 
                       int beta, lowfloat0 *C, int i0,
                       int n1, int m1, int n2, int m2)
{
/*-------------declarations --------------------------------------------*/

     int i,j,k,l,ci,cj;
 
/* ----------- implementation  ----------------------------------------*/ 

 
     for (i=0; i<n1; i++)
           for (j=0; j<m1; j++)
                 for (k=0; k<n2; k++)
                       for (l=0; l<m2; l++)
                       {
                             ci=n2*(i+i0)+k;
                             cj=m2*j+l;
                             C[ci*m1*m2+cj]=alpha*A[i*m1+j]*B[k*m2+l]+beta*C[ci*m1*m2+cj];
                       }

     return;

}


/****************************************************************************************/
/*											*/
/*      Newton_Step									*/
/*											*/
/****************************************************************************************/


int Newton_Step0   (ode_sys *system, solution *u, val_type tn,val_type h, 
                    gauss_method *method,solver_stat *thestatptr,
                    lowfloat0 *MM0, int *ipiv)
{

     int extern thread_count;

/* ---------- First initializations ------------------------------------*/

     int neq,ns;
     parameters *params;
     val_type0 *z0,*li0;
    
     neq=system->neq;
     ns=method->ns;
     params=&system->params;

     val_type0 fz[ns*neq];
     lowfloat0 fl[(ns+1)*neq];
     z0=thestatptr->z0;
     li0=thestatptr->li0;

/*------ declarations --------------------------------------------------*/

     int i,is,isn;
     int info;

     int n,lda,ldb,nrhs;
     char trans;
  
/* ----------- implementation  --------------------------------------*/


#ifdef PARALLEL
#      pragma omp parallel for num_threads(thread_count) private(isn)
#endif  
     for (is = 0; is<ns; is++)
     {
           isn=neq*is;
           params->ipar[0]=system->cod[0];
           thestatptr->fcn++;
           system->f0(neq,tn+method->hc[is],&z0[isn],&fz[isn],params); 
           for (i=0; i<neq; i++) fl[isn+i]=fz[isn+i]*method->hb0[is]-li0[isn+i];		
     }

     for (i=0; i<neq; i++) fl[neq*ns+i]=0.;

     trans = 'N';
     n = (ns+1)*neq; // number of rows of B
     nrhs=1;         // number of rigth sides
     lda = (ns+1)*neq;
     ldb = 1;

     info= GETRS0 (LAPACK_ROW_MAJOR,trans,n,nrhs,MM0,lda,ipiv,fl,ldb);

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
           for (i=0; i<neq; i++) li0[isn+i]=li0[isn+i]+fl[isn+i];
     }

     Yi_update0 (u,z0,system,method,thestatptr);
     return(0);

}	



