/*------------------------------------------------------------------------------*/
/*										*/
/*         GaussCoefficients.c							*/
/*     										*/
/*              GaussCoefficients ()	                                        */
/*		   Implicit Runge-Kutta (based Gauss quadrature)  coefficients. */
/*                 Method coefficients: m,b,c. 			                */
/*                 Interpolation coefficientes: nu.                             */
/*              			                                        */
/*              GaussCoefficientsNewton ()                                      */
/*                                                                              */      
/* -----------------------------------------------------------------------------*/

#include <GaussCoefficients.h>


/******************************************************************************/
/* 					   				      */
/* GaussCoefficients	         	  				      */
/*                                        				      */
/*									      */
/******************************************************************************/

void GaussCoefficients (char *path,gauss_method *method,toptions *options)
{
     int i,j;
     int info;
     int ns=method->ns;
     val_type sum;

     char mydir[20];
     FILE *fileM,*fileA,*fileB,*fileC,*fileNU,*fileMUNU;

     char filenameM[STRMAX];
     char filenameA[STRMAX];
     char filenameB[STRMAX];
     char filenameC[STRMAX];
     char filenameNU[STRMAX];
     char filenameMUNU[STRMAX];


     method->m =(val_type *) malloc(ns*ns*sizeof(val_type));
     method->a =(val_type *) malloc(ns*ns*sizeof(val_type));
     method->b = (val_type *) malloc(ns*sizeof(val_type));     
     method->hb = (val_type *) malloc(ns*sizeof(val_type));      
     method->c = (val_type *) malloc((ns)*sizeof(val_type));
     method->hc = (val_type *) malloc((ns)*sizeof(val_type));
     method->nu = (val_type *) malloc(ns*ns*sizeof(val_type));
     method->munu = (val_type *) malloc(ns*ns*sizeof(val_type));
     method->orderedindices=(int *) malloc(ns*sizeof(int));

  
     strcpy(mydir,path);

     switch (ns)
     { 
     case 6: 

           // orderindices (ascending order).
           method->orderedindices[0]=0;   
	   method->orderedindices[1]=5;
       	   method->orderedindices[2]=1;
	   method->orderedindices[3]=4;
	   method->orderedindices[4]=2;
	   method->orderedindices[5]=3;

           strcat(mydir,"S6/");


     break;

     case 8: 

           // orderindices (ascending order).
           method->orderedindices[0]=0;   
	   method->orderedindices[1]=7;
           method->orderedindices[2]=1;
	   method->orderedindices[3]=6;
	   method->orderedindices[4]=2;
	   method->orderedindices[5]=5;
	   method->orderedindices[6]=3;
	   method->orderedindices[7]=4;

           strcat(mydir,"S8/");

     break;

     case 9: 

           // orderindices (ascending order).
           method->orderedindices[0]=0;   
	   method->orderedindices[1]=8;
           method->orderedindices[2]=1;
	   method->orderedindices[3]=7;
	   method->orderedindices[4]=2;
	   method->orderedindices[5]=6;
	   method->orderedindices[6]=3;
	   method->orderedindices[7]=5;
	   method->orderedindices[8]=4;

           strcat(mydir,"S9/");

     break;

     case 16: 

           // orderindices (ascending order).
           method->orderedindices[0]=0;   
	   method->orderedindices[1]=15;
       	   method->orderedindices[2]=1;
	   method->orderedindices[3]=14;
	   method->orderedindices[4]=2;
	   method->orderedindices[5]=13;
	   method->orderedindices[6]=3;
	   method->orderedindices[7]=12;
	   method->orderedindices[8]=4;
	   method->orderedindices[9]=11;
	   method->orderedindices[10]=5;
	   method->orderedindices[11]=10;
	   method->orderedindices[12]=6;
           method->orderedindices[13]=9;
           method->orderedindices[14]=7;
           method->orderedindices[15]=8;

           strcat(mydir,"S16/");

     break;

     case 7: 

           // orderindices (ascending order).
           method->orderedindices[0]=0;   
	   method->orderedindices[1]=6;
           method->orderedindices[2]=1;
	   method->orderedindices[3]=5;
	   method->orderedindices[4]=2;
	   method->orderedindices[5]=4;
	   method->orderedindices[6]=3;

           strcat(mydir,"S7/");

     break;


     default:
          printf("Coefficients not defined (I) \n");
     break;

     }

     strcpy(filenameM,mydir);
     strcpy(filenameA,mydir);
     strcpy(filenameB,mydir);
     strcpy(filenameC,mydir);
     strcpy(filenameNU,mydir);
     strcpy(filenameMUNU,mydir);


     strcat(filenameM,"DMCoef.bin");
     strcat(filenameA,"DACoef.bin");
     strcat(filenameB,"DBCoef.bin");
     strcat(filenameC,"DCCoef.bin");
     strcat(filenameNU,"DNUCoef.bin");
     strcat(filenameMUNU,"DMUNUCoef.bin");

 
    fileM = fopen(filenameM,"rb");
    if (fileM == NULL) printf("File doesnt exists\n");
    fileA = fopen(filenameA,"rb");
    if (fileA == NULL) printf("File doesnt exists\n");
    fileB = fopen(filenameB,"rb");
    if (fileB == NULL) printf("File doesnt exists\n");
    fileC = fopen(filenameC,"rb");
    if (fileC == NULL) printf("File doesnt exists\n");
    fileNU = fopen(filenameNU,"rb");
    if (fileNU == NULL) printf("File doesnt exists\n");
    fileMUNU = fopen(filenameMUNU,"rb");
    if (fileMUNU == NULL) printf("File doesnt exists\n");



    info=fread(method->m, sizeof(val_type),ns*ns,fileM);
    if (info == -1) printf("Error fread command\n");
    info=fread(method->a, sizeof(val_type),ns*ns,fileA);
    if (info == -1) printf("Error fread command\n");
    info=fread(method->b, sizeof(val_type),ns,fileB);
    if (info == -1) printf("Error fread command\n");
    info=fread(method->c, sizeof(val_type),ns,fileC);
    if (info == -1) printf("Error fread command\n");
    info=fread(method->nu, sizeof(val_type),ns*ns,fileNU);
    if (info == -1) printf("Error fread command\n");
    info=fread(method->munu, sizeof(val_type),ns*ns,fileMUNU);
    if (info == -1) printf("Error fread command\n");



/*---- Verify symplectic condition -------------------------------------------*/

#    ifdef IOUT
     for (i=0; i<ns; i++)
     {
        printf("\n");
        for (j=0; j<ns; j++) printf ("%lg,", method->m[i*ns+j]+method->m[j*ns+i]-1.);
     }

     printf("\n");
#    endif

/*---- Calculate hb coefficients ---------------------------------------------*/

     sum=0.;
     for (i=1; i<ns-1; i++)
     {
        method->hb[i]=(options->h)*method->b[i];
        sum+=method->hb[i];
     }
     
     method->hb[0]=((options->h)-sum)/2;
     method->hb[ns-1]=((options->h)-sum)/2;

/*---- Calculate hc coefficients ---------------------------------------------*/

     for (i=0; i<ns; i++)
     {
        method->hc[i]=(options->h)*method->c[i];
     }


     fclose(fileM);
     fclose(fileA);
     fclose(fileB);
     fclose(fileC);
     fclose(fileNU);
     fclose(fileMUNU);

     return;

}




/******************************************************************************/
/* 					   				      */
/* GaussCoefficientsNewtonNEW            				      */
/*                                        				      */
/*									      */
/******************************************************************************/

void GaussCoefficientsNewton (char *path,ode_sys *system, 
                              gauss_method *method,toptions *options)
{

     int i,j;
     int info;
     int ns=method->ns;
     int neq=system->neq; 
     int mm=method->mm;
     int smm=method->smm;
     val_type h=options->h;
     val_type B[ns*ns],D[mm*smm],sigma[smm];
     val_type hB32[ns*ns];
     val_type m32[ns*ns];

     val_type aa,bb;
     int lda,ldb,ldc;    

     char mydir[20];
     FILE *fileQ1,*fileQ2,*fileSigma,*fileAlpha;

     char filenameQ1[STRMAX];
     char filenameQ2[STRMAX];
     char filenameSigma[STRMAX];
     char filenameAlpha[STRMAX];

     method->Q1 = (val_type *) malloc(ns*mm*sizeof(val_type));
     method->Q2 = (val_type *) malloc(ns*smm*sizeof(val_type));
     method->sigma2 = (val_type *) malloc(smm*sizeof(val_type));
     method->alpha = (val_type *) malloc(mm*sizeof(val_type));
     method->alpha2 = (val_type *) malloc(mm*sizeof(val_type));
     method->hDQ2T= (val_type *) malloc(mm*ns*sizeof(val_type));
     method->hDT= (val_type *) malloc(smm*mm*sizeof(val_type));
     method->BQ1= (val_type *) malloc(ns*mm*sizeof(val_type));
     method->BQ2= (val_type *) malloc(ns*smm*sizeof(val_type));
     method->hBAB= (val_type *) malloc((ns*ns)*sizeof(val_type));

  
     strcpy(mydir,path);

     switch (ns)
     { 
     case 6: 
           strcat(mydir,"S6/");
     break;

     case 8: 
           strcat(mydir,"S8/");
     break;

     case 9: 
           strcat(mydir,"S9/");
     break;

     case 16: 
           strcat(mydir,"S16/");
     break;

     case 7: 
           strcat(mydir,"S7/");
     break;


     default:
          printf("Coefficients not defined (II)\n");
     break;

     }

     strcpy(filenameQ1,mydir);
     strcpy(filenameQ2,mydir);
     strcpy(filenameSigma,mydir);
     strcpy(filenameAlpha,mydir);

     strcat(filenameQ1,"DQ1.bin");
     strcat(filenameQ2,"DQ2.bin");
     strcat(filenameSigma,"DSIGMA.bin");
     strcat(filenameAlpha,"DALPHA.bin");


    fileQ1 = fopen(filenameQ1,"rb");
    if (fileQ1 == NULL) printf("File doesnt exists\n");
    fileQ2 = fopen(filenameQ2,"rb");
    if (fileQ2 == NULL) printf("File doesnt exists\n");
    fileSigma = fopen(filenameSigma,"rb");
    if (fileSigma == NULL) printf("File doesnt exists\n");
    fileAlpha = fopen(filenameAlpha,"rb");
    if (fileAlpha == NULL) printf("File doesnt exists\n");

    info=fread(method->Q1, sizeof(val_type),ns*mm,fileQ1);
    if (info == -1) printf("Error fread command\n");
    info=fread(method->Q2, sizeof(val_type),ns*mm,fileQ2);
    if (info == -1) printf("Error fread command\n");
    info=fread(sigma, sizeof(val_type),mm,fileSigma);
    if (info == -1) printf("Error fread command\n");
    info=fread(method->alpha, sizeof(val_type),mm,fileAlpha);
    if (info == -1) printf("Error fread command\n");


/* ---- Precompute auxs matrices----------------------------------------------------------*/

     /*----  B-----------------------*/

     for (i=0; i<ns; i++)
     {
           for (j=0; j<ns; j++)
           {
                if (i==j) B[i*ns+j]=method->b[i]; 
                else B[i*ns+j]=0.;
           }
     }

     /*----  D-----------------------*/

     for (i=0; i<mm; i++)
     {
           for (j=0; j<smm; j++)
           {
                if (i==j) D[i*smm+j]=sigma[i]; 
                else D[i*smm+j]=0.;
           }
     }


     /* ---- sigma2= sigma**2  ---- */
     
     for (i=0; i<mm; i++) method->sigma2[i]=sigma[i]*sigma[i];

     /* ---- alpha2= alpha**2  ---- */
     
     for (i=0; i<mm; i++) method->alpha2[i]=method->alpha[i]*method->alpha[i];

     /* ---- hDT -------------- */

     for (i=0; i<smm; i++)
          for (j=0; j<mm; j++) 
          {
                 if (i==j) method->hDT[i*mm+j]=h*sigma[i];
                 else method->hDT[i*mm+j]=0.; 
          } 

     /* ----  hDQ2T ---------- */

     aa = h; bb = 0.;
     lda=smm; ldb=smm; ldc=ns;

     GEMM (CblasRowMajor,CblasNoTrans,CblasTrans,mm,ns,smm,aa,	
           D, lda,method->Q2, ldb, bb, method->hDQ2T, ldc);

     /* -------  BQ1 ---------- */

     aa = 1.; bb = 0.;
     lda=ns; ldb=mm; ldc=mm;

     GEMM (CblasRowMajor,CblasNoTrans,CblasNoTrans,ns,mm,ns,aa,	
           B, lda,method->Q1, ldb, bb, method->BQ1, ldc);


     /* -------  BQ2 ---------- */

     aa = 1.; bb = 0.;
     lda=ns; ldb=smm; ldc=smm;

     GEMM (CblasRowMajor,CblasNoTrans,CblasNoTrans,ns,smm,ns,aa,	
           B ,lda,method->Q2, ldb, bb, method->BQ2, ldc);


     /* -------- hBAB -----------------------*/


     for (i=0; i<ns; i++)
           for (j=0; j<ns; j++)
           {      
                 if (i==j) hB32[i*ns+j]=method->hb[i];
                     else  hB32[i*ns+j]=0.;                 
                 m32[i*ns+j]=method->m[i*ns+j];
           }

     aa = 1.; bb = 0.;
     lda=ns; ldb=ns; ldc=ns;

     GEMM (CblasRowMajor,CblasNoTrans,CblasNoTrans,ns,ns,ns,aa,	
           hB32, lda, m32, ldb, bb, method->hBAB, ldc);

/*  ------ Tests -------------------------- */


/*

     printf("\nQ1\n");
     for (i=0; i<ns; i++)
     {
           for (j=0; j<mm; j++) printf("%lg,",method->Q1[i*mm+j]);
           printf("\n");
     }

     printf("\nQ2\n");
     for (i=0; i<ns; i++)
     {
           for (j=0; j<smm; j++) printf("%lg,",method->Q2[i*smm+j]);
           printf("\n");
     }

     printf("\nsigma\n");
     for (i=0; i<smm; i++) printf("%lg,",sigma[i]);

     printf("\nalpha\n");
     for (i=0; i<mm; i++) printf("%lg,",method->alpha[i]);

     printf("\nhDQ2T\n");
     for (i=0; i<mm; i++)
     {
           for (j=0; j<ns; j++) printf("%lg,",method->hDQ2T[i*ns+j]);
           printf("\n");
     }

     printf("\nhDT\n");
     for (i=0; i<smm; i++)
     {
           for (j=0; j<mm; j++) printf("%lg,",method->hDT[i*mm+j]);
           printf("\n");
     }


     printf("\nBQ1\n");
     for (i=0; i<ns; i++)
     {
           for (j=0; j<mm; j++) printf("%lg,",method->BQ1[i*mm+j]);
           printf("\n");
     }


     printf("\nBQ2\n");
     for (i=0; i<ns; i++)
     {
           for (j=0; j<smm; j++) printf("%lg,",method->BQ2[i*smm+j]);
           printf("\n");
     }

     printf("\n Coeficcients: hBAB\n");
     for (i=0; i<ns; i++)
     {
         for (j=0; j<ns; j++) printf("%lg,",method->hBAB[i*ns+j]);
         printf("\n");
     }

*/

     fclose(fileQ1);
     fclose(fileQ2);
     fclose(fileSigma);
     fclose(fileAlpha);

     return;

}




