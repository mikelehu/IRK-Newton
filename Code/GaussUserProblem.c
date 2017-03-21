/*------------------------------------------------------------------------------*/
/*										*/
/*    GaussUserProblem.c							*/
/*	Functions:								*/
/*         test_jac()                                                           */
/*	   Ode1()= OdePendulumStiff():						*/
/*	   Ham1()= HamPendulumStiff():						*/
/*	   Jac1()= JacPendulumStiff():						*/
/*         Ode2()= OdeNbody():							*/
/*	   Ham2()= HamNBody():							*/
/*	   Jac2()= JacNBody():							*/
/*										*/
/*	   Ode3,...,Ode10 : empty functions					*/
/*	   Ham3,...,Ham10 : empty functions					*/
/*	   Jac3,...,Jac10 : empty functions					*/
/*										*/
/*										*/
/*------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <def.h>
#include <quadmath.h>


void test_jac (int neq, val_type t,val_type *u,void (*Ode)(), 
               val_type *Jac,parameters *params)
{

 /* ----------Test Jacobian -----------------------------------------*/


     int i,j;

     val_type epsilon;
     val_type Y1[neq],Y2[neq];
     val_type FY1[neq],FY2[neq];
     val_type APPROX[neq],V[neq];

     epsilon=pow(2,-32);


 /* ----------- implementation  ---------------------------------------*/


     printf("\nTest Jacobian: A=(f(y+epsilon)-f(y-epsilon))/(2*epsilon) ; B=Jac ; A-B \n"); 
  

     for (j=0; j<neq; j++)
     {
           printf("\nj=%i. df/duj\n",j);
           for (i=0; i<neq; i++) V[i]=0.;
           V[j]=1.;          

           for (i=0; i<neq; i++)
           {
               Y1[i]=u[i]+epsilon*V[i];
               Y2[i]=u[i]-epsilon*V[i];
           }

           Ode(neq,t,Y1,FY1,params); 
           Ode(neq,t,Y2,FY2,params);
           
           for (i=0; i<neq; i++)
           {          
               APPROX[i]= (FY1[i]-FY2[i])/(2*epsilon);
               printf("%lg,%lg,%lg\n",APPROX[i], Jac[i*neq+j], (APPROX[i]-Jac[i*neq+j]));
           }
     }
 

     return;
}




/*------------------------------------------------------------------------------*/
/*										*/
/*       Pendulum Stiff Problem:						*/
/*	    Ode1()=OdePendulum Stiff						*/
/*	    Ham1()=HamPendulum Stiff()						*/
/*	    Jac1()=JacPendulum Stiff()						*/
/*										*/
/*------------------------------------------------------------------------------*/

void Ode1 (int neq, val_type t,val_type *u,val_type *f,parameters *params)

{
    
/* ---------- First initializations -----------------------------------------*/

 /*------ declarations -------------------------------------------------*/ 

     val_type C1,C2,C3,C4,C5,C6,C7,C8;
     val_type Q1,Q2,P1,P2;
     val_type P2P1, cosQ2, cos2Q2;
     val_type sinQ1,sinQ2,sin2Q2,sinQ1Q2;
     val_type aux1,aux2;
     val_type daux1,daux2,daux3,daux4;

 /* ----------- implementation  ---------------------------------------*/

     C1=params->rpar[0];
     C2=params->rpar[1];
     C3=params->rpar[2];
     C4=params->rpar[3];
     C5=params->rpar[4];
     C6=params->rpar[5];
     C7=params->rpar[6];
     C8=params->rpar[7];

     Q1=u[0];
     Q2=u[1];
     P1=u[2];
     P2=u[3];
           

     P2P1=P2-P1;
     cosQ2= COS(Q2);
     cos2Q2= COS(2*Q2);
     sinQ1= SIN(Q1);
     sinQ2= SIN(Q2);
     sin2Q2= SIN(2*Q2);
     sinQ1Q2= SIN(Q1+Q2);
     
     aux1=C4-C5*cos2Q2;
     aux2=C3*cosQ2;

     daux1=(aux2*P2+2*C2*P2P1)/aux1;
     daux2=(2*C1*P2+aux2*P2P1)/aux1;
     daux3=(C3*P2*P2P1)/aux1;
     daux4=C5*(C1*P2*P2+aux2*P2*P2P1+C2*P2P1*P2P1)/(aux1*aux1);
    

     f[0]= -daux1;
     f[1]= daux1+daux2;    
     f[2]= -C7*sinQ1Q2-C6*sinQ1;
     f[3]= -2*C8*Q2+sinQ2*daux3+2*sin2Q2*daux4-C7*sinQ1Q2;


/*
     printf("\nOdefun6\n");
     printf("%.20lg,%.20lg,%.20lg,%.20lg\n", f[0],f[1],f[2],f[3]);           
*/
                        
     return;

}



val_type Ham1 (int neq,solution *u,parameters *params)
{

/* ---------- First initializations -----------------------------------------*/

     val_type *uu;
     uu=u->uu;

/*------ declarations -------------------------------------------------*/ 

     val_type C1,C2,C3,C4,C5,C6,C7,C8;
     val_type Q1,Q2,P1,P2;
     val_type P2P1, cosQ1,cosQ2, cos2Q2, cosQ1Q2;
     val_type H;


/* ----------- implementation  ---------------------------------------*/


     C1=params->rpar[0];
     C2=params->rpar[1];
     C3=params->rpar[2];
     C4=params->rpar[3];
     C5=params->rpar[4];
     C6=params->rpar[5];
     C7=params->rpar[6];
     C8=params->rpar[7];

     Q1=uu[0];
     Q2=uu[1];
     P1=uu[2];
     P2=uu[3];

     P2P1=P2-P1;
     cosQ1= COS(Q1);
     cosQ2= COS(Q2);
     cos2Q2= COS(2*Q2);
     cosQ1Q2= COS(Q1+Q2);

     H=((C1*P2*P2+C2*P2P1*P2P1+C3*P2*P2P1*cosQ2) / (C4-C5*cos2Q2))
        + C8*Q2*Q2- C6*cosQ1-C7*cosQ1Q2;

/*
     printf("\nHamiltonian\n");
     printf("H=%lg\n",H);
*/

     return(H);


}


void Jac1 (int neq, val_type t,val_type *u,val_type *Jac,parameters *params)

{
    

 /*------ declarations -------------------------------------------------*/ 

     int i,j;

     val_type c1,c2,c3,c4,c5,c6,c7,c8;
     val_type Q1,Q2,P1,P2;
     val_type P2P1,cosQ1,cosQ2,cos2Q2,cosQ1Q2,sinQ2,sin2Q2;
     val_type aux1,aux2,aux3,aux4,aux5,aux6;
     val_type aux11,aux12,aux13,aux21,aux22,aux23,aux24;
     val_type aux41,aux42,aux43,aux44;
     val_type aux51,aux52,aux53;
     val_type daux1,daux41,daux42;
     val_type aux101,aux102,aux103,aux104,aux105,aux106,aux107,aux108,aux109;
     val_type aux110,aux111,aux112,aux113,aux114;


 /* ----------- implementation  ---------------------------------------*/

     c1=params->rpar[0];
     c2=params->rpar[1];
     c3=params->rpar[2];
     c4=params->rpar[3];
     c5=params->rpar[4];
     c6=params->rpar[5];
     c7=params->rpar[6];
     c8=params->rpar[7];
    
     Q1=u[0];
     Q2=u[1];
     P1=u[2];
     P2=u[3];
           

     P2P1=P2-P1;
     cosQ1 = cos(Q1);
     cosQ2 = cos(Q2);
     cos2Q2 = cos(2*Q2);
     cosQ1Q2 = cos(Q1+Q2);
     sinQ2 = sin(Q2);
     sin2Q2 = sin(2*Q2);

     aux1 = c4-c5*cos2Q2;
     aux2 = c3*cosQ2;
     aux3 = c6*cosQ1;
     aux4 = c7*cosQ1Q2;
     aux5 = 2*c5*sin2Q2;
     aux6 = -c3*sinQ2;

     aux11 = aux2*P2;
     aux12 = aux2*P2P1;
     aux13 = aux2*P2*P2P1;

     aux21 = aux6*P2;
     aux22 = aux6*P2P1;
     aux23 = aux6*P2*P2P1;
     aux24 = -2*aux23*aux5;

     aux41 = aux1*aux1;
     aux42 = aux1*aux41;
     aux43 = P2*P2;
     aux44 = P2P1*P2P1;

     aux51 = c1*aux43+aux13+c2*aux44;
     aux52 = (aux11+2*c2*P2P1);
     aux53 = 2*c1*P2+aux12;

     daux1 = 1/aux1;
     daux41 = 1/aux41;
     daux42 = 1/aux42;

     aux101 = daux1*aux21-aux5*aux52*daux41;
     aux102 = -2*c2*daux1;
     aux103 = (aux2+2*c2)*daux1;
     aux104 = aux22*daux1-aux5*aux53*daux41;
     aux105 = 2*(aux2+c1+c2)*daux1;
     aux106 = -c3*P2*daux1*sinQ2;
     aux107 = (c3*P2+c3*P2P1)*sinQ2*daux1;
     aux108 = 2*c5*aux5*aux51*daux42;
     aux109 = aux13*daux1;
     aux110 = 4*c5*aux51*cos2Q2*daux41;
     aux111 = aux24*daux41;
     aux112 = -2*c5*aux52*sin2Q2*daux41;
     aux113 = 2*c5*daux41*(aux53+aux52)*sin2Q2;
     aux114 = 2*sin2Q2*aux108;
     
     Jac[0]=0.;   
     Jac[1]=-aux101;  
     Jac[2]=-aux102;  
     Jac[3]=-aux103;  
          
     Jac[4]=0.;  
     Jac[5]=aux101+aux104; 
     Jac[6]=-aux103;
     Jac[7]=aux105;

     Jac[8]=-aux3-aux4;
     Jac[9]=-aux4;
     Jac[10]=0.;
     Jac[11]=0.;

     Jac[12]=-aux4;
     Jac[13]=-2*c8+aux110+aux109-aux4-aux111- aux114;
     Jac[14]=aux106+aux112;                    
     Jac[15]=aux107+aux113;  


#    ifdef TESTJAC
     printf("\n ******** TEST JACOBIAN ******** \n");
     test_jac (neq,t,u,Ode1,Jac,params);
#    endif


/*
     printf("\nJacobian\n");
     for (i=0; i<4; i++)
     {
         for (j=0; j<4; j++) printf("%lg,",Jac[i*4+j]);
         printf("\n");     
     }

*/

     return;
}

/*------------------------------------------------------------------------------*/
/*										*/
/*       NBody Problem: 							*/
/*	    Ode2()=OdeNBody							*/
/*	    Ham2()=HamNBody():							*/
/*	    Jac2()=JacNBody():							*/
/*										*/
/*------------------------------------------------------------------------------*/

void Ode2 (int neq, val_type t,val_type *u,val_type *f,parameters *params)

{
/* ---------- First initializations ------------------------------------*/
  
     int dim;
     dim=3;
    
/*------ declarations -------------------------------------------------*/ 

     int i,id,i1,i2,j,j1,j2;
     int ix,jx;
     int nd,nbody;
     val_type d3,qij;
     val_type *Gm;

/* ----------- implementation  ---------------------------------------*/

/*   params->ipar=orderplanet (ascending order)*/

     nbody=neq/(2*dim);
     nd=neq/2;
     Gm=params->rpar;
    
     switch (params->eval)
     {
     case 1: /*OdePlanetsq*/

           for (ix=0; ix<nbody; ix++)
           {
                 i=params->ipar[ix];		 
                 i1=i*dim;
                 i2=nd+i1;
                 for (id=0; id<dim; id++) f[i1+id]=u[i2+id];
           }
              
     break;

     case 2: /*OdePlanetsv*/      

           for (ix=0; ix<nbody; ix++)
           {
                 i=params->ipar[ix];
                 i1=i*dim;
                 i2=nd+i1;
                 for (id=0; id<dim; id++) f[i2+id]=0.;
           }
        
           for (ix=0; ix<nbody; ix++)
           {
                 i=params->ipar[ix];
                 i1=i*dim;
                 i2=nd+i1;

                 for (jx=ix+1; jx<nbody; jx++)
                 {
                       j=params->ipar[jx];
                       j1=j*dim;
                       j2=nd+j1;
                       d3=0.;
                       for (id=0; id<dim; id++) d3+=(u[i1+id]-u[j1+id])*(u[i1+id]-u[j1+id]);
                       d3=SQRT(d3)*d3;

                       for (id=0; id<dim; id++) 
                       {
                             qij=(u[i1+id]-u[j1+id])/d3;
                             f[i2+id]-=Gm[j]*qij;
                             f[j2+id]+=Gm[i]*qij;
                       }   
                 }
           }
         
     break;
                                     
     default:     

        
           for (ix=0; ix<nbody; ix++)
           {
                 i=params->ipar[ix];
                 i=ix;
                 i1=i*dim;
                 i2=nd+i1;
                 for (id=0; id<dim; id++) 
                 {
                       f[i1+id]=u[i2+id];
                       f[i2+id]=0.;
                 }
           }

           for (ix=0; ix<nbody; ix++)
           {
                 i=params->ipar[ix];
                 i1=i*dim;
                 i2=nd+i1;
                 for (jx=ix+1; jx<nbody; jx++)
                 {
                       j=params->ipar[jx];
                       j1=j*dim;
                       j2=nd+j1;
                       d3=0.;
                       for (id=0; id<dim; id++) d3+=(u[i1+id]-u[j1+id])*(u[i1+id]-u[j1+id]);
                       d3=SQRT(d3)*d3;

                       for (id=0; id<dim; id++) 
                       {
                             qij=(u[i1+id]-u[j1+id])/d3;
                             f[i2+id]-=Gm[j]*qij;
                             f[j2+id]+=Gm[i]*qij;
                       }   
                 }

           }
 
     break;

}

/*
    printf("\nOdefun2\n");
    for (i=0; i<neq;i++) printf("%lg,", f[i]);
    printf("\n");
*/

    return ;

}


val_type Ham2 (int neq,solution *u,parameters *params)
{

/* ---------- First initializations ------------------------------------*/
  
     int dim;
     dim=3;

/*------ declarations -------------------------------------------------*/ 

     val_type *uu;  
     int i,j,i1,j1;
     int ix,jx;
     int nbody,node,nd;
     val_type *d,*Gm;
     val_type H,Pot;
 
/* ----------- implementation  ----------------------------------------*/
     uu=u->uu;
     Gm=params->rpar;
     nbody=neq/(2*dim);
     d = malloc(nbody*nbody*sizeof(val_type));
     node=neq;
     nd=node/2;

     H=0.;
 
     for (ix=0; ix<nbody; ix++)
     {
           i=params->ipar[ix];
           i1=dim*i;
           for (jx=ix+1; jx<nbody; jx++)
           {
                 j=params->ipar[jx];
                 j1=dim*j;
                 d[i*nbody+j]=SQRT(POW(uu[i1]-uu[j1],2)+
                                   POW(uu[i1+1]-uu[j1+1],2)+POW(uu[i1+2]-uu[j1+2],2));
                 d[j*nbody+i]=d[i*nbody+j];
           }
     }

     for (ix=0; ix<nbody; ix++)
     {
           i=params->ipar[ix];
           i1=nd+dim*i;
           H+=Gm[i]*(POW(uu[i1],2)+POW(uu[i1+1],2)+POW(uu[i1+2],2));
     }

     H=H/2.;
     Pot=0;

     for (ix=0; ix<nbody-1; ix++)
     {
          i=params->ipar[ix];
          for (jx=ix+1; jx<nbody; jx++)
          {
               j=params->ipar[jx];
               Pot+=Gm[i]*Gm[j]/d[i*nbody+j];
          }
     }

     H=H-Pot;
 
     free(d);
     return(H);

}


void Jac2 (int neq, val_type t,val_type *u,val_type *Jac,parameters *params)

{

/* ---------- First initializations ------------------------------------*/
  
     int dim,nbody;
     dim=3;

     nbody=neq/(2*dim);
    
/*------ declarations -------------------------------------------------*/ 

     int i,id,i1,i2,j,j1;
     int nd;
     val_type d3[nbody*nbody],d5[nbody*nbody];
     val_type *Gm;
     val_type Gmid3,Gmid5xx,Gmid5xy,Gmid5xz,Gmid5yy,Gmid5zz,Gmid5yz;
     val_type xx,yy,zz;
     val_type aux,aux1;

/* ----------- implementation  ---------------------------------------*/
     
     nd=neq/2;
     Gm=params->rpar;


     for (i=0; i<neq; i++)
         for (j=0; j<neq; j++) Jac[i*neq+j]=0.;           


     for (i=0; i<nbody; i++)
     {
           i1=dim*i;

           for (j=i+1; j<nbody; j++)
           {
                 j1=dim*j;
                
                 xx=(u[j1]-u[i1]);
                 yy=(u[j1+1]-u[i1+1]);
                 zz=(u[j1+2]-u[i1+2]);
                 aux=xx*xx+yy*yy+zz*zz;
                 d3[i*nbody+j]=POW(aux,3./2);
                 d5[i*nbody+j]=POW(aux,5./2);
                 d3[j*nbody+i]=d3[i*nbody+j];
                 d5[j*nbody+i]=d5[i*nbody+j];

           }
     }

    

     for (i=0; i<nbody; i++)
     {

      i1=i*dim;
      i2=nd+i1;

           for (id=0; id<dim; id++) Jac[(i1+id)*neq+(i2+id)]=1.;

           Gmid3=0.; 
           Gmid5xx=0.;
           Gmid5xy=0.;
           Gmid5xz=0.;
           Gmid5yz=0.;
           Gmid5yy=0.;
           Gmid5zz=0.;

           for (j=0; j<nbody; j++)
           {    

                 j1=j*dim;
                 if (i!=j)
                 {
                       Gmid3-=Gm[j]/d3[i*nbody+j];
                       xx=(u[j1]-u[i1]);
                       yy=(u[j1+1]-u[i1+1]);
                       zz=(u[j1+2]-u[i1+2]);
                       Gmid5xx+=Gm[j]*(xx*xx)/d5[i*nbody+j];
                       Gmid5xy+=Gm[j]*(xx*yy)/d5[i*nbody+j];
                       Gmid5xz+=Gm[j]*(xx*zz)/d5[i*nbody+j];
                       Gmid5yz+=Gm[j]*(yy*zz)/d5[i*nbody+j];
                       Gmid5yy+=Gm[j]*(yy*yy)/d5[i*nbody+j];
                       Gmid5zz+=Gm[j]*(zz*zz)/d5[i*nbody+j];
                       
                 }   
           }


           for (j=0; j<nbody; j++)
           {
                 
                 j1=j*dim;

                 if (i==j)
                 {
                       Jac[i2*neq+j1]= Gmid3+3*Gmid5xx;       				        //jac[0]
                       Jac[i2*neq+j1+1]= 3*Gmid5xy;            				        //jac[1]
                       Jac[i2*neq+j1+2]= 3*Gmid5xz;                                             //jac[2]
                       Jac[(i2+1)*neq+j1]=Jac[i2*neq+j1+1];      			        //jac[3]=jac[1]
                       Jac[(i2+1)*neq+j1+1]=Gmid3+3*Gmid5yy;    			        //jac[4]
                       Jac[(i2+1)*neq+j1+2]=3*Gmid5yz;            				//jac[5]
                       Jac[(i2+2)*neq+j1]=Jac[i2*neq+j1+2];        				//jac[6]=jac[2]
                       Jac[(i2+2)*neq+j1+1]=Jac[(i2+1)*neq+j1+2]; 				//jac[7]=jac[5]
                       Jac[(i2+2)*neq+j1+2]=Gmid3+3*Gmid5zz;      				//jac[8]

                 }
                 else
                 {
                       xx=(u[j1]-u[i1]);
                       yy=(u[j1+1]-u[i1+1]);
                       zz=(u[j1+2]-u[i1+2]);

                       aux1= Gm[j]/d3[i*nbody+j];

                       Jac[i2*neq+j1]= aux1-3*Gm[j]*xx*xx/d5[i*nbody+j];    			//jac[0]
                       Jac[i2*neq+j1+1]= -3*Gm[j]*xx*yy/d5[i*nbody+j];	     			//jac[1]
                       Jac[i2*neq+j1+2]= -3*Gm[j]*xx*zz/d5[i*nbody+j];	   	         	//jac[2]

                       Jac[(i2+1)*neq+j1]= Jac[i2*neq+j1+1];					//jac[3]=jac[1]
                       Jac[(i2+1)*neq+j1+1]= aux1-3*Gm[j]*yy*yy/d5[i*nbody+j];			//jac[4]
                       Jac[(i2+1)*neq+j1+2]= -3*Gm[j]*yy*zz/d5[i*nbody+j];			//jac[5]

                       Jac[(i2+2)*neq+j1]= Jac[i2*neq+j1+2];					//jac[6]=jac[2]
                       Jac[(i2+2)*neq+j1+1]= Jac[(i2+1)*neq+j1+2];				//jac[7]=jac[5]
                       Jac[(i2+2)*neq+j1+2]= aux1-3*Gm[j]*zz*zz/d5[i*nbody+j];	                //jac[8]


                 } 

           }

     }
    

#    ifdef TESTJAC
     printf("\n ******** TEST JACOBIAN ******** \n");
     test_jac (neq,t,u,Ode2,Jac,params);
#    endif

    return ;

}


/*------------------------------------------------------------------------------*/
/*										*/
/*        Problem-3:	 							*/
/*	    Ode3()=:								*/
/*	    Ham3()=:								*/
/*										*/
/*------------------------------------------------------------------------------*/

void Ode3 (int neq, val_type t,val_type *u,val_type *f,parameters *params)
{
                         
     return ;

}

val_type Ham3 (int neq,solution *u,parameters *params)
{

     return(0);

}


void Jac3 (int neq, val_type t,val_type *u,val_type *Jac,parameters *params)
{

     return ;

}


/*------------------------------------------------------------------------------*/
/*										*/
/*        Problem-4:	 							*/
/*	    Ode4()= 								*/
/*	    Ham4()= 								*/
/*	    Jac4()=								*/
/*										*/
/*------------------------------------------------------------------------------*/

void Ode4 (int neq, val_type t,val_type *u,val_type *f,parameters *params)
{                      
     return;

}

val_type Ham4 (int neq,solution *u,parameters *params)
{
     return(0);

}

void Jac4 (int neq, val_type t,val_type *u,val_type *Jac,parameters *params)
{         
     return;

}


/*------------------------------------------------------------------------------*/
/*										*/
/*        Problem-5:   Double PEndulum NEW 13.12.2016				*/
/*	    Ode5()=								*/
/*	    Ham5()=								*/
/*										*/
/*------------------------------------------------------------------------------*/

void Ode5 (const int neq,const val_type t,const val_type *u,
           val_type *f, const parameters *params)
{

                        
     return;

}

val_type Ham5 (const int neq, const solution *u, const parameters *params)
{


     return(0);


}

void Jac5 (int neq, val_type t,val_type *u,val_type *Jac,parameters *params)
{


     return;
}

/*------------------------------------------------------------------------------*/
/*										*/
/*        Problem-6:   								*/
/*	    Ode6(), Ham6(), Jac6()						*/
/*	    									*/
/*										*/
/*------------------------------------------------------------------------------*/

void Ode6 (const int neq,const val_type t,const val_type *u,
           val_type *f, const parameters *params)
{

                        
     return;

}

val_type Ham6 (const int neq, const solution *u, const parameters *params)
{

     return(0);


}

void Jac6 (int neq, val_type t,val_type *u,val_type *Jac,parameters *params)
{


     return;
}


/*------------------------------------------------------------------------------*/
/*										*/
/*        Problem-7:	 							*/
/*	    Ode7()=								*/
/*	    Ham7()=								*/
/*	    Jac7()=								*/
/*										*/
/*------------------------------------------------------------------------------*/

void Ode7 (int neq, val_type t,val_type *u,val_type *f,parameters *params)
{
     return;
}

val_type Ham7 (int neq,solution *u,parameters *params)
{
     return(0.);

}

void Jac7 (int neq, val_type t,val_type *u,val_type *Jac,parameters *params)
{
     return;
}

/*------------------------------------------------------------------------------*/
/*										*/
/*        Problem-8:	 							*/
/*	    Ode8()=								*/
/*	    Ham8()=								*/
/*	    Jac8()=								*/
/*										*/
/*------------------------------------------------------------------------------*/

void Ode8 (int neq, val_type t,val_type *u,val_type *f,parameters *params)
{
     return;
}

val_type Ham8 (int neq,solution *u,parameters *params)
{
     return(0.);

}

void Jac8 (int neq, val_type t,val_type *u,val_type *Jac,parameters *params)
{
     return;
}

/*------------------------------------------------------------------------------*/
/*										*/
/*        Problem-9:	 							*/
/*	    Ode9()=								*/
/*	    Ham9()=								*/
/*	    Jac9()=								*/
/*										*/
/*------------------------------------------------------------------------------*/

void Ode9 (int neq, val_type t,val_type *u,val_type *f,parameters *params)
{
     return;
}

val_type Ham9 (int neq,solution *u,parameters *params)
{
     return(0.);

}

void Jac9 (int neq, val_type t,val_type *u,val_type *Jac,parameters *params)
{
     return;
}

/*------------------------------------------------------------------------------*/
/*										*/
/*        Problem-10:	 							*/
/*	    Ode10()=								*/
/*	    Ham10()=								*/
/*	    Jac10()=								*/
/*										*/
/*------------------------------------------------------------------------------*/

void Ode10 (int neq, val_type t,val_type *u,val_type *f,parameters *params)
{
     return;
}

val_type Ham10 (int neq,solution *u,parameters *params)
{
     return(0.);

}

void Jac10 (int neq, val_type t,val_type *u,val_type *Jac,parameters *params)
{
     return;
}


