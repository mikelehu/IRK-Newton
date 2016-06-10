/*------------------------------------------------------------------------------*/
/*										*/
/*    GaussUserProblem0.c							*/
/*	Functions:								*/
/*	   Ode1()= OdePendulum():						*/
/*	   Ham1()= HamPendulum():						*/
/*	   Jac1()= JacPendulum():						*/
/*         Ode2()= OdeNbody():							*/
/*	   Ham2()= HamNBody():							*/
/*	   Jac2()= JacNBody():							*/
/*										*/
/*	   Ode3,...,Ode10 : empty functions					*/
/*	   Ham3,...,Ham10 : empty functions					*/
/*	   Jac3,...,Jac10 : empty functions					*/
/*										*/
/*         Ode11= OdePendulumX(): ideal integrator				*/
/*         Ode12= OdeNBodyX():    ideal integrator				*/
/*										*/
/*------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <def.h>
#include <quadmath.h>


/*------------------------------------------------------------------------------*/
/*										*/
/*       OdePendulum Problem: 							*/
/*	    Ode1()=OdePendulum							*/
/*	    Ham1()=HamPendulum():						*/
/*	    Jac1()=JacPendulum():						*/
/*										*/
/*------------------------------------------------------------------------------*/

void Ode1low (int neq, val_type t,val_type0 *u,val_type0 *f,parameters *params)

{
    
 /* ---------- First initializations -----------------------------------------*/

 /*------ declarations -------------------------------------------------*/ 

     val_type0 C1,C2,C3,C4,C5,C6,C7;
     val_type0 Q1,Q2,P1,P2,Q12;
     val_type0 cosQ1Q2,sinQ1Q2,sinQ1,sinQ2;
     val_type0 aux0,aux1,dQ1aux1,aux2,aux3,aux4;

 /* ----------- implementation  ---------------------------------------*/

     C1=params->rpar0[0];
     C2=params->rpar0[1];
     C3=params->rpar0[2];
     C4=params->rpar0[3];
     C5=params->rpar0[4];
     C6=params->rpar0[5];
     C7=params->rpar0[6];

     Q1=u[0];
     Q2=u[1];
     P1=u[2];
     P2=u[3];
     Q12=(u[0]-u[1]);
           
     cosQ1Q2= COS0(Q12);
     sinQ1Q2= SIN0(Q12);
     sinQ1= SIN0(Q1);
     sinQ2= SIN0(Q2);

     aux0= C5*sinQ1Q2*sinQ1Q2;
     aux1= C4+aux0;
     dQ1aux1= 2*C5*sinQ1Q2*cosQ1Q2;
     aux2= aux1*aux1;
     aux3= C3*cosQ1Q2;
     aux4= (-1/aux2)*(C1*P1*P1+C2*P2*P2+P1*P2*aux3)*dQ1aux1-(C3*P1*P2*sinQ1Q2)/aux1;

     f[0]= (2*C1*P1+aux3*P2)/aux1;
     f[1]= (2*C2*P2+aux3*P1)/aux1;
     f[2]=-(aux4+C6*sinQ1);
     f[3]=-(-aux4+C7*sinQ2);               
                        
     return ;

}



val_type0 Ham1low (int neq,solution *u,parameters *params)
{

/* ---------- First initializations -----------------------------------------*/

     val_type0 *uu;
     uu=u->uu;

/*------ declarations -------------------------------------------------*/ 

     val_type0 C1,C2,C3,C4,C5,C6,C7;
     val_type0 Q1,Q2,P1,P2,Q12;
     val_type0 cosQ1Q2,sinQ1Q2,cosQ1,cosQ2;
     val_type0 H;

/* ----------- implementation  ---------------------------------------*/

     C1=params->rpar0[0];
     C2=params->rpar0[1];
     C3=params->rpar0[2];
     C4=params->rpar0[3];
     C5=params->rpar0[4];
     C6=params->rpar0[5];
     C7=params->rpar0[6];

     Q1=uu[0];
     Q2=uu[1];
     P1=uu[2];
     P2=uu[3];
     Q12=Q1-Q2;
           
     cosQ1Q2= COS(Q12);
     sinQ1Q2= SIN(Q12);
     cosQ1= COS(Q1);
     cosQ2= COS(Q2);

     H= ((C1*P1*P1+C2*P2*P2+C3*P1*P2*cosQ1Q2) / (C4+C5*sinQ1Q2*sinQ1Q2))
        - C6*cosQ1-C7*cosQ2;

     return(H);

}


void Jac1low (int neq, val_type t,val_type0 *u,lowfloat0 *Jac,parameters *params)

{
    
 /* ---------- First initializations -----------------------------------------*/

 /*------ declarations -------------------------------------------------*/ 

     lowfloat0 C1,C2,C3,C4,C5,C6,C7;
     lowfloat0 Q1,Q2,P1,P2,Q12;
     lowfloat0 cosQ1Q2,sinQ1Q2,cosQ1,cosQ2;
     lowfloat0 aux0,aux1,dQ1aux1,aux2,aux3;
     lowfloat0 daux5,daux6,daux9,daux10,daux17,daux22;
     lowfloat0 daux33,daux34,daux37,daux38,daux41,daux42;
     lowfloat0 daux45,daux46,daux49,daux50,daux53,daux54,daux55,daux56;
     lowfloat0 daux66,daux69,daux80,daux83,daux89,daux104;

 //    int i,j;

 /* ----------- implementation  ---------------------------------------*/

     C1=params->rpar0[0];
     C2=params->rpar0[1];
     C3=params->rpar0[2];
     C4=params->rpar0[3];
     C5=params->rpar0[4];
     C6=params->rpar0[5];
     C7=params->rpar0[6];

     Q1=u[0];
     Q2=u[1];
     P1=u[2];
     P2=u[3];
     Q12=(u[0]-u[1]);
           
    
     cosQ1Q2= COS0(Q12);
     sinQ1Q2= SIN0(Q12);
     cosQ1= COS0(Q1);
     cosQ2=COS0(Q2);

     aux0= C5*sinQ1Q2*sinQ1Q2;
     aux1= C4+aux0;
     dQ1aux1= 2*C5*sinQ1Q2*cosQ1Q2;
     aux2= aux1*aux1;
     aux3= C3*cosQ1Q2;

     daux5=-sinQ1Q2;
     daux6=sinQ1Q2;
     daux9=cosQ1Q2;
     daux10=-cosQ1Q2;
     daux17=cosQ1;
     daux22=cosQ2;
     daux33=2*C5*daux9*sinQ1Q2;
     daux34=2*C5*daux10*sinQ1Q2;
     daux37=daux33;
     daux38=daux34;
     daux41=2*C5*cosQ1Q2*daux9+2*C5*daux5*sinQ1Q2;
     daux42=2*C5*cosQ1Q2*daux10+2*C5*daux6*sinQ1Q2;
     daux45=2*aux1*daux37;
     daux46=2*aux1*daux38;
     daux49=C3*daux5;     
     daux50=C3*daux6;
     daux53=-((C3*daux9*P1*P2)/aux1)-(daux49*dQ1aux1*P1*P2)/aux2-
              (daux41*(C1*P1*P1+aux3*P1*P2+C2*P2*P2))/aux2+
              (daux45*dQ1aux1*(C1*P1*P1+aux3*P1*P2+C2*P2*P2))/(aux2*aux2)+(C3*daux37*P1*P2*sinQ1Q2)/(aux1*aux1);

     daux54=-((C3*daux10*P1*P2)/aux1)-(daux50*dQ1aux1*P1*P2)/aux2-
              (daux42*(C1*P1*P1+aux3*P1*P2+C2*P2*P2))/aux2+
              (daux46*dQ1aux1*(C1*P1*P1+aux3*P1*P2+C2*P2*P2))/(aux2*aux2)+(C3*daux38*P1*P2*sinQ1Q2)/(aux1*aux1);

     daux55=-((dQ1aux1*(2*C1*P1+aux3*P2))/aux2)-(C3*P2*sinQ1Q2)/aux1;
     daux56=-((dQ1aux1*(aux3*P1+2*C2*P2))/aux2)-(C3*P1*sinQ1Q2)/aux1;
     
     daux66=-((2*C1*P1+aux3*P2)/(aux1*aux1));
     daux69=P2/aux1;
     daux80=-((aux3*P1+2*C2*P2)/(aux1*aux1));
     daux83=P1/aux1;
     daux89=-C6;
     daux104=-C7;

     Jac[0]=daux37*daux66+daux49*daux69;
     Jac[1]=daux38*daux66+daux50*daux69;
     Jac[2]=(2*C1)/aux1;
     Jac[3]=aux3/aux1;
     Jac[4]=daux37*daux80+daux49*daux83;
     Jac[5]=daux38*daux80+daux50*daux83;
     Jac[6]=aux3/aux1;
     Jac[7]=(2*C2)/aux1;
     Jac[8]=-daux53+daux17*daux89;
     Jac[9]=-daux54;
     Jac[10]=-daux55;
     Jac[11]=-daux56;
     Jac[12]=daux53;
     Jac[13]=daux104*daux22+daux54;
     Jac[14]=daux55;
     Jac[15]=daux56;      
                       
     return ;

}

/*------------------------------------------------------------------------------*/
/*										*/
/*       NBody Problem: 							*/
/*	    Ode2()=OdeNBody							*/
/*	    Ham2()=HamNBody():							*/
/*	    Jac2()=JacNBody():							*/
/*										*/
/*------------------------------------------------------------------------------*/

void Ode2low (int neq, val_type t,val_type0 *u,val_type0 *f,parameters *params)

{
/* ---------- First initializations ------------------------------------*/
  
     int dim;
     dim=3;
    
/*------ declarations -------------------------------------------------*/ 

     int i,id,i1,i2,j,j1,j2;
     int nd,nbody;
     val_type0 d3,qij;
     val_type0 *Gm;

/* ----------- implementation  ---------------------------------------*/

     nbody=neq/(2*dim);
     nd=neq/2;
     Gm=params->rpar0;
    
     switch (params->ipar[0])
     {
     case 1: /*OdePlanetsq*/

           for (i=0; i<nbody; i++)
           {
                 i1=i*dim;
                 i2=nd+i1;
                 for (id=0; id<dim; id++) f[i1+id]=u[i2+id];
           }
              
     break;

     case 2: /*OdePlanetsv*/      

           for (i=0; i<nbody; i++)
           {
                 i1=i*dim;
                 i2=nd+i1;
                 for (id=0; id<dim; id++) f[i2+id]=0.;
           }
        
           for (i=0; i<nbody; i++)
           {
                 i1=i*dim;
                 i2=nd+i1;

                 for (j=i+1; j<nbody; j++)
                 {
                       j1=j*dim;
                       j2=nd+j1;
                       d3=0.;
                       for (id=0; id<dim; id++) d3+=POW0(u[i1+id]-u[j1+id],2);
                       d3=POW0(d3,3./2);

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
        
           for (i=0; i<nbody; i++)
           {
                 i1=i*dim;
                 i2=nd+i1;
                 for (id=0; id<dim; id++) 
                 {
                       f[i1+id]=u[i2+id];
                       f[i2+id]=0.;
                 }
           }

           for (i=0; i<nbody; i++)
           {
                 i1=i*dim;
                 i2=nd+i1;
                 for (j=i+1; j<nbody; j++)
                 {
                       j1=j*dim;
                       j2=nd+j1;
                       d3=0.;
                       for (id=0; id<dim; id++) d3+=POW0(u[i1+id]-u[j1+id],2);
                       d3=POW0(d3,3./2);

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


    return ;

}


val_type0 Ham2low (int neq,solution *u,parameters *params)
{

/* ---------- First initializations ------------------------------------*/
  
     int dim;
     dim=3;

/*------ declarations -------------------------------------------------*/ 

     val_type0 *uu;  
     int i,j,i1,j1;
     int nbody,node,nd;
     val_type0 *d,*Gm;
     val_type0 H,Pot;
 
/* ----------- implementation  ----------------------------------------*/
     uu=u->uu;
     Gm=params->rpar0;
     nbody=neq/(2*dim);
     d = malloc(nbody*nbody*sizeof(val_type));
     node=neq;
     nd=node/2;

     H=0.;
 
     for (i=0; i<nbody; i++)
     {
           i1=dim*i;
           for (j=i+1; j<nbody; j++)
           {
                 j1=dim*j;
                 d[i*nbody+j]=SQRT0(POW0(uu[i1]-uu[j1],2)+
                                    POW0(uu[i1+1]-uu[j1+1],2)+POW0(uu[i1+2]-uu[j1+2],2));
                 d[j*nbody+i]=d[i*nbody+j];
           }
     }

     for (i=0; i<nbody; i++)
     {
           i1=nd+dim*i;
           H+=Gm[i]*(POW0(uu[i1],2)+POW0(uu[i1+1],2)+POW0(uu[i1+2],2));
     }

     H=H/2.;
     Pot=0;

     for (i=0; i<nbody-1; i++)
          for (j=i+1; j<nbody; j++)
               Pot+=Gm[i]*Gm[j]/d[i*nbody+j];


     H=H-Pot;
 
     free(d);
     return(H);

}


void Jac2low (int neq, val_type t,val_type0 *u,lowfloat0 *Jac,parameters *params)

{

/* ---------- First initializations ------------------------------------*/
  
     int dim,nbody;
     dim=3;

     nbody=neq/(2*dim);
    
/*------ declarations -------------------------------------------------*/ 

     int i,id,i1,i2,j,j1,j2;
     int nd;
     val_type0 d3[nbody*nbody],d5[nbody*nbody];
     val_type0 *Gm;
     val_type0 Gmid3,Gmid5xx,Gmid5xy,Gmid5xz,Gmid5yy,Gmid5zz,Gmid5yz;
     val_type0 xx,yy,zz;
     val_type0 aux,aux1;

/* ----------- implementation  ---------------------------------------*/
     
     nd=neq/2;
     Gm=params->rpar0;


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
                 d3[i*nbody+j]=POW0(aux,3./2);
                 d5[i*nbody+j]=POW0(aux,5./2);
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
    

    return ;

}



/*------------------------------------------------------------------------------*/
/*										*/
/*        Problem-3:	 							*/
/*	    Ode3()=								*/
/*	    Ham3()=								*/
/*	    Jac3()=								*/
/*										*/
/*------------------------------------------------------------------------------*/

void Ode3low (int neq, val_type t,val_type0 *u,val_type0 *f,parameters *params)
{
     return;
}

val_type0 Ham3low (int neq,solution *u,parameters *params)
{
     return(0.);

}

void Jac3low (int neq, val_type t,val_type0 *u,lowfloat0 *Jac,parameters *params)
{
     return;
}

/*------------------------------------------------------------------------------*/
/*										*/
/*        Problem-4:	 							*/
/*	    Ode4()=								*/
/*	    Ham4()=								*/
/*	    Jac4()=								*/
/*										*/
/*------------------------------------------------------------------------------*/

void Ode4low (int neq, val_type t,val_type0 *u,val_type0 *f,parameters *params)
{
     return;
}

val_type0 Ham4low (int neq,solution *u,parameters *params)
{
     return(0.);

}

void Jac4low (int neq, val_type t,val_type0 *u,lowfloat0 *Jac,parameters *params)
{
     return;
}

/*------------------------------------------------------------------------------*/
/*										*/
/*        Problem-5:	 							*/
/*	    Ode5()=								*/
/*	    Ham5()=								*/
/*	    Jac5()=								*/
/*										*/
/*------------------------------------------------------------------------------*/

void Ode5low (int neq, val_type t,val_type0 *u,val_type0 *f,parameters *params)
{
     return;
}

val_type0 Ham5low (int neq,solution *u,parameters *params)
{
     return(0.);

}

void Jac5low (int neq, val_type t,val_type0 *u,lowfloat0 *Jac,parameters *params)
{
     return;
}

/*------------------------------------------------------------------------------*/
/*										*/
/*        Problem-6:	 							*/
/*	    Ode6()=								*/
/*	    Ham6()=								*/
/*	    Jac6()=								*/
/*										*/
/*------------------------------------------------------------------------------*/

void Ode6low (int neq, val_type t,val_type0 *u,val_type0 *f,parameters *params)
{
     return;
}

val_type0 Ham6low (int neq,solution *u,parameters *params)
{
     return(0.);

}

void Jac6low (int neq, val_type t,val_type0 *u,lowfloat0 *Jac,parameters *params)
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

void Ode7low (int neq, val_type t,val_type0 *u,val_type0 *f,parameters *params)
{
     return;
}

val_type0 Ham7low (int neq,solution *u,parameters *params)
{
     return(0.);

}

void Jac7low (int neq, val_type t,val_type0 *u,lowfloat0 *Jac,parameters *params)
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

void Ode8low (int neq, val_type t,val_type0 *u,val_type0 *f,parameters *params)
{
     return;
}

val_type0 Ham8low (int neq,solution *u,parameters *params)
{
     return(0.);

}

void Jac8low (int neq, val_type t,val_type0 *u,lowfloat0 *Jac,parameters *params)
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

void Ode9low (int neq, val_type t,val_type0 *u,val_type0 *f,parameters *params)
{
     return;
}

val_type0 Ham9low (int neq,solution *u,parameters *params)
{
     return(0.);

}

void Jac9low (int neq, val_type t,val_type0 *u,lowfloat0 *Jac,parameters *params)
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

void Ode10low (int neq, val_type t,val_type0 *u,val_type0 *f,parameters *params)
{
     return;
}

val_type0 Ham10low (int neq,solution *u,parameters *params)
{
     return(0.);

}

void Jac10low (int neq, val_type t,val_type0 *u,lowfloat0 *Jac,parameters *params)
{
     return;
}


/*------------------------------------------------------------------------------*/
/*										*/
/*       OdePendulum Problem: 							*/
/*		  Ode11= OdePendulumX (Ideal integrator)			*/
/*										*/
/*										*/
/*------------------------------------------------------------------------------*/

void Ode11low (int neq, val_type t,val_type *u,val_type *f,parameters *params)

{
 /* --------------------------------------------------------------------------*/
 /* specific Odependulum() for "ideal integrator" esperiment                  */
 /*     - input:      quadruple---> double                                    */
 /*     - operations: double                                                  */
 /*     - output:     double ---> quadruple.                                  */
 /* --------------------------------------------------------------------------*/

 /* ---------- First initializations -----------------------------------------*/

 /*------ declarations -------------------------------------------------*/ 

     int i;

     double fdouble[neq];

     double C1,C2,C3,C4,C5,C6,C7;
     double Q1,Q2,P1,P2,Q12;
     double cosQ1Q2,sinQ1Q2,sinQ1,sinQ2;
     double aux0,aux1,dQ1aux1,aux2,aux3,aux4;

 /* ----------- implementation  ---------------------------------------*/

     C1=params->rpar0[0];
     C2=params->rpar0[1];
     C3=params->rpar0[2];
     C4=params->rpar0[3];
     C5=params->rpar0[4];
     C6=params->rpar0[5];
     C7=params->rpar0[6];

     Q1=u[0];
     Q2=u[1];
     P1=u[2];
     P2=u[3];
     Q12=(u[0]-u[1]);
            
     cosQ1Q2= cos(Q12);
     sinQ1Q2= sin(Q12);
     sinQ1= sin(Q1);
     sinQ2= sin(Q2);

     aux0= C5*sinQ1Q2*sinQ1Q2;
     aux1= C4+aux0;
     dQ1aux1= 2*C5*sinQ1Q2*cosQ1Q2;
     aux2= aux1*aux1;
     aux3= C3*cosQ1Q2;
     aux4= (-1/aux2)*(C1*P1*P1+C2*P2*P2+P1*P2*aux3)*dQ1aux1-(C3*P1*P2*sinQ1Q2)/aux1;

     fdouble[0]= (2*C1*P1+aux3*P2)/aux1;
     fdouble[1]= (2*C2*P2+aux3*P1)/aux1;
     fdouble[2]=-(aux4+C6*sinQ1);
     fdouble[3]=-(-aux4+C7*sinQ2);               
              
     for (i=0; i<neq; i++) f[i]=fdouble[i];
          
     return ;

}


/*------------------------------------------------------------------------------*/
/*										*/
/*       NBody Problem:								*/
/*	       Ode12=OdeNBodyX (Ideal integrator)				*/
/*										*/
/*------------------------------------------------------------------------------*/

void Ode12low (int neq, val_type t,val_type *u,val_type *f,parameters *params)

{

 /* --------------------------------------------------------------------------*/
 /* specific Odependulum() for "ideal integrator" esperiment                  */
 /*     - input:      quadruple---> double                                    */
 /*     - operations: double                                                  */
 /*     - output:     double ---> quadruple.                                  */
 /* --------------------------------------------------------------------------*/


/* ---------- First initializations ------------------------------------*/
  
     int dim,nd,nbody;
     dim=3;
     nbody=neq/(2*dim);
     nd=neq/2;
    
/*------ declarations -------------------------------------------------*/ 

     int i,id,i1,i2,j,j1,j2;
     double d3,qij;
     val_type *Gm;

     double fdouble[neq],Udouble[neq],Gmdouble[nbody];

/* ----------- implementation  ---------------------------------------*/

     Gm=params->rpar0;

     for (i=0;i<nbody;i++) Gmdouble[i]=Gm[i];
     for (i=0;i<neq;i++) Udouble[i]=u[i];
    
     switch (params->ipar[0])
     {
     case 1: /*OdePlanetsq*/

           for (i=0; i<nbody; i++)
           {
                 i1=i*dim;
                 i2=nd+i1;
                 for (id=0; id<dim; id++) fdouble[i1+id]=Udouble[i2+id];
           }
              
     break;

     case 2: /*OdePlanetsv*/      

           for (i=0; i<nbody; i++)
           {
                 i1=i*dim;
                 i2=nd+i1;
                 for (id=0; id<dim; id++) fdouble[i2+id]=0.;
           }
        
           for (i=0; i<nbody; i++)
           {
                 i1=i*dim;
                 i2=nd+i1;

                 for (j=i+1; j<nbody; j++)
                 {
                       j1=j*dim;
                       j2=nd+j1;
                       d3=0.;
                       for (id=0; id<dim; id++) d3+=POW(Udouble[i1+id]-Udouble[j1+id],2);
                       d3=POW(d3,3./2);

                       for (id=0; id<dim; id++) 
                       {
                             qij=(Udouble[i1+id]-Udouble[j1+id])/d3;
                             fdouble[i2+id]-=Gmdouble[j]*qij;
                             fdouble[j2+id]+=Gmdouble[i]*qij;
                       }   
                 }
           }
         
     break;
                                     
     default:     
        
           for (i=0; i<nbody; i++)
           {
                 i1=i*dim;
                 i2=nd+i1;
                 for (id=0; id<dim; id++) 
                 {
                       fdouble[i1+id]=Udouble[i2+id];
                       fdouble[i2+id]=0.;
                 }
           }

           for (i=0; i<nbody; i++)
           {
                 i1=i*dim;
                 i2=nd+i1;
                 for (j=i+1; j<nbody; j++)
                 {
                       j1=j*dim;
                       j2=nd+j1;
                       d3=0.;
                       for (id=0; id<dim; id++) d3+=POW(Udouble[i1+id]-Udouble[j1+id],2);
                       d3=POW(d3,3./2);

                       for (id=0; id<dim; id++) 
                       {
                             qij=(Udouble[i1+id]-Udouble[j1+id])/d3;
                             fdouble[i2+id]-=Gmdouble[j]*qij;
                             fdouble[j2+id]+=Gmdouble[i]*qij;
                       }   
                 }

           }
 
     break;

}

    for (i=0; i<neq; i++) f[i]=fdouble[i];

    return ;

}




