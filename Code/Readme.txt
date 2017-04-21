****************************************************************************
Readme.txt
****************************************************************************
     Gauss Implicit Runge-Kutta implementation. 
     Newton iteration.

     version: 1.1 (20-04-2017).

     Article: "Efficient implementation of symplectic implicit
               Runge-Kutta schemes with simplified Newton iterations" (2017)


********************************************************************************
CONTENTS:

   def.c:               Parameters and general definitions we use in the code.
                        You must specify math-functions of your Odefun. 

   Terminal-IRKNEWTON.c:     An example to show how to call the numerical integration.
   math-IRKNEWTON.c:	     An auxiliar file to call from mathematica (double precision). 

   GaussCoefficients.c:      Coefficients define Inplicit Runge-Kutta method.

   Common-IRKNEWTON.c:       Numerical integration method. We will find these functions:

	IRKNEWTON () : 	      IRK Newton integration (main function).
	NSS_Step():           Newton Simplified step.	
        NSS_Step_plus();      Newton Simplified step (efficient).
        NSS_MIX_Step();       Newton Simplified-Mix step.

        NSS_Solve():          Compute one iteration of Newton Simplified.
        NSS_Solve_plus():     Compute one iteration of Newton Simplified (efficient).

        Compute_R():	      Function for Newton Simplified efficient. 
        Compute_Z():	      Function for Newton Simplified efficient. 
        Compute_W1():	      Function for Newton Simplified efficient. 
        Compute_W2():	      Function for Newton Simplified efficient. 
        Compute_DL():	      Function for Newton Simplified efficient. 
        Compute_GG():	      Function for Newton Simplified efficient. 

        Default_Stage_init:      Initialization of Stages, L_{n,i}=0	
        Interpolated_Stage_init	 Initialization of stages, L_{n,i}=G(L_{n-1,i})		

        Yi_update():          Y_i= y + (e+ sum_j^s m_{ij} L_j) 
        Yi_update_Classic():  Y_i= y + sum_j^s m_{ij} L_j 
											
	MyOutput():	      An output function defined by the user.

	NormalizedDistance()     : Check of convergence of the iteration.
	StopCriterion ()         : Stopping criterion (double precision).
	StopCriterionFloat ():   : Stopping criterion (float precision).
	RemoveDigitsFcn()        : Rounding a floating point number with p-r significant binary digits.
							
   Problems.c:  Double pendulum Stiff and N-Body ode system:
														
	OdeDPS()=Ode system of Double Pendulum Stiff.
	HamDPS()=Hamiltonian of Double Pendulum Stiff.
        JacDPS()=Jacobian of Double Pendulum Stiff.	
	Ode2()=Ode system of NBody.
	Ham2()=Hamiltonian of NBody.
        Jac2()=Jacobian of NBody.

   math-IRKNEWTON.tm:	Mathlink file.
			

*********************************************************************************
USE OF THE CODE:

   The Function call to integrate an ODE system is:

     int thread_count = 1;   // the user can set the  number of threads for parallel computations. 
     IRKNEWTON (t0,t1,h, &gsmethod, &u, &system, &options, &thestat);

       t0,t1:     interval of numerical integration.
       h:         step size.
       method:    structure with information about the method to be used in the integration: 
                  Butcher tableau and the reformulation of the symplectic IRK scheme (see below)
       u:         initial values of the problem (structure with two parts: uu for initial values 
                  and ee for accumulated round-off error values) 
       system:    structure with the information of the problem to be solved (function that evaluates 
                  the differential equations, the Hamiltonian, list of parameters of the problem...)
       options:   options of the integration: tolerances, function to be called after each step, 
                  name of the file where save intermediate values, initialization of the stages at each step...
       thestat:   internal data of the integration: number of iterations, number of steps, number of 
                  calls to the ODE function...
      
      
    All these structures are defined in the file "def.h". Let's see them:
       method:    It defines the scheme used to solve the equations. Its members are:
                      int ns;        the number of stages 
                      val_type *c,*b,*a;	The Butcher tableau of the method.
                      val_type *m;   mij=aij/bj coefficients derived from the reformulation to maintain symplecticity mij+mji-1=0.
                      val_type *hc;  hc=h*c (to avoid recalculations).
                      val_type *hb;  hb=h*b (to avoid recalculations).
                      val_type *munu;  coefficients for interpolation.
                      int *orderedindices;   ascending order indices for bi coefficients                  
                  The function defined in the file GaussCoefficients.c loads all these values for 
                  Gauss collocation schemes. The informations needed by this function 
                  are the number of stages and the step size h.
       system:    This structure has the information of the ODE: 
                      int neq:       Number of equations or dimension of the problem. 
                      void (*f)():   the function that evaluates the differential equations, The user has to code it and set here its name.
                      val_type (*Ham)(): the function that evaluates the Hamiltonian (if the problem is Hamiltonian),
                      val_type (*Jac)(): the function that evaluates the Jacobian, 
                      parameters:    structure with the parameters that define the problem (list of real parameters and list of integer parameters),
                      int cod[2]:    if the problem is a partitioned problem there is the possibility to define which part must be evaluated first.
       options:   The user can set several options for the integration process:
                      val_type *rtol,*atol;    relative and absolute tolerances	
                      void (*TheOutput)();     The function that will be called after each step (or after "sampling" steps).	  
                      int sampling;            For long time integrations TheOutput function will be called once after "sampling" steps 
                      int rdigits,mrdigits;    When we want a computation with r digits less of accuracy we can set rdigits = r
                      IRKNEWTON_Step_Fn:       Newton method; NSS_Step, NSS_Step_plus, NSS_MIX_Step.
                      char filename[STRMAX];   Output filename.
                      void (*StageInitFn)();   The function that initializes the stages at the beginning of the step. We offer two functions, but the user can set its own function. 


*********************************************************************************
OUTPUT:


       the value at t1 is returned in the variable u.
       If the user uses MyOutput function then the file set in options.filename will contain the values at each sampling time.

*********************************************************************************
PARAMETERS (file: def.h):

   You can specify next parameters:

   PARALLEL: we want parallel execution.

   MAXIT :     maximum number of Newton-like iterations 
   RTOL,ATOL:  Newton iteration tolerance. 


   #define DIR_TERM :  // Path Coefficients for terminal executions  
   #define DIR_MATH :  // Path Coefficients for mathematica executions


********************************************************************************

INSTALATION (Ubuntu 16.04 lts):

   Two ways to execute :

     Terminal execution:
          make term-IRKNEWTON
          ./Terminal-IRKNEWTON.exe

      Mathematica (double)
          make math-IRKNEWTON
          Execution from mathematica notebook (see examples).

   The following libraries must be installed: CBLAS, LAPACKE.

********************************************************************************
EXAMPLES (diretories):

    NCDP-STIFF (Many K)  : Mathematica notebooks for the integration of different k values.
    NCDP-STIFF (K=0)     : Mathematica notebooks for the integration of "nstat" perturbed 
                          initial values of k=0 Double Pendulum Stiff problem.
    NCDP-STIFF (K=(32)^2 : Mathematica notebooks for the integration of "nstat" perturbed 
                          initial values of k=(32)^2 Double Pendulum Stiff problem.
    NCDP-STIFF (K=(64)^2 : Mathematica notebooks for the integration of "nstat" perturbed 
                          initial values of k=(64)^2 Double Pendulum Stiff problem.   


    
    



