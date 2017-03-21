****************************************************************************
Readme.txt
****************************************************************************
     Gauss Implicit Runge-Kutta implementation. 
     Newton iteration.

     version: 1 (20-03-2017).

     Article: "Efficient implementation of symplectic implicit
               Runge-Kutta schemes with simplified Newton iterations" (2017)


********************************************************************************
CONTENTS:

   prec.c:		Precision.
   def.c:               Parameters and general definitions we use in the code.
                        You must specify math-functions of your Odefun. 

   GaussTerminal.c:     An example to show how to call the numerical integration.
   math-Gauss.c:	An auxiliar file to call from mathematica (double precision). 

   GaussInitData.c:     Contain "InitialData" with the initial values for double pendulum stiff problem.

   GaussCoefficients.c: Coefficients mij,bi,ci that define Inplicit Runge-Kutta method.

   GaussCommon.c:       Numerical integration method. We will find these functions:

	RKG () : 	      IRK Newton iteration method.
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

	Li_init();	      Initialization function for Li stages.
        Yi_update():          Y_i= y + (e+ sum_j^s m_{ij} L_j) 
        Yi_update_Classic():  Y_i= y + sum_j^s m_{ij} L_j 
											
	TheOutput():	      Output function for RKG ().

	NormalizedDistance()     : check of convergence of the iteration.
	StopCriterion ()         : Stopping criterion (double precision).
	StopCriterionFloat ():   : Stopping criterion (float precision).
	RemoveDigitsFcn()        : Rounding a floating point number with p-r significant binary digits.
							
   GaussUserProblem.c:  Double pendulum Stiff and N-Body ode system:
														
	Ode1()=Ode system of Double Pendulum Stiff.
	Ham1()=Hamiltonian of Double Pendulum Stiff.
        Jac1()=Jacobian of Double Pendulum Stiff.	
	Ode2()=Ode system of NBody.
	Ham2()=Hamiltonian of NBody.
        Jac2()=Jacobian of NBody.

   math-Gauss.tm:	Mathlink file.
			

*********************************************************************************
OPTIONS:

   You will have to specify next options for the numerical integration:

   ns= number of stages of Inplicit Runge-Kutta method.
   eda= name of differential equation.

   algorith= You must specify one of the next implementations options.
	=  1 Newton Simplified method (LU decomposition of full sd x sd matrix). 
	=  2 Newton Simplified method (efficient).
	=  3 Newton Simplified-Mix method (efficient).
   
   h= stepsize.
   sampling: we sample the numerical results once every "sampling" steps.
   filename = output binary filename. 

   thread_count: number of threads for parallel computation. 

   note: neq, t0, tf and some others options of the problems are 
         initialized with initial values.

*********************************************************************************
PARAMETERS (file: def.h):

   You can specify next parameters:

   IOUT: default form (enable otuput binary file). 
   PARALLEL: we want parallel execution.

   MAXIT :     maximum number of Newton-like iterations 
   RTOL,ATOL:  Newton iteration tolerance. 
   PREC:
       1=DOUBLEPRECISION


*********************************************************************************
RESULTS:

   filename (output binary format).

********************************************************************************
INSTALATION (Ubuntu 16.04 lts):

   Two ways to execute :

     Terminal execution:
          make term-Gauss
          ./GaussTerminal.exe

      Mathematica (double)
          make math-Gauss
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


    
    



