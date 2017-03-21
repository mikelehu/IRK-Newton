#include <def.h>
#include <math.h>
#include <stdio.h> 
#include <stdlib.h>
#include <cblas.h>
//#include <string.h>   


void GaussCoefficients
(char *path,gauss_method *method,toptions *options
);

void GaussCoefficientsNewton   
(char *path,ode_sys *system,gauss_method *method,
 toptions *options
);
