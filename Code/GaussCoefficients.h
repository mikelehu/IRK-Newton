#include <def.h>
#include <math.h>
#include <stdio.h> 
#include <stdlib.h>
#include <string.h>   
#include <cblas.h>


void GaussCoefficients
(char *path,gauss_method *method,val_type h
);

void GaussCoefficientsNewton   
(char *path,ode_sys *system,gauss_method *method,
 val_type h
);
