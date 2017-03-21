#include <def.h>

void test_jac
(int neq, val_type t,val_type *u,void (*Ode)(), 
 val_type *Jac,parameters *params
);
void Ode1 (int neq, val_type t,val_type *u,val_type *f,parameters *params);
val_type Ham1 (int neq,solution *u,parameters *params);
void Jac1 (int neq, val_type t,val_type *u,val_type *Jac,parameters *params);
void Ode2 (int neq, val_type t,val_type *u,val_type *f,parameters *params);
val_type Ham2 (int neq,solution *u,parameters *params);
void Jac2 (int neq, val_type t,val_type *u,val_type *Jac,parameters *params);
void Ode3 (int neq, val_type t,val_type *u,val_type *f,parameters *params);
val_type Ham3 (int neq,solution *u,parameters *params);
void Jac3 (int neq, val_type t,val_type *u,val_type *Jac,parameters *params);
void Ode4 (int neq, val_type t,val_type *u,val_type *f,parameters *params);
val_type Ham4 (int neq,solution *u,parameters *params);
void Jac4 (int neq, val_type t,val_type *u,val_type *Jac,parameters *params);
void Ode5 (int neq, val_type t,val_type *u,val_type *f,parameters *params);
val_type Ham5 (int neq,solution *u,parameters *params);
void Jac5 (int neq, val_type t,val_type *u,val_type *Jac,parameters *params);
void Ode6 (int neq, val_type t,val_type *u,val_type *f,parameters *params);
val_type Ham6 (int neq,solution *u,parameters *params);
void Jac6 (int neq, val_type t,val_type *u,val_type *Jac,parameters *params);
void Ode7 (int neq, val_type t,val_type *u,val_type *f,parameters *params);
val_type Ham7 (int neq,solution *u,parameters *params);
void Jac7 (int neq, val_type t,val_type *u,val_type *Jac,parameters *params);
void Ode8 (int neq, val_type t,val_type *u,val_type *f,parameters *params);
val_type Ham8 (int neq,solution *u,parameters *params);
void Jac8 (int neq, val_type t,val_type *u,val_type *Jac,parameters *params);
void Ode9 (int neq, val_type t,val_type *u,val_type *f,parameters *params);
val_type Ham9 (int neq,solution *u,parameters *params);
void Jac9 (int neq, val_type t,val_type *u,val_type *Jac,parameters *params);
void Ode10 (int neq, val_type t,val_type *u,val_type *f,parameters *params);
val_type Ham10 (int neq,solution *u,parameters *params);
void Jac10 (int neq, val_type t,val_type *u,val_type *Jac,parameters *params);

