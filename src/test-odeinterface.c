#include <stdio.h>
#include "OpenSMOKE_Interface.h"

void batchReactor (const double * y, const double dt, double * dy) {
  dy[0] = 0.5;
  dy[1] = 1.;
}

int main (void) {

  int neq = 2;
  double dt = 1.;
  double y0[] = {0., 0.};

  OpenSMOKE_ODESolver (&batchReactor, neq, dt, y0);

  printf ("---------------------------------------------\n");
  printf ("ODE system results:\n");
  for (int i=0; i<neq; i++) {
    printf ("y[%d] = %f\n", i, y0[i]);
  }
  printf ("---------------------------------------------\n");

  return 0;
}
