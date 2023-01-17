#include <stdio.h>
#include "OpenSMOKE_Interface.h"

typedef struct {
  double k1, k2, k3;
} UserData;

void richardson (const double * y, const double dt, double * dy, void * args) {

  UserData data = *(UserData *)args;
  double k1 = data.k1;
  double k2 = data.k2;
  double k3 = data.k3;

  //double k1, k2, k3;
  //k1 = 0.04, k2 = 1e7, k3 = 1e4;

  double y1 = y[0], y2 = y[1], y3 = y[2];

  dy[0] = -k1*y1 + k2*y2*y3;
  dy[1] =  k1*y1 - k2*y2*y3 - k3*y2*y2;
  dy[2] =  k3*y2*y2;
}

int main (void) {
  int neq = 3;
  double y0[] = {1., 0., 0.};
  ///double dt = 4.e+10;
  double dt = 1.;

  UserData data;
  data.k1 = 0.04;
  data.k2 = 1e4;
  data.k3 = 3e7;

  OpenSMOKE_ODESolver (&richardson, neq, dt, y0, &data);

  for (int i=0; i<neq; i++) {
    printf ("y0[%d] = %.16f\n", i, y0[i]);
  }

  return 0;
}


