#include <stdio.h>
#include <math.h>
#include "OpenSMOKE_Interface.h"

double sq (double x) {
  return x*x;
}

double D0 = 1e-3;
double L0 = 10.*1e-3*0.3;
double rho1 = 250., rho2 = 25.;
double Dmix2 = 1.e-3;
double BY = 0.5;

void d2law (const double * y, const double dt, double * dy, void * args) {

  dy[0] = -8.*rho2*Dmix2/rho1*log (1. + BY)/log (L0/sqrt (y[0]));
}

int main (void) {

  OpenSMOKE_InitODESolver ();

  int neq = 1;
  double y0[] = {D0*D0};
  double dt = 0.0001;
  double tadend = 5.;
  double tend = tadend*sq (D0) / Dmix2;

  FILE * fp = fopen ("profiles", "w");

  for (double t=0; t<=tend; t += dt) {
    fprintf (fp, "%.16f %.16f %.16f %.16f\n", t, t*Dmix2/sq(D0), y0[0]/sq (D0), sqrt (y0[0])/D0 );
    fflush (fp);
    OpenSMOKE_ODESolver (&d2law, neq, dt, y0, NULL);
  }
  fclose (fp);

  OpenSMOKE_CleanODESolver ();
  return 0;
}


