#include <stdio.h>
#include <math.h>
#include "OpenSMOKE_Interface.h"

#define D0 1.e-3
#define NT 4

double sq (double x) {
  return x*x;
}

//double rho1 = 500., rho2 = 25.;
//double rho1 = 2500., rho2 = 25.;
//double rho1 = 5500., rho2 = 25.;
//double rho1 = 2500., rho2 = 25.;
double rho1 = 50., rho2 = 1.;
double Dmix2 = 1.e-1;
double BY = 0.5;

void d2law (const double * y, const double dt, double * dy, void * args) {
  double D = y[0];
  double L0 = NT*(D0*0.5);
  double Ldiag = L0*sq(2);
  //double De = 2.*Ldiag;
  double De = 4.*D0;
  //dy[0] = -8.*rho2*Dmix2/rho1/sq(D)*1./(2./D - 2./De)*log( (1. - 0.) / (1. - BY) );
  dy[0] = -8.*rho2*Dmix2/rho1*log (2) / log (De/sqrt(D));
}

int main (void) {

  OpenSMOKE_InitODESolver ();

  int neq = 1;
  double y0[] = {sq(D0)};
  //double y0[] = {D0};
  double dt = 0.0000001;
  double tadend = 60.;
  double tend = tadend*sq (D0) / Dmix2;

  FILE * fp = fopen ("profiles", "w");

  for (double t=0; t<=tend; t += dt) {
    //fprintf (fp, "%.16f %.16f %.16f %.16f\n", t, t*Dmix2/sq(D0), sq(y0[0]/D0), y0[0]/D0 );
    fprintf (fp, "%.16f %.16f %.16f %.16f\n", t, t*Dmix2/sq(D0), y0[0]/sq(D0), sqrt(y0[0])/D0 );
    fflush (fp);
    OpenSMOKE_ODESolver (&d2law, neq, dt, y0, NULL);

    //if (y0[0] < 0.3*D0)
    if (sqrt(y0[0]) < 0.3*D0)
      break;
  }
  fclose (fp);

  OpenSMOKE_CleanODESolver ();
  return 0;
}

