
#include <stdio.h>
#include <string.h>
#include "OpenSMOKE_Interface.h"

int main (void) {
  OpenSMOKE_Init();
  OpenSMOKE_ReadKinetics();
  OpenSMOKE_ReadLiquidKinetics();
  OpenSMOKE_ReadLiquidProperties();

  printf ("NumberOfSpecies   = %d\n", OpenSMOKE_NumberOfSpecies());
  printf ("NumberOfReactions = %d\n", OpenSMOKE_NumberOfReactions());
  printf ("PiGreco           = %f\n", OpenSMOKE_Printpi());

  for (int i=0; i<OpenSMOKE_NumberOfSpecies(); i++) {
    double MWi = OpenSMOKE_MW(i);
    printf ("MWi = %f\n", MWi);
  }

  for (int i=0; i<OpenSMOKE_NumberOfSpecies(); i++) {
    printf ("Species = %s\n", OpenSMOKE_NamesOfSpecies(i));
  }

  for (int i=0; i<OpenSMOKE_NumberOfSpecies(); i++) {
    printf ("IndexOfSpecies %s = %d\n", OpenSMOKE_NamesOfSpecies(i), OpenSMOKE_IndexOfSpecies(OpenSMOKE_NamesOfSpecies(i)));
  }

  double xgas[OpenSMOKE_NumberOfSpecies()];
  for (int i=0; i<OpenSMOKE_NumberOfSpecies(); i++)
    xgas[i] = 0.;
  xgas[OpenSMOKE_IndexOfSpecies("O2")] = 0.21;
  xgas[OpenSMOKE_IndexOfSpecies("N2")] = 0.79;
  for (int i=0; i<OpenSMOKE_NumberOfSpecies(); i++)
    printf ("x[%d] = %f\n", i, xgas[i]);

  OpenSMOKE_GasProp_SetTemperature(298.);
  OpenSMOKE_GasProp_SetPressure(101325.);
  printf ("Air Dynamic Viscosity = %f\n", OpenSMOKE_GasProp_DynamicViscosity(xgas));
  printf ("Air thermal conductivity = %f\n", OpenSMOKE_GasProp_ThermalConductivity(xgas));

  printf ("Air MolecularWeight_From_MoleFractions = %f\n", OpenSMOKE_MolecularWeight_From_MoleFractions(xgas));
  printf ("Air MolecularWeight_From_MassFractions = %f\n", OpenSMOKE_MolecularWeight_From_MassFractions(xgas));

  double wgas[OpenSMOKE_NumberOfSpecies()];
  double MWGasMix;
  OpenSMOKE_MassFractions_From_MoleFractions (wgas, &MWGasMix, xgas);

  printf ("Air mix molecular weight = %f\n", MWGasMix);
  for (int i=0; i<OpenSMOKE_NumberOfSpecies(); i++) {
    printf ("wgas[%d] = %f\n", i, wgas[i]);
  }

  OpenSMOKE_MoleFractions_From_MassFractions (xgas, &MWGasMix, wgas);
  printf ("Air mix molecular weight = %f\n", MWGasMix);
  for (int i=0; i<OpenSMOKE_NumberOfSpecies(); i++) {
    printf ("xgas[%d] = %f\n", i, xgas[i]);
  }

  printf ("Air density (298K) = %f\n", OpenSMOKE_GasProp_Density_IdealGas(298., 101325., MWGasMix));
  printf ("Air density (773K) = %f\n", OpenSMOKE_GasProp_Density_IdealGas(773., 101325., MWGasMix));
  printf ("CpMix              = %f\n", OpenSMOKE_GasProp_cpMolar_Mixture_From_MoleFractions(xgas));

  return 0;
}
