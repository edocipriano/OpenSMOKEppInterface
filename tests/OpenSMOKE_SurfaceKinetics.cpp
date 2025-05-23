#include <iostream>
#include "OpenSMOKE_Interface.h"

#define kinfolder "../kinetics/surface/char/kinetics"
#define rho 1

int main () {
  OpenSMOKE_Init();
  OpenSMOKE_ReadKinetics (kinfolder);
  OpenSMOKE_ReadSurfaceKinetics (kinfolder);

  unsigned int NGS = OpenSMOKE_NumberOfSpecies();
  unsigned int NSS = OpenSMOKE_NumberOfSurfaceSpecies();
  unsigned int NSR = OpenSMOKE_NumberOfSurfaceReactions();

  std::cerr<< "Number of gas species = " << NGS << std::endl
           << "Number of surface species = " << NSS << std::endl
           << "Number of surface reactions = " << NSR << std::endl;

  for (int i = 0; i < NGS; i++)
    std::cerr << "Species " << i << " is " << OpenSMOKE_NamesOfSpecies(i) << std::endl;

  double Temperature = 1000.;
  double Pressure = 101325.;

  // c = number of gas species
  // z = number of site species
  // a = number of bulk species
  // gamma = number of site species
  std::vector<double> c(NGS), z(NSS), a(0), gamma(NSS);
  std::vector<double> R(NGS), Rgas(NGS), Rsite(NSS), Rbulk(0), RsitePhases(NSS);

  for (int i = 0; i < NGS; i++)
    c[i] = 1;

  for (int i = 0; i < NSS; i++)
    z[i] = 1;

  for (int i = 0; i < NSS; i++)
    gamma[i] = 1; // 0 set to zero all surface reactions

  OpenSMOKE_GasProp_SetTemperature (Temperature);
  OpenSMOKE_GasProp_SetPressure (Pressure);
  OpenSMOKE_GasProp_ReactionRates (c.data());
  OpenSMOKE_GasProp_FormationRates (R.data());

  for (unsigned int i = 0; i < NGS; i++)
    std::cerr<< "R species " << i << ": " << R[i] << std::endl;

  double Qh = OpenSMOKE_GasProp_HeatRelease (R.data());
  std::cerr<< "Heat released by the gas phase chemical reactions: " << Qh << std::endl;

  OpenSMOKE_SurProp_SetTemperature (Temperature);
  OpenSMOKE_SurProp_SetPressure (Pressure);
  OpenSMOKE_SurProp_ReactionRates (c.data(), z.data(), a.data(), gamma.data());
  OpenSMOKE_SurProp_KineticConstants();
  //OpenSMOKE_SurProp_FormationRatesGasOnly (Rgas.data());
  OpenSMOKE_SurProp_FormationRates (Rgas.data(), Rsite.data(), Rbulk.data(), RsitePhases.data());

  for (unsigned int i = 0; i < NGS; i++)
    std::cerr<< "Rgas species " << i << ": " << Rgas[i] << std::endl;

  for (unsigned int i = 0; i < NSS; i++)
    std::cerr<< "Rsite species " << i << ": " << Rsite[i] << std::endl;

  for (unsigned int i = 0; i < 0; i++)
    std::cerr<< "Rbulk species " << i << ": " << Rbulk[i] << std::endl;

  for (unsigned int i = 0; i < NSS; i++)
    std::cerr<< "RsitePhases species " << i << ": " << RsitePhases[i] << std::endl;

  double Qr = OpenSMOKE_SurProp_HeatRelease (Rgas.data(), Rsite.data(), Rbulk.data());
  std::cerr<< "Heat released by the surface chemical reactions: " << Qr << std::endl;

  OpenSMOKE_Clean();
  return 0;
}
