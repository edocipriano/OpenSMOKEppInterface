#include <iostream>
#include <vector>
#include <math.h>
#include "OpenSMOKE_Interface.h"

#define kinfolder "../kinetics/surface/char/kinetics"
#define rho 1
#define R_J_mol_K 8.314

double Rexact (double T, double P, double xO2, double xCO2) {
  double A1 = 0.095;  // kg/m2/s/Pa
  double E1 = 108000; // J/mol

  double A2 = 7.55;   // kg/m2/s/Pa^0.45
  double E2 = 148500; // J/mol

  double P1 = P*xO2;
  double P2 = P*xCO2;

  double R1 = A1*pow (P1, 1.00)*exp (-E1/R_J_mol_K/T);
  double R2 = A2*pow (P2, 0.45)*exp (-E2/R_J_mol_K/T);

  return -(R1 + R2)/12.;
}

int main () {
  OpenSMOKE_Init();
  OpenSMOKE_ReadKinetics (kinfolder);
  OpenSMOKE_ReadSurfaceKinetics (kinfolder);

  unsigned int NGS = OpenSMOKE_NumberOfSpecies();
  unsigned int NSS = OpenSMOKE_NumberOfSurfaceSpecies();
  unsigned int NSR = OpenSMOKE_NumberOfSurfaceReactions();
  unsigned int NSP = 1;

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
  // gamma = number of surface phases
  std::vector<double> c(NGS), z(NSS), a(0), gamma(NSP);
  std::vector<double> R(NGS), Rgas(NGS), Rsite(NSS), Rbulk(0), RsitePhases(NSS);

  // Concentration of the species in gas phase
  double ctot = Pressure / R_J_mol_K / Temperature * 1e-3;
  std::vector<double> x(NGS);
  for (int i = 0; i < NGS; i++)
    x[i] = 1./(double)NGS;

  for (int i = 0; i < NGS; i++)
    c[i] = ctot*x[i];

  // Mole fraction of the surface species. The surface/char scheme assume a pure
  // solid surface, therefore z = 1.
  for (int i = 0; i < NSS; i++)
    z[i] = 1;

  // Site density (material property)
  // For the surface/char scheme, which does not depend on the Gamma value, this
  // parameter can be set to any number different from zero
  for (int i = 0; i < NSP; i++)
    gamma[i] = 1;

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

  std::cerr<< std::endl
           << "Total carbon formation rate from Hassan et al. = "
           << Rexact (Temperature, Pressure, x[OpenSMOKE_IndexOfSpecies ("O2")],
               x[OpenSMOKE_IndexOfSpecies ("CO2")])
           << " kmol/m2/s"
           << std::endl;

  std::cerr<< "Total carbon formation rate from OpenSMOKE++ = "
           << Rsite[0]
           << " kmol/m2/s"
           << std::endl;

  OpenSMOKE_Clean();
  return 0;
}

