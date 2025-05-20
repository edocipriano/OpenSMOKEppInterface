#include <iostream>
#include <vector>
#include "OpenSMOKE_Interface.h"

#define kinfolder "../kinetics/biomass/Solid-only-2003/kinetics"

int main () {
  OpenSMOKE_Init();
  OpenSMOKE_ReadKinetics (kinfolder);

  std::cerr<< "Number of gas phase species = "
    << OpenSMOKE_NumberOfSpecies() << std::endl
    << "Number of gas phase reactions = "
    << OpenSMOKE_NumberOfReactions() << std::endl;

  OpenSMOKE_ReadSolidKinetics (kinfolder);

  std::cerr<< "Number of solid phase species = "
    << OpenSMOKE_NumberOfSolidSpecies() << std::endl
    << "Number of solid phase reactions = "
    << OpenSMOKE_NumberOfSolidReactions() << std::endl;

  for (int jj=0; jj<OpenSMOKE_NumberOfSolidSpecies(); jj++) {
    std::cerr<< "Species[" << jj << "] = " << OpenSMOKE_NamesOfSolidSpecies (jj) << std::endl
      << "MW[" << jj << "] = " << OpenSMOKE_MW_Solid (jj) << std::endl
      << "Index[" << OpenSMOKE_NamesOfSolidSpecies (jj) << "] = " <<
      OpenSMOKE_IndexOfSolidSpecies (OpenSMOKE_NamesOfSolidSpecies (jj)) << std::endl;
  }

  std::vector<double> x (OpenSMOKE_NumberOfSolidSpecies());

  x[OpenSMOKE_IndexOfSolidSpecies ("CELL")] = 0.5;
  x[OpenSMOKE_IndexOfSolidSpecies ("LIGO")] = 0.5;

  OpenSMOKE_SolProp_SetTemperature(300.);
  OpenSMOKE_SolProp_SetPressure(101325.);
  double cp = OpenSMOKE_SolProp_HeatCapacity (x.data());
  std::cerr<< "Heat capacity = " << cp << std::endl;

  OpenSMOKE_Clean();
  return 0;
}
