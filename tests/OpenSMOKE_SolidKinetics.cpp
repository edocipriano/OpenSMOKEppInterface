#include <iostream>
#include "OpenSMOKE_Interface.h"

#define kinfolder "../kinetics/biomass/solid-only-2003/kinetics"

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

  OpenSMOKE_Clean();
  return 0;
}
