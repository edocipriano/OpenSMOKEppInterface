#include <iostream>
#include "OpenSMOKE_Interface.h"

#define kinfolder "../kinetics/skeletal/methanol/kinetics"

int main () {
  OpenSMOKE_Init();
  OpenSMOKE_ReadKinetics (kinfolder);

  std::cerr<< "Number of gas phase species = "
    << OpenSMOKE_NumberOfSpecies() << std::endl
    << "Number of gas phase reactions = "
    << OpenSMOKE_NumberOfReactions() << std::endl;

  for (int jj=0; jj<OpenSMOKE_NumberOfSpecies(); jj++) {
    std::cerr<< "Species[" << jj << "] = " << OpenSMOKE_NamesOfSpecies (jj) << std::endl
      << "MW[" << jj << "] = " << OpenSMOKE_MW (jj) << std::endl
      << "Index[" << OpenSMOKE_NamesOfSpecies (jj) << "] = " <<
      OpenSMOKE_IndexOfSpecies (OpenSMOKE_NamesOfSpecies (jj)) << std::endl;
  }

  OpenSMOKE_Clean();
  return 0;
}
