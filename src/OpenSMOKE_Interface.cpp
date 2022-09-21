
#include "OpenSMOKE_Interface.h"
#include "OpenSMOKEpp"
#include "math/PhysicalConstants.h"
#include "maps/Maps_CHEMKIN"
#include "maps/ThermodynamicsMap_Liquid_CHEMKIN.h"
#include "maps/KineticsMap_Liquid_CHEMKIN.h"
#include "thermodynamics/mixture/mixtureL/mixtureL.h"
#include "thermodynamics/mixture/mixtureG/mixtureG.h"
#include "dictionary/OpenSMOKE_Dictionary"

OpenSMOKE::ThermodynamicsMap_CHEMKIN*           thermodynamicsMapXML;
OpenSMOKE::KineticsMap_CHEMKIN*                 kineticsMapXML;
OpenSMOKE::TransportPropertiesMap_CHEMKIN*      transportMapXML;

OpenSMOKE::ThermodynamicsMap_Liquid_CHEMKIN* thermodynamicsLiquidMapXML;
OpenSMOKE::KineticsMap_Liquid_CHEMKIN*       kineticsLiquidMapXML;

speciesMap* species_map;

struct {
  double T, P;
} gasdata;

#ifdef __cplusplus
extern "C" {
#endif

void OpenSMOKE_Init (void) {

  thermodynamicsMapXML = NULL;
  kineticsMapXML = NULL;
  transportMapXML = NULL;
  thermodynamicsLiquidMapXML = NULL;
  kineticsLiquidMapXML = NULL;
  species_map = NULL;
}

void OpenSMOKE_ReadKinetics (void) {

  boost::filesystem::path path_kinetics = "kinetics";
  boost::property_tree::ptree ptree;
  boost::property_tree::read_xml( (path_kinetics / "kinetics.xml").c_str(), ptree );

  double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

  thermodynamicsMapXML = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(ptree);
  transportMapXML = new OpenSMOKE::TransportPropertiesMap_CHEMKIN(ptree);
  kineticsMapXML = new OpenSMOKE::KineticsMap_CHEMKIN(*thermodynamicsMapXML,ptree);

  double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
  std::cout<< "Time to read XML file: " << tEnd - tStart << endl; 
}

void OpenSMOKE_ReadLiquidKinetics (void) {

  boost::filesystem::path path_kinetics = "kinetics";
  boost::property_tree::ptree ptree;
  boost::property_tree::read_xml( (path_kinetics / "kinetics.liquid.xml").string(), ptree );

  double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

  thermodynamicsLiquidMapXML = new OpenSMOKE::ThermodynamicsMap_Liquid_CHEMKIN(ptree);
  kineticsLiquidMapXML = new OpenSMOKE::KineticsMap_Liquid_CHEMKIN(*thermodynamicsLiquidMapXML, ptree, 1); 

  double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
  std::cout<< "Time to read XML file: " << tEnd - tStart << endl; 

}

void OpenSMOKE_ReadLiquidProperties (void) {
  species_map = new speciesMap("LiquidProperties");
}

int OpenSMOKE_NumberOfSpecies (void) {
  return thermodynamicsMapXML->NumberOfSpecies();
}

int OpenSMOKE_NumberOfReactions (void) {
  return kineticsMapXML->NumberOfReactions();
}

double OpenSMOKE_Printpi (void) {
  return PhysicalConstants::pi;
}

void OpenSMOKE_GasProp_SetPressure (const double P) {
  gasdata.P = P;
  thermodynamicsMapXML->SetPressure(P);
  transportMapXML->SetPressure(P);
  kineticsMapXML->SetPressure(P);
}

void OpenSMOKE_GasProp_SetTemperature (const double T) {
  gasdata.T = T;
  thermodynamicsMapXML->SetTemperature(T);
  transportMapXML->SetTemperature(T);
  kineticsMapXML->SetTemperature(T);
}

double OpenSMOKE_MW (const int i) {
  return thermodynamicsMapXML->MW(i);
}

const char* OpenSMOKE_NamesOfSpecies (const int i) {
  return thermodynamicsMapXML->NamesOfSpecies()[i].c_str();
}

int OpenSMOKE_IndexOfSpecies (const char* s) {
  return thermodynamicsMapXML->IndexOfSpecies(s) - 1;
}

double OpenSMOKE_GasProp_DynamicViscosity (double* x) {
  return transportMapXML->DynamicViscosity(x);
}

double OpenSMOKE_GasProp_ThermalConductivity (double* x) {
  return transportMapXML->ThermalConductivity(x);
}

double OpenSMOKE_GasProp_cpMolar_Mixture_From_MoleFractions (const double* x) {
  double MW = OpenSMOKE_MolecularWeight_From_MoleFractions (x);
  double Cp = thermodynamicsMapXML->cpMolar_Mixture_From_MoleFractions(x);
  return Cp/MW;
}

double OpenSMOKE_GasProp_Concentrations (double T, double P) {
  return P/PhysicalConstants::R_J_kmol/T;
}

double OpenSMOKE_GasProp_Density_IdealGas (double T, double P, double MW) {
  return OpenSMOKE_GasProp_Concentrations(T,P)*MW;
}

double OpenSMOKE_MolecularWeight_From_MoleFractions (const double* x) {
  return thermodynamicsMapXML->MolecularWeight_From_MoleFractions(x);
}

double OpenSMOKE_MolecularWeight_From_MassFractions (const double* x) {
  return thermodynamicsMapXML->MolecularWeight_From_MassFractions(x);
}

void OpenSMOKE_MassFractions_From_MoleFractions (double* y, double* MW, const double* x) {
  return thermodynamicsMapXML->MassFractions_From_MoleFractions(y,*MW,x);
}

void OpenSMOKE_MoleFractions_From_MassFractions (double* x, double* MW, const double* y) {
  return thermodynamicsMapXML->MoleFractions_From_MassFractions(x,*MW,y);
}

#ifdef __cplusplus
}
#endif

