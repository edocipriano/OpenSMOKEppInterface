// OpenSMOKE Interface
#include "OpenSMOKE_Interface.h"
#include "OpenSMOKE_ODESystem_Interface.h"

// OpenSMOKE++ Definitions
#include "OpenSMOKEpp"
#include "math/PhysicalConstants.h"

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"
#include "maps/ThermodynamicsMap_Liquid_CHEMKIN.h"
#include "maps/KineticsMap_Liquid_CHEMKIN.h"
#include "maps/ThermodynamicsMap_Solid_CHEMKIN.h"
#include "maps/KineticsMap_Solid_CHEMKIN.h"
#include "maps/ThermodynamicsMap_Surface_CHEMKIN.h"
#include "maps/KineticsMap_Surface_CHEMKIN.h"

// OpenSMOKE++ Thermodynamics and liquid phase properties
#include "thermodynamics/mixture/mixtureL/mixtureL.h"
#include "thermodynamics/mixture/mixtureG/mixtureG.h"

// OpenSMOKE++ Dictionaries
#include "dictionary/OpenSMOKE_Dictionary"

// ODE solvers
#include "math/native-ode-solvers/MultiValueSolver"
#include "math/external-ode-solvers/ODE_Parameters.h"

// Global pointers
OpenSMOKE::ThermodynamicsMap_CHEMKIN*           thermodynamicsMapXML;
OpenSMOKE::KineticsMap_CHEMKIN*                 kineticsMapXML;
OpenSMOKE::TransportPropertiesMap_CHEMKIN*      transportMapXML;

OpenSMOKE::ThermodynamicsMap_Liquid_CHEMKIN* thermodynamicsLiquidMapXML;
OpenSMOKE::KineticsMap_Liquid_CHEMKIN*       kineticsLiquidMapXML;

OpenSMOKE::ThermodynamicsMap_Solid_CHEMKIN* thermodynamicsSolidMapXML;
OpenSMOKE::KineticsMap_Solid_CHEMKIN*       kineticsSolidMapXML;

OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN* thermodynamicsSurfaceMapXML;
OpenSMOKE::KineticsMap_Surface_CHEMKIN* kineticsSurfaceMapXML;

speciesMap* species_map;
mixtureL* liqmix;

// Create ODE parameters class
OpenSMOKE::ODE_Parameters* ode_parameters_;

// Helper variables
std::vector<std::string> liqspecies;

#ifdef __cplusplus
extern "C" {
#endif

#define cleanptr(PTR) \
  if (PTR != NULL) {  \
    delete PTR;       \
    PTR = NULL;       \
  }

void OpenSMOKE_Init (void) {

  thermodynamicsMapXML = NULL;
  kineticsMapXML = NULL;
  transportMapXML = NULL;
  thermodynamicsLiquidMapXML = NULL;
  kineticsLiquidMapXML = NULL;
  thermodynamicsSolidMapXML = NULL;
  thermodynamicsSurfaceMapXML = NULL;
  kineticsSurfaceMapXML = NULL;
  kineticsSolidMapXML = NULL;
  species_map = NULL;
  liqmix = NULL;
}

void OpenSMOKE_InitODESolver (void) {

  //ode_solver = new OdeSMOKE::MultiValueSolver<methodGear>;
  ode_parameters_ = new OpenSMOKE::ODE_Parameters;
}

void OpenSMOKE_Clean (void) {
  cleanptr (thermodynamicsMapXML);
  cleanptr (kineticsMapXML);
  cleanptr (transportMapXML);
  cleanptr (thermodynamicsLiquidMapXML);
  cleanptr (kineticsLiquidMapXML);
  cleanptr (thermodynamicsSolidMapXML);
  cleanptr (kineticsSolidMapXML);
  cleanptr (thermodynamicsSurfaceMapXML);
  cleanptr (kineticsSurfaceMapXML);
  //cleanptr (species_map); // FIXME: strange behavior
  //cleanptr (liqmix);      // FIXME: strange behavior
}

//void OpenSMOKE_Clean (void) {
//  delete thermodynamicsMapXML;
//  delete kineticsMapXML;
//  delete transportMapXML;
//  delete species_map;
//  //delete thermodynamicsLiquidMapXML;
//  //delete kineticsLiquidMapXML;
//  thermodynamicsMapXML = NULL;
//  kineticsMapXML = NULL;
//  transportMapXML = NULL;
//  species_map = NULL;
//
//  //delete species_map;
//}

void OpenSMOKE_CleanODESolver (void) {
  //delete ode_solver;
  delete ode_parameters_;
  ode_parameters_ = NULL;
}

void OpenSMOKE_ReadKinetics (const char* kinfolder) {

  boost::filesystem::path path_kinetics = boost::filesystem::exists(kinfolder) ? kinfolder : "kinetics";
  boost::property_tree::ptree ptree;
  boost::property_tree::read_xml( (path_kinetics / "kinetics.xml").c_str(), ptree );

  double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

  thermodynamicsMapXML = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(ptree);
  transportMapXML = new OpenSMOKE::TransportPropertiesMap_CHEMKIN(ptree);
  kineticsMapXML = new OpenSMOKE::KineticsMap_CHEMKIN(*thermodynamicsMapXML,ptree);

  double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
  std::cout<< "Time to read XML file: " << tEnd - tStart << endl; 
}

void OpenSMOKE_ReadLiquidKinetics (const char* kinfolder) {

  boost::filesystem::path path_kinetics = boost::filesystem::exists(kinfolder) ? kinfolder : "kinetics";
  boost::property_tree::ptree ptree;
  if (boost::filesystem::exists (path_kinetics / "kinetics.liquid.xml")) {
    boost::property_tree::read_xml( (path_kinetics / "kinetics.liquid.xml").string(), ptree );

    double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

    thermodynamicsLiquidMapXML = new OpenSMOKE::ThermodynamicsMap_Liquid_CHEMKIN(ptree);
    kineticsLiquidMapXML = new OpenSMOKE::KineticsMap_Liquid_CHEMKIN(*thermodynamicsLiquidMapXML, ptree, 1); 

    double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
    std::cout<< "Time to read XML file: " << tEnd - tStart << endl; 
  }
  else
    std::cout<< "Unable to read kinetics.liquid.xml" << std::endl;
}

void OpenSMOKE_ReadSolidKinetics (const char* kinfolder) {

  boost::filesystem::path path_kinetics = boost::filesystem::exists(kinfolder) ? kinfolder : "kinetics";
  boost::property_tree::ptree ptree;
  if (boost::filesystem::exists (path_kinetics / "kinetics.solid.xml")) {
    boost::property_tree::read_xml( (path_kinetics / "kinetics.solid.xml").string(), ptree );

    double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

    thermodynamicsSolidMapXML = new OpenSMOKE::ThermodynamicsMap_Solid_CHEMKIN(ptree);
    kineticsSolidMapXML = new OpenSMOKE::KineticsMap_Solid_CHEMKIN(*thermodynamicsSolidMapXML, ptree, 1);
     
    double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
    std::cout<< "Time to read XML file: " << tEnd - tStart << endl; 
  }
  else
    std::cout<< "Unable to read kinetics.solid.xml" << std::endl;
}

void OpenSMOKE_ReadSurfaceKinetics (const char* kinfolder) {

  boost::filesystem::path path_kinetics = boost::filesystem::exists(kinfolder) ? kinfolder : "kinetics";
  boost::property_tree::ptree ptree;
  if (boost::filesystem::exists (path_kinetics / "kinetics.surface.xml")) {
    boost::property_tree::read_xml( (path_kinetics / "kinetics.surface.xml").string(), ptree );

    double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

    thermodynamicsSurfaceMapXML = new OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN(ptree);
    kineticsSurfaceMapXML = new OpenSMOKE::KineticsMap_Surface_CHEMKIN(*thermodynamicsSurfaceMapXML, ptree);

    double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
    std::cout<< "Time to read XML file: " << tEnd - tStart << endl;
  }
  else
    std::cout<< "Unable to read kinetics.surface.xml" << std::endl;
}

void OpenSMOKE_ReadLiquidProperties (const char* liqpropfolder) {
  std::string opensmoke_interface_root = std::getenv ("OPENSMOKE_INTERFACE");
  boost::filesystem::path path_liqprop
    = boost::filesystem::exists(liqpropfolder)
    ? liqpropfolder
    : opensmoke_interface_root.append ("/kinetics/LiquidProperties/LiquidProperties");
  species_map = new speciesMap (path_liqprop);

  // Setup liquid species
  liqspecies.resize (OpenSMOKE_NumberOfLiquidSpecies());
  for (int jj=0; jj<liqspecies.size(); jj++) {
    liqspecies[jj] = OpenSMOKE_NamesOfLiquidSpecies(jj);
    liqspecies[jj] = liqspecies[jj].substr(0, liqspecies[jj].size()-3);
  }

  // Create liquid mixture
  liqmix = new mixtureL (liqspecies, *species_map);
}

int OpenSMOKE_NumberOfSpecies (void) {
  return thermodynamicsMapXML->NumberOfSpecies();
}

int OpenSMOKE_NumberOfLiquidSpecies (void) {
  return thermodynamicsLiquidMapXML->number_of_liquid_species();
}

int OpenSMOKE_NumberOfSolidSpecies (void) {
  return thermodynamicsSolidMapXML->number_of_solid_species();
}

int OpenSMOKE_NumberOfSurfaceSpecies (void) {
  return thermodynamicsSurfaceMapXML->number_of_site_species();
}

int OpenSMOKE_NumberOfReactions (void) {
  return kineticsMapXML->NumberOfReactions();
}

int OpenSMOKE_NumberOfSolidReactions (void) {
  return kineticsSolidMapXML->NumberOfReactions();
}

int OpenSMOKE_NumberOfSurfaceReactions (void) {
  return kineticsSurfaceMapXML->NumberOfReactions();
}

double OpenSMOKE_Printpi (void) {
  return PhysicalConstants::pi;
}

void OpenSMOKE_GasProp_SetPressure (const double P) {
  thermodynamicsMapXML->SetPressure(P);
  transportMapXML->SetPressure(P);
  kineticsMapXML->SetPressure(P);
}

void OpenSMOKE_GasProp_SetTemperature (const double T) {
  thermodynamicsMapXML->SetTemperature(T);
  transportMapXML->SetTemperature(T);
  kineticsMapXML->SetTemperature(T);
}

void OpenSMOKE_LiqProp_SetPressure (const double P) {
  thermodynamicsLiquidMapXML->SetPressure(P);
  kineticsLiquidMapXML->SetPressure(P);
}

void OpenSMOKE_LiqProp_SetTemperature (const double T) {
  thermodynamicsLiquidMapXML->SetTemperature(T);
  kineticsLiquidMapXML->SetTemperature(T);
}

void OpenSMOKE_SurProp_SetTemperature (const double T) {
  thermodynamicsSurfaceMapXML->SetTemperature (T);
  kineticsSurfaceMapXML->SetTemperature (T);
}

void OpenSMOKE_SurProp_SetPressure (const double P) {
  thermodynamicsSurfaceMapXML->SetPressure (P);
  kineticsSurfaceMapXML->SetPressure (P);
}

void OpenSMOKE_SolProp_SetPressure (const double P) {
  thermodynamicsSolidMapXML->SetPressure(P);
  kineticsSolidMapXML->SetPressure(P);
}

void OpenSMOKE_SolProp_SetTemperature (const double T) {
  thermodynamicsSolidMapXML->SetTemperature(T);
  kineticsSolidMapXML->SetTemperature(T);
} 

void OpenSMOKE_GasProp_KineticConstants (void) {
  kineticsMapXML->KineticConstants();
}

void OpenSMOKE_GasProp_ReactionRates (const double * c) {
  kineticsMapXML->ReactionRates(c);
}

void OpenSMOKE_GasProp_FormationRates (double * R) {
  kineticsMapXML->FormationRates(R);
}

void OpenSMOKE_SolProp_KineticConstants (void) {
  kineticsSolidMapXML->KineticConstants();
}

void OpenSMOKE_SolProp_ReactionRates (const double * cgas, const double * csolid) {
  kineticsSolidMapXML->ReactionRates(cgas, csolid);
}

void OpenSMOKE_SolProp_FormationRates (double * Rgas, double * Rsolid) {
  kineticsSolidMapXML->FormationRates(Rgas, Rsolid);
}

double OpenSMOKE_GasProp_HeatRelease (const double * R) {
  return kineticsMapXML->HeatRelease(R);
}

double OpenSMOKE_GasProp_kPlanckMix (const double * x) {
  return transportMapXML->kPlanckMix (x);
}

double OpenSMOKE_SolProp_HeatRelease (const double * Rgas, const double * Rsolid) {
  return kineticsSolidMapXML->HeatRelease(Rgas, Rsolid);
}

void OpenSMOKE_SurProp_KineticConstants (void) {
  kineticsSurfaceMapXML->KineticConstants();
}

void OpenSMOKE_SurProp_ReactionRates (const double* c, const double* z,
    const double* a, const double* gamma) {
  kineticsSurfaceMapXML->ReactionRates (c, z, a, gamma);
}

void OpenSMOKE_SurProp_FormationRatesGasOnly (double * Rgas) {
  kineticsSurfaceMapXML->FormationRates (Rgas);
}

void OpenSMOKE_SurProp_FormationRates (double * Rgas, double * Rsite,
    double * Rbulk, double * RsitePhases) {
  kineticsSurfaceMapXML->FormationRates (Rgas, Rsite, Rbulk, RsitePhases);
}

double OpenSMOKE_SurProp_HeatRelease (const double* Rgas, const double* Rsite,
    const double* Rbulk) {
  return kineticsSurfaceMapXML->HeatRelease (Rgas, Rsite, Rbulk);
}

double OpenSMOKE_MW (const int i) {
  return thermodynamicsMapXML->MW(i);
}

double OpenSMOKE_MW_Solid (const int i) {
  unsigned int ngs = thermodynamicsSolidMapXML->number_of_gas_species();
  return thermodynamicsSolidMapXML->MW(i + ngs);
}

const char* OpenSMOKE_NamesOfSpecies (const int i) {
  return thermodynamicsMapXML->NamesOfSpecies()[i].c_str();
}

const char* OpenSMOKE_NamesOfLiquidSpecies (const int i) {
  return thermodynamicsLiquidMapXML->vector_names_liquid_species()[i].c_str();
}

const char* OpenSMOKE_NamesOfSolidSpecies (const int i) {
  return thermodynamicsSolidMapXML->vector_names_solid_species()[i].c_str();
}

int OpenSMOKE_IndexOfSpecies (const char* s) {
  return thermodynamicsMapXML->IndexOfSpecies(s) - 1;
}

int OpenSMOKE_IndexOfSolidSpecies (const char* s) {
  unsigned int ngs = thermodynamicsSolidMapXML->number_of_gas_species();
  return thermodynamicsSolidMapXML->IndexOfSpecies(s) - 1 - ngs;
}

double OpenSMOKE_GasProp_DynamicViscosity (double* x) {
  return transportMapXML->DynamicViscosity(x);
}

double OpenSMOKE_GasProp_ThermalConductivity (double* x) {
  return transportMapXML->ThermalConductivity(x);
}

double OpenSMOKE_GasProp_HeatCapacity (const double* x) {
  double MW = OpenSMOKE_MolecularWeight_From_MoleFractions (x);
  double Cp = thermodynamicsMapXML->cpMolar_Mixture_From_MoleFractions(x);
  return Cp/MW;
}

double OpenSMOKE_SolProp_HeatCapacity (const double* x) {
  std::vector<double> Cp (thermodynamicsSolidMapXML->number_of_solid_species());
  thermodynamicsSolidMapXML->cpMolar_SolidSpecies(Cp.data());

  double CpMix = 0.;
  for (unsigned int i=0; i<thermodynamicsSolidMapXML->number_of_solid_species(); i++)
    CpMix += x[i]*Cp[i]/OpenSMOKE_MW_Solid(i);

  return CpMix;
}

void OpenSMOKE_GasProp_HeatCapacity_PureSpecies (double * cp) {
  thermodynamicsMapXML->cpMolar_Species (cp);
  for (unsigned int i=0; i<OpenSMOKE_NumberOfSpecies(); i++)
    cp[i] /= OpenSMOKE_MW (i);
}

double OpenSMOKE_GasProp_Dmix (const double* x, const int i) {
  const double baseline_diffusion = 1e-10;
  double Diffs[thermodynamicsMapXML->NumberOfSpecies()];
  transportMapXML->MassDiffusionCoefficients (Diffs, x, transportMapXML->is_species_bundling());
  return Diffs[i] + baseline_diffusion;
}

double OpenSMOKE_GasProp_Concentrations (double T, double P) {
  return P/PhysicalConstants::R_J_kmol/T;
}

double OpenSMOKE_GasProp_Density_IdealGas (double T, double P, double MW) {
  return OpenSMOKE_GasProp_Concentrations(T,P)*MW;
}

double OpenSMOKE_LiqProp_Density_PureSpecies (const char* s, double T, double P) {
  return species_map->rhoL(s, T, P);
}

double OpenSMOKE_LiqProp_Density_Mix_AddVol (double T, double P, const double* x) {
  double y[OpenSMOKE_NumberOfLiquidSpecies()];
  double MWmix = 0.;
  thermodynamicsLiquidMapXML->LiquidMassFractions_From_LiquidMoleFractions (y, MWmix, x);
  double urho = 0.;
  for (unsigned int i=0; i<OpenSMOKE_NumberOfLiquidSpecies(); i++) {
    std::string species_name = OpenSMOKE_NamesOfLiquidSpecies (i);
    species_name.erase (species_name.length() - 3);
    urho += y[i]/species_map->rhoL(species_name, T, P);
  }
  return 1./urho;
}

double OpenSMOKE_LiqProp_DynamicViscosity_PureSpecies (const char* s, double T) {
  return species_map->etaL(s, T);
}

double OpenSMOKE_LiqProp_DynamicViscosity_Mix (double T, const double* x) {
  double mumix = 0.;
  for (unsigned int i=0; i<OpenSMOKE_NumberOfLiquidSpecies(); i++) {
    std::string species_name = OpenSMOKE_NamesOfLiquidSpecies (i);
    species_name.erase (species_name.length() - 3);
    double etaLi = species_map->etaL (species_name, T);
    mumix += x[i] * std::log (etaLi);
  }
  return std::exp (mumix);
}

double OpenSMOKE_LiqProp_ThermalConductivity_PureSpecies (const char* s, double T) {
  return species_map->lambdaL(s, T);
}

double OpenSMOKE_LiqProp_ThermalConductivity_Mix (double T, const double* x) {
  // Vredeveld method (1973) - From Prausnitz page 10.60
  double y[OpenSMOKE_NumberOfLiquidSpecies()];
  double MWmix = 0.;
  thermodynamicsLiquidMapXML->LiquidMassFractions_From_LiquidMoleFractions (y, MWmix, x);
  double lambdaMix = 0.;
  for (unsigned int i=0; i<OpenSMOKE_NumberOfLiquidSpecies(); i++) {
    std::string species_name = OpenSMOKE_NamesOfLiquidSpecies (i);
    species_name.erase (species_name.length() - 3);
    double lambdaLi = species_map->lambdaL(species_name, T);
    lambdaMix += y[i] * std::pow(lambdaLi, -2.);
  }
  return std::pow (lambdaMix, -0.5);
}

double OpenSMOKE_LiqProp_HeatCapacity_PureSpecies (const char* s, double T) {
  return species_map->cpL(s, T);
}

double OpenSMOKE_LiqProp_HeatCapacity_Mix (double T, const double* x) {
  double y[OpenSMOKE_NumberOfLiquidSpecies()];
  double MWmix = 0.;
  thermodynamicsLiquidMapXML->LiquidMassFractions_From_LiquidMoleFractions (y, MWmix, x);
  double cpmix = 0.;
  for (unsigned int i=0; i<OpenSMOKE_NumberOfLiquidSpecies(); i++) {
    std::string species_name = OpenSMOKE_NamesOfLiquidSpecies (i);
    species_name.erase (species_name.length() - 3);
    double cpi = species_map->cpL (species_name, T);
    cpmix += y[i] * cpi;
  }
  return cpmix;
}

double OpenSMOKE_LiqProp_VaporPressure (const char* s, double T, double P) {
  return species_map->pVap(s, T, P);
}

double OpenSMOKE_LiqProp_VaporizationEnthalpy (const char* s, double T) {
  return species_map->deltaHv(s, T);
}

double OpenSMOKE_LiqProp_Sigma (const char* s, double T) {
  return species_map->sigma(s, T);
}

double OpenSMOKE_LiqProp_Dmix_PerkinsGeankopolis (double T, double P, const double* x, const int i) {
  if (OpenSMOKE_NumberOfLiquidSpecies() == 1) {
    return 0;
  }
  else {
    Eigen::MatrixXd Dinf = liqmix->Dinf (T,P);
    double Dmix = 0.;
    for (int j=0; j<OpenSMOKE_NumberOfLiquidSpecies(); j++) {
      if (i != j) {
        double etaLj = OpenSMOKE_LiqProp_DynamicViscosity_PureSpecies (liqspecies[j].data(), T);
        double xjDjmuj = x[j]*Dinf(i,j)*std::pow (etaLj, 0.8);
        Dmix += xjDjmuj;
      }
    }
    Dmix /= std::pow (OpenSMOKE_LiqProp_DynamicViscosity_Mix(T,x), 0.8);
    return Dmix;
  }
}

double OpenSMOKE_LiqProp_Dmix_Cullinan (double T, double P, const double* x, const int i) {
  if (OpenSMOKE_NumberOfLiquidSpecies() == 1) {
    return 0;
  }
  else {
    Eigen::MatrixXd Dinf = liqmix->Dinf (T,P);
    double Dmix = 1.;
    for (int j=0; j<OpenSMOKE_NumberOfLiquidSpecies(); j++) {
      Dmix *= std::pow (Dinf(i,j), x[j]);
    }
    return Dmix;
  }
}

double OpenSMOKE_LiqProp_Dmix_LefflerCullinan (double T, double P, const double* x, const int i) {
  if (OpenSMOKE_NumberOfLiquidSpecies() == 1) {
    return 0;
  }
  else {
    Eigen::MatrixXd Dinf = liqmix->Dinf (T,P);
    double Dmix = 1.;
    for (int j=0; j<OpenSMOKE_NumberOfLiquidSpecies(); j++) {
      double etaLj = OpenSMOKE_LiqProp_DynamicViscosity_PureSpecies (liqspecies[j].data(), T);
      Dmix *= std::pow (Dinf(i,j)*etaLj, x[j]);
    }
    Dmix /= OpenSMOKE_LiqProp_DynamicViscosity_Mix(T,x);
    return Dmix;
  }
}

double OpenSMOKE_MolecularWeight_From_MoleFractions (const double* x) {
  return thermodynamicsMapXML->MolecularWeight_From_MoleFractions(x);
}

double OpenSMOKE_MolecularWeight_From_MassFractions (const double* x) {
  return thermodynamicsMapXML->MolecularWeight_From_MassFractions(x);
}

double OpenSMOKE_LiquidMolecularWeight_From_LiquidMoleFractions (const double* x) {
  return thermodynamicsLiquidMapXML->LiquidMolecularWeight_From_LiquidMoleFractions (x);
}

double OpenSMOKE_LiquidMolecularWeight_From_LiquidMassFractions (const double* x) {
  return thermodynamicsLiquidMapXML->LiquidMolecularWeight_From_LiquidMassFractions (x);
}

void OpenSMOKE_MassFractions_From_MoleFractions (double* y, double* MW, const double* x) {
  thermodynamicsMapXML->MassFractions_From_MoleFractions (y, *MW, x);
}

void OpenSMOKE_MoleFractions_From_MassFractions (double* x, double* MW, const double* y) {
   thermodynamicsMapXML->MoleFractions_From_MassFractions (x, *MW ,y);
}

void OpenSMOKE_SolidMoleFractions_From_SolidMassFractions (double* x, double* MW, const double* y) {
  thermodynamicsSolidMapXML->SolidMoleFractions_From_SolidMassFractions (x, *MW, y);
}

void OpenSMOKE_LiquidMassFractions_From_LiquidMoleFractions (double* y, double* MW, const double* x) {
  thermodynamicsLiquidMapXML->LiquidMassFractions_From_LiquidMoleFractions (y, *MW, x);
}

void OpenSMOKE_LiquidMoleFractions_From_LiquidMassFractions (double* x, double* MW, const double* y) {
  thermodynamicsLiquidMapXML->LiquidMoleFractions_From_LiquidMassFractions (x, *MW, y);
}

void OpenSMOKE_CheckAndCorrectSumOfFractions (double* x, int n) {
  std::vector<double> v(n);
  for (int jj=0; jj<v.size(); jj++)
    v[jj] = (x[jj] > 0.) ? x[jj] : 0.;
  OpenSMOKE::CheckAndCorrectSumOfFractions (v);
  for (int jj=0; jj<v.size(); jj++)
    x[jj] = v[jj];
}

double OpenSMOKE_GetMixtureFractionFromMassFractions (const double* y, const double* yfuel, const double* yox) {
  std::vector<std::string> names = thermodynamicsMapXML->NamesOfSpecies();
  std::vector<double> yfuelstd (thermodynamicsMapXML->NumberOfSpecies());
  std::vector<double> yoxstd (thermodynamicsMapXML->NumberOfSpecies());
  for (int jj=0; jj<thermodynamicsMapXML->NumberOfSpecies(); jj++) {
    yfuelstd[jj] = yfuel[jj];
    yoxstd[jj] = yox[jj];
  }
  return thermodynamicsMapXML->GetMixtureFractionFromMassFractions (y, names, yfuelstd, names, yoxstd);
}

void OpenSMOKE_ODESolver
(
  odefunction ode,
  unsigned int neq,
  double dt,
  double * y,
  void * args
)
{
  // Set time step and initial values
  double t0 = 0., tf = t0 + dt;
  Eigen::VectorXd y0_eigen(neq);
  for (int i=0; i<y0_eigen.size(); i++)
    y0_eigen[i] = y[i];

  // Create the solver
  typedef OdeSMOKE::KernelDense<OpenSMOKE::ODESystem_Interface> denseOde;
  typedef OdeSMOKE::MethodGear<denseOde> methodGear;
  OdeSMOKE::MultiValueSolver<methodGear> ode_solver;

  //// Set the ODE system of equations function
  ode_solver.SetUserArgs(args);
  ode_solver.SetSystemOfEquations(ode);

  //// Set initial conditions
  ode_solver.SetInitialConditions(t0, y0_eigen);

  //// Set linear algebra options
  ode_solver.SetLinearAlgebraSolver(ode_parameters_->linear_algebra());
  ode_solver.SetFullPivoting(ode_parameters_->full_pivoting());

  //// Set relative and absolute tolerances
  ode_solver.SetAbsoluteTolerances(ode_parameters_->absolute_tolerance());
  ode_solver.SetRelativeTolerances(ode_parameters_->relative_tolerance());

  //Eigen::VectorXd ymins(neq), ymaxs(neq);
  //ymins.setConstant (0.);
  //ymaxs.setConstant (1.);
  //ymins[neq-1] = -1e-32;
  //ymaxs[neq-1] = +1e+32;

  //// Set minimum and maximum values
  //ode_solver.SetMinimumValues(ymins);
  //ode_solver.SetMaximumValues(ymaxs);

  // Solve the system
  double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
  OdeSMOKE::OdeStatus status = ode_solver.Solve(tf);
  double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

  // Check the solution
  if (status > 0)
  {
    Eigen::VectorXd yf_eigen;
    ode_solver.Solution(yf_eigen);
    for (int i=0; i<yf_eigen.size(); i++)
      y[i] = yf_eigen[i];
  }
  else
  {
    std::cout<< "Solution not found." << std::endl;
  }
}

#ifdef __cplusplus
}
#endif

