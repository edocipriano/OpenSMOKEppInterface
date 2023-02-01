
#include "OpenSMOKE_Interface.h"
#include "OpenSMOKE_ODESystem_Interface.h"

// OpenSMOKE++ Definitions
#include "OpenSMOKEpp"
#include "math/PhysicalConstants.h"

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"
#include "maps/ThermodynamicsMap_Liquid_CHEMKIN.h"
#include "maps/KineticsMap_Liquid_CHEMKIN.h"

// OpenSMOKE++ Thermodynamics and liquid phase properties
#include "thermodynamics/mixture/mixtureL/mixtureL.h"
#include "thermodynamics/mixture/mixtureG/mixtureG.h"

// OpenSMOKE++ Dictionaries
#include "dictionary/OpenSMOKE_Dictionary"

// ODE solvers
#include "math/native-ode-solvers/MultiValueSolver"
#include "math/external-ode-solvers/ODE_Parameters.h"

// Utilities
#include "idealreactors/utilities/Utilities"
#include "utilities/ropa/OnTheFlyROPA.h"
#include "utilities/cema/OnTheFlyCEMA.h"
#include "utilities/ontheflypostprocessing/OnTheFlyPostProcessing.h"
#include "utilities/kineticsmodifier/KineticsModifier.h"

// PolimiSoot Analyzer
#include "utilities/soot/polimi/OpenSMOKE_PolimiSoot_Analyzer.h"

// Batch reactor
#include "idealreactors/batch/BatchReactor"

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

void OpenSMOKE_GasProp_KineticConstants (void) {
  kineticsMapXML->KineticConstants();
}

void OpenSMOKE_GasProp_ReactionRates (const double * c) {
  kineticsMapXML->ReactionRates(c);
}

void OpenSMOKE_GasProp_FormationRates (double * R) {
  kineticsMapXML->FormationRates(R);
}

double OpenSMOKE_GasProp_HeatRelease (const double * R) {
  return kineticsMapXML->HeatRelease(R);
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

void OpenSMOKE_ODESolver
(
  odefunction ode,
  int neq,
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

  // Create ODE parameters class
  OpenSMOKE::ODE_Parameters ode_parameters_;

  // Create the solver
  typedef OdeSMOKE::KernelDense<OpenSMOKE::ODESystem_Interface> denseOde;
  typedef OdeSMOKE::MethodGear<denseOde> methodGear;
  OdeSMOKE::MultiValueSolver<methodGear> ode_solver;

  // Set the ODE system of equations function
  ode_solver.SetUserArgs(args);
  ode_solver.SetSystemOfEquations(ode);

  // Set initial conditions
  ode_solver.SetInitialConditions(t0, y0_eigen);

  // Set linear algebra options
  ode_solver.SetLinearAlgebraSolver(ode_parameters_.linear_algebra());
  ode_solver.SetFullPivoting(ode_parameters_.full_pivoting());

  // Set relative and absolute tolerances
  ode_solver.SetAbsoluteTolerances(ode_parameters_.absolute_tolerance());
  ode_solver.SetRelativeTolerances(ode_parameters_.relative_tolerance());

  // Set minimum and maximum values
  //ode_solver.SetMinimumValues(0.);
  //ode_solver.SetMaximumValues(1.);

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

