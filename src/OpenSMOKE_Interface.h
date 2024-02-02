/**
# OpenSMOKE Interface

C Interface for the OpenSMOKE++ library. This module allows
the OpenSMOKE++ library to be used with a high level
interface, and using a shared library, in order to save
compilation time, facilitate the use of the library and
with the possibility to call OpenSMOKE++ from C codes.
*/

#ifndef OPENSMOKE_INTERFACE_H_
#define OPENSMOKE_INTERFACE_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef void(*odefunction)(const double * y, const double t, double * dy, void * args);

/**
## Initialize, Finalize and Read Data

Setup and cleaning functions for the OpenSMOKE pointers
*/

/**
### *OpenSMOKE_Init()*: Initialize the OpenSMOKE pointers
*/
void
  OpenSMOKE_Init (void);

/**
### *OpenSMOKE_InitODESolver()*: Initialize ODE Solver pointers
*/
void
  OpenSMOKE_InitODESolver (void);

/**
### *OpenSMOKE_Clean()*: Free the OpenSMOKE pointers from memory
*/
void
  OpenSMOKE_Clean (void);

/**
### *OpenSMOKE_CleanODESolver()*: Free the ODE solver pointers from memory
*/
void
  OpenSMOKE_CleanODESolver (void);

/**
## Read Kinetics

Function for reading kinetics and LiquidProperties, and return
basic information.
*/

/**
### *OpenSMOKE_ReadKinetics()*: Read gas phase kinetics

* *kinfolder*: path to the kinetics folder containing kinetics.xml
*/
void
  OpenSMOKE_ReadKinetics (const char* kinfolder);

/**
### *OpenSMOKE_ReadLiquidKinetics()*: Read liquid phase kinetics

* *kinfolder*: path to the kinetics folder containing kinetics.liquid.xml
*/
void
  OpenSMOKE_ReadLiquidKinetics (const char* kinfolder);

/**
### *OpenSMOKE_ReadSolidKinetics()*: Read solid phase kinetics

* *kinfolder*: path to the kinetics folder containing kinetics.solid.xml
*/
void
  OpenSMOKE_ReadSolidKinetics (const char* kinfolder);

/**
### *OpenSMOKE_ReadLiquidProperties()*: Read liquid properties folder

* *liqpropfolder*: path to LiquidProperties
*/
void
  OpenSMOKE_ReadLiquidProperties (const char* liqpropfolder);

/**
### *OpenSMOKE_NumberOfSpecies()*: Number of species in gas phase
*/
int
  OpenSMOKE_NumberOfSpecies (void);

/**
### *OpenSMOKE_NumberOfLiquidSpecies()*: Number of species in liquid phase
*/
int
  OpenSMOKE_NumberOfLiquidSpecies (void);

/**
### *OpenSMOKE_NumberOfSolidSpecies()*: Number of species in solid phase
*/
int
  OpenSMOKE_NumberOfSolidSpecies (void);

/**
### *OpenSMOKE_NumberOfReactions()*: Number of reactions in gas phase
*/
int
  OpenSMOKE_NumberOfReactions (void);

/**
### *OpenSMOKE_NumberOfSolidReactions()*: Number of reactions in gas solid
*/
int
  OpenSMOKE_NumberOfSolidReactions (void);
/**
### *OpenSMOKE_Printpi()*: Print pi constant (debug)
*/
double
  OpenSMOKE_Printpi (void);

/**
### *OpenSMOKE_NamesOfSpecies()*: Names of the gas phase species

* *i*: index of the species based on the number of gas species
*/
const char*
  OpenSMOKE_NamesOfSpecies (const int i);

/**
### *OpenSMOKE_NamesOfLiquidSpecies()*: Names of the liquid phase species

* *i*: index of the species based on the number of liquid species
*/
const char*
  OpenSMOKE_NamesOfLiquidSpecies (const int i);

/**
### *OpenSMOKE_NamesOfSolidSpecies()*: Names of the solid phase species

* *i*: index of the species based on the number of ONLY solid species 
*/
const char*
  OpenSMOKE_NamesOfSolidSpecies (const int i);

/**
### *OpenSMOKE_IndexOfSpecies()*: Index of the specific species in gas phase

* *s*: name of the species under investigation
*/
int
  OpenSMOKE_IndexOfSpecies (const char* s);

/**
### *OpenSMOKE_IndexOfSolidSpecies()*: Index of the specific species in solid phase

* *s*: name of the species under investigation
*/
int
  OpenSMOKE_IndexOfSolidSpecies (const char* s);

/**
# Gas Phase Thermodynamics

Functions for the calculation of thermodynamic properties in gas phase.
*/

/**
### *OpenSMOKE_MW()*: Molecular weights of the gas phase species

* *i*: index of the species
*/
double
  OpenSMOKE_MW (const int i);

/**
### *OpenSMOKE_GasProp_SetPressure()*: Set Pressure in gas phase

* *P*: pressure in gas phase
*/
void
  OpenSMOKE_GasProp_SetPressure (const double P);

/**
### *OpenSMOKE_GasProp_SetTemperature()*: Set Temperature in gas phase

* *T*: temperature in gas phase
*/
void
  OpenSMOKE_GasProp_SetTemperature (const double T);

/**
### *OpenSMOKE_LiqProp_SetPressure()*: Set Pressure in gas phase

* *P*: pressure in gas phase
*/
void
  OpenSMOKE_LiqProp_SetPressure (const double P);

/**
### *OpenSMOKE_LiqProp_SetTemperature()*: Set Temperature in gas phase

* *T*: temperature in gas phase
*/
void
  OpenSMOKE_LiqProp_SetTemperature (const double T);

/**
# Solid Phase Thermodynamics

Functions for the calculation of thermodynamic properties in solid phase.
*/

/**
### *OpenSMOKE_MW()*: Molecular weights of the solid phase species

* *i*: index of the species
*/
double
  OpenSMOKE_MW_Solid (const int i);

/**
### *OpenSMOKE_SolProp_SetPressure()*: Set Pressure in gas solid

* *P*: pressure in solid phase
*/
void
  OpenSMOKE_SolProp_SetPressure (const double P);

/**
### *OpenSMOKE_SolProp_SetTemperature()*: Set Temperature in solid phase

* *T*: temperature in solid phase
*/
void
  OpenSMOKE_SolProp_SetTemperature (const double T);

/**
### *OpenSMOKE_GasProp_DynamicViscosity()*: Density of the gas phase mixture

* *x*: mole fractions in gas phase
*/
double
  OpenSMOKE_GasProp_DynamicViscosity (double* x);

/**
### *OpenSMOKE_GasProp_ThermalConductivity()*: Thermal conductivity of the gas phase mixture

* *x*: mole fractions in gas phase
*/
double
  OpenSMOKE_GasProp_ThermalConductivity (double* x);

/**
### *OpenSMOKE_GasProp_Concentrations()*: Concentration of the gas phase mixture

* *T*: temperature in gas phase
* *P*: pressure in gas phase
*/
double
  OpenSMOKE_GasProp_Concentrations (double T, double P);

/**
### *OpenSMOKE_GasProp_Density_IdealGas()*: Density of the gas phase mixture as an ideal gas

* *T*: temperature in gas phase
* *P*: pressure in gas phase
* *MW*: gas mixture molecular weight
*/
double
  OpenSMOKE_GasProp_Density_IdealGas (double T, double P, double MW);

/**
### *OpenSMOKE_GasProp_HeatCapacity()*: Specific heat capacity of the gas phase mixture

* *x*: mole fractions in gas phase
*/
double
  OpenSMOKE_GasProp_HeatCapacity (const double* x);

/**
### *OpenSMOKE_SolProp_HeatCapacity()*: Specific heat capacity of the solid phase mixture

* *x*: mole fractions in solid phase
*/
double
  OpenSMOKE_SolProp_HeatCapacity (const double* x);

/**
### *OpenSMOKE_GasProp_SpeciesHeatCapacity()*: Specific heat capacity of species i in the gas phase

* *i*: index of species
*/
void
  OpenSMOKE_GasProp_HeatCapacity_PureSpecies (double * cp);

/**
### *OpenSMOKE_GasProp_Dmix()* Mixture diffusion coefficients

* *x*: mole fractions in gas phase
* *i*: index of the species
*/
double
  OpenSMOKE_GasProp_Dmix (const double* x, const int i);

/**
### *OpenSMOKE_MolecularWeight_From_MoleFractions()*: Mixture molecular weight from mole fractions in gas phase

* *x*: mole fractions in gas phase
*/
double
  OpenSMOKE_MolecularWeight_From_MoleFractions (const double* x);

/**
### *OpenSMOKE_MolecularWeight_From_MassFractions()*: Mixture molecular weight from mass fractions in gas phase

* *x*: mass fractions in gas phase
*/
double
  OpenSMOKE_MolecularWeight_From_MassFractions (const double* x);

/**
### *OpenSMOKE_MassFractions_From_MoleFractions()*: Mass fractions and mixture molecular weight from mole fractions in gas phase

* *y*: mass fractions in gas phase (to be computed)
* *MW*: gas phase mixture molecular weight (to be computed)
* *x*: mole fractions in gas phase
*/
void
  OpenSMOKE_MassFractions_From_MoleFractions (double* y, double* MW, const double* x);

/**
### *OpenSMOKE_MoleFractions_From_MassFractions()*: Mole fractions and mixture molecular weight from mass fractions in gas phase

* *x*: mole fractions in gas phase (to be computed)
* *MW*: gas phase mixture molecular weight (to be computed)
* *x*: mass fractions in gas phase
*/
void
  OpenSMOKE_MoleFractions_From_MassFractions (double* x, double* MW, const double* y);

/**
## Gas Phase Kinetics

Functions for the calculation of kinetic properties in gas phase.
*/

/**
### *OpenSMOKE_GasProp_KineticConstants()*: Compute kinetic constants in gas phase
*/
void
  OpenSMOKE_GasProp_KineticConstants (void);

/**
### *OpenSMOKE_GasProp_ReactionRates()*: Compute the reaction rates in gas phase

* *c*: concentration of each chemical species in gas phase
*/

void
  OpenSMOKE_GasProp_ReactionRates (const double * c);

/**
### *OpenSMOKE_GasProp_FormationRates()*: Return the formation rates in gas phase

* *R*: formation rate [kmol/m3/s] of each chemical
species in gas phase (to be computed)
*/
void OpenSMOKE_GasProp_FormationRates (double * R);

/**
### *OpenSMOKE_GasProp_HeatRelease()*: Compute the heat released from the gas phase reactions

* *R*: formation rate [kmol/m3/s] of each chemical species in gas phase
*/
double
  OpenSMOKE_GasProp_HeatRelease (const double * R);

/**
## Solid Phase Kinetics

Functions for the calculation of kinetic properties in solid phase.
*/

/**
### *OpenSMOKE_SolProp_KineticConstants()*: Compute kinetic constants in solid phase
*/
void
  OpenSMOKE_SolProp_KineticConstants (void);

/**
### *OpenSMOKE_SolProp_ReactionRates()*: Compute the reaction rates in solid phase

* *cgas*: concentration of each chemical species in gas phase
* *csold*: concentration of each chemical species in solid phase
*/

void
  OpenSMOKE_SolProp_ReactionRates (const double * cgas, const double * csolid);

/**
### *OpenSMOKE_SolidProp_FormationRates()*: Return the formation rates in solid phase

* *Rgas*: formation rate [kmol/m3/s] of each chemical
species in gas phase (to be computed)
* *Rsolid*: formation rate [kmol/m3/s] of each chemical
species in gas phase (to be computed)
*/
void OpenSMOKE_SolProp_FormationRates (double * Rgas, double * Rsolid);

/**
### *OpenSMOKE_SolProp_HeatRelease()*: Compute the heat released from the solid phase reactions

* *Rgas*: formation rate [kmol/m3/s] of each chemical species in gas phase
* *Rsolid*: formation rate [kmol/m3/s] of each chemical species in solid phase
*/
double
  OpenSMOKE_SolProp_HeatRelease (double * Rgas, double * Rsolid);

/**
## Liquid Phase Thermodynamics

Functions for the calculation of thermodynamic properties in liquid phase.
*/

/**
### *OpenSMOKE_LiqProp_Density_PureSpecies()*: Density of a pure species in liquid phase

* *s*: name of the species in liquid phase
* *T*: temperature in liquid phase
* *P*: pressure in liquid phase
*/
double
  OpenSMOKE_LiqProp_Density_PureSpecies (const char* s, double T, double P);

/**
### *OpenSMOKE_LiqProp_Density_Mix_AddVol()*: Density of the liquid mixture with the added volume method

* *T*: temperature in liquid phase
* *P*: pressure in liquid phase
* *x*: mole fractions in liquid phase
*/
double
  OpenSMOKE_LiqProp_Density_Mix_AddVol (double T, double P, const double* x);

/**
### *OpenSMOKE_LiqProp_DynamicViscosity_PureSpecies()*: Dynamic viscosity of a pure species in liquid phase

* *s*: name of the species in liquid phase
* *T*: temperature in liquid phase
*/
double
  OpenSMOKE_LiqProp_DynamicViscosity_PureSpecies (const char* s, double T);

/**
### *OpenSMOKE_LiqProp_DynamicViscosity_Mix()*: Dynamic viscosity of the liquid mixture

* *T*: temperature in liquid phase
* *x*: mole fractions in liquid phase
*/
double
  OpenSMOKE_LiqProp_DynamicViscosity_Mix (double T, const double* x);

/**
### *OpenSMOKE_LiqProp_ThermalConductivity_PureSpecies()*: Thermal conductivity of a pure species in liquid phase

* *s*: name of the species in liquid phase
* *T*: temperature in liquid phase
*/
double
  OpenSMOKE_LiqProp_ThermalConductivity_PureSpecies (const char* s, double T);

/**
### *OpenSMOKE_LiqProp_ThermalConductivity_Mix()*: Thermal conductivity of the liquid mixture

* *T*: temperature in liquid phase
* *x*: mole fractions in liquid phase
*/
double
  OpenSMOKE_LiqProp_ThermalConductivity_Mix (double T, const double* x);

/**
### *OpenSMOKE_LiqProp_HeatCapacity_PureSpecies()*: Specific heat capacity of a pure species in liquid phase

* *s*: name of the species in liquid phase
* *T*: temperature in liquid phase
*/
double
  OpenSMOKE_LiqProp_HeatCapacity_PureSpecies (const char* s, double T);

/**
### *OpenSMOKE_LiqProp_HeatCapacity_Mix()*: Specific heat capacity of the liquid mixture

* *T*: temperature in liquid phase
* *x*: mole fractions in liquid phase
*/
double
  OpenSMOKE_LiqProp_HeatCapacity_Mix (double T, const double* x);

/**
### *OpenSMOKE_LiqProp_VaporPressure()*: Vapor pressure of a pure species in liquid phase

* *s*: name of the chemical species
* *T*: temperature in liquid phase
* *P*: pressure in liquid phase
*/
double
  OpenSMOKE_LiqProp_VaporPressure (const char* s, double T, double P);

/**
### *OpenSMOKE_LiqProp_VaporizationEnthalpy()*: Latent heat of evaporation

* *s*: name of the chemical species
* *T*: temperature in liquid phase
*/
double
  OpenSMOKE_LiqProp_VaporizationEnthalpy (const char* s, double T);

/**
### *OpenSMOKE_LiqProp_Dmix_PerkinsGeankopolis()* Mixture diffusion coefficients (Perkins Geankopolis model)

* *T*: temperature in liquid phase
* *P*: Pressure in liquid phase
* *x*: mole fractions in liquid phase
* *i*: index of the species
*/
double
  OpenSMOKE_LiqProp_Dmix_PerkinsGeankopolis (double T, double P, const double* x, const int i);

/**
### *OpenSMOKE_LiqProp_Dmix_Cullinan()* Mixture diffusion coefficients (Cullinan model)

* *T*: temperature in liquid phase
* *P*: Pressure in liquid phase
* *x*: mole fractions in liquid phase
* *i*: index of the species
*/
double
  OpenSMOKE_LiqProp_Dmix_Cullinan (double T, double P, const double* x, const int i);

/**
### *OpenSMOKE_LiqProp_Dmix_LefflerCullinan()* Mixture diffusion coefficients (Leffler Cullinan model)

* *T*: temperature in liquid phase
* *P*: Pressure in liquid phase
* *x*: mole fractions in liquid phase
* *i*: index of the species
*/
double
  OpenSMOKE_LiqProp_Dmix_LefflerCullinan (double T, double P, const double* x, const int i);

/**
### *OpenSMOKE_LiquidMolecularWeight_From_LiquidMoleFractions()*: Mixture molecular weight from mole fractions in liquid phase

* *x*: mole fractions in liquid phase
*/
double
  OpenSMOKE_LiquidMolecularWeight_From_LiquidMoleFractions (const double* x);

/**
### *OpenSMOKE_LiquidMolecularWeight_From_LiquidMassFractions()*: Mixture molecular weight from mass fractions in liquid phase

* *x*: mass fractions in liquid phase
*/
double
  OpenSMOKE_LiquidMolecularWeight_From_LiquidMassFractions (const double* x);

/**
### *OpenSMOKE_LiquidMassFractions_From_LiquidMoleFractions()*: Mass fractions and mixture molecular weight from mole fractions in liqid phase

* *y*: mass fractions in liquid phase (to be computed)
* *MW*: liquid phase mixture molecular weight (to be computed)
* *x*: mole fractions in liquid phase
*/
void
  OpenSMOKE_LiquidMassFractions_From_LiquidMoleFractions (double* y, double* MW, const double* x);

/**
### *OpenSMOKE_LiquidMoleFractions_From_LiquidMassFractions()*: Mole fractions and mixture molecular weight from mass fractions in liqid phase

* *x*: mole fractions in liquid phase (to be computed)
* *MW*: liquid phase mixture molecular weight (to be computed)
* *y*: mass fractions in liquid phase
*/
void
  OpenSMOKE_LiquidMoleFractions_From_LiquidMassFractions (double* x, double* MW, const double* y);

/**
## Utilities

Functions for post-processing and numerical algortihms.
*/

/**
### *OpenSMOKE_CheckAndCorrectSumOfFractions()*: Mixture fraction using the Bilger formula

* *x*: mass or mole fractions vector
* *n*: length of the vector
*/
void
  OpenSMOKE_CheckAndCorrectSumOfFractions (double* x, int n);

/**
### *OpenSMOKE_GetMixtureFractionFromMassFractions()*: Mixture fraction using the Bilger formula

* *y*: mole fractions in gas phase
* *yfuel*: stoichiometric mole fractions of the fuel
* *yox*: stoichiometric mole fractions of the oxidizer
*/
double
  OpenSMOKE_GetMixtureFractionFromMassFractions (const double* y, const double* yfuel, const double* yox);

/**
### *OpenSMOKE_ODESolver()*: Native OpenSMOKE ODE solver

* *ode*: pointer to the function implementing the ODE system
* *neq*: number of equations
* *dt*: time step
* *y*: unknowns of the problem (modified)
* *args*: additional struct with arguments required by the odefunction
*/
void
  OpenSMOKE_ODESolver (odefunction ode, unsigned int neq, double dt, double * y, void * args);

#ifdef __cplusplus
}
#endif

#endif

