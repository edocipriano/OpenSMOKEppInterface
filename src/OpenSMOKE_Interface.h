
#ifndef OPENSMOKE_INTERFACE_H_
#define OPENSMOKE_INTERFACE_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef void(*odefunction)(const double * y, const double t, double * dy, void * args);

void OpenSMOKE_Init (void);
void OpenSMOKE_InitODESolver (void);
void OpenSMOKE_Clean (void);
void OpenSMOKE_CleanODESolver (void);
void OpenSMOKE_ReadKinetics (const char* kinfolder);
void OpenSMOKE_ReadLiquidKinetics (const char* kinfolder);
void OpenSMOKE_ReadLiquidProperties (const char* liqpropfolder);
int OpenSMOKE_NumberOfSpecies (void);
int OpenSMOKE_NumberOfLiquidSpecies (void);
int OpenSMOKE_NumberOfReactions (void);
double OpenSMOKE_Printpi (void);
void OpenSMOKE_GasProp_SetPressure (const double P);
void OpenSMOKE_GasProp_SetTemperature (const double T);
void OpenSMOKE_GasProp_KineticConstants (void);
void OpenSMOKE_GasProp_ReactionRates (const double * c);
void OpenSMOKE_GasProp_FormationRates (double * R);
double OpenSMOKE_GasProp_HeatRelease (const double * R);
double OpenSMOKE_MW (const int i);
const char* OpenSMOKE_NamesOfSpecies (const int i);
const char* OpenSMOKE_NamesOfLiquidSpecies (const int i);
int OpenSMOKE_IndexOfSpecies (const char* s);
double OpenSMOKE_GasProp_DynamicViscosity (double* x);
double OpenSMOKE_GasProp_ThermalConductivity (double* x);
double OpenSMOKE_GasProp_Concentrations (double T, double P);
double OpenSMOKE_GasProp_Density_IdealGas (double T, double P, double MW);
double OpenSMOKE_GasProp_cpMolar_Mixture_From_MoleFractions (const double* x);
double OpenSMOKE_LiqProp_Density_PureSpecies (const char* s, double T, double P);
double OpenSMOKE_LiqProp_Density_Mix_AddVol (double T, double P, const double* x);
double OpenSMOKE_LiqProp_DynamicViscosity_PureSpecies (const char* s, double T);
double OpenSMOKE_LiqProp_DynamicViscosity_Mix (double T, const double* x);
double OpenSMOKE_LiqProp_ThermalConductivity_PureSpecies (const char* s, double T);
double OpenSMOKE_LiqProp_ThermalConductivity_Mix (double T, const double* x);
double OpenSMOKE_LiqProp_cp_PureSpecies (const char* s, double T);
double OpenSMOKE_LiqProp_cp_Mix (double T, const double* x);
double OpenSMOKE_MolecularWeight_From_MoleFractions (const double* x);
double OpenSMOKE_MolecularWeight_From_MassFractions (const double* x);
double OpenSMOKE_LiquidMolecularWeight_From_LiquidMoleFractions (const double* x);
double OpenSMOKE_LiquidMolecularWeight_From_LiquidMassFractions (const double* x);
void OpenSMOKE_MassFractions_From_MoleFractions (double* y, double* MW, const double* x);
void OpenSMOKE_MoleFractions_From_MassFractions (double* x, double* MW, const double* y);
void OpenSMOKE_LiquidMassFractions_From_LiquidMoleFractions (double* y, double* MW, const double* x);
void OpenSMOKE_LiquidMoleFractions_From_LiquidMassFractions (double* x, double* MW, const double* y);
double OpenSMOKE_GetMixtureFractionFromMassFractions (const double* y, const double* yfuel, const double* yox);
void OpenSMOKE_ODESolver (odefunction ode, unsigned int neq, double dt, double * y, void * args);

#ifdef __cplusplus
}
#endif

#endif

