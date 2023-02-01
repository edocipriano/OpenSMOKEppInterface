
#ifndef OPENSMOKE_INTERFACE_H_
#define OPENSMOKE_INTERFACE_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef void(*odefunction)(const double * y, const double t, double * dy, void * args);

void OpenSMOKE_Init (void);
void OpenSMOKE_ReadKinetics (void);
void OpenSMOKE_ReadLiquidKinetics (void);
void OpenSMOKE_ReadLiquidProperties (void);
int OpenSMOKE_NumberOfSpecies (void);
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
int OpenSMOKE_IndexOfSpecies (const char* s);
double OpenSMOKE_GasProp_DynamicViscosity (double* x);
double OpenSMOKE_GasProp_ThermalConductivity (double* x);
double OpenSMOKE_GasProp_Concentrations (double T, double P);
double OpenSMOKE_GasProp_Density_IdealGas (double T, double P, double MW);
double OpenSMOKE_GasProp_cpMolar_Mixture_From_MoleFractions (const double* x);
double OpenSMOKE_MolecularWeight_From_MoleFractions (const double* x);
double OpenSMOKE_MolecularWeight_From_MassFractions (const double* x);
void OpenSMOKE_MassFractions_From_MoleFractions (double* y, double* MW, const double* x);
void OpenSMOKE_MoleFractions_From_MassFractions (double* x, double* MW, const double* y);
void OpenSMOKE_ODESolver (odefunction ode, int neq, double dt, double * y, void * args);

#ifdef __cplusplus
}
#endif

#endif

