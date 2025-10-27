import platform
from ctypes import *

system = platform.system()

if system == 'Darwin':
    lib_name = 'libopensmoke.dylib'
elif system == 'Linux':
    lib_name = 'libopensmoke.so'
elif system == 'Windows':
    lib_name = 'opensmoke.dll'
else:
    raise NotImplementedError(f"Unsupported operating system: {system}")

OpenSMOKE = cdll.LoadLibrary ("../build/" + lib_name)

OpenSMOKE_Init = OpenSMOKE.OpenSMOKE_Init
OpenSMOKE_Init.argtypes = None
OpenSMOKE_Init.restype = None

OpenSMOKE_InitODESolver = OpenSMOKE.OpenSMOKE_InitODESolver
OpenSMOKE_InitODESolver.argtypes = None
OpenSMOKE_InitODESolver.restype = None

OpenSMOKE_Clean = OpenSMOKE.OpenSMOKE_Clean
OpenSMOKE_Clean.argtypes = None
OpenSMOKE_Clean.restype = None

OpenSMOKE_CleanODESolver = OpenSMOKE.OpenSMOKE_CleanODESolver
OpenSMOKE_CleanODESolver.argtypes = None
OpenSMOKE_CleanODESolver.restype = None

OpenSMOKE_ReadKinetics = OpenSMOKE.OpenSMOKE_ReadKinetics
OpenSMOKE_ReadKinetics.argtypes = [c_char_p]
OpenSMOKE_ReadKinetics.restype = None

OpenSMOKE_ReadLiquidKinetics = OpenSMOKE.OpenSMOKE_ReadLiquidKinetics
OpenSMOKE_ReadLiquidKinetics.argtypes = [c_char_p]
OpenSMOKE_ReadLiquidKinetics.restype = None

OpenSMOKE_ReadSolidKinetics = OpenSMOKE.OpenSMOKE_ReadSolidKinetics
OpenSMOKE_ReadSolidKinetics.argtypes = [c_char_p]
OpenSMOKE_ReadSolidKinetics.restype = None

OpenSMOKE_ReadLiquidProperties = OpenSMOKE.OpenSMOKE_ReadLiquidProperties
OpenSMOKE_ReadLiquidProperties.argtypes = [c_char_p]
OpenSMOKE_ReadLiquidProperties.restype = None

OpenSMOKE_NumberOfSpecies = OpenSMOKE.OpenSMOKE_NumberOfSpecies
OpenSMOKE_NumberOfSpecies.argtypes = None
OpenSMOKE_NumberOfSpecies.restype = c_int

OpenSMOKE_NumberOfLiquidSpecies = OpenSMOKE.OpenSMOKE_NumberOfLiquidSpecies
OpenSMOKE_NumberOfLiquidSpecies.argtypes = None
OpenSMOKE_NumberOfLiquidSpecies.restype = c_int

OpenSMOKE_NumberOfSolidSpecies = OpenSMOKE.OpenSMOKE_NumberOfSolidSpecies
OpenSMOKE_NumberOfSolidSpecies.argtypes = None
OpenSMOKE_NumberOfSolidSpecies.restype = c_int

OpenSMOKE_NumberOfReactions = OpenSMOKE.OpenSMOKE_NumberOfReactions
OpenSMOKE_NumberOfReactions.argtypes = None
OpenSMOKE_NumberOfReactions.restype = c_int

OpenSMOKE_NumberOfSolidReactions = OpenSMOKE.OpenSMOKE_NumberOfSolidReactions
OpenSMOKE_NumberOfSolidReactions.argtypes = None
OpenSMOKE_NumberOfSolidReactions.restype = c_int

OpenSMOKE_Printpi = OpenSMOKE.OpenSMOKE_Printpi
OpenSMOKE_Printpi.argtypes = None
OpenSMOKE_Printpi.restype = c_double

OpenSMOKE_NamesOfSpecies = OpenSMOKE.OpenSMOKE_NamesOfSpecies
OpenSMOKE_NamesOfSpecies.argtypes = [c_int]
OpenSMOKE_NamesOfSpecies.restype = c_char_p

OpenSMOKE_NamesOfLiquidSpecies = OpenSMOKE.OpenSMOKE_NamesOfLiquidSpecies
OpenSMOKE_NamesOfLiquidSpecies.argtypes = [c_int]
OpenSMOKE_NamesOfLiquidSpecies.restype = c_char_p

OpenSMOKE_NamesOfSolidSpecies = OpenSMOKE.OpenSMOKE_NamesOfSolidSpecies
OpenSMOKE_NamesOfSolidSpecies.argtypes = [c_int]
OpenSMOKE_NamesOfSolidSpecies.restype = c_char_p

OpenSMOKE_IndexOfSpecies = OpenSMOKE.OpenSMOKE_IndexOfSpecies
OpenSMOKE_IndexOfSpecies.argtypes = [c_char_p]
OpenSMOKE_IndexOfSpecies.restype = c_int

OpenSMOKE_IndexOfSolidSpecies = OpenSMOKE.OpenSMOKE_IndexOfSolidSpecies
OpenSMOKE_IndexOfSolidSpecies.argtypes = [c_char_p]
OpenSMOKE_IndexOfSolidSpecies.restype = c_int

OpenSMOKE_MW = OpenSMOKE.OpenSMOKE_MW
OpenSMOKE_MW.argtypes = [c_int]
OpenSMOKE_MW.restype = c_double

OpenSMOKE_GasProp_SetPressure = OpenSMOKE.OpenSMOKE_GasProp_SetPressure
OpenSMOKE_GasProp_SetPressure.argtypes = [c_double]
OpenSMOKE_GasProp_SetPressure.restype = None

OpenSMOKE_GasProp_SetTemperature = OpenSMOKE.OpenSMOKE_GasProp_SetTemperature
OpenSMOKE_GasProp_SetTemperature.argtypes = [c_double]
OpenSMOKE_GasProp_SetTemperature.restype = None

OpenSMOKE_LiqProp_SetPressure = OpenSMOKE.OpenSMOKE_LiqProp_SetPressure
OpenSMOKE_LiqProp_SetPressure.argtypes = [c_double]
OpenSMOKE_LiqProp_SetPressure.restype = None

OpenSMOKE_LiqProp_SetTemperature = OpenSMOKE.OpenSMOKE_LiqProp_SetTemperature
OpenSMOKE_LiqProp_SetTemperature.argtypes = [c_double]
OpenSMOKE_LiqProp_SetTemperature.restype = None

OpenSMOKE_MW_Solid = OpenSMOKE.OpenSMOKE_MW_Solid
OpenSMOKE_MW_Solid.argtypes = [c_int]
OpenSMOKE_MW_Solid.restype = c_double

OpenSMOKE_GasProp_DynamicViscosity = OpenSMOKE.OpenSMOKE_GasProp_DynamicViscosity
OpenSMOKE_GasProp_DynamicViscosity.argtypes = [POINTER(c_double)]
OpenSMOKE_GasProp_DynamicViscosity.restype = c_double

OpenSMOKE_GasProp_ThermalConductivity = OpenSMOKE.OpenSMOKE_GasProp_ThermalConductivity
OpenSMOKE_GasProp_ThermalConductivity.argtypes = [POINTER(c_double)]
OpenSMOKE_GasProp_ThermalConductivity.restype = c_double

OpenSMOKE_GasProp_Concentrations = OpenSMOKE.OpenSMOKE_GasProp_Concentrations
OpenSMOKE_GasProp_Concentrations.argtypes = [c_double, c_double]
OpenSMOKE_GasProp_Concentrations.restype = c_double

OpenSMOKE_GasProp_Density_IdealGas = OpenSMOKE.OpenSMOKE_GasProp_Density_IdealGas
OpenSMOKE_GasProp_Density_IdealGas.argtypes = [c_double, c_double, c_double]
OpenSMOKE_GasProp_Density_IdealGas.restype = c_double

OpenSMOKE_GasProp_HeatCapacity = OpenSMOKE.OpenSMOKE_GasProp_HeatCapacity
OpenSMOKE_GasProp_HeatCapacity.argtypes = [POINTER(c_double)]
OpenSMOKE_GasProp_HeatCapacity.restype = c_double

OpenSMOKE_GasProp_HeatCapacity_PureSpecies = OpenSMOKE.OpenSMOKE_GasProp_HeatCapacity_PureSpecies
OpenSMOKE_GasProp_HeatCapacity_PureSpecies.argtypes = [POINTER(c_double)]
OpenSMOKE_GasProp_HeatCapacity_PureSpecies.restype = None

OpenSMOKE_GasProp_Dmix = OpenSMOKE.OpenSMOKE_GasProp_Dmix
OpenSMOKE_GasProp_Dmix.argtypes = [POINTER(c_double), c_int]
OpenSMOKE_GasProp_Dmix.restype = c_double

OpenSMOKE_MolecularWeight_From_MoleFractions = OpenSMOKE.OpenSMOKE_MolecularWeight_From_MoleFractions
OpenSMOKE_MolecularWeight_From_MoleFractions.argtypes = [POINTER(c_double)]
OpenSMOKE_MolecularWeight_From_MoleFractions.restype = c_double

OpenSMOKE_MolecularWeight_From_MassFractions = OpenSMOKE.OpenSMOKE_MolecularWeight_From_MassFractions
OpenSMOKE_MolecularWeight_From_MassFractions.argtypes = [POINTER(c_double)]
OpenSMOKE_MolecularWeight_From_MassFractions.restype = c_double

OpenSMOKE_MassFractions_From_MoleFractions = OpenSMOKE.OpenSMOKE_MassFractions_From_MoleFractions
OpenSMOKE_MassFractions_From_MoleFractions.argtypes = [POINTER(c_double), POINTER(c_double), POINTER(c_double)]
OpenSMOKE_MassFractions_From_MoleFractions.restype = None

OpenSMOKE_MoleFractions_From_MassFractions = OpenSMOKE.OpenSMOKE_MoleFractions_From_MassFractions
OpenSMOKE_MoleFractions_From_MassFractions.argtypes = [POINTER(c_double), POINTER(c_double), POINTER(c_double)]
OpenSMOKE_MoleFractions_From_MassFractions.restype = None

OpenSMOKE_GasProp_KineticConstants = OpenSMOKE.OpenSMOKE_GasProp_KineticConstants
OpenSMOKE_GasProp_KineticConstants.argtypes = None
OpenSMOKE_GasProp_KineticConstants.restype = None

OpenSMOKE_GasProp_ReactionRates = OpenSMOKE.OpenSMOKE_GasProp_ReactionRates
OpenSMOKE_GasProp_ReactionRates.argtypes = [POINTER(c_double)]
OpenSMOKE_GasProp_ReactionRates.restype = None

OpenSMOKE_GasProp_FormationRates = OpenSMOKE.OpenSMOKE_GasProp_FormationRates
OpenSMOKE_GasProp_FormationRates.argtypes = [POINTER(c_double)]
OpenSMOKE_GasProp_FormationRates.restype = None

OpenSMOKE_GasProp_HeatRelease = OpenSMOKE.OpenSMOKE_GasProp_HeatRelease
OpenSMOKE_GasProp_HeatRelease.argtypes = [POINTER(c_double)]
OpenSMOKE_GasProp_HeatRelease.restype = c_double

OpenSMOKE_LiqProp_Density_PureSpecies = OpenSMOKE.OpenSMOKE_LiqProp_Density_PureSpecies
OpenSMOKE_LiqProp_Density_PureSpecies.argtypes = [c_char_p, c_double, c_double]
OpenSMOKE_LiqProp_Density_PureSpecies.restype = c_double

OpenSMOKE_LiqProp_Density_Mix_AddVol = OpenSMOKE.OpenSMOKE_LiqProp_Density_Mix_AddVol
OpenSMOKE_LiqProp_Density_Mix_AddVol.argtypes = [c_double, c_double, POINTER(c_double)]
OpenSMOKE_LiqProp_Density_Mix_AddVol.restype = c_double

OpenSMOKE_LiqProp_DynamicViscosity_PureSpecies = OpenSMOKE.OpenSMOKE_LiqProp_DynamicViscosity_PureSpecies
OpenSMOKE_LiqProp_DynamicViscosity_PureSpecies.argtypes = [c_char_p, c_double]
OpenSMOKE_LiqProp_DynamicViscosity_PureSpecies.restype = c_double

OpenSMOKE_LiqProp_DynamicViscosity_Mix = OpenSMOKE.OpenSMOKE_LiqProp_DynamicViscosity_Mix
OpenSMOKE_LiqProp_DynamicViscosity_Mix.argtypes = [c_double, c_char_p]
OpenSMOKE_LiqProp_DynamicViscosity_Mix.restype = c_double

OpenSMOKE_LiqProp_ThermalConductivity_PureSpecies = OpenSMOKE.OpenSMOKE_LiqProp_ThermalConductivity_PureSpecies
OpenSMOKE_LiqProp_ThermalConductivity_PureSpecies.argtypes = [c_char_p, c_double]
OpenSMOKE_LiqProp_ThermalConductivity_PureSpecies.restype = c_double

OpenSMOKE_LiqProp_ThermalConductivity_Mix = OpenSMOKE.OpenSMOKE_LiqProp_ThermalConductivity_Mix
OpenSMOKE_LiqProp_ThermalConductivity_Mix.argtypes = [c_double, POINTER(c_double)]
OpenSMOKE_LiqProp_ThermalConductivity_Mix.restype = c_double

OpenSMOKE_LiqProp_HeatCapacity_PureSpecies = OpenSMOKE.OpenSMOKE_LiqProp_HeatCapacity_PureSpecies
OpenSMOKE_LiqProp_HeatCapacity_PureSpecies.argtypes = [c_char_p, c_double]
OpenSMOKE_LiqProp_HeatCapacity_PureSpecies.restype = c_double

OpenSMOKE_LiqProp_HeatCapacity_Mix = OpenSMOKE.OpenSMOKE_LiqProp_HeatCapacity_Mix
OpenSMOKE_LiqProp_HeatCapacity_Mix.argtypes = [c_double, c_char_p]
OpenSMOKE_LiqProp_HeatCapacity_Mix.restype = c_double

OpenSMOKE_LiqProp_VaporPressure = OpenSMOKE.OpenSMOKE_LiqProp_VaporPressure
OpenSMOKE_LiqProp_VaporPressure.argtypes = [c_char_p, c_double, c_double]
OpenSMOKE_LiqProp_VaporPressure.restype = c_double

OpenSMOKE_LiqProp_VaporizationEnthalpy = OpenSMOKE.OpenSMOKE_LiqProp_VaporizationEnthalpy
OpenSMOKE_LiqProp_VaporizationEnthalpy.argtypes = [c_char_p, c_double]
OpenSMOKE_LiqProp_VaporizationEnthalpy.restype = c_double

OpenSMOKE_LiqProp_Dmix_PerkinsGeankopolis = OpenSMOKE.OpenSMOKE_LiqProp_Dmix_PerkinsGeankopolis
OpenSMOKE_LiqProp_Dmix_PerkinsGeankopolis.argtypes = [c_double, c_double, POINTER(c_double), c_int]
OpenSMOKE_LiqProp_Dmix_PerkinsGeankopolis.restype = c_double

OpenSMOKE_LiqProp_Dmix_Cullinan = OpenSMOKE.OpenSMOKE_LiqProp_Dmix_Cullinan
OpenSMOKE_LiqProp_Dmix_Cullinan.argtypes = [c_double, c_double, POINTER(c_double), c_int]
OpenSMOKE_LiqProp_Dmix_Cullinan.restype = c_double

OpenSMOKE_LiqProp_Dmix_LefflerCullinan = OpenSMOKE.OpenSMOKE_LiqProp_Dmix_LefflerCullinan
OpenSMOKE_LiqProp_Dmix_LefflerCullinan.argtypes = [c_double, c_double, POINTER(c_double), c_int]
OpenSMOKE_LiqProp_Dmix_LefflerCullinan.restype = c_double

OpenSMOKE_LiquidMolecularWeight_From_LiquidMoleFractions = OpenSMOKE.OpenSMOKE_LiquidMolecularWeight_From_LiquidMoleFractions
OpenSMOKE_LiquidMolecularWeight_From_LiquidMoleFractions.argtypes = [POINTER(c_double)]
OpenSMOKE_LiquidMolecularWeight_From_LiquidMoleFractions.restype = c_double

OpenSMOKE_LiquidMolecularWeight_From_LiquidMassFractions = OpenSMOKE.OpenSMOKE_LiquidMolecularWeight_From_LiquidMassFractions
OpenSMOKE_LiquidMolecularWeight_From_LiquidMassFractions.argtypes = [POINTER(c_double)]
OpenSMOKE_LiquidMolecularWeight_From_LiquidMassFractions.restype = c_double

OpenSMOKE_LiquidMassFractions_From_LiquidMoleFractions = OpenSMOKE.OpenSMOKE_LiquidMassFractions_From_LiquidMoleFractions
OpenSMOKE_LiquidMassFractions_From_LiquidMoleFractions.argtypes = [POINTER(c_double), POINTER(c_double), POINTER(c_double)]
OpenSMOKE_LiquidMassFractions_From_LiquidMoleFractions.restype = None

OpenSMOKE_LiquidMoleFractions_From_LiquidMassFractions = OpenSMOKE.OpenSMOKE_LiquidMoleFractions_From_LiquidMassFractions
OpenSMOKE_LiquidMoleFractions_From_LiquidMassFractions.argtypes = [POINTER(c_double), POINTER(c_double), POINTER(c_double)]
OpenSMOKE_LiquidMoleFractions_From_LiquidMassFractions.restype = None

OpenSMOKE_CheckAndCorrectSumOfFractions = OpenSMOKE.OpenSMOKE_CheckAndCorrectSumOfFractions
OpenSMOKE_CheckAndCorrectSumOfFractions.argtypes = [POINTER(c_double), c_int]
OpenSMOKE_CheckAndCorrectSumOfFractions.restype = None

OpenSMOKE_GetMixtureFractionFromMassFractions = OpenSMOKE.OpenSMOKE_GetMixtureFractionFromMassFractions
OpenSMOKE_GetMixtureFractionFromMassFractions.argtypes = [POINTER(c_double), POINTER(c_double), POINTER(c_double)]
OpenSMOKE_GetMixtureFractionFromMassFractions.restype = c_double

c_odefunction = CFUNCTYPE(None, POINTER(c_double), c_double, POINTER(c_double), c_void_p)
OpenSMOKE_ODESolver = OpenSMOKE.OpenSMOKE_ODESolver
OpenSMOKE_ODESolver.argtypes = [POINTER(c_odefunction), c_int, c_double, POINTER(c_double), c_void_p]

OpenSMOKE_ReadSolidKinetics = OpenSMOKE.OpenSMOKE_ReadSolidKinetics
OpenSMOKE_ReadSolidKinetics.argtypes = [c_char_p]
OpenSMOKE_ReadSolidKinetics.restype = None

OpenSMOKE_NumberOfSolidSpecies = OpenSMOKE.OpenSMOKE_NumberOfSolidSpecies
OpenSMOKE_NumberOfSolidSpecies.argtypes = None
OpenSMOKE_NumberOfSolidSpecies.restype = c_int

OpenSMOKE_NumberOfSolidReactions = OpenSMOKE.OpenSMOKE_NumberOfSolidReactions
OpenSMOKE_NumberOfSolidReactions.argtypes = None
OpenSMOKE_NumberOfSolidReactions.restype = c_int

OpenSMOKE_SolProp_SetPressure = OpenSMOKE.OpenSMOKE_SolProp_SetPressure
OpenSMOKE_SolProp_SetPressure.argtypes = [c_double]
OpenSMOKE_SolProp_SetPressure.restype = None

OpenSMOKE_SolProp_SetTemperature = OpenSMOKE.OpenSMOKE_SolProp_SetTemperature
OpenSMOKE_SolProp_SetTemperature.argtypes = [c_double]
OpenSMOKE_SolProp_SetTemperature.restype = None

OpenSMOKE_SolProp_KineticConstants = OpenSMOKE.OpenSMOKE_SolProp_KineticConstants
OpenSMOKE_SolProp_KineticConstants.argtypes = None
OpenSMOKE_SolProp_KineticConstants.restype = None

OpenSMOKE_SolProp_ReactionRates = OpenSMOKE.OpenSMOKE_SolProp_ReactionRates
OpenSMOKE_SolProp_ReactionRates.argtypes = [POINTER(c_double), POINTER(c_double)]
OpenSMOKE_SolProp_ReactionRates.restype = None

OpenSMOKE_SolProp_FormationRates = OpenSMOKE.OpenSMOKE_SolProp_FormationRates
OpenSMOKE_SolProp_FormationRates.argtypes = [POINTER(c_double), POINTER(c_double)]
OpenSMOKE_SolProp_FormationRates.restype = None

OpenSMOKE_SolProp_HeatRelease = OpenSMOKE.OpenSMOKE_SolProp_HeatRelease
OpenSMOKE_SolProp_HeatRelease.argtypes = [POINTER(c_double), POINTER(c_double)]
OpenSMOKE_SolProp_HeatRelease.restype = c_double

OpenSMOKE_MW_Solid = OpenSMOKE.OpenSMOKE_MW_Solid
OpenSMOKE_MW_Solid.argtypes = [c_int]
OpenSMOKE_MW_Solid.restype = c_double

OpenSMOKE_NamesOfSolidSpecies = OpenSMOKE.OpenSMOKE_NamesOfSolidSpecies
OpenSMOKE_NamesOfSolidSpecies.argtypes = [c_int]
OpenSMOKE_NamesOfSolidSpecies.restype = c_char_p

OpenSMOKE_IndexOfSolidSpecies = OpenSMOKE.OpenSMOKE_IndexOfSolidSpecies
OpenSMOKE_IndexOfSolidSpecies.argtypes = [c_char_p]
OpenSMOKE_IndexOfSolidSpecies.restype = c_int

OpenSMOKE_SolProp_HeatCapacity = OpenSMOKE.OpenSMOKE_SolProp_HeatCapacity
OpenSMOKE_SolProp_HeatCapacity.argtypes = [POINTER(c_double)]
OpenSMOKE_SolProp_HeatCapacity.restype = c_double
