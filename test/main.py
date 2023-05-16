
from ctypes import *

OpenSMOKE = cdll.LoadLibrary("../build/libopensmoke.dylib")

OpenSMOKE.OpenSMOKE_Init()
OpenSMOKE.OpenSMOKE_ReadKinetics("../kinetics/one-step/n-heptane/kinetics")
OpenSMOKE.OpenSMOKE_ReadLiquidKinetics("../kinetics/one-step/n-heptane/kinetics")
OpenSMOKE.OpenSMOKE_ReadLiquidProperties("../kinetics/LiquidProperties/LiquidProperties")

OpenSMOKE.OpenSMOKE_GasProp_SetTemperature(c_double(773.))
OpenSMOKE.OpenSMOKE_GasProp_SetPressure(c_double(101325.))

NGS = OpenSMOKE.OpenSMOKE_NumberOfSpecies()
NGR = OpenSMOKE.OpenSMOKE_NumberOfReactions()

print(type(NGS))

print("NGS = ", NGS)
print("NGR = ", NGR)

