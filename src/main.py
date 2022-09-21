
from ctypes import *

OpenSMOKE = cdll.LoadLibrary("libopensmoke.so")

OpenSMOKE.OpenSMOKE_ReadKinetics()
OpenSMOKE.OpenSMOKE_ReadLiquidKinetics()
OpenSMOKE.OpenSMOKE_ReadLiquidProperties()

OpenSMOKE.OpenSMOKE_GasProp_SetTemperature(c_double(773.))
OpenSMOKE.OpenSMOKE_GasProp_SetPressure(c_double(101325.))

NGS = OpenSMOKE.OpenSMOKE_NumberOfSpecies()
NGR = OpenSMOKE.OpenSMOKE_NumberOfReactions()

print("NGS = ", NGS)
print("NGR = ", NGR)

