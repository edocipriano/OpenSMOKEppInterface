import os, sys
sys.path.append(os.environ.get('OPENSMOKE_INTERFACE')+'/python_interface')
from pyopensmoke import *

kinfolder = "../kinetics/biomass/solid-only-2003/kinetics"

OpenSMOKE_Init()
OpenSMOKE_ReadSolidKinetics (kinfolder.encode('utf-8'))

print ("OpenSMOKE_NumberOfSolidSpecies = ", OpenSMOKE_NumberOfSolidSpecies(), file=sys.stderr)
print ("OpenSMOKE_NumberOfSolidReactions = ", OpenSMOKE_NumberOfSolidReactions(), file=sys.stderr)

for i in range (1, OpenSMOKE_NumberOfSolidSpecies()):
  print ("Index = ", i, " is ", OpenSMOKE_NamesOfSolidSpecies (i).decode ('utf-8'), file=sys.stderr)

for i in range (1, OpenSMOKE_NumberOfSolidSpecies()):
  species = OpenSMOKE_NamesOfSolidSpecies (i)
  print ("Species: ", species.decode ('utf-8'), " has index ", OpenSMOKE_IndexOfSolidSpecies (species), file=sys.stderr)

for i in range (1, OpenSMOKE_NumberOfSolidSpecies()):
  species = OpenSMOKE_NamesOfSolidSpecies (i)
  print ("MW[", species.decode ('utf-8'), "] = ", OpenSMOKE_MW_Solid (i), file=sys.stderr)

OpenSMOKE_Clean()
