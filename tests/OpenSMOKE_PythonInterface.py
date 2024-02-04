import os, sys
sys.path.append(os.environ.get('OPENSMOKE_INTERFACE')+'/python_interface')
from pyopensmoke import *

kinfolder = "../kinetics/skeletal/methanol/kinetics"

OpenSMOKE_Init()
OpenSMOKE_ReadKinetics (kinfolder.encode('utf-8'))
OpenSMOKE_Clean()
