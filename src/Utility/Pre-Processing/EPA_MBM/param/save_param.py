#!/usr/bin/env python3
from pylib import *
close("all")

#save setup parameters
C=zdata()
C.hydro=read_schism_param('param.nml',1)
C.sed=read_schism_param('sediment.nml',1)
C.wave=read_schism_param('wwminput.nml',1)
C.icm=read_schism_param('icm.nml',1)
C.hydro_sed=read_schism_param('old/param_sed.nml',1) #outdated for backward accommendation
savez('param',C)
