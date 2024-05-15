#!/usr/bin/env python3
#---------------------------------------------------------------------
#import system lib
#---------------------------------------------------------------------
import os,sys

Libs=['pylib','mylib','schism_file']
if not set(Libs).issubset(set(sys.modules.keys())):
   #---------------------------------------------------------------------
   #spyder pylab library
   #C:/Program Files/Python/Python3[57]/Lib/site-packages/matplotlib/pylab.py
   #---------------------------------------------------------------------
   pversion=sys.version.split(' ')[0]
   #print(pversion)

   #---------------------------------------------------------------------
   #libraries of packages
   #---------------------------------------------------------------------
   #matplotlib
   import matplotlib as mpl
   HNAME=str(os.getenv('HOSTNAME')); TNAME=str(os.getenv('TACC_SYSTEM'))
   # if ('frontera' in HNAME) or ('stampede2' in HNAME): mpl.use('tkagg')
   # if ('frontera' in TNAME) or ('stampede2' in TNAME): mpl.use('tkagg')
   mpl.use('agg')
   from matplotlib import pyplot as plt
   from matplotlib import cbook, mlab
   from matplotlib.dates import *
   if hasattr(mpl.dates,'set_epoch'):
      try:
         mpl.dates.set_epoch('0000-12-31')
      except:
         pass
   from matplotlib.pyplot import *

   import platform
   if platform.system().lower()=='windows':
      try:
         if get_ipython().__class__.__name__!='ZMQInteractiveShell': mpl.use('Qt5Agg')
      except:
         pass

   #numpy
   import numpy as np
   from numpy import *
   from numpy.random import *
   from numpy.linalg import *
   import numpy.ma as ma
   #from numpy.fft import *

   #scipy
   import scipy as sp
   from scipy import interpolate
   from scipy.fftpack import fft, ifft
   #from scipy import (optimize,interpolate,io,signal)

   #pandas
   import pandas as pd

   #misc
   import re
   import datetime
   #from io import StringIO
   #import imp
   #import importlib as imp

   #proj
   from pyproj import Transformer
   #from pyproj import Proj, transform

   #netcdf
   from netCDF4 import Dataset

   #excel
   try:
      import xlsxwriter as xw
   except:
      pass

  #pickle
   try:
       import pickle
       import copy
       from copy import copy as scopy
       from copy import deepcopy as dcopy
   except:
       pass

   #mpi4py
   try:
      from mpi4py import MPI
   except:
       pass

   #url download
   try:
      import urllib
      from urllib.request import urlretrieve as urlsave
      import ssl
      try:
          _create_unverified_https_context = ssl._create_unverified_context
      except AttributeError:
          # Legacy Python that doesn't verify HTTPS certificates by default
          pass
      else:
          # Handle target environment that doesn't support HTTPS verification
          ssl._create_default_https_context = _create_unverified_https_context
   except:
      pass

   #reload
   try:
     from importlib import reload
   except:
     pass

   #sympy
   #try:
   #  from sympy import init_session as sym_init
   #except:
   #  pass

   #---------------------------------------------------------------------
   #libraries of self-defined modules
   #---------------------------------------------------------------------
   #import mylib as mylib
   import mylib    

   sys.modules['mylib'] = mylib 
   from mylib import (get_xtick,close_data_loop,datenum,
        loadz,zdata,savez,find_cs,convert_matfile,
        smooth,daytime_length,move_figure,lpfilt,mdivide,signa,
        inside_polygon,command_outputs,near_pts,proj,proj_pts,rewrite,rewrite_input,
        get_prj_file,mfft,read_shapefile_data,write_shapefile_data,
        ReadNC,WriteNC,harmonic_fit,harmonic_analysis,get_hycom,compute_contour,
        get_stat,get_subplot_position,get_subplot_position2,load_bathymetry,plot_taylor_diagram,
        convert_dem_format,get_hpc_command,least_square_fit,read_yaml,read_excel,
        write_excel,rtext)

   import schism_file as schism_file
   sys.modules['schism_file'] = schism_file
   from schism_file import (read_schism_hgrid, read_schism_bpfile,getglob,
        schism_grid,schism_vgrid,schism_bpfile,sms2grd,read_schism_vgrid,save_schism_grid,
        compute_zcor,read_schism_param,write_schism_param,read_schism_local_to_global,
        create_schism_vgrid,srank,grd2sms,scatter_to_schism_grid,delete_schism_grid_element,
        read_schism_prop,read_schism_reg,interp_schism_3d)

   if os.getenv('HOME')!=None:
       sys.path.append(os.getenv('HOME'))

   #sys.modules['loadz'] = mylib #in case oldmodule name used
   #sys.modules['read_schism_file'] = schism_file #in case oldmodule name used

   #import mpas_file
   #from mpas_file import (read_mpas_grid)

   #---------------------------------------------------------------------
   #alias
   #---------------------------------------------------------------------
   from os.path import exists as fexist
   from mylib import savez as save_npz; mylib.save_npz=savez
   from mylib import zdata as npz_data; mylib.npz_data=zdata
   from mylib import least_square_fit as lsq; mylib.least_square_fit=lsq
   from mylib import move_figure as mvfig
   from mylib import find_cs as find_continuous_sections; mylib.find_continuous_sections=find_cs
