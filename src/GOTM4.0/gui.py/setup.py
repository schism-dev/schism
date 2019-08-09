import sys

try:
    import modulefinder
    #import win32com
    #for p in win32com.__path__[1:]:
    #    modulefinder.AddPackagePath("win32com", p)
    for extra in ["win32com","win32com.shell"]: #,"win32com.mapi"
        __import__(extra)
        m = sys.modules[extra]
        for p in m.__path__[1:]:
            modulefinder.AddPackagePath(extra, p)
except ImportError:
    # no build path setup, no worries.
    pass

from distutils.core import setup
import py2exe

from distutils.filelist import findall
import os, os.path
import matplotlib

def adddir(path,localtarget=None):
    if localtarget==None: localtarget = path
    for f in findall(path):
        localname = os.path.join(localtarget, f[len(path)+1:])
        if 'CVS' in localname: continue
        own_data_files.append((os.path.dirname(localname),[f]))

own_data_files = []

own_data_files.append(('',['logo.png']))
own_data_files.append(('',['C:\Program Files\Python24\MSVCP71.dll']))

adddir(matplotlib.get_data_path(),'matplotlibdata')
adddir('defaultscenarios')
adddir('reporttemplates')
adddir('schemas')

setup(
    console=['gotm.py'],
    options={'py2exe': {
                'packages' : ['matplotlib', 'pytz'],
                'includes' : ['sip'],
                'excludes' : ['_gtkagg', '_tkagg', '_wxagg','Tkconstants','Tkinter','tcl','wx'],
                'dll_excludes': ['libgdk-win32-2.0-0.dll', 'libgobject-2.0-0.dll', 'libgdk_pixbuf-2.0-0.dll','wxmsw26uh_vc.dll','tcl84.dll','tk84.dll'],
            }},
    data_files=own_data_files
)

