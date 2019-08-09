import os, re, datetime, xml.dom.minidom, tempfile, shutil, StringIO

import common, xmlstore, scenario

# Import NetCDF file format support
#import pycdf
from pynetcdf import NetCDFFile
import matplotlib.numerix

# Abstract class that contains one or more variables that can be plotted.
# Classes deriving from it must support the virtual methods below.
class PlotVariableStore:

    def __init__(self):
        pass

    def getVariableNames(self):
        return []

    def getVariableLongNames(self):
        return dict([(name,self.getVariable(name).getLongName()) for name in self.getVariableNames()])

    def getVariable(self,varname):
        return None

    def getVariableTree(self,path,otherstores={}):
        xmlschema = xml.dom.minidom.parse(path)
        vardict = self.getVariableLongNames()
        found = set(self.filterNodes(xmlschema.documentElement,vardict))
        remaining = [nodename for nodename in vardict if nodename not in found]
        other = None
        for ch in xmlschema.getElementsByTagName('element'):
            if ch.getAttribute('id')=='other':
                other = ch
                break
        if other!=None:
            if len(remaining)==0:
                other.parentNode.removeChild(other)
            else:
                for varid in sorted(remaining,cmp=lambda x,y: cmp(vardict[x].lower(), vardict[y].lower())):
                    el = xmlschema.createElement('element')
                    el.setAttribute('id',varid)
                    el.setAttribute('label',vardict[varid])
                    el.setAttribute('type','bool')
                    other.appendChild(el)
        return xmlstore.TypedStore(xmlschema,otherstores=otherstores)

    def filterNodes(self,node,vardict):
        nodeid = node.getAttribute('id')
        assert nodeid!='', 'Node lacks "id" attribute.'
        nodeids = []
        if nodeid in vardict:
            if not node.hasAttribute('label'):
                node.setAttribute('label',vardict[nodeid])
            node.setAttribute('type','bool')
            nodeids.append(nodeid)
        for ch in common.findDescendantNodes(node,['element']):
            nodeids += self.filterNodes(ch,vardict)
        if len(nodeids)==0 and nodeid!='other':
            node.parentNode.removeChild(node)
        return nodeids

# Abstract class that represents a variable that can be plotted.
# Classes deriving from it must support the virtual methods below.
class PlotVariable:

    def __init__(self):
        pass

    # getName()
    #   Return type: string
    #   Returns the short name (or identifier) of the variable.
    def getName(self):
        return ''

    # getLongName()
    #   Return type: string
    #   Returns the pretty name for the variable.
    def getLongName(self):
        return ''

    # getUnit()
    #   Return type: string
    #   Returns the unit of the variable.
    def getUnit(self):
        return ''

    # getDimensions()
    #   Return type: tuple of strings
    #   Returns the names of the dimensions of the variable; currently supported dimensions: "time", "z".
    def getDimensions(self):
        return ()

    # getValues(bounds, staggered=False)
    #   Return type: tuple of (numpy/numeric) arrays
    #   Returns the arrays with coordinates (in order of dimnesions returned by getDimensions), and the value array.
    #   Coordinates may be given as 1D arrays (if the coordinates are constant across all other dimensions), or as arrays with
    #   the same numbers of dimensions as the value array (i.e., for every value a coordinate is supplied). If staggered is True,
    #   coordinates must be given at the interfaces and values at the centers; then coordinate arrays must have one extra value in
    #   every dimension compared to the value array.
    def getValues(self,bounds,staggered=False,coordinatesonly=False):
        return ()

class LinkedFileVariableStore(PlotVariableStore):

    class LinkedFileVariable(PlotVariable):

        def __init__(self,store,data,index):
            self.store = store
            self.data = data
            self.index = index

        def getName(self):
            return self.data[0]

        def getLongName(self):
            return self.data[1]

        def getUnit(self):
            return self.data[2]

        def getDimensions(self):
            filetype = self.store.type
            if filetype=='profilesintime':
                return ('time','z')
            elif filetype=='pointsintime':
                return ('time',)
            else:
                assert False, 'Cannot plot variables from file of unknown type "%s".' % filetype

        def getValues(self,bounds,staggered=False,coordinatesonly=False):
            data = self.store.getData()
            res = []

            timebounds = common.findindices(bounds[0],data[0])
            res.append(data[0][timebounds[0]:timebounds[1]+1])
            if len(data)==2:
                res.append(data[1][timebounds[0]:timebounds[1]+1,self.index])
            elif len(data)==3:
                res.append(data[1])
                res.append(data[2][timebounds[0]:timebounds[1]+1,:,self.index])
            else:
                assert False, 'Cannot handle variables with %i dimensions; only know how to deal with 2 or 3 dimensions.' % len(data)
                
            if staggered:
                for idim in range(len(res)-1):
                    delta = (res[idim][1:]-res[idim][:-1])/2
                    res[idim] = matplotlib.numerix.concatenate(([res[idim][0]-delta[0]],res[idim][:-1]+delta,[res[idim][-1]+delta[-1]]),0)
                    
            return res

    def __init__(self,node,datafile=None):
        self.vardata = []

        if datafile==None: datafile = node.getValue()
        self.datafile = datafile

        finfo = common.findDescendantNode(node.templatenode,['fileinfo'])
        self.type = finfo.getAttribute('type')
        
        fvars = common.findDescendantNode(finfo,['filevariables'])
        if fvars!=None:
            for ch in fvars.childNodes:
                if ch.nodeType==ch.ELEMENT_NODE and ch.localName=='filevariable':
                    longname = ch.getAttribute('label')
                    unit = ch.getAttribute('unit')
                    name = longname
                    self.vardata.append([name,longname,unit])

        self.data = None

    def getVariableNames(self):
        if not self.datafile.isValid(): return []
        return [data[0] for data in self.vardata]

    def getVariableLongNames(self):
        if not self.datafile.isValid(): return {}
        return dict([(data[0],data[1]) for data in self.vardata])

    def getVariable(self,varname):
        if not self.datafile.isValid(): return None
        for (index,data) in enumerate(self.vardata):
            if data[0]==varname:
                return self.LinkedFileVariable(self,data,index)
        assert False, 'Variable with name "%s" not found in store.' % varname

    def getData(self,callback=None):
        if self.data!=None: return self.data
        if not self.datafile.isValid(): return None
        
        varcount = len(self.vardata)
        
        # Get the size of the file (in bytes, may be None if the size is not known)
        # This will be used in combination with the position of the file pointer to report progress.
        filesize = self.datafile.getSize()
        
        # Access the data through some read-only file-like object.
        f = self.datafile.getAsReadOnlyFile()

        # Compile regular expression for reading dates.
        datetimere = re.compile('(\d\d\d\d).(\d\d).(\d\d) (\d\d).(\d\d).(\d\d)')

        filetype = self.type
        if filetype=='pointsintime':
            times = []
            values = []
            iline = 0
            while True:
                # Read a line (stop if end-of-file was reached)
                line = f.readline()
                if line=='': break
                iline += 1
                
                # Read the date + time
                datematch = datetimere.match(line)
                if datematch==None:
                    raise Exception('Line %i does not start with time (yyyy-mm-dd hh:mm:ss). Line contents: %s' % (iline,line))
                refvals = map(int,datematch.group(1,2,3,4,5,6)) # Convert matched strings into integers
                curdate = datetime.datetime(*refvals)
                
                # Read values.
                data = line[datematch.end()+1:].split()
                if len(data)<varcount:
                    raise Exception('Line %i contains only %i observations, where %i are expected.' % (iline,len(data),varcount))
                data = map(float,data)
                
                # Store time and values.
                times.append(curdate)
                values.append(data)
                
                # Inform caller about progress
                if callback!=None and iline%1000==0:
                    progress = None
                    if filesize!=None:
                        try:
                            progress = float(f.tell())/filesize
                        except:
                            progress = None
                    callback('processed %i lines.' % iline,progress=progress)
            
            # Convert sequences to numpy arrays and store these.
            times = matplotlib.numerix.array(times,matplotlib.numerix.PyObject)
            values = matplotlib.numerix.array(values,matplotlib.numerix.Float32)
            self.data = (times,values)
        elif filetype=='profilesintime':
            times = []
            depths = []
            values = []
            uniquedepths = {}
            iline = 0
            while True:
                # Read a line (stop if end-of-file was reached)
                line = f.readline()
                if line=='': break
                iline += 1
                
                # Read date & time
                datematch = datetimere.match(line)
                if datematch==None:
                    raise Exception('Line %i does not start with time (yyyy-mm-dd hh:mm:ss). Line contents: %s' % (iline,line))
                refvals = map(int,datematch.group(1,2,3,4,5,6)) # Convert matched strings into integers
                curdate = datetime.datetime(*refvals)

                # Get the number of observations and the depth direction.
                (depthcount,updown) = map(int, line[datematch.end()+1:].split())

                # Create arrays that will contains depths and observed values.
                curdepths = matplotlib.numerix.empty((depthcount,),matplotlib.numerix.Float32)
                curvalues = matplotlib.numerix.empty((depthcount,varcount),matplotlib.numerix.Float32)
                
                # Depths can be increasing (updown==1) or decreasing (updown!=1)
                if updown==1:
                    depthindices = range(0,depthcount,1)
                else:
                    depthindices = range(depthcount-1,-1,-1)
                
                # Now parse the specified number of observations to cretae the profiles.
                prevdepth = None
                for idepthline in depthindices:
                    if callback!=None and iline%1000==0:
                        pos = f.tell()
                        callback('processed %i lines.' % iline,progress=0.5*pos/filesize)
                        
                    # Read line
                    line = f.readline()
                    if line=='':
                        raise Exception('Premature end-of-file after line %i; expected %i more observations.' % (iline,depthcount-depthindices.index(idepthline)))
                    iline += 1
                    
                    # Read values (depth followed by data) and check.
                    linedata = map(float,line.split())
                    if len(linedata)<varcount+1:
                        raise Exception('Line %i contains only %i value(s), where %i (1 time and %i observations) are expected.' % (iline,len(linedata),varcount+1,varcount))
                    if prevdepth!=None:
                        if linedata[0]==prevdepth:
                            raise Exception('Found duplicate observation for depth %.4f at line %i.' % (linedata[0],iline))
                        if updown==1:
                            if linedata[0]<prevdepth:
                                raise Exception('Observation depth decreases from %.4f to %.4f at line %i, but the profile depth was set to increase from first to last observation.' % (prevdepth,linedata[0],iline))
                        elif linedata[0]>prevdepth:
                            raise Exception('Observation depth increases from %.4f to %.4f at line %i, but the profile depth was set to decrease from first to last observation.' % (prevdepth,linedata[0],iline))
                    prevdepth = linedata[0]
                    
                    # Store current observation
                    uniquedepths[linedata[0]] = True
                    curdepths[idepthline] = linedata[0]
                    curvalues[idepthline,:] = linedata[1:varcount+1]
                    
                # Append the profiles for the current time to the list.
                times.append(curdate)
                depths.append(curdepths)
                values.append(curvalues)
                
                # Inform caller about progress.
                if callback!=None and iline%1000==0:
                    pos = f.tell()
                    callback('processed %i lines.' % iline,progress=0.5*pos/filesize)
                    
            # Convert sequence with times to numpy array.
            times = matplotlib.numerix.array(times,matplotlib.numerix.PyObject)
            
            # Create depth grid to interpolate on to. Use the observation depths if less than 200,
            # otherwise create a equidistant 200-point grid between the minimum and maximum depth.
            uniquedepths = uniquedepths.keys()
            uniquedepths.sort()
            if len(uniquedepths)<200:
                depthgrid = matplotlib.numerix.array(uniquedepths,matplotlib.numerix.Float32)
            else:
                depthstep = (uniquedepths[-1]-uniquedepths[0])/200
                depthgrid = matplotlib.numerix.arange(uniquedepths[0],uniquedepths[-1]+depthstep,depthstep)
                
            # Grid observed profiles to depth grid.
            griddedvalues = matplotlib.numerix.empty((times.shape[0],depthgrid.shape[0],varcount),matplotlib.numerix.Float32)
            for it in range(len(times)):
                griddedvalues[it,:,:] = common.interp1(depths[it],values[it],depthgrid)
                if callback!=None and (it+1)%20==0:
                    callback('gridded %i profiles.' % (it+1),progress=.5+float(it+1)/len(times)/2)
                
            # Store time grid, depth grid and observations.
            self.data = (times,depthgrid,griddedvalues)
        else:
            assert False, 'Cannot plot variables from file of unknown type "%s".' % filetype
            
        # Close data file
        f.close()
        
        # Return data
        return self.data

# Class that represents a GOTM result.
#   Inherits from PlotVariableStore, as it contains variables that can be plotted.
#   Contains a link to the scenario from which the result was created (if available)
class Result(PlotVariableStore):
    reportdirname = 'reporttemplates'
    reportname2path = None

    @staticmethod
    def getReportTemplates():
        if Result.reportname2path==None:
            Result.reportname2path = {}
            #sourcedir = os.path.join(os.path.dirname(__file__),Result.reportdirname)
            sourcedir = Result.reportdirname
            if os.path.isdir(sourcedir):
                for filename in os.listdir(sourcedir):
                    if filename=='CVS': continue
                    fullpath = os.path.join(sourcedir,filename)
                    if os.path.isdir(fullpath):
                        if os.path.isfile(os.path.join(fullpath,'index.xml')):
                            Result.reportname2path[filename] = fullpath
                        else:
                            print 'WARNING: template directory "%s" does not contain "index.xml"; it will be ignored.' % fullpath
                    else:
                        print 'WARNING: template directory "%s" contains "%s" which is not a directory; the latter will be ignored.' % (sourcedir,filename)
            else:
                print 'WARNING: no report templates will be available, because subdirectory "%s" is not present!' % Result.reportdirname
        return Result.reportname2path

    class ResultVariable(PlotVariable):
        def __init__(self,result,varname):
            PlotVariable.__init__(self)
            self.result = result
            self.varname = str(varname)

        def getName(self):
            return self.varname

        def getLongName(self):
            nc = self.result.getcdf()
            return nc.variables[self.varname].long_name

        def getUnit(self):
            nc = self.result.getcdf()
            return nc.variables[self.varname].units

        def getDimensions(self):
          nc = self.result.getcdf()

          dimnames = nc.variables[self.varname].dimensions
          dimcount = len(dimnames)
          if   dimcount==3:
              if dimnames==('time','lat','lon'):
                  return ('time',)
          elif dimcount==4:
              if (dimnames==('time','z','lat','lon')) or (dimnames==('time','z1','lat','lon')):
                  return ('time','z')
          raise Exception('This variable has dimensions %s; I do not know how to handle such variables.' % str(dimnames))

        def getValues(self,bounds,staggered=False,coordinatesonly=False):
          nc = self.result.getcdf()
            
          v = nc.variables[self.varname]
          dims = v.dimensions
          dimcount = len(dims)
          res = []

          # Get time coordinates and time bounds.
          (t,t_stag) = self.result.getTime()
          timebounds = common.findindices(bounds[0],t)
          if not staggered:
              res.append(t[timebounds[0]:timebounds[1]+1])
          else:
              res.append(t_stag[timebounds[0]:timebounds[1]+2])

          # Add depth dimension (if available)
          if dimcount==4:
            (z,z1,z_stag,z1_stag) = self.result.getDepth()
            depthbounds = (0,z.shape[1])
            if dims[1]=='z':
                if staggered:
                    res.append(z_stag[timebounds[0]:timebounds[1]+2,depthbounds[0]:depthbounds[1]+2])
                else:
                    res.append(z[timebounds[0]:timebounds[1]+1,depthbounds[0]:depthbounds[1]+1])
            elif dims[1]=='z1':
                if staggered:
                    res.append(z1_stag[timebounds[0]:timebounds[1]+2,depthbounds[0]:depthbounds[1]+2])
                else:
                    res.append(z1[timebounds[0]:timebounds[1]+1,depthbounds[0]+1:depthbounds[1]+2])

          # Add value data, if needed.
          if not coordinatesonly:
            if dimcount==4:
                # Four-dimensional variable: longitude, latitude, depth, time
                try:
                    dat = v[timebounds[0]:timebounds[1]+1,depthbounds[0]:depthbounds[1]+1,0,0]
                except Exception, e:
                    raise Exception('Unable to read values for NetCDF variable "%s". Error: %s' % (self.varname,str(e)))
                res.append(dat)
            elif dimcount==3:
                # Three-dimensional variable: longitude, latitude, time
                try:
                    dat = v[timebounds[0]:timebounds[1]+1,0,0]
                except Exception, e:
                    raise Exception('Unable to read values for NetCDF variable "%s". Error: %s' % (self.varname,str(e)))
                res.append(dat)
            else:
                raise Exception('This variable has dimensions %s; I do not know how to handle such variables.' % unicode(dimensions))

            # Mask missing data.
            if hasattr(v,'missing_value'):
                res[-1] = matplotlib.numerix.ma.masked_array(res[-1],res[-1]==v.missing_value)
              
          return res

    def __init__(self):
        PlotVariableStore.__init__(self)
        
        self.scenario = None
        self.tempdir = None
        self.datafile = None
        self.nc = None
        self.changed = False
        self.gotmoutput = None

        # Cached coordinates
        self.t = None
        self.t_stag = None
        self.z = None
        self.z1 = None
        self.z_stag = None
        self.z1_stag = None
        
        self.store = xmlstore.TypedStore('schemas/result/gotmgui.xml')
        
        self.path = None

    def hasChanged(self):
        return self.changed or self.store.changed

    def getTempDir(self,empty=False):
        if self.tempdir!=None:
            if empty:
                for f in os.listdir(self.tempdir): 
                    os.remove(os.path.join(self.tempdir,f))
        else:
            self.tempdir = tempfile.mkdtemp('','gotm-')
            print 'Created temporary result directory "%s".' % self.tempdir
        return self.tempdir
        
    def saveNetCDF(self,path):
        shutil.copyfile(self.datafile,path)

    def save(self,path,addfiguresettings=True):
        assert self.datafile!=None, 'The result object was not yet attached to a result file (NetCDF).'

        # Create a ZIP container to hold the result.
        container = xmlstore.DataContainerZip(path,'w')

        if not addfiguresettings:
            # First clear all figure settings.
            self.store.root.getLocation(['FigureSettings']).clearValue(recursive=True)

        # Add the XML file describing result properties.            
        df = xmlstore.DataFileXmlNode(self.store.root.valuenode)
        df_added = container.addItem(df,'result.xml')
        df_added.release()
        df.release()
        self.store.resetChanged()

        # If we have a link to the scenario, add it to the result file.
        if self.scenario!=None:
            fscen = StringIO.StringIO()
            self.scenario.saveAll(fscen,claim=False)
            df = xmlstore.DataFileMemory(fscen.getvalue(),'scenario.gotmscenario')
            fscen.close()
            container.addItem(df)
            df.release()
        
        # Add the result data (NetCDF)
        container.addFile(self.datafile,'result.nc')
        
        # Make changes to container persistent (this closes the ZIP file), and release it.
        container.persistChanges()
        container.release()

        # Saved all changes; reset "changed" state
        self.changed = False
        
        # Store the path we saved to.
        self.path = path

    def load(self,path):
        # Basic check: does the specified file exist?
        if not os.path.exists(path): raise Exception('File "%s" does not exist.' % path)

        container = xmlstore.DataContainerZip(path,'r')
        files = container.listFiles()
        if 'scenario.gotmscenario' not in files:
            raise Exception('The archive "%s" does not contain "scenario.gotmscenario"; it cannot be a GOTM result.' % path)
        if 'result.nc' not in files:
            raise Exception('The archive "%s" does not contain "result.nc"; it cannot be a GOTM result.' % path)

        # Create a temporary directory to which we can unpack the archive.
        tempdir = self.getTempDir()

        df = container.getItem('scenario.gotmscenario')
        self.scenario = scenario.Scenario.fromSchemaName(scenario.guiscenarioversion)
        self.scenario.loadAll(df)
        df.release()

        df = container.getItem('result.nc')
        resultpath = os.path.join(tempdir,'result.nc')
        df.saveToFile(resultpath)
        df.release()

        df = container.getItem('result.xml')
        if df!=None:
            f = df.getAsReadOnlyFile()
            valuedom = xml.dom.minidom.parse(f)
            self.store.setStore(valuedom)
            f.close()
            df.release()

        # Close the archive
        container.release()

        # Attach the result, try to open the CDF file
        self.datafile = resultpath
        self.getcdf()

        # Reset "changed" status.
        self.changed = False
        self.path = path

    def setFigure(self,source):
        name = source.root.getLocation(['Name']).getValue()
        assert name!=None, 'Name attribute of figure is not set; cannot store figure.'

        setroot = self.store.root.getLocation(['FigureSettings'])
        for fig in setroot.getLocationMultiple(['Figure']):
            if fig.getLocation(['Name']).getValue()==name: break
        else:
            fignodename = source.root.templatenode.getAttribute('id')
            fig = setroot.addChild(fignodename)
            assert fig!=None, 'Unable to add new node "%s".' % fignodename
        fig.copyFrom(source.root,replace=True)
        fig.getLocation(['Name']).setValue(name)

    def getFigure(self,name,target):
        setroot = self.store.root.getLocation(['FigureSettings'])
        for fig in setroot.getLocationMultiple(['Figure']):
            if fig.getLocation(['Name']).getValue()==name:
                target.root.copyFrom(fig,replace=True)
                return True
        return False

    def unlink(self):
        if self.nc!=None:
            # Close NetCDF result file.
            self.nc.close()
            self.nc = None
        if self.tempdir!=None:
            # Delete temporary directory.
            print 'Deleting temporary result directory "%s".' % self.tempdir
            shutil.rmtree(self.tempdir)
            self.tempdir = None

    def attach(self,srcpath,scenario=None,copy=True):
        self.scenario = scenario
        
        if copy:
            # Create a copy of the result file.
            tempdir = self.getTempDir(empty=True)
            datafile = os.path.join(tempdir,'result.nc')
            shutil.copyfile(srcpath,datafile)
        else:
            datafile = srcpath

        # Store link to result file, and try to open the CDF file
        self.datafile = datafile
        self.getcdf()

        # Attached to an existing result: we consider it unchanged.
        self.changed = False
        
        if self.scenario!=None and self.scenario.path!=None and self.scenario.path.endswith('.gotmscenario'):
            self.path = self.scenario.path[:-12]+'gotmresult'
        else:
            self.path = None

    def getcdf(self):
        if self.nc!=None: return self.nc
        assert self.datafile!=None, 'The result object has not yet been attached to an actual result.'
        try:
          #pycdf self.nc = pycdf.CDF(str(self.datafile))
          self.nc = NetCDFFile(self.datafile)
        #pycdf except pycdf.CDFError, e:
        except Exception, e:
            raise Exception('An error occured while opening the NetCDF file "%s": %s' % (self.datafile,str(e)))
        return self.nc

    def getVariableNames(self,plotableonly=True):
        nc = self.getcdf()

        # Get names of NetCDF variables
        try:
          #pycdf vars = nc.variables()
          vars = nc.variables
          if plotableonly:
              # Only take variables with 3 or 4 dimensions
              varNames = []
              for v in vars.keys():
                  #pycdf dimnames = vars[v][0]
                  dimnames = vars[v].dimensions
                  dimcount = len(dimnames)
                  if   dimcount==3:
                      if dimnames==('time','lat','lon'):
                          varNames += [v]
                  elif dimcount==4:
                      if (dimnames==('time','z','lat','lon')) | (dimnames==('time','z1','lat','lon')):
                          varNames += [v]
          else:
              # Take all variables
              varNames = vars.keys()

        #pycdf except pycdf.CDFError, msg:
        except Exception, e:
            raise Exception('CDFError: '+str(e))

        return varNames

    def getVariableLongNames(self):
      varnames = self.getVariableNames()
      nc = self.getcdf()
      vardict = {}
      for varname in varnames:
          #pycdf varname_str = str(varname)
          #pycdf vardict[varname] = nc.var(varname_str).long_name
          vardict[varname] = nc.variables[varname].long_name
      return vardict

    def getVariable(self,varname,check=True):
        varname = str(varname)
        if check:
            nc = self.getcdf()
            #pycdf vars = nc.variables()
            vars = nc.variables
            if not (varname in vars): return None
        return self.ResultVariable(self,varname)

    def getTime(self):
        if self.t==None:
            nc = self.getcdf()

            # Get time coordinate (in seconds since reference date)
            secs = nc.variables['time'][:]
            
            # Convert time-in-seconds to Python datetime objects.
            dateref = self.getReferenceDate()
            t = matplotlib.numerix.zeros((secs.shape[0],),matplotlib.numerix.PyObject)
            for it in range(t.shape[0]):
                t[it] = dateref + datetime.timedelta(secs[it]/3600/24)

            # Get time step (of output, not simulation!)
            if secs.shape[0]==1:
                # Only one time step saved: try to get output time step from scenario. If this
                # does not work, assume the simulation started with MinN == 1, allowing us to
                # use the first time (in seconds) as time step. If that does not work, default
                # to one day.
                dt = None
                if self.scenario!=None:
                    dt = self.scenario.getProperty(['output','dtsave'],usedefault=True)
                if dt==None:
                    if secs[0]>0:
                        dt=secs[0]
                    else:
                        dt=3600*24
            else:
                dt = secs[1]-secs[0]

            # Create staggered time grid.
            t_stag = matplotlib.numerix.zeros((secs.shape[0]+1,),matplotlib.numerix.PyObject)
            halfdt = datetime.timedelta(seconds=dt/2.)
            t_stag[0]  = t[0]-halfdt
            t_stag[1:] = t[:]+halfdt
            
            # Cache time grids.
            self.t = t
            self.t_stag = t_stag
            
        return (self.t,self.t_stag)

    def getDepth(self):
        if self.z==None:
            nc = self.getcdf()

            # Get layer heights
            h = nc.variables['h'][:,:,0,0]
            
            # Get depths of interfaces
            z1 = h.cumsum(1)
            z1 = matplotlib.numerix.concatenate((matplotlib.numerix.zeros((z1.shape[0],1),matplotlib.numerix.typecode(z1)),z1),1)
            bottomdepth = z1[0,-1]-nc.variables['zeta'][0,0,0]
            z1 -= bottomdepth

            # Get depth of layer centers
            z = z1[:,1:z1.shape[1]]-0.5*h

            # Interpolate in time to create staggered grid
            z1_med = matplotlib.numerix.concatenate((matplotlib.numerix.take(z1,(0,),0),z1,matplotlib.numerix.take(z1,(-1,),0)),0)
            z_stag = 0.5 * (z1_med[0:z1_med.shape[0]-1,:] + z1_med[1:z1_med.shape[0],:])
            
            z_med = matplotlib.numerix.concatenate((z,matplotlib.numerix.take(z1,(-1,),1)),1)
            z_med = matplotlib.numerix.concatenate((matplotlib.numerix.take(z_med,(0,),0),z_med,matplotlib.numerix.take(z_med,(-1,),0)),0)
            z1_stag = 0.5 * (z_med[0:z_med.shape[0]-1,:] + z_med[1:z_med.shape[0],:])

            self.z = z
            self.z1 = z1
            self.z_stag = z_stag
            self.z1_stag = z1_stag

        return (self.z,self.z1,self.z_stag,self.z1_stag)

    def getReferenceDate(self):
      # Retrieve reference date/time.
      nc = self.getcdf()
      #pycdf timeunit = nc.var('time').units
      timeunit = nc.variables['time'].units
      datematch = re.compile('(\d\d\d\d)[-\/](\d\d)-(\d\d) (\d\d):(\d\d):(\d\d)').search(timeunit, 1)
      assert datematch!=None, 'Unable to parse "units" attribute of "time" variable in NetCDF file!'
      refvals = map(int,datematch.group(1,2,3,4,5,6)) # Convert matched strings into integers
      dateref = datetime.datetime(*refvals)
      return dateref

    def getVariableTree(self,path):
        otherstores = {}
        if self.scenario!=None: otherstores['scenario'] = self.scenario
        return PlotVariableStore.getVariableTree(self,path,otherstores=otherstores)
    
