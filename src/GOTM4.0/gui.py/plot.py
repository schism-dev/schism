import datetime

import matplotlib
matplotlib.use('Qt4Agg')
matplotlib.rcParams['numerix'] = 'numpy'
import matplotlib.numerix,matplotlib.numerix.ma,matplotlib.colors
import matplotlib.dates
import matplotlib.pylab
import matplotlib.backends.backend_agg

import common,xmlstore,data

class MonthFormatter(matplotlib.dates.DateFormatter):
    def __init__(self):
        matplotlib.dates.DateFormatter.__init__(self,'%b')

    def __call__(self, x, pos=None):
        return matplotlib.dates.DateFormatter.__call__(self,x,pos)[0]
        
class VariableTransform(data.PlotVariable):
    def __init__(self,sourcevar):
        data.PlotVariable.__init__(self)
        self.sourcevar = sourcevar

    def getName(self):
        return self.sourcevar.getName()

    def getLongName(self):
        return self.sourcevar.getLongName()

    def getUnit(self):
        return self.sourcevar.getLongName()

    def getDimensions(self):
        return self.sourcevar.getDimensions()
        
class VariableSlice(VariableTransform):
    def __init__(self,variable,slicedimension,slicecoordinate):
        VariableTransform.__init__(self,variable)
        self.slicedim = slicedimension
        self.sliceval = slicecoordinate

        dims = self.sourcevar.getDimensions()
        for (i,d) in enumerate(dims):
            if d==self.slicedim: break
        else:
            assert False, 'Slice dimension "%s" is not present for this variable.' % self.slicedim
        self.islicedim = i

    def getDimensions(self):
        dims = self.sourcevar.getDimensions()
        return [d for d in dims if d!=self.slicedim]

    def getValues(self,bounds,staggered=False,coordinatesonly=False):
        bounds.insert(self.islicedim,(self.sliceval,self.sliceval))
        dims = self.getDimensions()
        data = self.sourcevar.getValues(bounds,staggered=staggered,coordinatesonly=coordinatesonly)
        if data==None or 0 in data[-1].shape: return None
        assert data[self.islicedim].ndim==1, 'Slicing is not (yet) supported for dimensions that have coordinates that depend on other dimensions.'
        ipos = data[self.islicedim].searchsorted(self.sliceval)
        if ipos==0 or ipos>=data[self.islicedim].shape[0]: return None
        leftx  = data[self.islicedim][ipos-1]
        rightx = data[self.islicedim][ipos]
        deltax = rightx-leftx
        stepx = self.sliceval-leftx
        if isinstance(deltax,datetime.timedelta):
            relstep = common.timedelta2float(stepx)/common.timedelta2float(deltax)
        else:
            relstep = stepx/deltax
        if len(dims)==1:
            data.pop(self.islicedim)
            for idat in range(len(data)):
                if data[idat].ndim==2:
                    if ipos>0 and ipos<len(data[self.islicedim]):
                        # centered: left and right bound available
                        left  = data[idat].take((ipos-1,),self.islicedim).squeeze()
                        right = data[idat].take((ipos,  ),self.islicedim).squeeze()
                        data[idat] = left + relstep*(right-left)
                    elif ipos==0:
                        # left-aligned (only right bound available)
                        data[idat]=data[idat].take((ipos,),self.islicedim).squeeze()
                    else:
                        # right-aligned (only left bound available)
                        data[idat]=data[idat].take((ipos-1,),self.islicedim).squeeze()
        else:
            assert False,'Cannot take slice because the result does not have 1 coordinate dimension (instead it has %i: %s).' % (len(dims),dims)
        return data

class Figure:

    def __init__(self,figure):
        self.figure = figure
        self.canvas = figure.canvas

        # Create store for the explicitly set properties
        self.properties = xmlstore.TypedStore('schemas/figure/gotmgui.xml')
        self.propertiesinterface = self.properties.getInterface()
        self.propertiesinterface.notifyOnDefaultChange = False
        self.propertiesinterface.addChangeHandler(self.onPropertyChanged)
        self.propertiesinterface.addStoreChangedHandler(self.onPropertyStoreChanged)
        
        # Create store for property defaults
        self.defaultproperties = xmlstore.TypedStore('schemas/figure/gotmgui.xml')

        # Set some default properties.
        self.defaultproperties.setProperty(['FontScaling'],100)
        self.defaultproperties.setProperty(['TimeAxis',  'Label'],'time')
        self.defaultproperties.setProperty(['DepthAxis', 'Label'],'depth (m)')

        self.properties.setDefaultStore(self.defaultproperties)

        self.sources = {}
        self.defaultsource = None
        self.updating = True
        self.dirty = False
        self.haschanged = False
        
        self.callbacks = {'completeStateChange':[]}
        
    def registerCallback(self,eventname,callback):
        assert eventname in self.callbacks, 'Event "%s" is unknown.' % eventname
        self.callbacks[eventname].append(callback)

    def setUpdating(self,allowupdates):
        if self.updating != allowupdates:
            self.updating = allowupdates
            if allowupdates and self.dirty: self.update()

    def onPropertyChanged(self,node):
        self.onPropertyStoreChanged()

    def onPropertyStoreChanged(self):
        self.haschanged = True
        self.update()

    def clearSources(self):
        self.sources = {}
        self.defaultsource = None

    def addDataSource(self,name,obj):
        self.sources[name] = obj
        if self.defaultsource==None: self.defaultsource = name

    def clearProperties(self):
        self.properties.root.clearValue(recursive=True)

    def setProperties(self,props):
        self.properties.setStore(props)
        self.update()

    def getPropertiesCopy(self):
        return self.properties.toXmlDom()

    def clearVariables(self):
        self.properties.root.getLocation(['Data']).removeChildren('Series')

    def addVariable(self,varname,source=None):
        datanode = self.properties.root.getLocation(['Data'])
        series = datanode.addChild('Series')
        series.getLocation(['Variable']).setValue(varname)
        if source!=None:
            series.getLocation(['Source']).setValue(source)
        self.update()

    def hasChanged(self):
        return self.haschanged

    def resetChanged(self):
        self.haschanged = False

    def update(self):
        if not self.updating:
            self.dirty = True
            return

        self.figure.clear()

        axes = self.figure.add_subplot(111)
        
        textscaling = self.properties.getProperty(['FontScaling'],usedefault=True)/100.
        
        # First scale the default font size; this takes care of all relative font sizes (e.g. "small")
        matplotlib.font_manager.fontManager.set_default_size(textscaling*matplotlib.rcParams['font.size'])
        
        # Now get some relevant font sizes.
        # Scale font sizes with text scaling parameter if they are absolute sizes.
        # (if they are strings, they are relative sizes already)
        titlesize = matplotlib.rcParams['axes.titlesize']
        axeslabelsize = matplotlib.rcParams['axes.labelsize']
        if not isinstance(titlesize,     basestring): titlesize     *=textscaling
        if not isinstance(axeslabelsize, basestring): axeslabelsize *=textscaling

        # Get the default line width
        deflinewidth = matplotlib.rcParams['lines.linewidth']
        deflinecolor = matplotlib.rcParams['lines.color']
        deflinecolor = matplotlib.colors.colorConverter.to_rgb(deflinecolor)
        deflinecolor = xmlstore.StoreColor.fromNormalized(*deflinecolor)

        # Get forced axes boundaries (will be None if not set; then we autoscale)
        tmin = self.properties.getProperty(['TimeAxis','Minimum'])
        tmax = self.properties.getProperty(['TimeAxis','Maximum'])
        zmin = self.properties.getProperty(['DepthAxis','Minimum'])
        zmax = self.properties.getProperty(['DepthAxis','Maximum'])

        # Link between dimension name (e.g., "time", "z") and axis (e.g., "x", "y")
        dim2data = {'time':{'forcedrange':[tmin,tmax],'datarange':[None,None]},
                    'z'   :{'forcedrange':[zmin,zmax],'datarange':[None,None]}}

        # Shortcuts to the nodes specifying the variables to plot.
        forceddatanode = self.properties.root.getLocation(['Data'])
        forcedseries = forceddatanode.getLocationMultiple(['Series'])

        # Shortcut to the node that will hold defaults for the plotted variables.
        defaultdatanode = self.defaultproperties.root.getLocation(['Data'])

        # This variable will hold all long names of the plotted variables.
        # It will be used to create the plot title.
        longnames = []

        for (iseries,seriesnode) in enumerate(forcedseries):
            # Get the name and data source of the variable to plot.
            varname   = seriesnode.getLocation(['Variable']).getValue()
            varsource = seriesnode.getLocation(['Source']).getValue()
            if varsource==None:
                # No data source specified; take default.
                assert self.defaultsource!=None, 'No data source set for variable "%s", but no default source available either.' % varname
                varsource = self.defaultsource
                
            # Get variable object.
            var = self.sources[varsource].getVariable(varname)
            assert var!=None, 'Source "%s" does not contain variable with name "%s".' % (varsource,varname)
            
            # Create default series information
            defaultseriesnode = defaultdatanode.getNumberedChild('Series',iseries,create=True)
            defaultseriesnode.getLocation(['Variable']).setValue(varname)
            defaultseriesnode.getLocation(['Source']).setValue(varsource)
            defaultseriesnode.getLocation(['PlotType2D']).setValue(0)
            defaultseriesnode.getLocation(['PlotType3D']).setValue(0)
            defaultseriesnode.getLocation(['LineWidth']).setValue(deflinewidth)
            defaultseriesnode.getLocation(['LineColor']).setValue(deflinecolor)
            defaultseriesnode.getLocation(['LogScale']).setValue(False)

            # Store the variable long name (to be used for building title)
            longnames.append(var.getLongName())

            # Build list of dimension boundaries for current variable.
            # For dimensions that have equal lower and upper bound, take a slice.
            dimbounds = []
            for dimname in var.getDimensions():
                if dimname not in dim2data: dim2data[dimname] = {}
                dimdata = dim2data[dimname]
                dimdata['used'] = True
                if 'forcedrange' in dimdata:
                    if dimdata['forcedrange'][0]==dimdata['forcedrange'][1] and (dimdata['forcedrange'][0]!=None):
                        var = VariableSlice(var,dimname,dimdata['forcedrange'][0])
                    else:
                        dimbounds.append(dimdata['forcedrange'])
                else:
                    dimbounds.append((None,None))
                    
            # Get final list of dimensions (after taking slices).
            dims = var.getDimensions()

            # Get the data (centered coordinates and values)
            data = var.getValues(dimbounds,staggered=False)
            
            # Skip this variable if no data are available.
            if data==None or 0 in data[-1].shape: continue
            
            # Find non-singleton dimensions (singleton dimension: dimension with length one)
            # Store singleton dimensions as fixed extra coordinates.
            gooddims = []
            newdims = []
            fixedcoords = []
            for idim,dimname in enumerate(dims):
                if data[-1].shape[idim]>1:
                    # Normal dimension (more than one coordinate)
                    gooddims.append(idim)
                    newdims.append(dimname)
                elif data[-1].shape[idim]==1:
                    # Singleton dimension
                    fixedcoords.append((dimname,data[idim][0]))
                    
            dims = newdims
            values = data[-1].squeeze()

            # Get effective number of independent dimensions (singleton dimensions removed)
            dimcount = len(dims)
            seriesnode.getLocation(['DimensionCount']).setValue(dimcount)

            # Get the plot type, based on the number of dimensions
            if dimcount==0:
                plottypenodename = 'PlotType2D'
            elif dimcount==1:
                plottypenodename = 'PlotType2D'
            elif dimcount==2:
                plottypenodename = 'PlotType3D'
            else:
                raise Exception('This variable has %i independent dimensions. Can only plot variables with 0, 1 or 2 independent dimensions.' % dimcount)
            plottype = seriesnode.getLocation([plottypenodename]).getValueOrDefault()

            # We use a staggered grid (coordinates at interfaces, values at centers) for certain 3D plot types.
            staggered = (plottypenodename=='PlotType3D' and plottype==0)

            # Get coordinate data (now that we know whether to use a staggered grid)
            data = var.getValues(dimbounds,staggered=staggered,coordinatesonly=True)
            data = [data[idim].squeeze() for idim in gooddims] + [values]

            # Get the minimum and maximum values; store these as default.
            defaultseriesnode.getLocation(['Minimum']).setValue(data[-1].min())
            defaultseriesnode.getLocation(['Maximum']).setValue(data[-1].max())

            # Mask values that are not within (minimum, maximum) range.
            minimum = seriesnode.getLocation(['Minimum']).getValue()
            maximum = seriesnode.getLocation(['Maximum']).getValue()
            if minimum!=None and maximum!=None:
                data[-1] = matplotlib.numerix.ma.masked_array(data[-1],matplotlib.numerix.logical_or(data[-1]<minimum, data[-1]>maximum))
            elif minimum!=None:
                data[-1] = matplotlib.numerix.ma.masked_array(data[-1],data[-1]<minimum)
            elif maximum!=None:
                data[-1] = matplotlib.numerix.ma.masked_array(data[-1],data[-1]>maximum)

            # Transform to log-scale if needed (first mask values <= zero)
            logscale = seriesnode.getLocation(['LogScale']).getValueOrDefault()
            if logscale:
                data[-1] = matplotlib.numerix.ma.masked_array(data[-1],data[-1]<=0.)
                data[-1] = matplotlib.numerix.ma.log10(data[-1])

            # Get label
            defaultlabel = '%s (%s)' % (var.getLongName(),var.getUnit())
            if logscale: defaultlabel = 'log10 '+defaultlabel
            defaultseriesnode.getLocation(['Label']).setValue(defaultlabel)
            label = seriesnode.getLocation(['Label']).getValueOrDefault()

            # Enumerate over the dimension of the variable.
            for idim,dimname in enumerate(dims):
                # Get minimum and maximum coordinates.
                if data[idim].ndim==1:
                    # Coordinates provided as vector (1D) valid over whole domain.
                    datamin = data[idim][0]
                    datamax = data[idim][-1]
                else:
                    # Coordinates are provided as multidimensional array, with a value for every
                    # coordinate (data point) in the domain. We assume that for a given point
                    # in the space of the other coordinates, the current cordinate increases
                    # monotonously (i.e., position 0 holds the lowest value and position -1 the
                    # highest)
                    datamin = data[idim].take((0, ),idim).min()
                    datamax = data[idim].take((-1,),idim).max()

                # Update effective dimension bounds                    
                if 'datarange' not in dim2data[dimname]: dim2data[dimname]['datarange'] = [None,None]
                effrange = dim2data[dimname]['datarange']
                if effrange[0]==None or datamin<effrange[0]: effrange[0] = datamin
                if effrange[1]==None or datamax>effrange[1]: effrange[1] = datamax

                # Convert time (datetime objects) to time unit used by MatPlotLib
                if dimname=='time': data[idim] = matplotlib.dates.date2num(data[idim])
            
            # Plot the data series
            if len(dims)==0:
                # Zero-dimensional coordinate space (i.e., only a single data value is available)
                # No plotting of coordinate-less data (yet)
                pass
            if len(dims)==1:
                # One-dimensional coordinate space (x). Use x-axis for coordinates, unless the
                # coordinate name is "z" (i.e., depth).
                xdim = 0
                datadim = 1
                if dims[xdim]=='z':
                    xdim = 1
                    datadim = 0
                linewidth = seriesnode.getLocation(['LineWidth']).getValueOrDefault()
                linecolor = seriesnode.getLocation(['LineColor']).getValueOrDefault()
                lines = axes.plot(data[xdim],data[datadim],'-',linewidth=linewidth,color=linecolor.getNormalized())
                if xdim==0:
                    dim2data[dims[0]]['axis'] = 'x'
                    axes.set_ylabel(label,size=axeslabelsize)
                else:
                    dim2data[dims[0]]['axis'] = 'y'
                    axes.set_xlabel(label,size=axeslabelsize)
            elif len(dims)==2:
                # Two-dimensional coordinate space (x,y). Use x-axis for first coordinate dimension,
                # and y-axis for second coordinate dimension.
                xdim = 0
                ydim = 1

                dim2data[dims[xdim]]['axis'] = 'x'
                dim2data[dims[ydim]]['axis'] = 'y'

                X = data[xdim]
                Y = data[ydim]
                Z = data[-1]
                
                # Get length of coordinate dimensions. Coordinates can be provided as vectors
                # valid over the whole domain, or as n-D array that match the shape of the values.
                if X.ndim==1:
                    xlength = X.shape[0]
                else:
                    xlength = X.shape[xdim]
                if Y.ndim==1:
                    ylength = Y.shape[0]
                else:
                    ylength = Y.shape[ydim]
                    
                # Adjust X dimension.
                if X.ndim==1:
                    X = X.reshape((1,-1)).repeat(ylength, 0)
                elif xdim<ydim:
                    X = X.transpose()
                    
                # Adjust Y dimension.
                if Y.ndim==1:
                    Y = Y.reshape((-1,1)).repeat(xlength, 1)
                elif xdim<ydim:
                    Y = Y.transpose()
                    
                # Adjust Z dimension.
                if xdim<ydim:
                    Z = Z.transpose()
                    
                if plottype==1:
                  cc = seriesnode.getLocation(['ContourCount']).getValue()
                  if cc!=None:
                    pc = axes.contourf(X,Y,Z,cc)
                  else:
                    pc = axes.contourf(X,Y,Z)
                  if cc==None:
                      defaultseriesnode.getLocation(['ContourCount']).setValue(len(pc.levels)-2)
                else:
                  #pc = axes.pcolor(X,Y,Z,shading='flat', cmap=matplotlib.pylab.cm.jet)
                  pc = axes.pcolormesh(X,Y,Z,shading='flat', cmap=matplotlib.pylab.cm.jet)
                  
                # Create colorbar
                if (Z.ravel()==Z[0,0]).all():
                    # Explicitly set color range; MatPlotLib 0.90.0 chokes on identical min and max.
                    pc.set_clim((Z[0,0]-1,Z[0,0]+1))
                    cb = self.figure.colorbar(pc)
                else:
                    cb = self.figure.colorbar(pc)
                    
                # Text for colorbar
                if label!='': cb.set_label(label,size=axeslabelsize)
                for l in cb.ax.xaxis.get_ticklabels():
                    l.set_size(l.get_size()*textscaling)
                for l in cb.ax.yaxis.get_ticklabels():
                    l.set_size(l.get_size()*textscaling)

            # Hold all plot properties so we can plot additional data series.
            axes.hold(True)

        # Remove unused series (remaining from previous plots that had more data series)
        defaultdatanode.removeChildren('Series',first=len(forcedseries))

        # Create and store title
        self.defaultproperties.setProperty(['Title'],', '.join(longnames))
        title = self.properties.getProperty(['Title'],usedefault=True)
        assert title!=None, 'Title must be available, either explicitly set or as default.'
        if title!='': axes.set_title(title,size=titlesize)

        # Store natural axes bounds (based on data ranges).
        self.defaultproperties.setProperty(['TimeAxis',  'Minimum'],dim2data['time']['datarange'][0])
        self.defaultproperties.setProperty(['TimeAxis',  'Maximum'],dim2data['time']['datarange'][1])
        self.defaultproperties.setProperty(['DepthAxis', 'Minimum'],dim2data['z']   ['datarange'][0])
        self.defaultproperties.setProperty(['DepthAxis', 'Maximum'],dim2data['z']   ['datarange'][1])
        
        #if tmin_eff==tmax_eff and tmin_eff!=None:
        #    tmin_eff -= datetime.timedelta(1)
        #    tmax_eff += datetime.timedelta(1)
        
        # Get effective ranges for each dimension (based on forced limits and natural data ranges)
        for dim,dat in dim2data.iteritems():
            effrange = dat['forcedrange'][:]
            if effrange[0]==None: effrange[0] = dat['datarange'][0]
            if effrange[1]==None: effrange[1] = dat['datarange'][1]
            dat['range'] = effrange

        # Configure time axis (if any).
        if 'time' in dim2data and 'axis' in dim2data['time']:
            timedata = dim2data['time']
            timeaxis = timedata['axis']
            
            # Obtain label for time axis.
            tlabel = self.properties.getProperty(['TimeAxis','Label'],usedefault=True)
            assert tlabel!=None, 'Time axis label must be available, either explicitly set or as default.'

            # Configure limits and label of time axis.
            if timeaxis=='x':
                taxis = axes.xaxis
                if tlabel!='': axes.set_xlabel(tlabel,size=axeslabelsize)
                axes.set_xlim(matplotlib.dates.date2num(timedata['range'][0]),matplotlib.dates.date2num(timedata['range'][1]))
            elif timeaxis=='y':
                taxis = axes.yaxis
                if tlabel!='': axes.set_ylabel(tlabel,size=axeslabelsize)
                axes.set_ylim(matplotlib.dates.date2num(timedata['range'][0]),matplotlib.dates.date2num(timedata['range'][1]))

            # Select tick type and spacing based on the time span to show.
            dayspan = (timedata['range'][1]-timedata['range'][0]).days
            if dayspan/365>10:
              # more than 10 years
              taxis.set_major_locator(matplotlib.dates.YearLocator(base=5))
              taxis.set_major_formatter(matplotlib.dates.DateFormatter('%Y'))
            elif dayspan/365>1:
              # less than 10 but more than 1 year
              taxis.set_major_locator(matplotlib.dates.YearLocator(base=1))
              taxis.set_major_formatter(matplotlib.dates.DateFormatter('%Y'))
            elif dayspan>61:
              # less than 1 year but more than 2 months
              taxis.set_major_locator(matplotlib.dates.MonthLocator(interval=1))
              taxis.set_major_formatter(MonthFormatter())
            elif dayspan>7:
              # less than 2 months but more than 1 day
              taxis.set_major_locator(matplotlib.dates.DayLocator(interval=15))
              taxis.set_major_formatter(matplotlib.dates.DateFormatter('%d %b'))
            elif dayspan>1:
              # less than 1 week but more than 1 day
              taxis.set_major_locator(matplotlib.dates.DayLocator(interval=1))
              taxis.set_major_formatter(matplotlib.dates.DateFormatter('%d %b'))
            else:
              # less than 1 day
              taxis.set_major_locator(matplotlib.dates.HourLocator(interval=1))
              taxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M'))
                      
        self.properties.setProperty(['HasTimeAxis'],'used' in dim2data['time'])

        # Configure depth axis (y-axis), if any.
        if 'z' in dim2data and 'axis' in dim2data['z']:
            zdata = dim2data['z']
            zaxis = zdata['axis']

            # Obtain label for depth axis.
            zlabel = self.properties.getProperty(['DepthAxis','Label'],usedefault=True)
            assert zlabel!=None, 'Depth axis label must be available, either explicitly set or as default.'

            # Configure limits and label of depth axis.
            if zaxis=='x':
                axes.set_xlim(zdata['range'][0],zdata['range'][1])
                if zlabel!='': axes.set_xlabel(zlabel,size=axeslabelsize)
            elif zaxis=='y':
                axes.set_ylim(zdata['range'][0],zdata['range'][1])
                if zlabel!='': axes.set_ylabel(zlabel,size=axeslabelsize)
        self.properties.setProperty(['HasDepthAxis'],'used' in dim2data['z'])
        
        # Scale the text labels
        for l in axes.get_xaxis().get_ticklabels():
            l.set_size(l.get_size()*textscaling)
        for l in axes.get_yaxis().get_ticklabels():
            l.set_size(l.get_size()*textscaling)

        # Draw the plot to screen.
        self.canvas.draw()
        
        for cb in self.callbacks['completeStateChange']: cb(len(forcedseries)>0)

        self.dirty = False
