import os, xml.dom.minidom, shutil

import xmlstore, plot, matplotlib

class Report:
    def __init__(self):
        self.store = xmlstore.TypedStore('schemas/report/gotmgui.xml')

        self.defaultstore = xmlstore.TypedStore('schemas/report/gotmgui.xml')

        # Set some default properties.
        self.defaultstore.setProperty(['Figures','Width'],10)
        self.defaultstore.setProperty(['Figures','Height'],8)
        self.defaultstore.setProperty(['Figures','Resolution'],96)
        self.defaultstore.setProperty(['Figures','FontScaling'],100)

        self.store.setDefaultStore(self.defaultstore)
        
    def generate(self,result,outputpath,templatepath,columncount=2,callback=None):
        xmldocument = xml.dom.minidom.parse(os.path.join(templatepath,'index.xml'))
        scenario = result.scenario
        
        # Get report settings
        figuresize = (self.store.getProperty(['Figures','Width'],usedefault=True),self.store.getProperty(['Figures','Height'],usedefault=True))
        dpi = self.store.getProperty(['Figures','Resolution'],usedefault=True)
        fontscaling = self.store.getProperty(['Figures','FontScaling'],usedefault=True)

        # Get list of variables to plot
        selroot = self.store.root.getLocation(['Figures','Selection'])
        plotvariables = [node.getValue() for node in selroot.children]

        steps = float(2+len(plotvariables))
        istep = 0

        # Create output directory if it does not exist yet.
        if not os.path.isdir(outputpath): os.mkdir(outputpath)

        for f in os.listdir(templatepath):
            fullpath = os.path.join(templatepath,f)
            if f.lower()!='index.xml' and os.path.isfile(fullpath): shutil.copy(fullpath,os.path.join(outputpath,f))

        for node in xmldocument.getElementsByTagName('gotm:scenarioproperty'):
            variablepath = node.getAttribute('variable')
            assert variablepath!='', 'gotm:scenarioproperty node in report template lacks "variable" attribute pointing to a location in the scenario.'
            variablenode = scenario.root.getLocation(variablepath.split('/'))
            assert variablenode!=None, 'Unable to locate "%s" in the scenario.' % variablepath
            val = variablenode.getValueAsString()
            node.parentNode.replaceChild(xmldocument.createTextNode(unicode(val)),node)
            node.unlink()

        scenarionodes = xmldocument.getElementsByTagName('gotm:scenario')
        assert len(scenarionodes)<=1, 'Found more than one "gotm:scenario" node in the report template.'
        if len(scenarionodes)>0:
            if callback!=None: callback(istep/steps,'Creating scenario description...')
            scenarionode = scenarionodes[0]

            sceninterface = xmlstore.TypedStoreInterface(scenario,showhidden=False,omitgroupers=True)

            scentable = xmldocument.createElement('table')
            scentable.setAttribute('id','tableScenario')

            totaldepth = sceninterface.getDepth(scenario.root)

            # Create columns.
            for i in range(totaldepth-2):
                col = xmldocument.createElement('col')
                col.setAttribute('width','25')
                scentable.appendChild(col)
            col = xmldocument.createElement('col')
            scentable.appendChild(col)
            col = xmldocument.createElement('col')
            scentable.appendChild(col)

            # Create rows
            for tr in sceninterface.toHtml(scenario.root,xmldocument,totaldepth-1,level=-1,hidedefaults=True):
                scentable.appendChild(tr)

            scenarionode.parentNode.replaceChild(scentable,scenarionode)

        istep += 1

        figuresnodes = xmldocument.getElementsByTagName('gotm:figures')
        assert len(figuresnodes)<=1, 'Found more than one "gotm:figures" node in the report template.'
        if len(figuresnodes)>0:
            figuresnode = figuresnodes[0]
        else:
            figuresnode = None
        if len(plotvariables)>0 and figuresnode!=None:
            mplfigure = matplotlib.figure.Figure(figsize=(figuresize[0]/2.54,figuresize[1]/2.54))
            canvas = matplotlib.backends.backend_agg.FigureCanvasAgg(mplfigure)
            fig = plot.Figure(mplfigure)
            fig.addDataSource('result',result)
            figurestable = xmldocument.createElement('table')
            icurvar = 0
            tr = None
            for varpath in plotvariables:
                varid = varpath.split('/')[-1]
                
                longname = result.getVariable(varid).getLongName()
                if callback!=None: callback(istep/steps,'Creating figure for %s...' % longname)
                
                if icurvar % columncount == 0:
                    tr = xmldocument.createElement('tr')
                    figurestable.appendChild(tr)
                fig.setUpdating(False)
                if not result.getFigure('result/'+varpath,fig.properties):
                    fig.clearProperties()
                    fig.addVariable(varid)
                fig.properties.setProperty(['FontScaling'],fontscaling)
                fig.setUpdating(True)
                filename = varid+'.png'
                outputfile = os.path.join(outputpath,filename)
                canvas.print_figure(outputfile,dpi=dpi)

                img = xmldocument.createElement('img')
                img.setAttribute('src',filename)
                img.setAttribute('alt',longname)
                img.setAttribute('style','width:%.2fcm' % figuresize[0])
                td = xmldocument.createElement('td')
                td.appendChild(img)
                tr.appendChild(td)

                icurvar = icurvar+1
                istep += 1
            for i in range(columncount - len(tr.childNodes)):
                tr.appendChild(xmldocument.createElement('td'))
            figuresnode.parentNode.replaceChild(figurestable,figuresnode)
        elif figuresnode!=None:
            figuresnode.parentNode.removeChild(figuresnode)
        
        if callback!=None: callback(istep/steps,'Writing HTML...')

        if outputpath!='':
            import codecs
            f = codecs.open(os.path.join(outputpath,'index.html'),'w','utf-8')
            xmldocument.writexml(f,encoding='utf-8')
            f.close()
        else:
            print xmldocument.toxml('utf-8')
        istep += 1

        if callback!=None: callback(istep/steps,'Done.')

