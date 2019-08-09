#!/usr/bin/python

#$Id: visualizer.py,v 1.17 2007-04-13 12:40:08 jorn Exp $

from PyQt4 import QtGui,QtCore

import commonqt, data, report

import matplotlib
matplotlib.use('Qt4Agg')
matplotlib.rcParams['numerix'] = 'numeric'

import sys,datetime
import xml.sax
import os.path

def loadResult(path):
    result = data.Result()

    try:
        if path.endswith('.gotmresult'):
            result.load(path)
        elif path.endswith('.nc'):
            result.attach(path)
        else:
            # We do not recognize this file type; try both GOTM result and NetCDF
            done = True
            try:
                result.attach(path)
            except Exception,e:
                done = False
            if not done:
                done = True
                try:
                    result.load(path)
                except Exception,e:
                    done = False
            if (not done):
                raise Exception('The file "'+path+'" is not a GOTM result or a NetCDF file.')
    except:
        result.unlink()
        result = None
        raise

    return result

class OpenWidget(QtGui.QWidget):
    def __init__(self,parent=None,mrupaths=[]):
        QtGui.QWidget.__init__(self,parent)

        self.pathOpen = commonqt.PathEditor(self,header='File to open: ',mrupaths=mrupaths)
        self.pathOpen.filter = 'GOTM result files (*.gotmresult);;NetCDF files (*.nc);;All files (*.*)'

        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.pathOpen)
        layout.setMargin(0)
        self.setLayout(layout)
        self.connect(self.pathOpen, QtCore.SIGNAL("onChanged()"), self.completeStateChanged)
        
    def completeStateChanged(self):
        self.emit(QtCore.SIGNAL('onCompleteStateChanged()'))

    def isComplete(self):
        return self.pathOpen.hasPath()

    def getResult(self):
        return loadResult(self.pathOpen.path())

class PageOpen(commonqt.WizardPage):

    def __init__(self,parent=None):
        commonqt.WizardPage.__init__(self, parent)

        self.label = QtGui.QLabel('Specify the location of the result you want to view.',self)
        self.openwidget = OpenWidget(self)

        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.label)
        layout.addWidget(self.openwidget)
        layout.addStretch()
        self.setLayout(layout)
        self.connect(self.openwidget, QtCore.SIGNAL('onCompleteStateChanged()'), self.completeStateChanged)

    def isComplete(self):
        return self.openwidget.isComplete()

    def saveData(self,mustbevalid):
        if not mustbevalid: return True
        try:
            result = self.openwidget.getResult()
        except Exception,e:
            QtGui.QMessageBox.critical(self, 'Unable to load result', str(e), QtGui.QMessageBox.Ok, QtGui.QMessageBox.NoButton)
            return False
        self.parent().shared['result'] = result
        self.parent().shared['scenario'] = result.scenario
        return True

class ConfigureReportWidget(QtGui.QWidget):
    def __init__(self,parent,result,rep):
        QtGui.QWidget.__init__(self,parent)
        
        self.result = result
        self.report = rep

        self.factory = commonqt.PropertyEditorFactory(self.report.store)

        reportname2path = data.Result.getReportTemplates()

        self.labTemplates = QtGui.QLabel('Report template:',self)
        self.comboTemplates = QtGui.QComboBox(parent)
        for (name,path) in reportname2path.items():
            self.comboTemplates.addItem(name,QtCore.QVariant(path))
        
        self.labOutput = QtGui.QLabel('Directory to save to:',self)
        self.pathOutput = commonqt.PathEditor(self,getdirectory=True)
        
        # Default report directory: result or scenario directory
        if self.result.path!=None:
            self.pathOutput.defaultpath = os.path.dirname(self.result.path)
        elif self.result.scenario!=None and self.result.scenario.path!=None:
            self.pathOutput.defaultpath = os.path.dirname(self.result.scenario.path)

        self.labVariables = QtGui.QLabel('Included variables:',self)
        self.treestore = self.result.getVariableTree('schemas/outputtree.xml')
        
        # Prepare selection based on report settings
        selroot = self.report.store.root.getLocation(['Figures','Selection'])
        for node in selroot.children:
            targetnode = self.treestore.root.getLocation(node.getValue().split('/'))
            if targetnode!=None: targetnode.setValue(True)
        
        self.model = commonqt.PropertyStoreModel(self.treestore,nohide=False,novalues=True,checkboxes=True)
        self.treeVariables = commonqt.ExtendedTreeView(self)
        self.treeVariables.header().hide()
        self.treeVariables.setModel(self.model)
        self.treeVariables.setSizePolicy(QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding)

        # Create labels+editors for figure settings
        self.editWidth = self.factory.createEditor(['Figures','Width'],self)
        self.labWidth = self.editWidth.createLabel()
        self.editHeight = self.factory.createEditor(['Figures','Height'],self)
        self.labHeight = self.editHeight.createLabel()
        self.editDpi = self.factory.createEditor(['Figures','Resolution'],self)
        self.labDpi = self.editDpi.createLabel()
        self.editFontScaling = self.factory.createEditor(['Figures','FontScaling'],self)
        self.labFontScaling = self.editFontScaling.createLabel()
        
        layout = QtGui.QGridLayout()
        layout.addWidget(self.labOutput,     0,0)
        layout.addWidget(self.pathOutput,    0,1)
        layout.addWidget(self.labTemplates,  1,0)
        layout.addWidget(self.comboTemplates,1,1)
        layout.addWidget(self.labVariables,  2,0,QtCore.Qt.AlignTop)
        layout.addWidget(self.treeVariables, 2,1)

        self.figbox = QtGui.QGroupBox('Figure settings',self)
        figlayout = QtGui.QGridLayout()
        figlayout.addWidget(self.labWidth,              3,0)
        figlayout.addWidget(self.editWidth.editor,      3,1)
        figlayout.addWidget(self.labHeight,             4,0)
        figlayout.addWidget(self.editHeight.editor,     4,1)
        figlayout.addWidget(self.labDpi,                5,0)
        figlayout.addWidget(self.editDpi.editor,        5,1)
        figlayout.addWidget(self.labFontScaling,        6,0)
        figlayout.addWidget(self.editFontScaling.editor,6,1)
        figlayout.setColumnStretch(1,1)
        self.figbox.setLayout(figlayout)
        layout.addWidget(self.figbox,3,0,1,2)
        
        layout.setMargin(0)
        self.setLayout(layout)

        self.connect(self.pathOutput, QtCore.SIGNAL('onChanged()'), self.completeStateChanged)

    def completeStateChanged(self):
        self.emit(QtCore.SIGNAL('onCompleteStateChanged()'))

    def isComplete(self):
        return self.pathOutput.hasPath()

    def generate(self):
        # Get path of target directory and template.
        templateindex = self.comboTemplates.currentIndex()
        templatepath = unicode(self.comboTemplates.itemData(templateindex).toString())
        outputpath = self.pathOutput.path()

        # Warn if the target directory is not empty.
        if os.path.isdir(outputpath) and len(os.listdir(outputpath))>0:
            ret = QtGui.QMessageBox.warning(self,'Directory is not empty','The specified target directory ("%s") contains one or more files, which may be overwritten. Do you want to continue?' % outputpath,QtGui.QMessageBox.Yes,QtGui.QMessageBox.No)
            if ret==QtGui.QMessageBox.No: return False

        # Update the list of selected variables.
        selroot = self.report.store.root.getLocation(['Figures','Selection'])
        selroot.removeAllChildren()
        for node in self.model.getCheckedNodes():
            if node.canHaveValue():
                ch = selroot.addChild('VariablePath')
                ch.setValue('/'.join(node.location))

        # Make changed report settings persistent
        self.factory.updateStore()

        # Generate the report and display the wait cursor while doing so.
        QtGui.QApplication.setOverrideCursor(QtGui.QCursor(QtCore.Qt.WaitCursor))
        self.report.generate(self.result,outputpath,templatepath,callback=self.onReportProgressed)
        QtGui.QApplication.restoreOverrideCursor()

        return True

    def onReportProgressed(self,progressed,status):
        self.emit(QtCore.SIGNAL('onReportProgressed'),progressed,status)
        
class PageReportGenerator(commonqt.WizardPage):
    def __init__(self,parent=None):
        commonqt.WizardPage.__init__(self, parent)

        self.result = parent.shared['result']
        self.report = report.Report()
        
        # Copy report settings from result.
        self.report.store.root.copyFrom(self.result.store.root.getLocation(['ReportSettings']),replace=True)

        self.label = QtGui.QLabel('You can generate a report that describes the scenario and the simulation results. A report consists of an HTML file, associated files (CSS, javascript) and image files for all figures.',self)
        self.label.setWordWrap(True)
        self.checkReport = QtGui.QCheckBox('Yes, I want to generate a report.', parent)
        self.reportwidget = ConfigureReportWidget(self,self.result,self.report)

        self.progressbar = QtGui.QProgressBar(self)
        self.progressbar.setRange(0,100)
        self.labStatus = QtGui.QLabel(self)
        self.progressbar.hide()
        self.labStatus.hide()

        layout = QtGui.QGridLayout()
        layout.addWidget(self.label,0,0,1,2)
        layout.addWidget(self.checkReport,1,0,1,2)
        layout.addWidget(self.reportwidget,2,1,1,1)

        layout.addWidget(self.progressbar,3,0,1,2)
        layout.addWidget(self.labStatus,4,0,1,2)

        layout.setRowStretch(5,1)
        layout.setColumnStretch(1,1)
        radiowidth = QtGui.QCheckBox().sizeHint().width()
        layout.setColumnMinimumWidth(0,radiowidth)

        self.setLayout(layout)

        self.connect(self.checkReport, QtCore.SIGNAL('stateChanged(int)'),        self.onCheckChange)
        self.connect(self.reportwidget,QtCore.SIGNAL('onCompleteStateChanged()'), self.completeStateChanged)
        self.connect(self.reportwidget,QtCore.SIGNAL('onReportProgressed'),       self.reportProgressed)
        self.onCheckChange()

    def onCheckChange(self):
        self.reportwidget.setVisible(self.checkReport.isChecked())
        self.completeStateChanged()

    def isComplete(self):
        if not self.checkReport.isChecked(): return True
        return self.reportwidget.isComplete()

    def saveData(self,mustbevalid):
        if mustbevalid and self.checkReport.isChecked():
            ret = self.reportwidget.generate()
            if ret:
                self.result.store.root.getLocation(['ReportSettings']).copyFrom(self.report.store.root,replace=True)
            return ret
        return True

    def reportProgressed(self,progressed,status):
        if self.progressbar.isHidden():
            self.label.setText('Please wait while the report is created...')
            self.checkReport.hide()
            self.reportwidget.hide()
            self.progressbar.show()
            self.labStatus.show()
            self.repaint()
            
        self.progressbar.setValue(round(progressed*100))
        self.labStatus.setText(status)
        QtGui.qApp.processEvents()

class PageSave(commonqt.WizardPage):

    def __init__(self,parent=None):
        commonqt.WizardPage.__init__(self, parent)

        self.result = parent.shared['result']

        self.label = QtGui.QLabel('Do you want to save the result of your simulation?',self)
        self.bngroup     = QtGui.QButtonGroup()
        self.radioNoSave = QtGui.QRadioButton('No, I do not want to save the result.', parent)
        self.radioSave   = QtGui.QRadioButton('Yes, I want to save the result to file.', parent)

        self.pathSave = commonqt.PathEditor(self,header='File to save to: ',save=True)
        self.pathSave.filter = 'GOTM result files (*.gotmresult);;NetCDF files (*.nc);;All files (*.*)'
        if self.result.path!=None: self.pathSave.setPath(self.result.path)

        self.checkboxAddFigures = QtGui.QCheckBox('Also save my figure settings.',self)
        self.checkboxAddFigures.setChecked(True)

        self.bngroup.addButton(self.radioNoSave, 0)
        self.bngroup.addButton(self.radioSave,   1)

        layout = QtGui.QGridLayout()
        layout.addWidget(self.label,      0,0,1,2)
        layout.addWidget(self.radioNoSave,1,0,1,2)
        layout.addWidget(self.radioSave,  2,0,1,2)
        layout.addWidget(self.pathSave,   3,1,1,1)
        layout.addWidget(self.checkboxAddFigures,4,1,1,1)

        layout.setRowStretch(5,1)
        layout.setColumnStretch(1,1)
        radiowidth = QtGui.QRadioButton().sizeHint().width()
        layout.setColumnMinimumWidth(0,radiowidth)

        self.setLayout(layout)

        self.connect(self.bngroup,  QtCore.SIGNAL('buttonClicked(int)'), self.onSourceChange)
        self.connect(self.pathSave, QtCore.SIGNAL('onChanged()'),        self.completeStateChanged)

        self.radioSave.setChecked(True)
        self.onSourceChange()

    def onSourceChange(self):
        checkedid = self.bngroup.checkedId()
        self.pathSave.setVisible(checkedid==1)
        self.checkboxAddFigures.setVisible(checkedid==1)
        self.completeStateChanged()

    def isComplete(self):
        checkedid = self.bngroup.checkedId()
        if   checkedid==0:
            return True
        elif checkedid==1:
            return self.pathSave.hasPath()

    def saveData(self,mustbevalid):
        if not mustbevalid: return True
        checkedid = self.bngroup.checkedId()
        if checkedid==1:
            targetpath = self.pathSave.path()
            if os.path.isfile(targetpath):
                ret = QtGui.QMessageBox.warning(self, 'Overwrite existing file?', 'There already exists a file at the specified location. Overwrite it?', QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
                if ret==QtGui.QMessageBox.No:
                    return False
            try:
                if targetpath.endswith('.nc'):
                    self.result.saveNetCDF(targetpath)
                else:
                    self.result.save(targetpath,addfiguresettings=self.checkboxAddFigures.isChecked())
                    self.owner.settings.addUniqueValue(('Paths','RecentResults'),'Path',targetpath)
            except Exception,e:
                print e
                QtGui.QMessageBox.critical(self, 'Unable to save result', str(e), QtGui.QMessageBox.Ok, QtGui.QMessageBox.NoButton)
                return False
        return True

    def doNotShow(self):
        return (not self.result.hasChanged())

class PageFinal(commonqt.WizardPage):
    
    def __init__(self,parent=None):
        commonqt.WizardPage.__init__(self, parent)

        self.label = QtGui.QLabel('You are now done.',self)

        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.label)
        layout.addStretch()
        self.setLayout(layout)

    def isComplete(self):
        return True

class PageVisualize(commonqt.WizardPage):
    
    def __init__(self,parent=None):
        commonqt.WizardPage.__init__(self, parent)

        self.varname = None

        self.result = parent.shared['result']
        self.treestore = self.result.getVariableTree('schemas/outputtree.xml')
        self.model = commonqt.PropertyStoreModel(self.treestore,nohide=False,novalues=True)

        self.treeVariables = commonqt.ExtendedTreeView(self)
        self.treeVariables.header().hide()
        self.connect(self.treeVariables, QtCore.SIGNAL('onSelectionChanged()'), self.OnVarSelected)
        self.treeVariables.setSizePolicy(QtGui.QSizePolicy.Minimum,QtGui.QSizePolicy.Expanding)
        self.treeVariables.setMaximumWidth(250)
        self.treeVariables.setModel(self.model)

        self.figurepanel = commonqt.FigurePanel(self)

        self.label = QtGui.QLabel('Here you can view the results of the simulation. Please choose a variable to be plotted from the menu.',self)

        layout = QtGui.QGridLayout()
        layout.addWidget(self.label,0,0,1,2)
        layout.addWidget(self.treeVariables,1,0)
        layout.addWidget(self.figurepanel,1,1)
        self.setLayout(layout)

        self.figurepanel.figure.addDataSource('result',self.result)

    def OnVarSelected(self):
        selected = self.treeVariables.selectedIndexes()
        if len(selected)==0: return
        node = selected[0].internalPointer()
        if node.hasChildren(): return
        
        varname = node.getId()
        varpath = '/'.join(node.location)
        QtGui.QApplication.setOverrideCursor(QtGui.QCursor(QtCore.Qt.WaitCursor))
        self.saveFigureSettings()
        self.varname = varname
        
        self.figurepanel.figure.setUpdating(False)
        
        props = self.figurepanel.figure.properties
        
        # Plot; first try stored figures, otherwise plot anew.
        if self.result.getFigure('result/'+varpath,props):
            # Set key properties to their required values.
            series = props.root.getLocation(['Data','Series'])
            series.getLocation(['Source']).setValue('result')
            series.getLocation(['Variable']).setValue(varname)
        else:
            self.figurepanel.plot(varname,'result')

        # Name the plot (later used as index for list of stored plot properties)
        props.root.getLocation(['Name']).setValue('result/'+varpath)
        
        self.figurepanel.figure.setUpdating(True)
        QtGui.QApplication.restoreOverrideCursor()

    def isComplete(self):
        return True

    def saveData(self,mustbevalid):
        self.saveFigureSettings()
        return True

    def saveFigureSettings(self):
        if self.varname!=None and self.figurepanel.figure.hasChanged():
            self.result.setFigure(self.figurepanel.figure.properties)

def main():
    # Debug info
    print 'Python version: '+str(sys.version_info)
    print 'PyQt4 version: '+QtCore.PYQT_VERSION_STR
    print 'Qt version: '+QtCore.qVersion()
    print 'xml version: '+xml.__version__

    # Create the application and enter the main message loop.
    createQApp = QtGui.QApplication.startingUp()
    if createQApp:
        app = QtGui.QApplication([' '])
    else:
        app = QtGui.qApp

    # Create wizard dialog
    wiz = commonqt.Wizard()
    wiz.setWindowTitle('Result visualizer')
    wiz.resize(800, 600)

    seq = [PageOpen,PageChooseAction,PageVisualize,PageReportGenerator,PageSave,PageFinal]

    # Get NetCDF file to open from command line or from FileOpen dialog.
    if len(sys.argv)>1:
        result = None
        try:
            result = loadResult(sys.argv[1])
        except Exception,e:
            QtGui.QMessageBox.critical(self, 'Unable to load result', unicode(e), QtGui.QMessageBox.Ok, QtGui.QMessageBox.NoButton)
        if result!=None:
            seq.pop(0)
            wiz.shared['result'] = result
            wiz.shared['scenario'] = result.scenario

    seq = commonqt.WizardSequence(seq)
    wiz.setSequence(seq)
    wiz.show()

    ret = app.exec_()
    page = None

    wiz.unlink()

    sys.exit(ret)

# If the script has been run (as opposed to imported), enter the main loop.
if (__name__=='__main__'): main()
