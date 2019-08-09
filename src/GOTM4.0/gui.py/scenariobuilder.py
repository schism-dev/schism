#!/usr/bin/python

#$Id: scenariobuilder.py,v 1.14 2007-04-18 09:31:09 jorn Exp $

from PyQt4 import QtGui,QtCore

import scenario, data, commonqt
import sys,xml, os.path

class ScenarioWidget(QtGui.QWidget):

    def __init__(self,parent=None,mrupaths=[]):
        QtGui.QWidget.__init__(self,parent)

        self.bngroup      = QtGui.QButtonGroup()
        self.radioNew     = QtGui.QRadioButton('Create a new scenario from a template.',self)
        self.radioOpen    = QtGui.QRadioButton('Open an existing scenario.',self)
        self.radioImport1 = QtGui.QRadioButton('Import a namelist-based scenario from an existing directory.',self)
        self.radioImport2 = QtGui.QRadioButton('Import a namelist-based scenario from a tar/gz archive.',self)

        self.labTemplate = QtGui.QLabel('Template:',self)
        default2path = scenario.Scenario.defaultname2path()
        self.comboTemplates = QtGui.QComboBox(parent)
        for (name,path) in default2path.items():
            self.comboTemplates.addItem(name,QtCore.QVariant(name))
        #self.comboTemplates.setSizePolicy(QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed)
        self.templatelayout = QtGui.QHBoxLayout()
        self.templatelayout.addWidget(self.labTemplate)
        self.templatelayout.addWidget(self.comboTemplates,1)
        self.templatelayout.addStretch()

        self.pathOpen    = commonqt.PathEditor(self,header='File to open: ',mrupaths=mrupaths)
        self.pathImport1 = commonqt.PathEditor(self,header='Directory to import: ',getdirectory=True)
        self.pathImport2 = commonqt.PathEditor(self,header='Archive to import: ')

        self.pathOpen.filter    = 'GOTM scenario files (*.gotmscenario);;GOTM result files (*.gotmresult);;dataless GOTM scenario files (*.xml);;All files (*.*)'
        self.pathImport2.filter = 'tar/gz files (*.tar.gz);;All files (*.*)'
        
        self.bngroup.addButton(self.radioNew,    0)
        self.bngroup.addButton(self.radioOpen,   1)
        self.bngroup.addButton(self.radioImport1,2)
        self.bngroup.addButton(self.radioImport2,3)

        layout = QtGui.QGridLayout()
        layout.addWidget(self.radioNew,       0,0,1,2)
        layout.addLayout(self.templatelayout, 1,1)
        layout.addWidget(self.radioOpen,      2,0,1,2)
        layout.addWidget(self.pathOpen,       3,1)
        layout.addWidget(self.radioImport1,   4,0,1,2)
        layout.addWidget(self.pathImport1,    5,1)
        layout.addWidget(self.radioImport2,   6,0,1,2)
        layout.addWidget(self.pathImport2,    7,1)

        layout.setColumnStretch(1,1)

        radiowidth = QtGui.QRadioButton().sizeHint().width()
        layout.setColumnMinimumWidth(0,radiowidth)

        layout.setMargin(0)
        
        self.setLayout(layout)

        self.connect(self.bngroup,     QtCore.SIGNAL('buttonClicked(int)'), self.onSourceChange)
        self.connect(self.pathOpen,    QtCore.SIGNAL('onChanged()'),        self.completeStateChanged)
        self.connect(self.pathImport1, QtCore.SIGNAL('onChanged()'),        self.completeStateChanged)
        self.connect(self.pathImport2, QtCore.SIGNAL('onChanged()'),        self.completeStateChanged)

        self.radioNew.setChecked(True)
        self.onSourceChange()

    def onSourceChange(self):
        self.setUpdatesEnabled(False)
        checkedid = self.bngroup.checkedId()
        self.labTemplate.setVisible(checkedid==0)
        self.comboTemplates.setVisible(checkedid==0)
        self.pathOpen.setVisible(checkedid==1)
        self.pathImport1.setVisible(checkedid==2)
        self.pathImport2.setVisible(checkedid==3)
        self.completeStateChanged()
        self.setUpdatesEnabled(True)

    def isComplete(self):
        checkedid = self.bngroup.checkedId()
        if   checkedid==0:
            return True
        elif checkedid==1:
            return self.pathOpen.hasPath()
        elif checkedid==2:
            return self.pathImport1.hasPath()
        elif checkedid==3:
            return self.pathImport2.hasPath()

    def getScenario(self):
        if not self.isComplete(): return None
        QtGui.QApplication.setOverrideCursor(QtGui.QCursor(QtCore.Qt.WaitCursor))
        try:
            checkedid = self.bngroup.checkedId()
            if   checkedid==0:
                index = self.comboTemplates.currentIndex()
                defname = unicode(self.comboTemplates.itemData(index).toString())
                defscenario = scenario.Scenario.getDefault(defname,scenario.guiscenarioversion)
                xmldom = defscenario.toXmlDom()
                scen = scenario.Scenario.fromSchemaName(scenario.guiscenarioversion)
                scen.setStore(xmldom)
            elif checkedid==1:
                path = self.pathOpen.path()
                if path.endswith('.gotmresult'):
                    try:
                        result = data.Result()
                        result.load(path)
                    except Exception,e:
                        raise Exception('An error occurred while loading the result: '+str(e))
                    scen = result.scenario
                    result.unlink()
                elif path.endswith('.xml'):
                    try:
                        scen = scenario.Scenario.fromSchemaName(scenario.guiscenarioversion)
                        scen.load(path)
                    except Exception,e:
                        raise Exception('An error occurred while loading the scenario: '+str(e))
                else:
                    try:
                        scen = scenario.Scenario.fromSchemaName(scenario.guiscenarioversion)
                        scen.loadAll(path)
                    except Exception,e:
                        raise Exception('An error occurred while loading the scenario: '+str(e))
            elif checkedid==2:
                try:
                    scen = scenario.Scenario.fromNamelists(self.pathImport1.path(),strict = False)
                except Exception,e:
                    raise Exception('Cannot parse namelist files. Error: '+str(e))
            elif checkedid==3:
                try:
                    scen = scenario.Scenario.fromNamelists(self.pathImport2.path(),strict = False)
                except Exception,e:
                    raise Exception('Cannot parse namelist files. Error: '+str(e))

        finally:
            QtGui.QApplication.restoreOverrideCursor()

        if checkedid!=0:
            # We have loaded a scenario from file. Look for empty nodes and reset these to their defaults.
            emptynodes = scen.root.getEmptyNodes()
            #emptynodes = [n for n in scen.root.getEmptyNodes() if not n.isHidden()]
            emptycount = len(emptynodes)
            if emptycount>0:
                QtGui.QMessageBox.information(self,'Scenario is incomplete','In this scenario %i variables do not have a value. These will be set to their default value.' % emptycount,QtGui.QMessageBox.Ok)
                scen.changed = True
            
        return scen

    def completeStateChanged(self):
        self.emit(QtCore.SIGNAL('onCompleteStateChanged()'))

class PageOpen(commonqt.WizardPage):

    def __init__(self,parent=None):
        commonqt.WizardPage.__init__(self, parent)

        self.label = QtGui.QLabel('How do you want to obtain a scenario?',self)
        self.scenariowidget = ScenarioWidget(self)
        self.connect(self.scenariowidget, QtCore.SIGNAL("onCompleteStateChanged()"),self.completeStateChanged)

        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.label)
        layout.addWidget(self.scenariowidget)
        layout.addStretch()
        self.setLayout(layout)

    def isComplete(self):
        return self.scenariowidget.isComplete()

    def saveData(self,mustbevalid):
        try:
            newscen = self.scenariowidget.getScenario()
        except Exception,e:
            QtGui.QMessageBox.critical(self, 'Unable to obtain scenario', str(e), QtGui.QMessageBox.Ok, QtGui.QMessageBox.NoButton)
            return False
        if 'result' in self.parent().shared:
            result = self.parent().shared.pop('result')
            result.unlink()
        if 'scenario' in self.parent().shared:
            oldscen = self.parent().shared['scenario']
            if oldscen!=None: oldscen.unlink()
        self.parent().shared['scenario'] = newscen
        return True

class PageLocation(commonqt.WizardPage):
    
    def __init__(self,parent=None):
        commonqt.WizardPage.__init__(self, parent)

        self.scenario = parent.shared['scenario']
        if self.scenario==None: raise Exception('No scenario available; this page should not have been available.')

        groupbox1 = QtGui.QGroupBox('Geographic location',self)

        self.factory = commonqt.PropertyEditorFactory(self.scenario)

        loclayout = QtGui.QGridLayout()
        self.lineName = self.factory.createEditor(['station','name'],self)
        self.labName = self.lineName.createLabel()
        self.lineLongitude = self.factory.createEditor(['station','longitude'],self)
        self.labLongitude = self.lineLongitude.createLabel()
        self.lineLatitude = self.factory.createEditor(['station','latitude'],self)
        self.labLatitude = self.lineLatitude.createLabel()
        self.lineDepth = self.factory.createEditor(['station','depth'],self)
        self.labDepth = self.lineDepth.createLabel()
        loclayout.addWidget(self.labName, 0,0)
        loclayout.addWidget(self.lineName.editor,0,1)
        loclayout.addWidget(self.labLongitude, 1,0)
        loclayout.addWidget(self.lineLongitude.editor,1,1)
        loclayout.addWidget(self.labLatitude, 2,0)
        loclayout.addWidget(self.lineLatitude.editor,2,1)
        loclayout.addWidget(self.labDepth, 3,0)
        loclayout.addWidget(self.lineDepth.editor,3,1)
        groupbox1.setLayout(loclayout)

        groupbox2 = QtGui.QGroupBox('Simulated period',self)

        periodlayout = QtGui.QGridLayout()
        self.lineStart = self.factory.createEditor(['time','start'],self)
        self.labStart = self.lineStart.createLabel()
        self.lineStop  = self.factory.createEditor(['time','stop'] ,self)
        self.labStop = self.lineStop.createLabel()
        periodlayout.addWidget(self.labStart, 0,0)
        periodlayout.addWidget(self.lineStart.editor,0,1)
        periodlayout.addWidget(self.labStop, 1,0)
        periodlayout.addWidget(self.lineStop.editor,1,1)
        groupbox2.setLayout(periodlayout)

        layout = QtGui.QVBoxLayout()
        layout.addWidget(groupbox1)
        layout.addWidget(groupbox2)
        layout.addStretch()
        self.setLayout(layout)

    def saveData(self,mustbevalid):
        if not self.factory.hasChanged(): return True
        
        if not mustbevalid:
            res = QtGui.QMessageBox.question(self,'Your scenario has changed','Do you want to preserve your changes?',QtGui.QMessageBox.Yes,QtGui.QMessageBox.No,QtGui.QMessageBox.NoButton)
            if res==QtGui.QMessageBox.No: return True

        self.factory.updateStore()
        return True

    def isComplete(self):
        return True

class PageDiscretization(commonqt.WizardPage):
    
    def __init__(self,parent=None):
        commonqt.WizardPage.__init__(self, parent)

        self.scenario = parent.shared['scenario']
        assert self.scenario!=None, 'No scenario available; this page should not have been available.'

        groupbox1 = QtGui.QGroupBox('Column structure',self)

        self.factory = commonqt.PropertyEditorFactory(self.scenario)

        gridlayout = QtGui.QGridLayout()
        self.bngroup = QtGui.QButtonGroup()
        self.radioEqual     = QtGui.QRadioButton('equal layers heights, optional zooming.',self)
        self.radioSigma     = QtGui.QRadioButton('custom sigma grid (relative depth fractions).',self)
        self.radioCartesian = QtGui.QRadioButton('custom cartesian grid (absolute layer heights).',self)

        self.editZoomSurface = self.factory.createEditor(['grid','ddu'],self)      
        self.editZoomBottom  = self.factory.createEditor(['grid','ddl'],self)      
        
        gridlayout.addWidget(self.radioEqual, 0,0,1,4)
        gridlayout.addWidget(self.editZoomSurface.createLabel(), 1,1)
        gridlayout.addWidget(self.editZoomSurface.editor, 1,2)
        gridlayout.addWidget(self.editZoomBottom.createLabel(), 2,1)
        gridlayout.addWidget(self.editZoomBottom.editor, 2,2)
        gridlayout.addWidget(self.radioSigma,3,0,1,4)
        gridlayout.addWidget(self.radioCartesian, 4,0,1,4)
        self.bngroup.addButton(self.radioEqual,     0)
        self.bngroup.addButton(self.radioSigma,     1)
        self.bngroup.addButton(self.radioCartesian, 2)
        groupbox1.setLayout(gridlayout)

        gridlayout.setColumnStretch(3,1)
        radiowidth = QtGui.QRadioButton().sizeHint().width()
        gridlayout.setColumnMinimumWidth(0,radiowidth)
        gridlayout.setMargin(0)

        layout = QtGui.QVBoxLayout()
        layout.addWidget(groupbox1)
        layout.addStretch()
        self.setLayout(layout)

    def saveData(self,mustbevalid):
        #checkedid = self.bngroup.checkedId()

        #if not self.factory.hasChanged(): return True
        
        #if not mustbevalid:
        #    res = QtGui.QMessageBox.question(self,'Your scenario has changed','Do you want to preserve your changes?',QtGui.QMessageBox.Yes,QtGui.QMessageBox.No,QtGui.QMessageBox.NoButton)
        #    if res==QtGui.QMessageBox.No: return True

        #self.factory.updateStore()
        return True

    def isComplete(self):
        return True

class PageAdvanced(commonqt.WizardPage):
    
    def __init__(self,parent=None):
        commonqt.WizardPage.__init__(self, parent)

        self.scenario = parent.shared['scenario']
        if self.scenario==None: raise Exception('No scenario available; this page should not have been available.')
        
        self.model = commonqt.PropertyStoreModel(self.scenario,nohide=False)

        self.tree = commonqt.ExtendedTreeView(self)
        #self.tree.header().hide()
        self.delegate = commonqt.PropertyDelegate()
        self.tree.setItemDelegate(self.delegate)
        self.tree.setModel(self.model)
        self.tree.setExpandedAll(maxdepth=1)
        self.tree.expandNonDefaults()

        lab = QtGui.QLabel('Here you can change all properties of your scenario. Values that differ from the default are shown in bold; they can be reset to their default through the context menu that opens after clicking the right mousebutton.',self)
        lab.setWordWrap(True)

        layout = QtGui.QVBoxLayout()
        layout.addWidget(lab)
        layout.addWidget(self.tree)
        self.setLayout(layout)

        self.connect(self.model, QtCore.SIGNAL('dataChanged(const QModelIndex&,const QModelIndex&)'),self.completeStateChanged)

    def showEvent(self,event):
        self.tree.header().resizeSection(0,.65*self.tree.width())

    def isComplete(self):
        return True

    def saveData(self,mustbevalid):
        if mustbevalid:
            # Find used file nodes that have not been supplied with data.
            filenodes = []
            for fn in self.scenario.root.getNodesByType('file'):
                if fn.isHidden(): continue
                value = fn.getValueOrDefault()
                if value==None or not value.isValid(): filenodes.append(fn)
            if len(filenodes)>0:
                vartext = '\n'.join([fn.getText(2) for fn in filenodes])
                QtGui.QMessageBox.critical(self,'Scenario is incomplete','The following variables will be used in the simulation, but have not been set:\n\n%s\n\nEither configure the scenario to not use these variables, or supply them with data.' % vartext,QtGui.QMessageBox.Ok,QtGui.QMessageBox.NoButton)
                return False

            # Find used nodes that have not been set, and lack a default value.
            errornodes = []
            for node in self.scenario.root.getEmptyNodes():
                if node.isHidden(): continue
                if node.getDefaultValue()==None: errornodes.append(node)
            if len(errornodes)>0:
                vartext = '\n'.join([node.getText(2) for node in errornodes])
                QtGui.QMessageBox.critical(self,'Scenario is incomplete','The following variables will be used in the simulation, but have not been set to a value:\n\n%s\n\nPlease set these variables to a value first.' % vartext,QtGui.QMessageBox.Ok,QtGui.QMessageBox.NoButton)
                return False
        return True
    
class PageSave(commonqt.WizardPage):

    def __init__(self,parent=None):
        commonqt.WizardPage.__init__(self, parent)

        self.scenario = parent.shared['scenario']

        self.label = QtGui.QLabel('The scenario has been modified. Do you want to save it?',self)
        self.bngroup     = QtGui.QButtonGroup()
        self.radioNoSave = QtGui.QRadioButton('No, I do not want to save the modified scenario.', parent)
        self.radioSave   = QtGui.QRadioButton('Yes, I want to save the scenario to file.', parent)

        self.pathSave = commonqt.PathEditor(self,header='File to save to: ',save=True)
        self.pathSave.filter = 'GOTM scenario files (*.gotmscenario);;All files (*.*)'

        if self.scenario.path!=None:
            self.pathSave.setPath(self.scenario.path)

        self.bngroup.addButton(self.radioNoSave, 0)
        self.bngroup.addButton(self.radioSave,   1)

        layout = QtGui.QGridLayout()
        layout.addWidget(self.label,0,0,1,2)
        layout.addWidget(self.radioNoSave,1,0,1,2)
        layout.addWidget(self.radioSave,2,0,1,2)
        layout.addWidget(self.pathSave,3,1,1,1)

        layout.setRowStretch(4,1)
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
                self.scenario.saveAll(targetpath)
            except Exception,e:
                QtGui.QMessageBox.critical(self, 'Unable to save scenario', str(e), QtGui.QMessageBox.Ok, QtGui.QMessageBox.NoButton)
                return False
            self.owner.settings.addUniqueValue(('Paths','RecentScenarios'),'Path',targetpath)
        return True

    def doNotShow(self):
        return (not self.scenario.hasChanged())

class PageFinal(commonqt.WizardPage):
    
    def __init__(self,parent=None):
        commonqt.WizardPage.__init__(self, parent)

        self.label = QtGui.QLabel('Your scenario is now complete.',self)

        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.label)
        layout.addStretch()
        self.setLayout(layout)

    def isComplete(self):
        return True

class SequenceEditScenario(commonqt.WizardSequence):
    def __init__(self):
        #commonqt.WizardSequence.__init__(self,[PageLocation,PageDiscretization,PageAdvanced,PageSave])
        commonqt.WizardSequence.__init__(self,[PageLocation,PageAdvanced,PageSave])

def main():
    # Debug info
    print 'Python version: '+str(sys.version_info)
    print 'PyQt4 version: '+QtCore.PYQT_VERSION_STR
    print 'Qt version: '+QtCore.qVersion()
    print 'xml version: '+xml.__version__

    # Create the application and enter the main message loop.
    createQApp = QtGui.QApplication.startingUp()
    if createQApp:
        app = QtGui.QApplication([" "])
    else:
        app = QtGui.qApp

    # Create wizard dialog
    wiz = commonqt.Wizard()
    wiz.setWindowTitle('Scenario builder')
    wiz.resize(800, 600)

    seq = commonqt.WizardSequence([PageOpen,SequenceEditScenario(),PageFinal])
    wiz.setSequence(seq)
    wiz.show()

    ret = app.exec_()
    page = None

    wiz.unlink()

    sys.exit(ret)

# If the script has been run (as opposed to imported), enter the main loop.
if (__name__=='__main__'): main()
