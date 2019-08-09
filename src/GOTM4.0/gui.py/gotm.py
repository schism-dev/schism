#!/usr/bin/python

from PyQt4 import QtGui,QtCore

import commonqt, scenario
import sys,xml,os

import scenariobuilder,simulator,visualizer

class PageIntroduction(commonqt.WizardPage):
    
    def __init__(self,parent=None):
        commonqt.WizardPage.__init__(self, parent)

        # For version only:
        import matplotlib,numpy,gotm,pynetcdf
        #import Numeric,pycdf

        versions = []
        versions.append(('Python','%i.%i.%i %s %i' % sys.version_info))
        versions.append(('Qt4',QtCore.qVersion()))
        versions.append(('PyQt4',QtCore.PYQT_VERSION_STR))
        versions.append(('numpy',numpy.__version__))
        versions.append(('matplotlib',matplotlib.__version__))
        #versions.append(('Numeric',Numeric.__version__))
        #versions.append(('pycdf',pycdf.pycdfVersion()))
        versions.append(('gotm',gotm.gui_util.getversion().rstrip()))
        
        strversions = '<table cellspacing="0" cellpadding="0">'
        for v in versions:
            strversions += '<tr><td>%s</td><td>&nbsp;</td><td>%s</td></tr>' % v
        strversions += '</table>'

        layout = QtGui.QVBoxLayout()

        self.label = QtGui.QLabel( \
            """<p>This is the Graphical User Interface to the <a href="http://www.gotm.net">General Ocean Turbulence Model (GOTM)</a>.</p>

<p>GOTM is a one-dimensional water column model for natural (marine and limnic) waters based on the Reynolds-averaged Navier-Stokes equations. Vertical mixing is  included through an extensive library of state-of-the-art turbulence closure models. The hydrodynamics may be forced by wind stresses, surface heat and buoyancy fluxes, solar radiation and prescribed external and internal pressure gradients.</p>

<p>GOTM includes also a library of ecosystem models, ranging from simple settling of suspended matter to low-, medium- and high-complexity biogeochemical formulations.</p>

<p>There is a number of ready-to-use scenarios available for GOTM, with coastal, shelf sea, open ocean and limnic applications, a few of them including ecosystem modelling. These can be downloaded from <a href="http://www.gotm.net/index.php?go=software&page=testcases">the GOTM web site</a>.</p>

<p>This program offers a user-friendly interface to all options supported by GOTM. It allows you to run existing test cases, or to create and configure a custom scenario. The program will guide you step by step through the process of setting up a scenario, doing the calculations and displaying the results.</p>

<p>For any questions, please consult <a href="http://www.gotm.net">www.gotm.net</a> or write an email to <a href="mailto:gotm-users@gotm.net">gotm-users@gotm.net</a>.<br></p>

<p>This program was developed by <a href="mailto:jorn.bruggeman@xs4all.nl">Jorn Bruggeman</a> from funding by <a href="http://www.bolding-burchard.com">Bolding & Burchard Hydrodynamics</a>.<br></p>
""",self)
        self.label.setWordWrap(True)
        self.label.setOpenExternalLinks(True)
        layout.addWidget(self.label)

        layout.addStretch()

        #strversions = '\n'.join(['%s %s' % v for v in versions])
        #self.labelVersions = QtGui.QLabel('Module versions:\n'+strversions,self)
        #layout.addWidget(self.labelVersions)

        layout.addStretch(1)

        self.labelVersions = QtGui.QLabel('Module versions:',self)
        layout.addWidget(self.labelVersions)
        
        self.textVersions = QtGui.QTextEdit(strversions,self)
        self.textVersions.setMaximumHeight(100)
        self.textVersions.setReadOnly(True)
        layout.addWidget(self.textVersions)

        self.setLayout(layout)

    def isComplete(self):
        return True

class PageChooseAction(commonqt.WizardPage):
    
    def __init__(self,parent=None):
        commonqt.WizardPage.__init__(self, parent)

        pathnodes = self.parent().settings.root.getLocationMultiple(['Paths','RecentScenarios','Path'])
        mruscenarios = [p.getValue() for p in pathnodes]

        pathnodes = self.parent().settings.root.getLocationMultiple(['Paths','RecentResults','Path'])
        mruresults = [p.getValue() for p in pathnodes]

        self.label = QtGui.QLabel('What would you like to do?',self)
        self.radioScenario = QtGui.QRadioButton('I want to create, view or edit a scenario.',self)
        self.radioResult = QtGui.QRadioButton('I want to view or process the result of a previous simulation.',self)
        self.scenariowidget = scenariobuilder.ScenarioWidget(self,mrupaths=mruscenarios)
        self.connect(self.scenariowidget, QtCore.SIGNAL('onCompleteStateChanged()'),self.completeStateChanged)
        self.resultwidget = visualizer.OpenWidget(self,mrupaths=mruresults)
        self.connect(self.resultwidget, QtCore.SIGNAL('onCompleteStateChanged()'),self.completeStateChanged)

        self.bngroup     = QtGui.QButtonGroup()
        self.bngroup.addButton(self.radioScenario,0)
        self.bngroup.addButton(self.radioResult,1)
        self.connect(self.bngroup, QtCore.SIGNAL('buttonClicked(int)'), self.onSourceChange)

        layout = QtGui.QGridLayout()
        layout.addWidget(self.label,0,0,1,2)
        layout.addWidget(self.radioScenario,1,0,1,2)
        layout.addWidget(self.scenariowidget,2,1,1,1)
        layout.addWidget(self.radioResult,3,0,1,2)
        layout.addWidget(self.resultwidget,4,1,1,1)
        
        radiowidth = QtGui.QRadioButton().sizeHint().width()
        layout.setColumnMinimumWidth(0,radiowidth)

        layout.setRowStretch(5,1)
        layout.setColumnStretch(1,1)
        
        self.setLayout(layout)

        self.radioScenario.setChecked(True)
        self.onSourceChange()

    def onSourceChange(self):
        checkedid = self.bngroup.checkedId()
        self.scenariowidget.setVisible(checkedid==0)
        self.resultwidget.setVisible(checkedid==1)
        self.completeStateChanged()

    def isComplete(self):
        checkedid = self.bngroup.checkedId()
        if checkedid==0:
            return self.scenariowidget.isComplete()
        elif checkedid==1:
            return self.resultwidget.isComplete()
        return False

    def saveData(self,mustbevalid):
        if not mustbevalid: return True
        checkedid = self.bngroup.checkedId()
        if checkedid==0:
            try:
                newscen = self.scenariowidget.getScenario()
            except Exception,e:
                QtGui.QMessageBox.critical(self, 'Unable to obtain scenario', str(e), QtGui.QMessageBox.Ok, QtGui.QMessageBox.NoButton)
                return False
            self.parent().shared['mainaction'] = 'scenario'
            self.owner.setProperty('result', None)
            self.owner.setProperty('scenario', newscen)

            # Add to list of most-recently-used scenarios
            if newscen.path!=None:
                self.owner.settings.addUniqueValue(('Paths','RecentScenarios'),'Path',newscen.path)
            
            return True
        if checkedid==1:
            self.parent().shared['mainaction'] = 'result'
            try:
                newresult = self.resultwidget.getResult()
            except Exception,e:
                QtGui.QMessageBox.critical(self, 'Unable to load result', str(e), QtGui.QMessageBox.Ok, QtGui.QMessageBox.NoButton)
                return False
            self.owner.setProperty('result', newresult)
            self.owner.setProperty('scenario', newresult.scenario)

            # Add to list of most-recently-used results
            if newresult.path!=None:
                self.owner.settings.addUniqueValue(('Paths','RecentResults'),'Path',newresult.path)

            return True
        return False

class ForkOnAction(commonqt.WizardFork):
    def getSequence(self):
        if self.wizard.shared['mainaction']=='scenario':
            return commonqt.WizardSequence([scenariobuilder.SequenceEditScenario(),simulator.PageProgress])
        else:
            return commonqt.WizardSequence([commonqt.WizardDummyPage])
def main():
    # Debug info
    print 'Python version: '+str(sys.version_info)
    print 'PyQt4 version: '+QtCore.PYQT_VERSION_STR
    print 'Qt version: '+QtCore.qVersion()
	
    # Create the application and enter the main message loop.
    createQApp = QtGui.QApplication.startingUp()
    if createQApp:
        app = QtGui.QApplication([' '])
    else:
        app = QtGui.qApp

    # Create wizard dialog
    wiz = commonqt.Wizard(closebutton = sys.platform!='win32')
    seq = commonqt.WizardSequence([PageIntroduction,PageChooseAction,ForkOnAction(wiz),visualizer.PageVisualize,visualizer.PageReportGenerator,visualizer.PageSave,visualizer.PageFinal])
    wiz.setSequence(seq)
    wiz.setWindowTitle('GOTM-GUI')
    wiz.resize(800, 600)
    wiz.show()

    ret = app.exec_()
    page = None

    wiz.unlink()

    sys.exit(ret)

# If the script has been run (as opposed to imported), enter the main loop.
if (__name__=='__main__'): main()
