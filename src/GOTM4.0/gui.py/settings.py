import xmlstore
import sys, os.path

class SettingsStore(xmlstore.TypedStore):
    def __init__(self):
        settingspath = self.getSettingsPath()
        if not os.path.isfile(settingspath): settingspath = None
        xmlstore.TypedStore.__init__(self,'schemas/settings/gotmgui.xml',settingspath)
        
        self.removeNonExistent(('Paths','RecentScenarios'),'Path')
        self.removeNonExistent(('Paths','RecentResults'),  'Path')

    @staticmethod    
    def getSettingsPath():
        if sys.platform == 'win32':
            from win32com.shell import shellcon, shell
            appdata = shell.SHGetFolderPath(0, shellcon.CSIDL_APPDATA, 0, 0)
            return os.path.join(appdata,'GOTM','settings.xml')
        else:
            return os.path.expanduser('~/.gotm.gui')

    def save(self):
        if not self.changed: return
        settingspath = self.getSettingsPath()
        settingsdir = os.path.dirname(settingspath)
        if not os.path.isdir(settingsdir): os.mkdir(settingsdir)
        xmlstore.TypedStore.save(self,settingspath)
        
    def removeNonExistent(self,parentlocation,nodename):
        """Removes nodes below specified location if their value is not
        a path to an existing file. Used to filter defunct most-recently-used
        files.
        """
        parent = self.root.getLocation(parentlocation)
        currentnodes = parent.getLocationMultiple([nodename])
        for i in range(len(currentnodes)-1,-1,-1):
            path = currentnodes[i].getValue()
            if not os.path.isfile(path):
                parent.removeChildren(nodename,first=i,last=i)

    def addUniqueValue(self,parentlocation,nodename,nodevalue):
        parent = self.root.getLocation(parentlocation)
        currentnodes = parent.getLocationMultiple([nodename])
        if len(currentnodes)>0:
            maxcount = int(currentnodes[0].templatenode.getAttribute('maxoccurs'))
            for i in range(len(currentnodes)-1,-1,-1):
                if currentnodes[i].getValue()==nodevalue:
                    parent.removeChildren(nodename,first=i,last=i)
                    currentnodes.pop(i)
            parent.removeChildren(nodename,first=maxcount-1)
        newnode = parent.addChild(nodename,position=0)
        newnode.setValue(nodevalue)