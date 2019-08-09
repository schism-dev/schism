# This module contains two classes for storing data in an XML tree:
#
# - Store: this class stores variables at locations in the XML tree specified by strings,
#   (e.g. ["root","foo","foo"] for /root/foo/foo) and can natively handle strings, integers,
#   floats, booleans and datetime variables. Conversion between types and from/to the XML
#   storage format is supported if a value type is specified during get/set calls.
#
# - TypedStore: this class combines a XML schema file (XSD-like) with an XML tree in
#   which values are stored. The latter is in fact an encapsulated "Store" described above.
#   Value types therefore do not need to specified in each call; they are obtained from
#   the schema definition. Additionally, the TypedStore supports (conditional) hiding of
#   nodes, notifications before/after changing of node values and node visiblity,
#   a set of default values, arbitrary data streams that are stored aside the XML value
#   tree, encapsulating containers such as ZIP, and many other features.

import datetime
import xml.dom.minidom, os
import zipfile, tarfile, tempfile, shutil, StringIO

import common

# dateformat: date format used for storing datetime objects in XML.
#   Used in conversion of (XML) string to datetime, and vice versa.
dateformat = '%Y-%m-%d %H:%M:%S'

# ------------------------------------------------------------------------------------------
# Store
# ------------------------------------------------------------------------------------------

# Store: class for storing 'properties' (i.e name,value pairs) in
#   hierarchical structure, using in-memory XML DOM. All values are stored as
#   strings, since XML is text-based; strings are converted to and from other
#   types (date, int, float, bool) whenever necessary.
class Store:

    class DataType:
        def __init__(self):
            self.refcount = 1
        
        @staticmethod
        def load(node,context):
            assert False, 'This virtual method MUST be overwritten by the inheriting class.'

        def save(self,node,context):
            assert False, 'This virtual method MUST be overwritten by the inheriting class.'

        def preparePersist(self,node,context):
            pass

        def persist(self,node,context):
            pass

        def release(self):
            assert self.refcount>0
            self.refcount -= 1
            if self.refcount==0: self.unlink()

        def addref(self):
            assert self.refcount>0
            self.refcount += 1

        def unlink(self):
            pass
    
    # =========================================================================================
    # PROTECTED
    # =========================================================================================
    # __init__: constructor
    def __init__(self,xmldocument=None,xmlroot=None):
        if isinstance(xmldocument,basestring):
            assert xmlroot==None,'Path to XML file specified, but also a (already parsed!) root node was supplied. This combination is invalid'
            xmldocument = xml.dom.minidom.parse(xmldocument)

        self.xmldocument = xmldocument
        if xmlroot==None: xmlroot = xmldocument.documentElement
        self.xmlroot = xmlroot
        self.xmlnamespace = self.xmldocument.namespaceURI

        self.context = {}

        # Links data type names (strings) to their respective Python classes.
        #   This is particularly relevant for the later TypedStore, which
        #   uses these data type names in its template XML files. Note that this is extensible:
        #   one can simple add other data types to the self.filetypes dictionary; just make sure
        #   their classes are derived from Store.DataType.
        self.filetypes = {'string'  :unicode,
                          'int'     :int,
                          'float'   :float,
                          'bool'    :bool,
                          'datetime':datetime.datetime}

    # =========================================================================================
    # PROTECTED
    # =========================================================================================
    # getText: gets all text directly below an XML element; may consist of multiple text nodes.
    def getText(self,node):
        rc = ''
        for ch in node.childNodes:
            if ch.nodeType == ch.TEXT_NODE: rc += ch.data
        return rc

    # =========================================================================================
    # PROTECTED
    # =========================================================================================
    # getText: sets text directly below an XML element, using one text node
    #   replaces any existing child text nodes.
    def setText(self,node,text):
        for ch in node.childNodes:
            if ch.nodeType == ch.TEXT_NODE:
                node.removeChild(ch)
                ch.unlink()
        val = self.xmldocument.createTextNode(text)
        node.insertBefore(val,node.firstChild)

    # =========================================================================================
    # PUBLIC
    # =========================================================================================
    # setProperty: sets specified location (list of ancestor names) to specified value.
    #   autoconverts specified value to string format.
    def setProperty(self,location,value,valuetype=None):
        node = common.findDescendantNode(self.xmlroot,location,create=True)
        assert node!=None, 'Unable to create new child node at "%s".' % location
        return self.setNodeProperty(node,value,valuetype)

    # =========================================================================================
    # PUBLIC
    # =========================================================================================
    # addProperty: adds a node at specified location (list of ancestor names) with specified
    #   value. Autoconverts specified value to string format.
    def addProperty(self,location,value):
        parentloc = location[:]
        name = parentloc.pop()
        parent = common.findDescendantNode(self.xmlroot,parentloc,create=True)
        assert parent!=None, 'Unable to locate or create parent node for "%s".' % location
        node = self.xmldocument.createElementNS(parent.namespaceURI,name)
        parent.appendChild(node)
        self.setNodeProperty(node,value)

    # =========================================================================================
    # PUBLIC
    # =========================================================================================
    # setNodeProperty: sets specified node to specified value.
    #   autoconverts specified value to string format.
    def setNodeProperty(self,node,value,valuetype=None):
        if valuetype!=None:
            if isinstance(valuetype,basestring):
                assert valuetype in self.filetypes, 'unknown type "%s" requested.' % valuetype
                valuetype = self.filetypes[valuetype]
            if not isinstance(value,valuetype): value = valuetype(value)
        if isinstance(value,Store.DataType):
            return value.save(node,self.context)
        else:
            value = self.packvalue(value,valuetype)
            if self.getText(node)==value: return False
            self.setText(node,value)
            return True

    # =========================================================================================
    # PUBLIC
    # =========================================================================================
    # getProperty: gets value at specified location (list of ancestor names).
    #   autoconverts value to the type requested (otherwise value = string).
    def getProperty(self,location,valuetype=str):
        node = common.findDescendantNode(self.xmlroot,location)
        if node==None: return None
        return self.getNodeProperty(node,valuetype=valuetype)

    # =========================================================================================
    # PUBLIC
    # =========================================================================================
    # getNodeProperty: gets value at node.
    #   autoconverts value to the type requested (otherwise value = string).
    def getNodeProperty(self,node,valuetype=str):
        if isinstance(valuetype,basestring):
            assert valuetype in self.filetypes, 'unknown type "%s" requested.' % valuetype
            valuetype = self.filetypes[valuetype]
        if issubclass(valuetype,Store.DataType):
            return valuetype.load(node,self.context)
        else:
            return self.unpackvalue(self.getText(node),valuetype=valuetype)

    def preparePersistNode(self,node,valuetype):
        value = self.getNodeProperty(node,valuetype)
        assert isinstance(value,Store.DataType)
        value.preparePersist(node,self.context)

    def persistNode(self,node,valuetype):
        value = self.getNodeProperty(node,valuetype)
        assert isinstance(value,Store.DataType)
        value.persist(node,self.context)

    # =========================================================================================
    # PUBLIC
    # =========================================================================================
    # clearProperty: removes all nodes with specified location (list of ancestor names).
    def clearProperty(self,location):
        parentloc = location[:]
        name = parentloc.pop()
        parent = common.findDescendantNode(self.xmlroot,parentloc,create=False)
        if parent==None: return
        for ch in parent.childNodes:
            if ch.nodeType==ch.ELEMENT_NODE and ch.localName==name:
                parent.removeChild(ch)

    # =========================================================================================
    # PUBLIC
    # =========================================================================================
    # clearNodeProperty: removes specified node.
    def clearNodeProperty(self,node):
        node.parentNode.removeChild(node)

    # =========================================================================================
    # PUBLIC
    # =========================================================================================
    # save: saves the current property tree to an XML document.
    def save(self,path):
        self.xmldocument.writexml(file(path,'w'),encoding='utf-8')

    # =========================================================================================
    # PUBLIC
    # =========================================================================================
    # packvalue: converts a value to a string representation suitable for storing in XML.
    def packvalue(self,value,valuetype=None):
        if valuetype!=None:
            if isinstance(valuetype,basestring):
                assert valuetype in self.filetypes, 'unknown type "%s" requested.' % valuetype
                valuetype = self.filetypes[valuetype]
            if not isinstance(value,valuetype): value = valuetype(value)
        assert not isinstance(value,Store.DataType)
        if isinstance(value,datetime.datetime):
            return value.strftime(dateformat)
        elif isinstance(value,bool):
            if value: return 'True'
            else:     return 'False'
        else:
            return unicode(value)

    # =========================================================================================
    # PUBLIC
    # =========================================================================================
    # unpackvalue: converts string representation of a value to the desired type.
    def unpackvalue(self,value,valuetype=str):
        if isinstance(valuetype,basestring):
            assert valuetype in self.filetypes, 'unknown type "%s" requested.' % valuetype
            valuetype = self.filetypes[valuetype]
        assert not issubclass(valuetype,Store.DataType)
        if valuetype==datetime.datetime:
            return common.parsedatetime(value,dateformat)
        elif valuetype==bool:
            return (value=='True')
        else:
            return valuetype(value)

class StoreColor(Store.DataType):
    def __init__(self,red=None,green=None,blue=None):
        Store.DataType.__init__(self)
        self.red = red
        self.green = green
        self.blue = blue

    @staticmethod
    def load(node,context):
        strcolor = common.getNodeText(node)
        if len(strcolor)>0:
            assert len(strcolor)==7 and strcolor[0]=='#', 'Colors must have exactly 7 characters and start with #.'
            strcolor = strcolor[1:]
            return StoreColor(int(strcolor[0:2],16),int(strcolor[2:4],16),int(strcolor[4:6],16))
        else:
            return StoreColor()
            
    @staticmethod
    def fromNormalized(r,g,b):
        return StoreColor(int(r*255),int(g*255),int(b*255))

    def save(self,node,context):
        assert self.red  ==None or self.red  <=255, 'Color channel red values should not exceed 255.'
        assert self.green==None or self.green<=255, 'Color channel green values should not exceed 255.'
        assert self.blue ==None or self.blue <=255, 'Color channel blue values should not exceed 255.'
        if self.isValid():
            common.setNodeText(node,'#%02x%02x%02x' % (self.red,self.green,self.blue))
        else:
            common.removeNodeChildren(node)
        
    def isValid(self):
        return self.red!=None and self.green!=None and self.blue!=None
        
    def getNormalized(self):
        assert self.isValid(), 'Cannot convert color to normalized tuple because the color object is not valid.'
        return (self.red/255.,self.green/255.,self.blue/255.)
        
    def __str__(self):
        if self.isValid():
            return '(%i, %i, %i)' % (self.red,self.green,self.blue)
        else:
            return 'None'
            
    def __eq__(self,other):
        if not isinstance(other,StoreColor): return False
        return self.red==other.red and self.green==other.green and self.blue==other.blue
        
    def __ne__(self,other):
        return not self.__eq__(other)
            

# ------------------------------------------------------------------------------------------
# DataFile
# ------------------------------------------------------------------------------------------

class DataContainer:
    def __init__(self):
        self.refcount = 1
    
    def listFiles(self):
        return []

    def getItem(self,name):
        return None

    def addItem(self,datafile,newname=None):
        return None

    def addFile(self,path,newname=None):
        df = DataContainerDirectory.DataFileFile(path)
        df_added = self.addItem(df,newname)
        df_added.release()
        df.release()

    def persistChanges(self):
        pass

    def addref(self):
        assert self.refcount>0
        self.refcount += 1
        return self

    def release(self):
        assert self.refcount>0
        self.refcount -= 1
        if self.refcount==0: self.unlink()

    def unlink(self):
        pass

# DataFile: class that encapsulates a data block, which can be a file on disk, an item
#    in a zip or tar/gz archive, or a memory block. It can be used as data type in the
#    xml stores.
class DataFile(Store.DataType):
    def __init__(self,name=''):
        Store.DataType.__init__(self)
        self.name = None
    
    @staticmethod
    def load(node,context):
        assert 'container' in context
        if 'cache' in context:
            uniquename = DataFile.getUniqueNodeName(node)
            cache = context['cache']
            if uniquename in cache: return cache[uniquename]
        container = context['container']
        if container==None: return DataFile()
        name = common.getNodeText(node)
        if name not in container.listFiles(): return DataFile()
        return container.getItem(name)

    def save(self,node,context):
        if 'cache' not in context: context['cache'] = {}
        cache = context['cache']
        uniquename = DataFile.getUniqueNodeName(node)
        if uniquename in cache: cache[uniquename].release()
        cache[uniquename] = self
        if self.name!=None:
            common.setNodeText(node,self.name)
        else:
            common.setNodeText(node,'')

    def preparePersist(self,node,context):
        assert 'targetcontainerpath' in context, 'preparePersist: "targetcontainerpath" not set in XML store context.'
        if not self.isValid(): return
        targetpath = context['targetcontainerpath']
        if self.isBelowPath(targetpath):
            print 'Reading "%s" into memory to prevent it from being overwritten.' % self.name
            memdf = DataFileMemory.fromDataFile(self)
            memdf.save(node,context)

    def persist(self,node,context):
        assert 'targetcontainer' in context, 'persist: "targetcontainer" not set in XML store context.'
        if not self.isValid(): return
        targetcontainer = context['targetcontainer']
        newname = node.localName+'.dat'
        df = targetcontainer.addItem(self,newname)
        if ('donotclaimtarget' in context) and context['donotclaimtarget']:
            common.setNodeText(node,df.name)
            df.release()
        else:
            df.save(node,context)

    @staticmethod
    def getUniqueNodeName(node):
        path = []
        while node.parentNode.parentNode!=None:
            path.append(node.localName)
            node = node.parentNode
        return '/'.join(path)

    def isValid(self):
        return self.name!=None

    # Default implementation of getAsReadOnlyFile: uses getData
    def getAsReadOnlyFile(self,textmode=True):
        data = self.getData(textmode=textmode,readonly=True)
        return StringIO.StringIO(data)

    # Default implementation of getData: uses getAsReadOnlyFile
    def getData(self,textmode=True,readonly=False):
        f = self.getAsReadOnlyFile(textmode=textmode)
        data = f.read()
        f.close()
        return data

    # Default implementation of saveToFile: uses the read-only file object returned by getAsReadOnlyFile
    def saveToFile(self,targetpath):
        assert self.isValid(), 'saveToFile: DataFile is not valid (i.e., empty).'
        fsource = self.getAsReadOnlyFile(textmode=False)
        ftarget = open(targetpath,'wb')
        shutil.copyfileobj(fsource,ftarget)
        ftarget.close()
        fsource.close()

    # Default implementation of addToZip: uses the read-only file object returned by getAsReadOnlyFile
    def addToZip(self,zfile,filename):
        assert self.isValid(), 'addToZip: DataFile is not valid (i.e., empty).'
        zfile.writestr(str(filename),self.getData(textmode=False,readonly=True))

    # Default implementation of addToTar: uses the read-only file object returned by getAsReadOnlyFile
    def addToTar(self,tfile,filename):
        assert self.isValid(), 'addToTar: DataFile is not valid (i.e., empty).'
        data = self.getData(textmode=False,readonly=True)
        tarinfo = tarfile.TarInfo(filename)
        tarinfo.size = len(data)
        import time
        tarinfo.mtime = time.time()
        tfile.addfile(tarinfo,StringIO.StringIO(data))

    # Default implementation of isBelowPath: returns False.
    def isBelowPath(self,path):
        return False

    def getSize(self):
        return None
        # Below an expensive way to get the size. Disabled: if the user really wants this,
        # he should do it himself.
        return len(self.getData(textmode=False,readonly=True))

class DataContainerDirectory(DataContainer):
    class DataFileFile(DataFile):
        def __init__(self,path):
            DataFile.__init__(self)
            assert os.path.isfile(path), 'Specified path "%s" is not a file.' % path
            self.path = path
            self.name = os.path.basename(self.path)

            # Create a lock on the file; while we live, its contents should persist.
            self.file = open(self.path,'r')

        def __del__(self):
            self.unlink()

        def getAsReadOnlyFile(self,textmode=True):
            if textmode:
                return open(self.path,'rU')
            else:
                return open(self.path,'rb')

        def saveToFile(self,targetpath):
            if self.path==targetpath: return
            shutil.copyfile(self.path,targetpath)

        def addToZip(self,zfile,filename):
            assert self.isValid()
            zfile.write(self.path,str(filename))

        def addToTar(self,tfile,filename):
            assert self.isValid()
            tfile.add(self.path,filename,recursive=False)

        def isBelowPath(self,path):
            owndir = os.path.normcase(os.path.dirname(self.path))
            return owndir.startswith(os.path.normcase(path))

        def getSize(self):
            return os.path.getsize(self.path)

        def unlink(self):
            if self.file==None: return
            self.file.close()
            self.file = None
   
    def __init__(self,path,create=False):
        DataContainer.__init__(self)
        assert os.path.isdir(path) or create
        if not os.path.isdir(path):
            try:
                os.mkdir(path)
            except Exception,e:
                raise Exception('Unable to create directory "%s". Error: %s' % (path,str(e)))
        self.path = path

    def getItem(self,name):
        sourcepath = os.path.join(self.path,name)
        if not os.path.isfile(sourcepath): return None
        return self.DataFileFile(sourcepath)

    def addItem(self,datafile,newname=None):
        assert datafile.isValid()
        if newname==None: newname = datafile.name
        targetpath = os.path.join(self.path,newname)
        datafile.saveToFile(targetpath)
        return self.DataFileFile(targetpath)

    def listFiles(self):
        res = []
        for fn in os.listdir(self.path):
            if os.path.isfile(os.path.join(self.path,fn)): res.append(fn)
        return res

class DataContainerZip(DataContainer):
    class DataFileZip(DataFile):
        def __init__(self,zipcontainer,name):
            DataFile.__init__(self)
            self.zipcontainer = zipcontainer
            self.zipcontainer.addref()
            self.name = name

        def __del__(self):
            self.unlink()

        def getData(self,textmode=True,readonly=False):
            assert self.zipcontainer!=None, 'DataFileZip.getData failed; ZIP file has been closed.'
            self.zipcontainer.setMode('r')
            return self.zipcontainer.zfile.read(self.name)

        def isBelowPath(self,path):
            assert self.zipcontainer!=None, 'DataFileZip.isBelowPath failed; ZIP file has been closed.'
            owndir = os.path.normcase(self.zipcontainer.path)
            return owndir.startswith(os.path.normcase(path))

        def getSize(self):
            return self.zipcontainer.zfile.getinfo(self.name).file_size

        def unlink(self):
            if self.zipcontainer==None: return
            self.zipcontainer.release()
            self.zipcontainer = None
    
    def __init__(self,path,mode='r'):
        DataContainer.__init__(self)
        if isinstance(path,basestring):
            assert os.path.isfile(path) or mode=='w', 'Cannot initialize DataContainerZip with supplied path; it does not point to an existing file, but is also not opened for writing.'
        elif isinstance(path,StringIO.StringIO):
            assert mode=='w', 'Can initialize DataContainerZip with StringIO object only in write-only mode.'
        elif isinstance(path,DataFile):
            assert mode=='r', 'Can initialize DataContainerZip with file-like object only in read-only mode.'
        else:
            assert False, 'Cannot initialize DataContainerZip with %s.' % path
        self.mode = None
        self.zfile = None
        self.path = path
        if isinstance(self.path,DataFile): self.path.addref()
        self.setMode(mode)

    def __del__(self):
        self.unlink()

    def unlink(self):
        if self.zfile==None: return
        self.zfile.close()
        self.zfile = None
        self.mode = None
        if isinstance(self.path,DataFile): self.path.release()
        self.path = None
    
    def getItem(self,name):
        if name not in self.listFiles(): return None
        return self.DataFileZip(self,name)

    def setMode(self,mode):
        assert mode in ('r','w','a'), 'DataContainerZip.setMode: mode must be "r", "w", or "a" (not "%s").' % mode
        if self.zfile!=None:
            if mode==self.mode: return
            self.zfile.close()
        self.mode = mode
        if isinstance(self.path,StringIO.StringIO):
            # Writing to in-memory data block.
            assert self.mode=='w', 'In-memory data blocks can only be written to, not read from.'
            self.zfile = zipfile.ZipFile(self.path,self.mode,zipfile.ZIP_DEFLATED)
        if isinstance(self.path,DataFile):
            # Reading from generic DataFile object.
            assert self.mode=='r', 'Data file objects can only be accessed as read-only zip file.'
            f = self.path.getAsReadOnlyFile()
            self.zfile = zipfile.ZipFile(f,self.mode,zipfile.ZIP_DEFLATED)
        else:
            # Reading from/writing to file.
            self.zfile = zipfile.ZipFile(self.path,self.mode,zipfile.ZIP_DEFLATED)

    def addItem(self,datafile,newname=None):
        if newname==None: newname = datafile.name
        if self.mode=='r': self.setMode('a')
        if isinstance(self.path,StringIO.StringIO):
            print 'Adding "%s" to in-memory archive...' % (newname,)
        else:
            print 'Adding "%s" to archive "%s"...' % (newname,self.path)
        datafile.addToZip(self.zfile,newname)
        return self.DataFileZip(self,newname)

    def listFiles(self):
        return self.zfile.namelist()

    def persistChanges(self):
        if self.zfile==None: return
        if isinstance(self.path,basestring):
            # Immediately re-open the ZIP file so we keep a lock on it.
            self.setMode('r')
        else:
            self.zfile.close()
            self.zfile = None

class DataContainerTar(DataContainer):
    class DataFileTar(DataFile):
        def __init__(self,tarcontainer,name):
            DataFile.__init__(self)
            self.tarcontainer = tarcontainer
            self.tarcontainer.addref()
            self.name = name

        def __del__(self):
            self.unlink()

        def getAsReadOnlyFile(self,textmode=True):
            assert self.tarcontainer!=None, 'DataFileTar.getAsReadOnlyFile failed; TAR file has been closed.'
            self.tarcontainer.setMode('r')
            return self.tarcontainer.tfile.extractfile(self.name)

        def isBelowPath(self,path):
            assert self.tarcontainer!=None, 'DataFileTar.isBelowPath failed; TAR file has been closed.'
            owndir = os.path.normcase(self.tarcontainer.path)
            return owndir.startswith(os.path.normcase(path))

        def getSize(self):
            return self.tarcontainer.tfile.getmember(self.name).size

        def unlink(self):
            if self.tarcontainer == None: return
            self.tarcontainer.release()
            self.tarcontainer = None

    def __init__(self,path,mode='r'):
        DataContainer.__init__(self)
        assert os.path.isfile(path) or mode=='w', 'Cannot initialize DataContainerTar with supplied path; it does not point to an existing file, but is also not opened for writing.'
        self.mode = None
        self.tfile = None
        self.path = path
        self.setMode(mode)

    def __del__(self):
        self.unlink()

    def unlink(self):
        if self.tfile==None: return
        self.tfile.close()
        self.mode = None
        self.tfile = None
    
    def getItem(self,name):
        if name not in self.listFiles(): return None
        return self.DataFileTar(self,name)

    def setMode(self,mode):
        assert mode in ('r','w'), 'DataContainerZip.setMode: mode must be "r" or "w" (not "%s").' % mode
        if mode=='w': mode='w:gz'
        if self.tfile!=None:
            if mode==self.mode: return
            self.tfile.close()
        self.mode = mode
        self.tfile = tarfile.open(str(self.path),self.mode)

    def addItem(self,datafile,newname=None):
        if newname==None: newname = datafile.name
        if self.mode=='r': self.setMode('a')
        print 'Adding "%s" to archive "%s"...' % (newname,self.path)
        datafile.addToTar(self.tfile,newname)
        return self.DataFileTar(self,newname)

    def listFiles(self):
        return self.tfile.getnames()

    def persistChanges(self):
        if self.tfile==None: return
        self.setMode('r')

class DataFileXmlNode(DataFile):
    def __init__(self,xmlnode,name=''):
        DataFile.__init__(self)
        self.xmlnode = xmlnode
        self.name = name

    def __del__(self):
        self.unlink()

    def getData(self,textmode=True,readonly=False):
        return self.xmlnode.toxml('utf-8')

    def saveToFile(self,targetpath):
        import codecs
        f = codecs.open(targetpath,'w','utf-8')
        self.xmlnode.writexml(f,encoding='utf-8')
        f.close()

    def unlink(self):
        self.xmlnode = None
        self.name = None

class DataFileMemory(DataFile):
    def __init__(self,data,name):
        DataFile.__init__(self)
        self.data = data
        self.name = name

    def __del__(self):
        self.unlink()

    @staticmethod
    def fromDataFile(df):
        return DataFileMemory(df.getData(),df.name)

    def getAsReadOnlyFile(self,textmode=True):
        return StringIO.StringIO(self.data)

    def getData(self,textmode=True,readonly=False):
        if readonly:
            return self.data
        else:
            return self.data[:]

    def unlink(self):
        self.data = None
        self.name = None

    def getSize(self):
        return len(self.data)

class TypedStoreInterface:
    def __init__(self,store,showhidden=True,omitgroupers=False):
        self.showhidden = showhidden
        self.omitgroupers = omitgroupers
        self.blockNotifyOfHiddenNodes = True
        self.notifyOnDefaultChange = True

        self.visibilityhandlers = []
        self.changehandlers = []
        self.beforechangehandlers = []
        self.storechangedhandlers = []

        store.connectInterface(self)

    def getChildCount(self,node):
        assert isinstance(node,TypedStore.Node), 'Supplied object is not of type "TypedStore.Node" (but "%s").' % node
        assert node.isValid(), 'Supplied node %s is invalid (has already been destroyed).' % node
        childcount = 0
        for child in node.children:
            if child.visible or self.showhidden:
                if self.omitgroupers and child.grouponly:
                    childcount += self.getChildCount(child)
                else:
                    childcount += 1
        return childcount

    def getChildren(self,node,showhidden=None):
        assert isinstance(node,TypedStore.Node), 'Supplied object is not of type "TypedStore.Node" (but "%s").' % node
        assert node.isValid(), 'Supplied node %s is invalid (has already been destroyed).' % node
        if showhidden==None: showhidden=self.showhidden
        res = []
        for child in node.children:
            if child.visible or showhidden:
                if self.omitgroupers and child.grouponly:
                    res += self.getChildren(child,showhidden=showhidden)
                else:
                    res.append(child)
        return res

    def getParent(self,node):
        assert isinstance(node,TypedStore.Node), 'Supplied object is not of type "TypedStore.Node" (but "%s").' % node
        assert node.isValid(), 'Supplied node %s is invalid (has already been destroyed).' % node
        if not self.omitgroupers: return node.parent
        par = node.parent
        while par.grouponly: par = par.parent
        return par

    def getChildByIndex(self,node,index,returnindex=False):
        assert isinstance(node,TypedStore.Node), 'Supplied object is not of type "TypedStore.Node" (but "%s").' % node
        assert node.isValid(), 'Supplied node %s is invalid (has already been destroyed).' % node
        for child in node.children:
            if child.visible or self.showhidden:
                if self.omitgroupers and child.grouponly:
                    index = self.getChildByIndex(child,index,returnindex=True)
                    if not isinstance(index,int): return index
                else:
                    if index==0: return child
                    index -= 1
        if returnindex:
            return index
        else:
            return None

    def getOwnIndex(self,node):
        assert isinstance(node,TypedStore.Node), 'Supplied object is not of type "TypedStore.Node" (but "%s").' % node
        assert node.isValid(), 'Supplied node %s is invalid (has already been destroyed).' % node
        ind = 0
        par = node.parent
        if self.omitgroupers and par.grouponly: ind = self.getOwnIndex(par)
        for (isib,sib) in enumerate(par.children):
            if sib is node or isib==node.futureindex: break
            if sib.visible or self.showhidden:
                if self.omitgroupers and sib.grouponly:
                    ind += self.getChildCount(sib)
                else:
                    ind += 1
        else:
            assert node.futureindex!=None, 'Could not find node "%s" in children of supposed parent, but future index was also not set. Data: %s' % (node,node.valuenode.toxml('utf-8'))
            assert node.futureindex==len(par.children), 'Could not find node "%s" in children of supposed parent, but future index (%i) was also not set to tailing position (%i).' % (node,node.futureindex,len(par.children))
        return ind

    def getDepth(self,node):
        assert isinstance(node,TypedStore.Node), 'Supplied object is not of type "TypedStore.Node" (but "%s").' % node
        assert node.isValid(), 'Supplied node %s is invalid (has already been destroyed).' % node
        childmax = 0
        for child in self.getChildren(node):
            curchilddepth = self.getDepth(child)
            if curchilddepth>childmax: childmax = curchilddepth
        return childmax+1

    def toHtml(self,node,xmldocument,totaldepth,level=0,hidedefaults=False):
        assert isinstance(node,TypedStore.Node), 'Supplied object is not of type "TypedStore.Node" (but "%s").' % node
        assert node.isValid(), 'Supplied node %s is invalid (has already been destroyed).' % node
        res = []

        tr = None
        if level>=0:
            tr = xmldocument.createElement('tr')

            for i in range(level):
                td = xmldocument.createElement('td')
                tr.appendChild(td)

            td1 = xmldocument.createElement('td')
            templatenode = node.templatenode
            td1.appendChild(xmldocument.createTextNode(node.getText(detail=1)))
            if level+1<totaldepth:
                td1.setAttribute('colspan',unicode(totaldepth-level))
            tr.appendChild(td1)

            td2 = xmldocument.createElement('td')
            if node.canHaveValue():
                val = node.getValueAsString(usedefault=True)
            else:
                val = ' '
            td2.appendChild(xmldocument.createTextNode(val))
            tr.appendChild(td2)

            res.append(tr)

        childtrs = []
        for child in self.getChildren(node):
            childnodes = self.toHtml(child,xmldocument,totaldepth,level+1,hidedefaults=hidedefaults)
            childtrs += childnodes
        res += childtrs

        if tr!=None and hidedefaults:
            isdefault = True
            if node.canHaveValue() and node.getValue()!=node.getDefaultValue():
                isdefault = False
            else:
                for childtr in childtrs:
                    if not childtr.hasAttribute('default'):
                        isdefault = False
                        break
            if isdefault:
                tr.setAttribute('style','display:none')
                tr.setAttribute('default','')

        return res

    # ---------------------------------------------------------------------------
    # Functions for connecting to events
    # ---------------------------------------------------------------------------

    def addVisibilityChangeHandler(self,beforecallback,aftercallback):
        self.visibilityhandlers.append((beforecallback,aftercallback))

    def addStoreChangedHandler(self,callback):
        self.storechangedhandlers.append(callback)

    def addChangeHandler(self,callback):
        self.changehandlers.append(callback)

    def addBeforeChangeHandler(self,callback):
        self.beforechangehandlers.append(callback)

    # ---------------------------------------------------------------------------
    # Functions called by store when events occur
    # ---------------------------------------------------------------------------

    def beforeVisibilityChange(self,node,shownew,showhide):
        assert isinstance(node,TypedStore.Node), 'Supplied object is not of type "TypedStore.Node" (but "%s").' % node
        assert node.isValid(), 'Supplied node %s is invalid (has already been destroyed).' % node
        if (self.getParent(node).isHidden() and self.blockNotifyOfHiddenNodes): return
        if self.blockNotifyOfHiddenNodes and node.isHidden() and not showhide: return
        if self.omitgroupers and node.grouponly:
            for ch in self.getChildren(node): self.emitBeforeVisibilityChange(ch,shownew,showhide)
        else:
            self.emitBeforeVisibilityChange(node,shownew,showhide)

    def afterVisibilityChange(self,node,shownew,showhide):
        assert isinstance(node,TypedStore.Node), 'Supplied object is not of type "TypedStore.Node" (but "%s").' % node
        assert node.isValid(), 'Supplied node %s is invalid (has already been destroyed).' % node
        if (self.getParent(node).isHidden() and self.blockNotifyOfHiddenNodes): return
        if self.blockNotifyOfHiddenNodes and node.isHidden() and not showhide: return
        if self.omitgroupers and node.grouponly:
            for ch in self.getChildren(node): self.emitAfterVisibilityChange(ch,shownew,showhide)
        else:
            self.emitAfterVisibilityChange(node,shownew,showhide)

    def afterStoreChange(self):
        for callback in self.storechangedhandlers: callback()

    def onBeforeChange(self,node,newvalue):
        assert isinstance(node,TypedStore.Node), 'Supplied object is not of type "TypedStore.Node" (but "%s").' % node
        assert node.isValid(), 'Supplied node %s is invalid (has already been destroyed).' % node
        if not (node.isHidden() and self.blockNotifyOfHiddenNodes):
            for callback in self.beforechangehandlers:
                if not callback(node,newvalue): return False
        return True

    def onChange(self,node):
        assert isinstance(node,TypedStore.Node), 'Supplied object is not of type "TypedStore.Node" (but "%s").' % node
        assert node.isValid(), 'Supplied node %s is invalid (has already been destroyed).' % node
        if node.isHidden() and self.blockNotifyOfHiddenNodes: return
        for callback in self.changehandlers: callback(node)

    def onDefaultChange(self,node):
        assert isinstance(node,TypedStore.Node), 'Supplied object is not of type "TypedStore.Node" (but "%s").' % node
        assert node.isValid(), 'Supplied node %s is invalid (has already been destroyed).' % node
        if self.notifyOnDefaultChange: self.onChange(node)

    # ---------------------------------------------------------------------------
    # Helper functions for emitting events
    # ---------------------------------------------------------------------------

    def emitBeforeVisibilityChange(self,node,visible,showhide):
        for callback in self.visibilityhandlers:
            if callback[0]!=None: callback[0](node,visible,showhide)

    def emitAfterVisibilityChange(self,node,visible,showhide):
        for callback in self.visibilityhandlers:
            if callback[1]!=None: callback[1](node,visible,showhide)

# ------------------------------------------------------------------------------------------
# TypedStore
# ------------------------------------------------------------------------------------------
            
# TypedStore: encapsulates the above store.
#   Adds the use of a second XML document (template) that describes the data types
#   of the nodes of the first DOM, and that describes dependencies between nodes.
#   Any node in the original document for which conditions are not met is hidden.
#   Nodes that are not described by the template are not allowed in the property store.
#   Node are obtained by traversing the tree (start: TypedStore.root).
class TypedStore:

    class Node:
        def __init__(self,controller,templatenode,valuenode,location,parent):
            assert templatenode.hasAttribute('id'),'Schema node %s lacks "id" attribute.' % location

            self.controller = controller
            self.store = controller.store
            self.templatenode = templatenode
            self.valuenode = valuenode
            self.location = location
            self.parent = parent
            self.children = []
            self.futureindex = None
            self.visible = (not self.templatenode.hasAttribute('hidden'))
            self.grouponly = self.templatenode.hasAttribute('grouponly')

            ids = []
            for templatechild in self.templatenode.childNodes:
                if templatechild.nodeType==templatechild.ELEMENT_NODE and templatechild.localName=='element':
                    childloc = self.location[:] + [templatechild.getAttribute('id')]
                    ids.append(childloc[-1])
                    
                    start = len(self.location)
                    root = self
                    while root.valuenode==None:
                        start -= 1
                        root = root.parent
                    
                    if templatechild.hasAttribute('maxoccurs'):
                        maxoccurs = int(templatechild.getAttribute('maxoccurs'))
                        valuechildren = common.findDescendantNodes(root.valuenode,childloc[start:])
                        assert len(valuechildren)<=maxoccurs, 'Number of children is greater than the imposed maximum (%i).' % maxoccurs
                        for valuechild in valuechildren:
                            self.children.append(TypedStore.Node(self.controller,templatechild,valuechild,childloc,parent=self))
                    else:
                        valuechild = common.findDescendantNode(root.valuenode,childloc[start:])                            
                        self.children.append(TypedStore.Node(self.controller,templatechild,valuechild,childloc,parent=self))

            if self.valuenode!=None:
                for ch in [ch for ch in self.valuenode.childNodes if ch.nodeType==ch.ELEMENT_NODE and ch.localName not in ids]:
                    print 'WARNING! Value "%s" below "%s" was unexpected and will be ignored.' % (ch.localName,self.location)
                    self.valuenode.removeChild(ch)

        def __str__(self):
            return str(self.location)

        def destroy(self):
            for ch in self.children:
                if ch!=None: ch.destroy()
            self.location = []
            self.children = []
            self.parent = None
            self.templatenode = None
            self.valuenode = None
            self.store = None
            
        def isValid(self):
            return self.store != None

        def getValue(self):
            if self.valuenode==None: return None
            valuetype = self.templatenode.getAttribute('type')
            assert valuetype!='', 'getValue was used on node without type (%s); canHaveValue should have showed that this node does not have a type.' % self
            return self.store.getNodeProperty(self.valuenode,valuetype=valuetype)

        def getDefaultValue(self):
            defaultstore = self.controller.defaultstore
            if defaultstore==None: return None
            defaultnode = defaultstore.mapForeignNode(self)
            if defaultnode==None: return None
            return defaultnode.getValue()

        def getValueOrDefault(self):
            val = self.getValue()
            if val==None: val = self.getDefaultValue()
            return val

        def setValue(self,value):
            if value==None:
                self.clearValue()
                return

            curval = self.getValue()
            if curval!=value:
                if self.controller.onBeforeChange(self,value):
                    valuetype = self.templatenode.getAttribute('type')
                    if self.valuenode==None: self.createValueNode()
                    changed = self.store.setNodeProperty(self.valuenode,value,valuetype)
                    self.controller.onChange(self)
                    return changed
            return False

        def clearValue(self,recursive=False,skipreadonly=False,deleteclones=True):
            # First clear children.
            if recursive:
                if deleteclones: self.removeAllChildren()
                for ch in self.children:
                    ch.clearValue(recursive=True,skipreadonly=skipreadonly,deleteclones=deleteclones)
                
            # Do not clear if (1) it is already cleared, (2) it is read-only and the user wants to respect that,
            # or (3) it is the root node.
            if self.valuenode==None or (skipreadonly and self.isReadOnly()) or self.parent==None: return
            
            # Clear if (1) this node can have no value - it must occur, and (2) the attached interfaces approve.
            if (not self.canHaveClones()) and self.controller.onBeforeChange(self,None):
                self.store.clearNodeProperty(self.valuenode)
                self.valuenode = None
                self.controller.onChange(self)

        def preparePersist(self,recursive=True):
            valuetype = self.templatenode.getAttribute('type')
            if valuetype!='' and self.valuenode!=None:
                valueclass = self.store.filetypes[valuetype]
                if issubclass(valueclass,Store.DataType):
                    self.store.preparePersistNode(self.valuenode,valueclass)
            if recursive:
                for ch in self.children: ch.preparePersist(recursive=True)

        def persist(self,recursive=True):
            valuetype = self.templatenode.getAttribute('type')
            if valuetype!='' and self.valuenode!=None:
                valueclass = self.store.filetypes[valuetype]
                if issubclass(valueclass,Store.DataType):
                    self.store.persistNode(self.valuenode,valueclass)
            if recursive:
                for ch in self.children: ch.persist(recursive=True)

        def getValueAsString(self,addunit = True,usedefault = False):
            templatenode = self.templatenode
            fieldtype = templatenode.getAttribute('type')
            if usedefault:
                value = self.getValueOrDefault()
            else:
                value = self.getValue()
            if value==None: return ''
            if fieldtype=='datetime':
                value = value.strftime(common.datetime_displayformat)
            if fieldtype=='bool':
                if value:
                    value = 'Yes'
                else:
                    value = 'No'
            elif fieldtype=='file':
                # Return filename only (not the path)
                if not value.isValid():
                    value = ''
                else:
                    value = value.name
            elif fieldtype=='select':
                # Get label of currently selected option
                optionsroot = common.findDescendantNode(templatenode,['options'])
                if optionsroot==None: raise Exception('Variable with "select" type lacks "options" element below.')
                for ch in optionsroot.childNodes:
                    if ch.nodeType==ch.ELEMENT_NODE and ch.localName=='option':
                        if value==int(ch.getAttribute('value')):
                            # We found the currently selected option; its label will serve as displayed value.
                            value = ch.getAttribute('label')
                            break
            else:
                value = unicode(value)

            # Append unit specifier (if available)
            if addunit and templatenode.hasAttribute('unit'):
                value = value + ' ' + templatenode.getAttribute('unit')

            return value

        def addChild(self,childname,position=None):
            index = -1
            templatenode = None

            # First see of already one instance of this child is in the tree; that makes finding the position easy.
            existingcount = 0
            for curindex,child in enumerate(self.children):
                if child.location[-1]==childname:
                    index = curindex
                    templatenode = child.templatenode
                    existingcount += 1
                elif index!=-1:
                    break
                    
            # If no insert position was specified, append at the end
            if position==None: position = existingcount
            
            if index!=-1:
                # Found an existing node with this name
                assert position>=0, 'Position must be positive, but is %i. Use position=None to append to the end.' % position
                assert position<=existingcount, 'Cannot insert child "%s" at position %i, because only %i nodes exist so far.' % (childname,position,existingcount)
                index = index+1-existingcount+position
            else:
                # Node with this name not yet in tree.
                assert position==0, 'Cannot insert child "%s" at position %i, because no node wih this name exists so far.' % (childname,position)

                # Enumerate over all template children of the parent we want to insert below.
                # Store a list of names of children that precede the node to be inserted.
                predecessors = []
                for templatenode in self.templatenode.childNodes:
                    if templatenode.nodeType==templatenode.ELEMENT_NODE and templatenode.localName=='element':
                        childid = templatenode.getAttribute('id')
                        if childid==childname: break
                        predecessors.append(childid)
                else:
                    # Could not find the specified child in the template.
                    return None

                # Enumerate over all actual children until we reach the point where the child should be inserted.
                index = 0
                for child in self.children:
                    curname = child.location[-1]
                    while len(predecessors)>0 and curname!=predecessors[0]:
                        predecessors.pop(0)
                    if len(predecessors)==0: break
                    index += 1
                    
            # Ensure the parent to insert below has a value node
            # (we need to insert the value node below it to give the child life)
            self.createValueNode()
            
            # Find the XML document
            doc = self.valuenode
            while doc.parentNode!=None: doc=doc.parentNode
            
            # Create the value node for the current child
            node = doc.createElementNS(self.valuenode.namespaceURI,childname)
            
            # Insert the value node
            if position>=existingcount:
                valuenode = self.valuenode.appendChild(node)
            else:
                valuenode = self.valuenode.insertBefore(node,self.children[index].valuenode)
                
            # Create the child (template + value)
            child = TypedStore.Node(self.controller,templatenode,valuenode,self.location+[childname],parent=self)
            assert child.canHaveClones(), 'Cannot add another child "%s" because there can exist only one child with this name.' % childname
            
            # Insert the child, and notify attached interfaces.
            child.futureindex = index
            self.controller.beforeVisibilityChange(child,True,False)
            self.children.insert(index,child)
            self.controller.afterVisibilityChange(child,True,False)
            child.futureindex = None
            
            # Return the newly inserted child.
            return child
            
        def createValueNode(self):
            if self.valuenode!=None: return
            parents = []
            root = self
            while root.valuenode==None:
                parents.insert(0,root)
                root = root.parent
            doc = root.valuenode
            while doc.parentNode!=None: doc=doc.parentNode
            valueroot = root.valuenode
            for par in parents:
                valueroot = common.findDescendantNode(valueroot,[par.getId()],create=True)
            self.valuenode = valueroot

        def getNumberedChild(self,childname,index,create=False):
            children = self.getLocationMultiple([childname])
            if index<len(children): return children[index]
            if not create: return None
            for ichild in range(index-len(children)+1):
                child = self.addChild(childname)
            return child

        def removeChildren(self,childname,first=0,last=None):
            ipos = 0
            ichildpos = -1
            while ipos<len(self.children):
                child = self.children[ipos]
                if child.location[-1]==childname:
                    assert child.canHaveClones(),'Cannot remove child "%s" because it must occur exactly one time.' % childname
                    ichildpos += 1
                    if last!=None and ichildpos>last: return
                    if ichildpos>=first:
                        self.controller.beforeVisibilityChange(child,False,False)
                        child = self.children.pop(ipos)
                        self.store.clearNodeProperty(child.valuenode)
                        self.controller.afterVisibilityChange(child,False,False)
                        child.destroy()
                        ipos -= 1
                ipos += 1
                
        def removeAllChildren(self,optionalonly=True):
            for ipos in range(len(self.children)-1,-1,-1):
                child = self.children[ipos]
                if (not optionalonly) or child.canHaveClones():
                    child.removeAllChildren(optionalonly=False)
                    self.controller.beforeVisibilityChange(child,False,False)
                    child = self.children.pop(ipos)
                    if child.valuenode!=None: self.store.clearNodeProperty(child.valuenode)
                    self.controller.afterVisibilityChange(child,False,False)
                    child.destroy()

        def getId(self):
            return self.location[-1]

        def getValueType(self):
            return self.templatenode.getAttribute('type')

        def getText(self,detail,minimumdetail = 0,capitalize=False):
            templatenode = self.templatenode
            ret = None
            if detail==2 and templatenode.hasAttribute('description'):
                ret = templatenode.getAttribute('description')
            elif detail>=1 and minimumdetail<=1 and templatenode.hasAttribute('label'):
                ret = templatenode.getAttribute('label')
            elif minimumdetail==0:
                ret = self.getId()
            if ret!=None and capitalize: ret = ret[0].upper() + ret[1:]
            return ret

        def getLocation(self,location):
            node = self
            for childname in location:
                if childname=='..':
                    assert self.parent!=None,'Cannot go up one level because we are at the root.'
                    node = node.parent
                elif childname!='' and childname!='.':
                    for node in node.children:
                        if node.location[-1]==childname: break
                    else:
                        return None
            return node

        def getLocationMultiple(self,location):
            # Get the first non-empty path term.
            path = location[:]
            target = ''
            while target=='' and len(path)>0: target = path.pop(0)
            if target=='': return [self]

            res = []
            for child in self.children:
                if child.location[-1]==target:
                    if len(path)==0:
                        res.append(child)
                    else:
                        res += child.getLocationMultiple(path)
            return res

        def isHidden(self):
            node = self
            while node!=None:
                if not node.visible: return True
                node = node.parent
            return False

        def isReadOnly(self):
            return self.templatenode.hasAttribute('readonly')

        def hasChildren(self):
            return len(self.children)>0
    
        def canHaveValue(self):
            return self.templatenode.hasAttribute('type')

        def canHaveClones(self):
            return self.templatenode.hasAttribute('maxoccurs')

        def getNodesByType(self,valuetype):
            res = []
            if self.getValueType()==valuetype:
                res.append(self)
            for ch in self.children:
                res += ch.getNodesByType(valuetype)
            return res

        def getEmptyNodes(self):
            res = []
            if self.canHaveValue() and self.getValue()==None:
                res.append(self)
            for ch in self.children:
                res += ch.getEmptyNodes()
            return res

        def updateVisibility(self,recursive=False,notify=True):
            templatenode = self.templatenode
            cond = common.findDescendantNode(templatenode,['condition'])
            if cond!=None:
                shownew = self.controller.checkCondition(cond,templatenode)
                if shownew!=self.visible:
                    if notify: self.controller.beforeVisibilityChange(self,shownew)
                    self.visible = shownew
                    if notify: self.controller.afterVisibilityChange(self,shownew)
            if recursive:
                for child in self.children: child.updateVisibility(recursive=True,notify=notify)

        def copyFrom(self,sourcenode,replace=True):
            # Copy node value (if both source and target can have a value)
            if self.canHaveValue() and sourcenode.canHaveValue():
                if replace or self.getValue()==None:
                    curval = sourcenode.getValue()
                    if isinstance(curval,Store.DataType): curval.addref()
                    self.setValue(curval)

            # If replacing previous contents, remove optional nodes (with minoccurs=0)
            if replace: self.removeAllChildren()
            prevchildname = None
            index = 0
            for sourcechild in sourcenode.children:
                childname = sourcechild.location[-1]
                if childname!=prevchildname:
                    index = 0
                    prevchildname = childname
                if sourcechild.canHaveClones():
                    child = self.getNumberedChild(childname,index,create=True)
                else:
                    child = self.getLocation([childname])
                if child==None: continue
                child.copyFrom(sourcechild,replace=replace)
                index += 1
            
    schemas = None
    defaults = None
    defaultname2scenarios = None
    
    @staticmethod
    def getDefaultSchemas():
        return {}

    @staticmethod
    def getDefaultValues():
        return {}

    @classmethod
    def schemaname2path(cls,schemaname=None):
        if cls.schemas == None:
            cls.schemas = cls.getDefaultSchemas()
        if schemaname==None: return cls.schemas
        if schemaname in cls.schemas:
            return cls.schemas[schemaname]
        else:
            return None

    @classmethod
    def defaultname2path(cls,defaultname=None):
        if cls.defaults == None:
            cls.defaults = cls.getDefaultValues()
        if defaultname==None: return cls.defaults
        if defaultname in cls.defaults:
            return cls.defaults[defaultname]
        else:
            return None

    @classmethod
    def getDefault(cls,name,version):
        if cls==TypedStore: return None
        if name==None: name = 'default'
        if cls.defaultname2scenarios==None: cls.defaultname2scenarios = {}
        if name not in cls.defaultname2scenarios: cls.defaultname2scenarios[name] = {}
        version2store = cls.defaultname2scenarios[name]
        if version in version2store:
            # We have the requested default with the requested version in our cache; return it.
            return version2store[version]
        elif 'source' not in version2store:
            # We do not have any version of the requested default; first obtain the source version.
            path = cls.defaultname2path(name)
            if path==None: return None
            sourcestore = cls.fromXmlFile(path,adddefault=False)
            version2store['source'] = sourcestore
            version2store[sourcestore.version] = sourcestore
            if sourcestore.version==version: return sourcestore
        # We now have the source version of the requested default, but we need another version. Convert.
        sourcestore = version2store['source']
        defstore = cls.fromSchemaName(version,adddefault=False)
        sourcestore.convert(defstore)
        version2store[version] = defstore
        return defstore

    @classmethod
    def fromSchemaName(cls,schemaname,valueroot=None,adddefault=True):
        assert cls!=TypedStore, 'fromSchemaName cannot be called on base class "TypedStore", only on derived classes. You probably need to create a derived class with versioning support.'
        path = cls.schemaname2path(schemaname)
        if path==None:
            raise Exception('Unable to locate XML schema file for "%s".' % schemaname)
        store = cls(path, valueroot, adddefault = adddefault)
        cls.schemas[schemaname] = store.schemadom
        return store

    @classmethod
    def fromXmlFile(cls,path,adddefault=True):
        assert cls!=TypedStore, 'fromXmlFile cannot be called on base class "TypedStore", only on derived classes. Use setStore if you do not require versioning.'
        if not os.path.isfile(path):
            raise Exception('Specified path "%s" does not exist, or is not a file.' % path)
        valuedom = xml.dom.minidom.parse(path)
        version = valuedom.documentElement.getAttribute('version')
        return cls.fromSchemaName(version,valuedom,adddefault=adddefault)

    def __init__(self,schemadom,valueroot=None,otherstores={},adddefault=True):

        # The template can be specified as a DOM object, or as string (i.e. path to XML file)
        schemapath = '.'
        if isinstance(schemadom,basestring):
            schemapath = schemadom
            if not os.path.isfile(schemapath):
                raise Exception('XML schema file "%s" does not exist.' % schemapath)
            schemadom = xml.dom.minidom.parse(schemapath)
        self.schemadom = schemadom
        
        # Resolve links to external documents
        links = self.schemadom.getElementsByTagName('link')
        for link in links:
            assert link.hasAttribute('path'), 'Link node does not have "path" attribute.'
            refpath = os.path.abspath(os.path.join(os.path.dirname(schemapath),link.getAttribute('path')))
            if not os.path.isfile(refpath):
                raise Exception('Linked XML schema file "%s" does not exist.' % refpath)
            childdom = xml.dom.minidom.parse(refpath)
            linkparent = link.parentNode
            newnode = common.copyNode(childdom.documentElement,linkparent,targetdoc=schemadom)
            for key in link.attributes.keys():
                if key!='path': newnode.setAttribute(key,link.getAttribute(key))
            linkparent.removeChild(link)
        
        # Get schema version
        self.version = self.schemadom.documentElement.getAttribute('version')
        self.originalversion = None

        # Events
        self.interfaces = []

        self.otherstores = otherstores

        # For every variable: build a list of variables/folders that depend on its value.
        if not self.schemadom.documentElement.hasAttribute('dependenciesbuilt'):
            self.buildDependencies()
            self.schemadom.documentElement.setAttribute('dependenciesbuilt','True')

        # Link to original source (if any)
        self.path = None

        # Clear store variables
        self.store = None
        self.defaultstore = None
        self.defaultinterface = None
        self.root = None
        
        # Add store with default values if requested and available.
        if adddefault:
            defscenario = self.getDefault(None,self.version)
            if defscenario!=None: self.setDefaultStore(defscenario,updatevisibility=False)

        # Now set current values in the store
        # NB: this must be done after default values are set, so that the default
        # values can be taken into account when checking conditions (via setStore)
        self.setStore(valueroot)

    def unlink(self):
        if self.root!=None: self.root.destroy()
        self.root = None

        self.setContainer(None)
        self.store = None
        self.interfaces = []

    def getInterface(self,showhidden=True,omitgroupers=False):
        return TypedStoreInterface(self,showhidden=showhidden,omitgroupers=omitgroupers)

    def setContainer(self,container):
        if 'cache' in self.store.context:
            for v in self.store.context['cache'].itervalues(): v.release()
            del self.store.context['cache']
        if 'container' in self.store.context and self.store.context['container']!=None:
            self.store.context['container'].release()
        self.store.context['container'] = container

    def setStore(self,valueroot):
        if self.root!=None: self.root.destroy()

        templateroot = self.schemadom.documentElement

        assert valueroot==None or isinstance(valueroot,basestring) or isinstance(valueroot,xml.dom.Node), 'Supplied value root must None, a path to an XML file, or an XML node, but is %s.' % valueroot

        valuedom = None
        if valueroot==None:
            impl = xml.dom.minidom.getDOMImplementation()
            assert templateroot.hasAttribute('id'), 'Root of the schema does not have attribute "id".'
            valuedom = impl.createDocument(None, templateroot.getAttribute('id'), None)
            valueroot = valuedom.documentElement
            valueroot.setAttribute('version',self.version)
        elif isinstance(valueroot,basestring):
            valuedom = xml.dom.minidom.parse(valueroot)
            valueroot = valuedom.documentElement
        elif valueroot.nodeType==valueroot.DOCUMENT_NODE:
            valuedom = valueroot
            valueroot = valuedom.documentElement
        else:
            valuedom = valueroot
            while valuedom.parentNode!=None: valuedom = valuedom.parentNode

        storeversion = valueroot.getAttribute('version')
        assert storeversion==self.version, 'Versions of the xml schema ("%s") and and the xml values ("%s") do not match.' % (self.version,storeversion)
                    
        self.store = Store(valuedom,xmlroot=valueroot)
        self.store.filetypes['select'] = int
        self.store.filetypes['file'] = DataFile
        self.store.filetypes['color'] = StoreColor
        self.root = TypedStore.Node(self,templateroot,self.store.xmlroot,[],None)
        self.changed = False
        self.setContainer(None)
        
        # Update the visibility of all nodes - based on conditions
        # Disable individual notifications because the single "storechanged" event emitted
        # below replaces them)
        self.root.updateVisibility(recursive=True,notify=False)
        
        # Notify attached interface about the store change.
        self.afterStoreChange()

    def setDefaultStore(self,store,updatevisibility=True):
        assert self.version==store.version
        self.defaultstore = store
        self.defaultinterface = self.defaultstore.getInterface()
        self.defaultinterface.addChangeHandler(self.onDefaultChange)
        
        # Default nodes are used in condition checking, so changing the default store
        # requires updating the visibility of all nodes. Do so, unless explicitly said not to.
        if updatevisibility: self.root.updateVisibility(recursive=True)

    def hasChanged(self):
        return self.changed

    def resetChanged(self):
        self.changed = False

    def setProperty(self,location,value):
        node = self.root.getLocation(location)
        if node==None: raise Exception('Cannot locate node at %s' % location)
        return node.setValue(value)
    
    def getProperty(self,location,usedefault=False):
        node = self.root.getLocation(location)
        if node==None: raise Exception('Cannot locate node at %s' % location)
        if usedefault:
            return node.getValueOrDefault()
        else:
            return node.getValue()

    def mapForeignNode(self,foreignnode):
        indices = []
        currentnode = foreignnode
        # First we walk up the tree from the supplied foreign node, in order to find the indices
        # of all involved ancestors.
        for name in reversed(foreignnode.location):
            if not currentnode.canHaveClones():
                # This node must appear once; its index can only be one.
                indices.insert(0,0)
            else:
                # This node can appear zero or more times. Find its index in the foreign store.
                siblings = currentnode.parent.getLocationMultiple([name])
                for (index,sib) in enumerate(siblings):
                    if sib is currentnode: break
                else:
                    assert False, 'Cannot find foreign node "%s" in list of its own siblings.' % name
                indices.insert(0,index)
            currentnode = currentnode.parent
        assert currentnode.parent==None, 'Location does not describe complete path to root. Currently at %s.' % currentnode
        
        # Now find the same location in our own store.
        currentnode = self.root
        for (name,index) in zip(foreignnode.location,indices):
            currentnode = currentnode.getNumberedChild(name,index)
            if currentnode==None: break
        return currentnode

    # buildDependencies: for every variable node, this creates lists of dependent nodes
    # (i.e. folders and variables that have one or more conditions that depend on the
    # variable under investigation). Essentially we convert lists of dependencies ('servant'-centric)
    # into lists of dependent nodes ('controller'-centric). We need the latter in order to selectively
    # re-check conditions (and hide/show corresponding nodes) after the value of
    # a dependency ('controller') changes.
    def buildDependencies(self,root=None,curpath='',curowner=None):
        if root==None: root=self.schemadom.documentElement
        for ch in root.childNodes:
            if ch.nodeType==ch.ELEMENT_NODE:
                if ch.localName=='element':
                    self.buildDependencies(root=ch,curpath=curpath+'/'+ch.getAttribute('id'),curowner=ch)
                elif ch.localName=='condition':
                    if ch.hasAttribute('source'): continue
                    if ch.hasAttribute('variable'):
                        deppath = ch.getAttribute('variable')
                        refnode = None
                        if deppath[0]!='/': refnode = curowner.parentNode
                        dep = self.getTemplateNode(deppath.split('/'),root=refnode)
                        assert dep!=None, 'Cannot locate dependency "%s" for node "%s".' % (ch.getAttribute('variable'),curpath)
                        deplist = common.findDescendantNode(dep,['dependentvariables'],create=True)
                        node = self.schemadom.createElementNS(deplist.namespaceURI,'dependentvariable')
                        node.setAttribute('path',curpath)
                        deplist.appendChild(node)
                    self.buildDependencies(root=ch,curpath=curpath,curowner=curowner)

    # checkCondition: checks whether then given condition (an XML node in the template) is currently met.
    #   "nodeCondition" is the "condition" XML node to check
    #   "ownernode" is the "variable" or "folder" XML node that 'owns' the condition
    #       (= the first ancestor that is not a condition itself)
    def checkCondition(self,nodeCondition,ownernode,ownstorename=None):
        assert nodeCondition.hasAttribute('type'), 'condition lacks "type" attribute in XML schema file.'
        src = nodeCondition.getAttribute('source')
        if src!='' and src!=ownstorename:
            if src not in self.otherstores: return True
            return self.otherstores[src].checkCondition(nodeCondition,ownernode,ownstorename=src)
        condtype = nodeCondition.getAttribute('type')
        if condtype=='eq' or condtype=='ne':
            # Check for required XML attributes
            assert nodeCondition.hasAttribute('variable'), 'condition lacks "variable" attribute in XML schema file.'
            assert nodeCondition.hasAttribute('value'), 'condition lacks "value" attribute in XML schema file.'

            valuepath = nodeCondition.getAttribute('variable')
            refnode = self.root
            if valuepath[0]!='/': refnode = self.root.getLocation(self.getTemplateNodePath(ownernode)).parent
            node = refnode.getLocation(valuepath.split('/'))
            assert node!=None, 'Cannot locate dependency "%s" for node "%s".' % (nodeCondition.getAttribute('variable'),self.getTemplateNodePath(ownernode))

            # Get the type and current value of the variable we depend on
            valuetype = node.getValueType()
            curvalue = node.getValueOrDefault()

            # If the node in question currently does not have a value, we cannot check the condition; just return 'valid'.
            if curvalue==None: return True

            # Get the reference value we will compare against
            refvalue = self.store.unpackvalue(nodeCondition.getAttribute('value'),valuetype)

            # Compare
            if condtype=='eq': return (curvalue==refvalue)
            if condtype=='ne': return (curvalue!=refvalue)
            
        elif condtype=='and' or condtype=='or':
            # Check every child condition.
            for ch in nodeCondition.childNodes:
                if ch.nodeType==ch.ELEMENT_NODE and ch.localName=='condition':
                    if self.checkCondition(ch,ownernode):
                        # OR query: one True guarantees success 
                        if condtype=='or': return True
                    else:
                        # AND query: one False guarantees failure 
                        if condtype=='and': return False
                        
            # We evaluated all children. If we are doing an OR, that means all
            # children returned False: we failed, if we are doing an AND, all
            # children returned True: we succeeded.
            if condtype=='and': return True
            return False
        else:
            raise Exception('unknown condition type "%s" in XML schema file.' % condtype)

    # getTemplateNode: obtains template node at given path
    # (path specification consists of array of node ids)
    def getTemplateNode(self,path,root=None):
        if root==None: root=self.schemadom.documentElement
        for childname in path:
            if childname=='..':
                assert not root.isSameNode(self.schemadom.documentElement)
                root = root.parentNode
            elif childname!='' and childname!='.':
                children = [ch for ch in root.childNodes if ch.nodeType==ch.ELEMENT_NODE and ch.localName=='element']
                for root in children:
                    if root.getAttribute('id')==childname:
                        break
                else:
                    assert False, 'Could not find child node "%s" while locating "%s". Children: %s' % (childname,path,','.join([ch.getAttribute('id') for ch in children]))
        return root

    # getNodePath: obtains path specification for given template node
    # (path specification consists of node ids with slash separators)
    def getTemplateNodePath(self,node):
        path = []
        while node.parentNode.parentNode!=None:
            path.insert(0,node.getAttribute('id'))
            node = node.parentNode
        return path

    def convert(self,target):        
        if isinstance(target,basestring):
            target = self.fromSchemaName(target)
        
        convertor = self.getConvertor(self.version,target.version)
        if convertor==None:
            raise Exception('No convertor available to convert version "%s" to "%s".' % (self.version,target.version))
        convertor.convert(self,target)

        return target

    convertorsfrom = {}

    # getConvertor(sourceid,targetid,directonly=False)
    #   Returns a convertor object, capable of converting between the specified versions.
    #   Use directonly=True to retrieve only direct conversion routes.
    #   Return None if no convertor is available that meets the specified criteria.
    @classmethod
    def getConvertor(cls,sourceid,targetid,directonly=False):
        # Try direct route first.
        if (sourceid in cls.convertorsfrom) and (targetid in cls.convertorsfrom[sourceid]):
            convertorclass = cls.convertorsfrom[sourceid][targetid]
            return convertorclass()

        # Direct route not available, try indirect routes
        if not directonly:
            indirectroutes = cls.findIndirectConversion(sourceid,targetid,depth='  ')
            if len(indirectroutes)>0:
                indirectroutes.sort(key=len)
                route = indirectroutes[0]
                chain = []
                for istep in range(len(route)-1):
                    convertor = cls.getConvertor(route[istep],route[istep+1],directonly=True)
                    chain.append(convertor)
                return ConvertorChain(chain)

        # No route available.
        return None

    # findIndirectConversion(sourceid,targetid)
    #   [internal use only!]
    #   Return conversion routes between sourceid and targetid, avoiding versions in disallowed.
    #   The depth argument is used for debug output only.
    @classmethod
    def findIndirectConversion(cls,sourceid,targetid,disallowed=[],depth=''):
        next = cls.convertorsfrom.get(sourceid,{}).keys()
        routes = []
        curdisallowed = disallowed[:]+[sourceid]
        for curnext in next:
            if curnext in curdisallowed: continue
            if curnext==targetid:
                routes.append([sourceid,curnext])
            else:
                childroutes = cls.findIndirectConversion(curnext,targetid,curdisallowed,depth=depth+'  ')
                for cr in childroutes:
                    routes.append([sourceid]+cr)
        return routes

    # addConvertor(convertorclass)
    #   [internal use only!]
    #   Register a convertor class.
    @classmethod
    def addConvertor(cls,convertorclass):
        sourceid = convertorclass.fixedsourceid
        targetid = convertorclass.fixedtargetid
        assert sourceid!=None, 'Error! Specified convertor class lacks a source identifier.'
        assert targetid!=None, 'Error! Specified convertor class lacks a target identifier.'
        if sourceid not in cls.convertorsfrom: cls.convertorsfrom[sourceid] = {}
        assert targetid not in cls.convertorsfrom[sourceid], 'Error! A class for converting from "%s" to "%s" was already specified previously.' % (sourceid,targetid)
        cls.convertorsfrom[sourceid][targetid] = convertorclass

    # hasConvertor(sourceid,targetid)
    #   Checks if a conversion route between the specified versions is available.
    #   Both direct and indirect (via another version) routes are ok.
    @classmethod
    def hasConvertor(cls,sourceid,targetid):
        # Try direct conversion
        if cls.getConvertor(sourceid,targetid)!=None:
            return True

        print 'Searching for indirect conversion routes from '+sourceid+' to '+targetid+'.'
        indirectroutes = cls.findIndirectConversion(sourceid,targetid,depth='  ')
        for indirect in indirectroutes:
            print indirect
        print 'Found '+str(len(indirectroutes))+' indirect routes.'

        return len(indirectroutes)>0

    # rankSources(targetid,sourceids,requireplatform=None)
    #   Rank a set of supplied versions/identifiers according to platform (i.e. gotmgui, gotm) and version
    #   Rank criterion is 'closeness' (in version and platform) to the reference targetid.
    @classmethod
    def rankSources(cls,targetid,sourceids,requireplatform=None):
        (targetplatform,targetversion) = targetid.split('-')
        targetversion = versionStringToInt(targetversion)

        # Decompose source ids into name and (integer) version, but only take source we can actually convert to the target version.
        sourceinfo = []
        for sid in sourceids:
            if sid==targetid or cls.hasConvertor(sid,targetid):
                (platform,version) = sid.split('-')
                if requireplatform==None or requireplatform==platform:
                    version = versionStringToInt(version)
                    sourceinfo.append((platform,version,sid))

        # Sort by platform (because we want the target platform first)
        sourceinfoclasses = {}
        for sinf in sourceinfo:
            if sinf[0] not in sourceinfoclasses: sourceinfoclasses[sinf[0]] = []
            sourceinfoclasses[sinf[0]].append(sinf)

        # Now sort per platform according to version (higher versions first)
        result = []
        for sourceplatform in sourceinfoclasses.keys():
            infos = sourceinfoclasses[sourceplatform]
            infos.sort(cmp=lambda x,y: cmp(y[1],x[1]))
            if sourceplatform==targetplatform:
                result = infos+result
            else:
                result += infos

        resultids = []
        for res in result: resultids.append(res[2])

        return resultids

    def load(self,path):
        if not os.path.isfile(path):
            raise Exception('Specified path "%s" does not exist, or is not a file.' % path)
        valuedom = xml.dom.minidom.parse(path)

        version = valuedom.documentElement.getAttribute('version')
        if self.version!=version:
            # The version of the saved scenario does not match the version of this scenario object; convert it.
            print 'Value file "%s" has version "%s"; starting conversion to "%s".' % (path,version,self.version)
            tempstore = self.fromSchemaName(version)
            tempstore.setStore(valuedom)
            tempstore.convert(self)
            tempstore.unlink()
            self.originalversion = version
        else:
            self.setStore(valuedom)

    def loadAll(self,path):
        if isinstance(path,basestring):
            # Basic check: does the specified source path exist?
            if not os.path.exists(path):
                raise Exception('Source path "%s" does not exist.' % path)

            # The specified source file may be a ZIP archive, or a directory that contains the extracted
            # contents of such an archive; decide which.
            if os.path.isdir(path):
                container = DataContainerDirectory(path)
            else:
                container = DataContainerZip(path)
        elif isinstance(path,DataFile):
            container = DataContainerZip(path)
        else:
            assert False,'Supplied source must be a path to file or directory, or a data file object.'

        # Get list of files in source container.
        files = container.listFiles()

        # Check for existence of main file.
        if self.storefilename not in files:
            raise Exception('The specified source does not contain "%s" and can therefore not be a %s.' % (self.storefilename,self.storetitle))

        # Read and parse the store XML file.
        datafile = container.getItem(self.storefilename)
        f = datafile.getAsReadOnlyFile()
        storedom = xml.dom.minidom.parse(f)
        f.close()
        datafile.release()
        
        # Get the schema version that the store values file matches.
        version = storedom.documentElement.getAttribute('version')
        if self.version!=version:
            # The version of the saved scenario does not match the version of this scenario object; convert it.
            print '%s "%s" has version "%s"; starting conversion to "%s".' % (self.storetitle,path,version,self.version)
            container.release()
            tempstore = self.fromSchemaName(version)
            tempstore.loadAll(path)
            tempstore.convert(self)
            tempstore.unlink()
            self.originalversion = version
        else:
            # Attach the parsed scenario (XML DOM).
            self.setStore(storedom)
            self.setContainer(container)

        # Store source path.
        if isinstance(path,basestring):
            self.path = path
        else:
            self.path = None

    def save(self,path):
        return self.store.save(path)

    def saveAll(self,path,targetversion=None,targetisdir = False,claim=True):
        if targetversion!=None and self.version!=targetversion:
            # First convert to the target version
            tempstore = self.convert(targetversion)

            # Now save the result of the conversion.
            tempstore.saveAll(path, targetversion = targetversion,targetisdir = targetisdir)

            # Convert back: by doing this the original store will be able to reference nodes
            # (of type "file") in the newly saved file.
            tempstore.convert(self)

            # Release the conversion result.
            tempstore.unlink()
        else:
            # First: fill nodes that are not set with the default value.
            self.root.copyFrom(self.defaultstore.root,replace=False)

            # Before opening the target container, allow nodes to prepare for saving to the specified path.
            # Specifically, nodes will read all files that might be overwritten into memory.
            if isinstance(path,basestring):
                self.store.context['targetcontainerpath'] = path
                self.root.preparePersist()
                del self.store.context['targetcontainerpath']

            # Open target container
            if isinstance(path,basestring):
                if targetisdir:
                    container = DataContainerDirectory(path,create=True)
                else:
                    container = DataContainerZip(path,mode='w')
            elif isinstance(path,StringIO.StringIO):
                container = DataContainerZip(path,mode='w')
                claim = False
            else:
                assert False,'Supplied target must be a path to file or directory, or a StringIO object.'

            # Allow all nodes to add custom data to the target container. This can change the values
            # in the XML store, and must therefore be done before the store is added to the container.
            self.store.context['targetcontainer'] = container
            self.store.context['donotclaimtarget'] = (not claim)
            self.root.persist()
            del self.store.context['donotclaimtarget']

            # Add XML store to the container
            df = DataFileXmlNode(self.store.xmldocument)
            df_added = container.addItem(df,self.storefilename)
            df_added.release()
            df.release()

            # Make the container save all changes and then release it.
            # Note if claim=True: even though we release it, many nodes (of type "file") may now hold
            # references to data in the saved container; the container will likely not be completely
            # released. On the other hand, the original sources that were locked before saving now
            # probably will be released (the nodes do not lock them any longer).
            container.persistChanges()
            container.release()

        if isinstance(path,basestring):
            self.path = path
        else:
            self.path = None
        
        self.resetChanged()

    def toxml(self,enc='utf-8'):
        return self.store.xmldocument.toxml(enc)

    def toXmlDom(self,target=None):
        return common.copyNode(self.store.xmlroot,target)

    # ----------------------------------------------------------------------------------------
    # Event handling
    # ----------------------------------------------------------------------------------------

    # connectInterface: connects an interface to the store. Interfaces provide events and
    #   can hide nodes with the hidden attribute from view, amongst others.
    def connectInterface(self,interface):
        self.interfaces.append(interface)

    # onChange: called internally after the default value of a node has changed.
    #   note that its argument will be a node in the DEFAULT store, not in the current store!
    def onDefaultChange(self,defaultnode):
        # Map node in default store to node in our own store.
        ownnode = self.mapForeignNode(defaultnode)
        if ownnode==None: return

        # Emit change event
        for i in self.interfaces: i.onDefaultChange(ownnode)

        # If the default is being used: update (visibility of) nodes that depend on the changed node.
        if ownnode.getValue()==None: self.updateDependantNodes(ownnode)

    # onChange: called internally after the value of a node has changed.
    def onChange(self,node):
        # Register that we changed.
        self.changed = True

        # Emit change event
        for i in self.interfaces: i.onChange(node)

        # Update (visibility of) nodes that depend on the changed node.
        self.updateDependantNodes(node)

    def updateDependantNodes(self,node):
        # Get nodes that depend on the changed node; if there are none, exit.
        deps = common.findDescendantNode(node.templatenode,['dependentvariables'])
        if deps==None: return

        # Now build a list of the dependant nodes; currently hidden nodes first, currently visible
        # nodes last, so that when we iterate over the list and switch visibilities first extra nodes
        # will appear, and only later some are removed (this prevents unnecessary automatic scrolling in GUI)
        depnodes = []
        for d in common.findDescendantNodes(deps,['dependentvariable']):
            varpath = d.getAttribute('path').split('/')
            varnode = self.root.getLocation(varpath)
            if varnode.visible:
                depnodes.append(varnode)
            else:
                depnodes.insert(0,varnode)
        for varnode in depnodes: varnode.updateVisibility()

    # onBeforeChange: called internally just before the value of a node changes.
    #   The return value decides if the change is allowed (return True) or denied (return False)
    def onBeforeChange(self,node,newvalue):
        for i in self.interfaces:
            if not i.onBeforeChange(node,newvalue): return False
        return True

    # onBeforeChange: called internally after the store (i.e., all values) have changed.
    def afterStoreChange(self):
        for i in self.interfaces: i.afterStoreChange()

    # onBeforeChange: called internally before a node is hidden (showhide=True) or deleted (showhide=False).
    def beforeVisibilityChange(self,node,visible,showhide=True):
        for i in self.interfaces: i.beforeVisibilityChange(node,visible,showhide)

    # onBeforeChange: called internally after a node is hidden (showhide=True) or deleted (showhide=False).
    def afterVisibilityChange(self,node,visible,showhide=True):
        for i in self.interfaces: i.afterVisibilityChange(node,visible,showhide)


# versionStringToInt(versionstring)
#   Converts a "major.minor.build" version string to a representative integer.
def versionStringToInt(versionstring):
    (major,minor,build) = versionstring.split('.')
    return int(major)*256*256 + int(minor)*256 + int(build)

# Convertor
#   Base class for conversion; derive custom convertors from this class.
class Convertor:
    fixedsourceid = None
    fixedtargetid = None

    def __init__(self,sourceid=None,targetid=None):
        if sourceid==None:
            if self.fixedsourceid==None:
                raise Exception('Convertor class was created without explicit version identifiers, but also lacks default identifiers.')
            sourceid = self.fixedsourceid
            targetid = self.fixedtargetid
        else:
            if self.fixedsourceid!=None:
                raise Exception('Convertor class was created with explicit version identifiers, but also has default identifiers.')
        
        self.sourceid = sourceid
        self.targetid = targetid

        (self.sourcename,self.sourceversion) = sourceid.split('-')
        (self.targetname,self.targetversion) = targetid.split('-')

        self.sourceversionint = versionStringToInt(self.sourceversion)
        self.targetversionint = versionStringToInt(self.targetversion)

        self.links = []
        self.defaults = []
        self.registerLinks()

    def registerLinks(self):
        pass

    def convert(self,source,target):
        # Try simple deep copy: nodes with the same name and location in both
        # source and target store will have their value copied.
        target.root.copyFrom(source.root)

        # Handle one-to-one links between source nodes and target nodes.
        for (sourcepath,targetpath) in self.links:
            if isinstance(sourcepath,basestring): sourcepath = sourcepath.split('/')
            if isinstance(targetpath,basestring): targetpath = targetpath.split('/')
            sourcenode = source.root.getLocation(sourcepath)
            if sourcenode==None:
                raise Exception('Cannot locate node "%s" in source.' % '/'.join(sourcepath))
            targetnode = target.root.getLocation(targetpath)
            if targetnode==None:
                raise Exception('Cannot locate node "%s" in target.' % '/'.join(targetpath))
            targetnode.copyFrom(sourcenode)

        # Reset target nodes to defaults where applicable.
        if len(self.defaults)>0:
            defscen = target.getDefault(None,target.version)
            for path in self.defaults:
                if isinstance(path,basestring): path = path.split('/')
                sourcenode = defscen.root.getLocation(path)
                if sourcenode==None:
                    raise Exception('Cannot locate node "%s" in default.')
                targetnode = target.root.getLocation(path)
                targetnode.copyFrom(sourcenode)

    def reverseLinks(self):
        newlinks = []
        for (sourcepath,targetpath) in self.links:
            newlinks.append((targetpath,sourcepath))
        return newlinks

# ConvertorChain
#   Generic class for multiple-step conversions
#   Conversion steps are specified at initialization as a list of convertors.
class ConvertorChain(Convertor):
    def __init__(self,chain):
        Convertor.__init__(self,chain[0].sourceid,chain[-1].targetid)
        self.chain = chain

    def convert(self,source,target):
        temptargets = []
        for istep in range(len(self.chain)-1):
            convertor = self.chain[istep]
            temptargetid = convertor.targetid
            print 'Converting to temporary target "%s".' % temptargetid
            temptarget = source.fromSchemaName(temptargetid)
            temptargets.append(temptarget)
            convertor.convert(source,temptarget)
            source = temptarget
        convertor = self.chain[-1]
        print 'Converting to final target "%s".' % target.version
        convertor.convert(source,target)
        for temptarget in temptargets: temptarget.unlink()

