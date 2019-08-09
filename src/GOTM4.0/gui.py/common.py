#$Id: common.py,v 1.25 2007-04-18 09:31:09 jorn Exp $

import datetime,time,sys,xml.dom.minidom
import matplotlib.numerix

# ------------------------------------------------------------------------------------------
# Date-time parsing variables and functions
# ------------------------------------------------------------------------------------------

# datetime_displayformat: date format used to display datetime objects in the GUI.
datetime_displayformat = '%Y-%m-%d %H:%M:%S'

# parsedatetime: Convert string to Python datetime object, using specified format.
#   Counterpart of datetime.strftime.
def parsedatetime(str,fmt):
    t1tmp = time.strptime(str,fmt) 
    return datetime.datetime(*t1tmp[0:6])
    
def timedelta2float(td):
    return td.days*3600*24 + td.seconds + td.microseconds/1e6

# ------------------------------------------------------------------------------------------
# Command line argument utility functions
# ------------------------------------------------------------------------------------------

# getNamedArgument: Get the value of a named command line argument, and removes both name
#   and value from the global list of command line arguments. Returns None if the command
#   line argument was not specified. If the script was called with 'script.py -d hello',
#   getNamedArgument('-d') will return 'hello'.
def getNamedArgument(name):
    try:
        iarg = sys.argv.index(name)
    except ValueError:
        return None
    val = sys.argv[iarg+1]
    del sys.argv[iarg+1]
    del sys.argv[iarg]
    return val

# ------------------------------------------------------------------------------------------
# XML helper functions
# ------------------------------------------------------------------------------------------

# findDescendantNode: Return the first child XML DOM node with the specified location
#   (location = array of path components) below the specified XML DOM node (root).
#   If create = True, the node will be created if it does not exist yet.
def findDescendantNode(root,location,create=False):
    assert root!=None,'findDescendantNode called on non-existent parent node (parent = None).'
    node = root
    for childname in location:
        if childname=='': continue
        foundchild = None
        for ch in node.childNodes:
            if ch.nodeType==ch.ELEMENT_NODE and ch.localName==childname:
                foundchild = ch
                break
        else:
            if create:
                doc = root
                while doc.parentNode!=None: doc=doc.parentNode
                foundchild = doc.createElementNS(node.namespaceURI,childname)
                node.appendChild(foundchild)
            else:
                return None
        node = foundchild
    return node

# findDescendantNodes: Return a list of all child XML DOM nodes with the specified location
#   (location = array of path components) below the specified XML DOM node (root).
def findDescendantNodes(root,location):
    parentloc = location[:]
    name = parentloc.pop()
    parent = findDescendantNode(root,parentloc,create=False)
    children = []
    if parent!=None:
        for ch in parent.childNodes:
            if ch.nodeType==ch.ELEMENT_NODE and ch.localName==name:
                children.append(ch)
    return children

def addDescendantNode(parent,location):
    doc = parent
    while doc.parentNode!=None: doc=doc.parentNode
    for name in location:
        node = doc.createElementNS(parent.namespaceURI,name)
        parent.appendChild(node)
        parent = node
    return parent

def getNodeText(node):
    rc = ''
    for ch in node.childNodes:
        if ch.nodeType == ch.TEXT_NODE: rc += ch.data
    return rc

def setNodeText(node,text,xmldocument=None):
    if xmldocument==None:
        xmldocument = node
        while xmldocument.parentNode!=None: xmldocument=xmldocument.parentNode
    for ch in node.childNodes:
        if ch.nodeType == ch.TEXT_NODE:
            node.removeChild(ch)
            ch.unlink()
    val = xmldocument.createTextNode(text)
    node.insertBefore(val,node.firstChild)
    
def removeNodeChildren(node):
    for ch in node.childNodes:
        node.removeChild(ch)
        ch.unlink()

def copyNode(sourcenode,newparent,targetdoc=None):
    if newparent==None:
        if targetdoc==None:
            impl = xml.dom.minidom.getDOMImplementation()
            targetdoc = impl.createDocument(None, None, None)
        newparent = targetdoc
    elif targetdoc==None:
        targetdoc = newparent
        while targetdoc.parentNode!=None: targetdoc = targetdoc.parentNode
    cpy = None
    if sourcenode.nodeType==sourcenode.ELEMENT_NODE:
        cpy = targetdoc.createElement(sourcenode.localName)
        cpy = newparent.appendChild(cpy)
        for key in sourcenode.attributes.keys():
            cpy.setAttribute(key,sourcenode.getAttribute(key))
    elif sourcenode.nodeType==sourcenode.TEXT_NODE:
        cpy = targetdoc.createTextNode(sourcenode.data)
        cpy = newparent.appendChild(cpy)
    else:
        print 'WARNING: do not know how to copy node with type %s. Skipping...' % sourcenode.nodeType
    if cpy!=None:
        for ch in sourcenode.childNodes: copyNode(ch,cpy,targetdoc)
    return cpy

# ------------------------------------------------------------------------------------------
# Numerical helper utilities
# ------------------------------------------------------------------------------------------

# Look for boundary indices of array based on desired value range.
def findindices(bounds,data):
    # Zero-based indices!
    start = 0
    stop = len(data)-1
    if bounds!=None:
        if bounds[0]!=None:
            while start<len(data) and data[start]<bounds[0]: start+=1
        if bounds[1]!=None:
            while stop>=0         and data[stop] >bounds[1]: stop-=1

        # Greedy: we want take the interval that fully encompasses the specified range.
        if start>0 and start<len(data):  start-=1
        if stop<len(data)-1 and stop>=0: stop +=1

        # Note: a start beyond the available range, or a stop before it, will now have resulted
        # in start index > stop index, i.e., an invalid range. The calling function must be able
        # to handle this scenario.
        
    return (start,stop)

# 1D linear inter- and extrapolation.
def interp1(x,y,X):
    assert x.ndim==1, 'Original coordinates must be supplied as 1D array.'
    assert X.ndim==1, 'New coordinates must be supplied as 1D array.'
    
    # Transpose because it is easiest in numpy to operate on last axis. (this because of
    # upcasting rules)
    y = y.transpose()
    
    # Create array to hold interpolated values
    Y = matplotlib.numerix.empty(y.shape[0:-1]+(X.shape[0],),matplotlib.numerix.typecode(y))
    
    # Find indices of interpolated X in original x.
    iX = x.searchsorted(X)
    
    # Get the bounds of valid indices (index<=0 means point before x-data, index>=x.shape[0]
    # means point beyond x-data)
    bounds = iX.searchsorted((0.5,x.shape[0]-0.5))
    
    # Shortcuts to indices for left and right bounds, and the left bound values.
    iX_high = iX[bounds[0]:bounds[1]]
    iX_low = iX_high-1
    x_low  = x[iX_low]
    y_low = y[...,iX_low]
    
    # Linear interpolation
    Y[...,bounds[0]:bounds[1]] = y_low + ((X[bounds[0]:bounds[1]]-x_low)/(x[iX_high]-x_low))*(y[...,iX_high]-y_low)
    
    # Set points beyond bounds to extreme values.
    Y[...,0:bounds[0]] = y.take((0,),-1)
    Y[...,bounds[1]:] = y.take((-1,),-1)
    
    # Undo the original transpose and return
    return Y.transpose()
    
    # Below: previous code for Python-based linear interpolation. Slow!!
    
    #newshape = [X.shape[0]]
    #for i in y.shape[1:]: newshape.append(i)
    #Y = matplotlib.numerix.zeros(newshape,matplotlib.numerix.typecode(y))
    #icurx = 0
    #for i in range(X.shape[0]):
    #    while icurx<x.shape[0] and x[icurx]<X[i]: icurx+=1
    #    if icurx==0:
    #        Y[i,:] = y[0,:]
    #    elif icurx>=x.shape[0]:
    #        Y[i,:] = y[-1,:]
    #    else:
    #        rc = (y[icurx,:]-y[icurx-1,:])/(x[icurx]-x[icurx-1])
    #        Y[i,:] = y[icurx-1,:] + rc*(X[i]-x[icurx-1])
    #return Y
