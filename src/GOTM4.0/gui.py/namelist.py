import re

# ------------------------------------------------------------------------------------------
# Namelist parsing utilities
# ------------------------------------------------------------------------------------------

class NamelistParseException(Exception):
    def __init__(self,error,filename=None,namelistname=None,variablename=None):
        Exception.__init__(self,error)
        self.filename     = filename
        self.namelistname = namelistname
        self.variablename = variablename

    def __str__(self):
        return Exception.__str__(self)+'.\nFile: '+str(self.filename)+', namelist: '+str(self.namelistname)+', variable: '+str(self.variablename)

class NamelistSubstitutions:
    # Commonly used regular expression (for parsing substitions in the .values file)
    subs_re = re.compile('\s*s/(\w+)/(.+)/')

    def __init__(self,valuesfile):
        self.subs = []
        
        # If the supplied argument is a string, it should be a path to a file.
        # Try to open it. Otherwise the supplied argument should be a file-like object.
        ownfile = False
        if isinstance(valuesfile,basestring):
            try:
                valuesfile = open(valuesfile,'rU')
            except Exception,e:
                raise NamelistParseException('Cannot open .values file "%s". Error: %s' % (path,str(e)))
            ownfile = True
            
        line = valuesfile.readline()
        while line!='':
            m = self.subs_re.match(line)
            if m!=None:
                #pat = re.compile(m.group(1),re.IGNORECASE)
                #self.subs.append((pat,m.group(2)))
                self.subs.append((m.group(1),m.group(2)))
            line = valuesfile.readline()
            
        # If we opened the file ourselves, close it.
        if ownfile:
            valuesfile.close()

    def substitute(self,text):
        for (old,new) in self.subs:
            #text = old.sub(new,text)
            text = text.replace(old,new)
        return text

class Namelist:

    varassign_re = re.compile('\s*(\w+)\s*=\s*')
    varstopchar_re = re.compile('[/,\n"\']')

    def __init__(self,name,data,filepath=None):
        self.name = name
        self.data = data
        self.filepath = filepath

    def __iter__(self):
        return self

    def next(self):
        ret = self.getNextVariable()
        if ret==None: raise StopIteration
        return ret

    def getNextVariable(self):
        if self.isEmpty(): return None
        
        match = self.varassign_re.match(self.data);
        if match==None:
            raise NamelistParseException('Cannot find a variable assignment (variable_name = ...). Current namelist data: "%s"' % (self.data,),self.filepath,self.name,None)
        foundvarname = match.group(1)
        
        self.data = self.data[match.end():]
        ipos = 0
        while True:
            match = self.varstopchar_re.search(self.data,pos=ipos)
            if match==None:
                raise NamelistParseException('End of variable value not found. Current namelist data: "%s"' % (self.data,),self.filepath,self.name,foundvarname)
            ch = match.group(0)
            if ch=='\'' or ch=='"':
                ipos = match.end(0)
                inextquote = self.data.find(ch,ipos)
                if inextquote==-1:
                    raise NamelistParseException('Opening quote %s was not matched by end quote.' % (ch,),self.filepath,self.name,foundvarname)
                ipos = inextquote+1
            else:
                # Found end of variable value.
                vardata = self.data[0:match.start()].rstrip()
                self.data = self.data[match.end():].lstrip()
                break

        return (foundvarname,vardata)

    def isEmpty(self):
        return len(self.data)==0

class NamelistFile:
    
    # Commonly used regular expressions, for:
    #   - locating the start of comments, or opening quotes.
    #   - locating the start of namelist (ampersand followed by list name).
    #   - locating the end of a namelist, i.e. a slash, or opening quotes.
    commentchar_re = re.compile('[!#"\']')
    namelistname_re = re.compile('\s*&\s*(\w+)\s*')
    stopchar_re = re.compile('[/"\']')

    def __init__(self,nmlfile,subs=[]):
        ownfile = False
        if isinstance(nmlfile,basestring):
            ownfile = True
            
            # Attempt to open namelist file and read all data
            try:
                nmlfile = open(nmlfile,'rU')
            except Exception,e:
                raise NamelistParseException('Cannot open namelist file. Error: %s' % (str(e),),path)
            self.path = nmlfile
        else:
            self.path = ''

        try:
            
            self.data = ''
            line = nmlfile.readline()
            iline = 1
            while line!='':

                # Strip comments, i.e. on every line, remove everything after (and including) the first exclamation
                # mark; ignore text between single and double quotes.
                ipos = 0
                match = self.commentchar_re.search(line,pos=ipos)
                while match!=None:
                    ch = match.group(0)
                    if ch=='\'' or ch=='"':
                        ipos = match.end(0)
                        inextquote = line.find(ch,ipos)
                        if inextquote==-1:
                            raise NamelistParseException('Line %i: opening quote %s was not matched by end quote.' % (iline,ch),path)
                        ipos = inextquote+1
                    else:
                        # Found start of comment; only keep everything preceding the start position.
                        line = line[:match.start(0)]
                        break
                    match = self.commentchar_re.search(line,pos=ipos)

                self.data += line
                line = nmlfile.readline()
                iline += 1
        finally:
            if ownfile: nmlfile.close()

        # Make substitutions (if any).
        for sub in subs:
            self.data = sub.substitute(self.data)

    def parseNextNamelist(self,expectedlist=None):
        match = self.namelistname_re.match(self.data)
        if match==None:
            raise NamelistParseException('No namelist found; expected ampersand followed by namelist name.',self.path)
        name = match.group(1)
        if expectedlist!=None and name!=expectedlist:
            raise NamelistParseException('Expected namelist "%s", but found "%s".' % (expectedlist,name),self.path,expectedlist)
        istart = match.end(0)
        ipos = istart
        while True:
            match = self.stopchar_re.search(self.data,pos=ipos)
            if match==None:
                raise NamelistParseException('End of namelist (slash) not found.',self.path,name)
            ch = match.group(0)
            ipos = match.end(0)
            if ch=='/':
                break
            else:
                inextquote = self.data.find(ch,ipos)
                if inextquote==-1:
                    raise NamelistParseException('Opening quote %s was not matched by end quote.' % (ch,),self.path,name)
                ipos = inextquote+1
        namelistdata = self.data[istart:ipos-1]
        self.data = self.data[ipos:]
        return Namelist(name,namelistdata,filepath=self.path)
