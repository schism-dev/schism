#!/usr/bin/env python

from __future__ import print_function
import os
import subprocess
import re
import sys

def get_git_revision_hash():
    try:
        ver = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD']).strip().decode('utf-8')
    except:
        ver = 'none'
    return ver

def gen_version_git():
    version_template     = "      character(LEN=32),parameter :: schism_version = '@{VERSION_SCHISM}', git_version = '@{VERSION_GIT}' " 
    version_path    = os.path.split( __file__)[0]
    query_path = os.path.join(version_path,"..")

    versionscratch_path = os.path.join(version_path,"_version")
    print(__file__)
    print(version_path)
    print(versionscratch_path)
    describe_re = re.compile(r"(?P<version>v?\d+\.\d+\.\d+)(?P<update>-\d+)?(?P<hash>-g[a-z0-9]+)?")

    try:
        with open(versionscratch_path, "w") as versionscratch:
            ok = subprocess.check_call(["git","describe","--always","--dirty"],stdout=versionscratch)
        with open(versionscratch_path,"r") as versionscratch:
            version_raw = versionscratch.readlines()[0].strip()

        
        is_dirty = "-dirty" in version_raw
        version_raw=version_raw.replace("-dirty","")
        m = describe_re.match(version_raw)
        if m is None:
            git_version = get_git_revision_hash()
            return git_version,None
        if m is None and len(version_raw) > 5: 
            schism_version = "semantic version not determined"
            git_version = get_git_revision_hash()
        else:
            schism_version = m.group("version")
            print(m.group("update"))
            if m.lastindex == 1:
                # This is an untouched tag with no commits
                git_version = get_git_revision_hash()
                print(git_version)
            if m.lastindex >= 2: 
                git_version = m.group("hash").replace("-g","")      
                if is_dirty or len(m.group("update"))>2:
                    nmod = m.group("update").replace("-","")
                    schism_version = schism_version+"mod"
                    git_version = git_version + " ({} commits since semantic tag, edits={})".format(nmod,is_dirty)          
        if git_version is None: git_version = get_git_revision_hash()
        return git_version,schism_version
    
    except Exception as inst:
        print(inst)
        if os.path.exists(versionscratch_path): 
            os.remove(versionscratch_path) 
        print("Error querying \"git describe\", offline from Git utlities?\nYou can create a file called schism_version_user.txt"
              " in this directory with your own version label on the first line.\n")
        return None,None

def gen_version_user(default_version = "develop"):
    schism_user_version_file = "schism_version_user.txt"
    print("Attempting to get version text manually from first line of \nsrc/Core/%s if file exists" % schism_user_version_file)
    user_version_path    = os.path.join(os.path.split( __file__)[0],schism_user_version_file) 
    if os.path.exists(user_version_path):
        with open(user_version_path,"r") as defaultfile:
            line = defaultfile.readline().strip()
            user_version = line if len(line) >=3 else default_version
    else:
        user_version = default_version
    assert len(user_version) > 3
    return user_version       


def gen_version(versionfile_path=None):
    version_template     = "      character(LEN=32),parameter :: schism_version = '@{VERSION_SCHISM}', git_version = '@{VERSION_GIT}' " 
    version_path    = os.path.split( __file__)[0]
    template_path = os.path.join(version_path,"schism_version.F90.template")
    query_path = os.path.join(version_path,"..")
    if versionfile_path is None:
        scriptpath=os.path.dirname(os.path.realpath(__file__))
        versionfile_path=os.path.join(scriptpath,"schism_version.F90")

    git_version,schism_version = gen_version_git()
    if schism_version is None:
        print("SCHISM version not available, searching for src/schism_user_version.txt or default")
        schism_version = gen_version_user()

    if git_version is None:
       print("Git hash not inferred from git describe, using rev-parse")
       git_version = get_git_revision_hash()
       if git_version is None: git_version = "unavailable"
  
    print(git_version)
    print(' SCHISM version:  {0:s}'.format(schism_version))
    print(' GIT commit       {0:s}'.format(git_version))
    with open(template_path,"r") as template:
        templatetxt = template.read()
    versiontxt = templatetxt.replace("@{VERSION_GIT}", git_version).replace("@{VERSION_SCHISM}",schism_version)
    with open(versionfile_path,"w") as versionfile:
        versionfile.write(versiontxt)   



if __name__=="__main__":
    if len(sys.argv) > 1:   
        outputfile = sys.argv[1]
    else:
        outputfile = None
    gen_version(outputfile)
