#!/usr/bin/env python

from __future__ import print_function
import os
import subprocess
import re
import sys

def gen_version_svn():
    version_template     = "      character(LEN=32),parameter :: schism_version = '@{VERSION_SCHISM}', svn_build = '@{VERSION_SVN}' " 
    version_path    = os.path.split( __file__)[0]
    query_path = os.path.join(version_path,"..")

    versionscratch_path = os.path.join(version_path,"_version")

    try:
        with open(versionscratch_path, "w") as versionscratch:
            ok = subprocess.check_call(["svnversion", query_path],stdout=versionscratch)
        with open(versionscratch_path,"r") as versionscratch:
            version_raw = versionscratch.readlines()[0].strip()

        with open(versionscratch_path, "w") as versionscratch:
            ok = subprocess.check_call(["svn","info",query_path],stdout=versionscratch)
        with open(versionscratch_path,"r") as versionscratch:
            url_raw = [x.strip()[4:] for x in versionscratch.readlines() if "URL:" in x][0]

        url = url_raw.split("/")

        release_branch_re = r"v?[0-9]{1,2}\.[0-9]{1,3}"
        official_tag_re = r"v?[0-9]{1,2}\.[0-9]{1,3}\.[0-9]{1,3}"
        if (url[-3] == "tags" and re.match(official_tag_re,url[-2])):
            schism_version = url[-2]
        elif (url[-3] == "branches" and re.match(release_branch_re,url[-2])):
            schism_version = url[-2] + " (untagged release)"
        elif url[-2] == "trunk":
            schism_version = "5.x.dev (trunk)"
        else:
            schism_version = "5.x.dev (branch %s)" % url[-2]    
        return version_raw,schism_version
    
    except Exception as inst:
        print(inst)
        if os.path.exists(versionscratch_path): 
            os.remove(versionscratch_path) 
        print("Offline from SVN utlities")
        return None,None

def gen_version_user(default_version = "5.6.1 or later"):
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
    version_template     = "      character(LEN=32),parameter :: schism_version = '@{VERSION_SCHISM}', svn_build = '@{VERSION_SVN}' " 
    version_path    = os.path.split( __file__)[0]
    template_path = os.path.join(version_path,"schism_version.F90.template")
    query_path = os.path.join(version_path,"..")
    if versionfile_path is None:
        scriptpath=os.path.dirname(os.path.realpath(__file__))
        versionfile_path=os.path.join(scriptpath,"schism_version.F90")

    version_raw,schism_version = gen_version_svn()
    if schism_version is None:
        print("SVN data not available")
        schism_version = gen_version_user()
        version_raw = "No SVN Revision Info"
  
    print(' SCHISM version:  %s;' % schism_version)
    print(' SVN version of schism:        '+ version_raw)
    with open(template_path,"r") as template:
        templatetxt = template.read()
    versiontxt = templatetxt.replace("@{VERSION_SVN}", version_raw).replace("@{VERSION_SCHISM}",schism_version)
    with open(versionfile_path,"w") as versionfile:
        versionfile.write(versiontxt)   



if __name__=="__main__":
    if len(sys.argv) > 1:   
        outputfile = sys.argv[1]
    else:
        outputfile = None
    gen_version(outputfile)
