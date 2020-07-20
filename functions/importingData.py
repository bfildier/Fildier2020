"""Module importingData

Contains
- functions to import variables between specific dates from SAM
"""

#---- Modules ----#

import glob
import sys,os
import string
import operator
import numpy as np
from dataFormat import *
from datetime import date, datetime, timedelta
from netCDF4 import Dataset

## Own functions
# currentpath = os.path.dirname(os.path.realpath(__file__))

from setEnv import *

#---- Functions ----#
def loadVar(varid,archivedir,simname,datatype='STAT',steprange=None):
    
    filepattern = os.path.join(archivedir,simname,"OUT_%s"%datatype,"%s_*.nc"%simname)
    files_all = glob.glob(filepattern)
    files_all.sort()
    files = files_all

    isSeq = varid.__class__ in [list,tuple]
    if isSeq:
        nvar = len(varid)
        varlist = [[] for j in range(nvar)]
    else:
        varlist = []
    
    for file in files:
        
        fh = Dataset(file,'r')
        if isSeq:
            for j in range(nvar):
                varlist[j].append(fh.variables[varid[j]][:])
        else:
            varlist.append(fh.variables[varid][:])
        fh.close()

    out = []
    if isSeq:
        for j in range(nvar):
            out.append(np.vstack(varlist[j]))
        out = tuple(out)
    else:
        if len(varlist[0].shape) > 1:
            out = np.vstack(varlist)
        else:
            out = np.hstack(varlist)
        
    return out

def get3DFilesBetweenSteps(archivedir,simname,steprange=None):

    filepattern = os.path.join(archivedir,simname,'OUT_3D',"%s_*.nc"%simname)
    files_all = glob.glob(filepattern)

    if steprange is None:
        files = files_all
    else:
        stepmin = steprange[0]
        stepmax = steprange[1]
        files = []
        for file in files_all:
            step = int(file.split('_')[-1].split('.')[0])
            if stepmin < step <= stepmax:
                files.append(file)

    files.sort()

    return files

def getFilesBetweenSteps(archivedir,simname,datatype='STAT',steprange=None):

    filepattern = os.path.join(archivedir,simname,'OUT_%s'%datatype,"%s_*.nc"%simname)
    files_all = glob.glob(filepattern)

    if steprange is None:
        files = files_all
    else:
        stepmin = steprange[0]
        stepmax = steprange[1]
        files = []
        for file in files_all:
            step = int(file.split('_')[-1].split('.')[0])
            if stepmin < step <= stepmax:
                files.append(file)

    files.sort()

    return files


