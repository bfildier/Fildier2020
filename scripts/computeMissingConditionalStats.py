#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute thermodynamic and dynamic contributions to extreme rainfall
@author: bfildier
"""

# load modules
import sys,os,glob
import re
import matplotlib
import matplotlib.pyplot as plt
import xarray as xr
import argparse
import pickle
from scipy.optimize import curve_fit

## Add own library to path
workdir = os.path.dirname(os.path.realpath(__file__))
#workdir = '/Users/bfildier/Code/analyses/Fildier2020/scripts'
thismodule = sys.modules[__name__]
moduledir = os.path.join(os.path.dirname(workdir),'src')
functionsdir = os.path.join(os.path.dirname(workdir),'functions')
sys.path.insert(0,moduledir)
sys.path.insert(0,functionsdir)
for includedir in [moduledir,functionsdir]:
    print("Own modules available:", [os.path.splitext(os.path.basename(x))[0]
                                     for x in glob.glob(os.path.join(includedir,'*.py'))])

from conditionalstats import *
from plot1DInvLog import *
#from plot2D import *
from importingData import *
#from savingResults import *
from scalingApproximations import *


def getSimulationInfo(simname):
    
    chunks = simname.split('_')
    case = chunks[0]
    Nxyz = chunks[3]
    Nx,Ny,Nz = np.array(Nxyz.split('x'),dtype=int)
    Nproc = np.max(np.array(Nxyz.split('x'),dtype=int))
    simroot = '_'.join(chunks[:4])
    expname = chunks[4]
    expchunks = expname.split('-')
    SST = int([expchunk[-3:] for expchunk in expchunks if re.match("SST*",expchunk)][0])
    
    return simroot, expname, SST, Nproc


#%%
if __name__ == '__main__':

    # Command-line arguments
    parser = argparse.ArgumentParser(description="Calculate rain statistics for given SAM simulation")
    parser.add_argument("-s","--simname", required=True, help="simulation name, RCE_*")
    parser.add_argument("-n","--ndays", default=50, help="number of days to analyze, starting from end")
    args = parser.parse_args()
    simname = args.simname
    ndays = args.ndays
    # simname = 'RCE_MPDATAxTKExCAMxSAM1MOM_4000x4000x15_256x256x64_TKE-SST308-r2'
    # ndays = 1 # days
    
    simroot, expname, SST, Nproc = getSimulationInfo(simname)
    
    print()
    print('* Analyzing run %s'%simname)
    print()
    
    # Directories
    archivedir = getArchivedir(machine='coriknl')
    resultdir = os.path.join(os.path.dirname(moduledir),'results',simroot,expname,
                             'stats_%ddays'%ndays)
    os.makedirs(resultdir,exist_ok=True)
    figuredir = os.path.join(os.path.dirname(moduledir),'figures')
    

    # Load 2D data
    filepattern_2D = os.path.join(archivedir,simname,"OUT_2D","%s_%s.2Dcom_*.nc"%(simname,Nproc))
    varids = 'Prec',
#    data2D = xr.open_mfdataset(filepattern_2D,decode_cf=False,data_vars=varids,combine='by_coords')
    vars2D2drop = ['PBLH', 'SHF', 'LHF', 'LWNS', 'LWNSC', 'LWNT', 'LWNTC',\
            'SOLIN', 'SWNS', 'SWNSC', 'SWNT', 'SWNTC', 'CWP', 'IWP', 'CLD', 'USFC', 'U200',\
            'VSFC', 'V200', 'W730', 'PSFC', 'SWVP', 'U850', 'V850', 'ZC', 'TB', 'ZE']
    data2D = xr.open_mfdataset(filepattern_2D,decode_cf=False,data_vars=varids,drop_variables=vars2D2drop)

    # Load 3D data
    if 'radhomo' in simname:
        NDmax = 100
    else:
        NDmax = int(data2D.time.values[-1])
    stepmax = NDmax*24*60*60//15
    stepmin = stepmax - ndays*24*60*60//15
    offset = 1*60*60//15 # = 1 hr lag
    
    steprange = (stepmin-offset,stepmax-offset)
    files_in_steprange = get3DFilesBetweenSteps(archivedir,simname,steprange)
#    data3D = xr.open_mfdataset(files_in_steprange,decode_cf=False,combine='by_coords')
    vars3D2drop = ['W','TABS','U', 'V', 'PP', 'QRAD', 'LWU', 'LWD', 'LWUS', 'LWDS', 'QN', 'QP']
    data3D = xr.open_mfdataset(files_in_steprange,decode_cf=False,drop_variables=vars3D2drop)

    ##-- Loading rain statistics
    
    print('- loading rain statistics')
    print()
    
    #- mean
    file_mean_pr = os.path.join(resultdir,'mean_pr.pickle')
    mean_pr = pickle.load(open(file_mean_pr,'rb'))
    #- rain stats
    file_dist_pr = os.path.join(resultdir,'dist_pr_IL.pickle')
    dist_pr_IL = pickle.load(open(file_dist_pr,'rb'))
    

    ##-- Compute conditional statistics

    # slice to analyze
    NTmax = floor(len(data2D.time)/10/24)*10*24
    s_end = slice(NTmax-24*ndays,NTmax)    
    
    # data
    qv = data3D.QV.values
    
    #- conditional statistics
    # define
    cdist_qv_on_pr_IL = ConditionalDistribution(name='QV',is3D=True,isTime=True,on=dist_pr_IL)
    # compute
    cdist_qv_on_pr_IL.computeConditionalMeanAndVariance(qv)


    ##-- Save results
    print('- saving')
    print()
    
    #- conditional statistics
    file_cdist_qv = os.path.join(resultdir,'cdist_qv_on_pr_IL.pickle')
    pickle.dump(cdist_qv_on_pr_IL,open(file_cdist_qv,'wb'))

    print('Success :)')

    sys.exit(0)
