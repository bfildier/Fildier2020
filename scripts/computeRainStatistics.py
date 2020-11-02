# -*- coding: utf-8 -*-
"""
Compute rainfall distributions.
"""

# load modules
import sys,os,glob
import re
import matplotlib
import matplotlib.pyplot as plt
import xarray as xr
import argparse
import pickle

## Add own library to path
workdir = os.path.dirname(os.path.realpath(__file__))
thismodule = sys.modules[__name__]
moduledir = os.path.join(os.path.dirname(workdir),'src')
functionsdir = os.path.join(os.path.dirname(workdir),'functions')
sys.path.insert(0,moduledir)
sys.path.insert(0,functionsdir)
for includedir in [moduledir,functionsdir]:
    print("Own modules available:", [os.path.splitext(os.path.basename(x))[0]
                                     for x in glob.glob(os.path.join(includedir,'*.py'))])

from conditionalstats import *
#from plot1DInvLog import *
#from plot2D import *
from importingData import *
#from savingResults import *


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

#%% Main script

if __name__ == "__main__":

#%% 
    
    # Command-line arguments
    parser = argparse.ArgumentParser(description="Calculate rain statistics for given SAM simulation")
    parser.add_argument("-s","--simname", required=True, help="simulation name, RCE_*")
    parser.add_argument("-n","--ndays", default=50, help="number of days to analyze, starting from end")
    args = parser.parse_args()
    simname = args.simname
    ndays = args.ndays
    # simname = 'RCE_MPDATAxTKExCAMxSAM1MOM_4000x4000x15_256x256x64_TKE-SST300-r1'
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

    # Load data
    filepattern_2D = os.path.join(archivedir,simname,"OUT_2D","%s_%s.2Dcom_*.nc"%(simname,Nproc))
    varids = 'Prec',
    data = xr.open_mfdataset(filepattern_2D,decode_cf=False,data_vars=varids)
    #data = xr.open_mfdataset(filepattern_2D,decode_cf=False,data_vars=varids,combine='by_coords')
    
    
    ##-- Compute rain statistics
    
    print('- computing statistics')
    print()

    # slice to analyze
    NTmax = floor(len(data.time)/10/24)*10*24
    s_end = slice(NTmax-24*ndays,NTmax)
    pr = data.Prec.values[s_end,:,:]
    
    #- Compute rain distribution    
    # Initialize
    dist_pr_IL = Distribution(name='pr',bintype='invlogQ',fill_last_decade=True)
    # Compute
    dist_pr_IL.computeDistribution(sample=pr)
    # Store locations of reference points in each bin
    dist_pr_IL.storeSamplePoints(sample=pr,verbose=True)
    # Compute inverse CDF on IL bins (fraction of rain mass above percentile)
    dist_pr_IL.computeInvCDF(sample=pr)
    # Compute bootstrapping to estimate uncertainty on percentiles
    nd_resample = 10*24 # time slices for 10 days
    dist_pr_IL.bootstrapPercentiles(sample=pr,nd_resample=nd_resample)
    
    # #- show
    # fig = plt.figure(figsize=(6,5))
    # ax = fig.add_subplot(111)
    # sQref = slice(0,49)
    # # subplotRanksILog(ax,dist_pr_IL.ranks[sQref],dist_pr_IL.invCDF[sQref]*100,\
    # #                       transformX=False)
    # subplotRanksILog(ax,dist_pr_IL.ranks[sQref],dist_pr_IL.percentiles[sQref]*100,\
    #                       transformX=False)
    # subplotRanksILog(ax,dist_pr_IL.ranks[sQref],dist_pr_IL.percentiles_Q1[sQref]*100,\
    #                       transformX=False)
    # subplotRanksILog(ax,dist_pr_IL.ranks[sQref],dist_pr_IL.percentiles_Q2[sQref]*100,\
    #                       transformX=False)
    # subplotRanksILog(ax,dist_pr_IL.ranks[sQref],dist_pr_IL.percentiles_Q3[sQref]*100,\
    #                       transformX=False)
    # x = np.flipud(1./(1-dist_pr_IL.ranks[sQref]/100.))
    # transformXaxisIL(ax,x)
    # # iQ_min = 4
    # # ax.set_xlim((x[iQ_min-1],0.8))
    # ax.set_xlabel('Percentile rank Q (%)')
    # ax.set_ylabel('Mass fraction above percentile (%)')
    # plt.show()

    #- Compute mean
    mean_pr = np.mean(pr)
    
    ##-- Save
    print('- saving')
    print()
    
    #- mean
    file_mean_pr = os.path.join(resultdir,'mean_pr.pickle')
    pickle.dump(mean_pr,open(file_mean_pr,'wb'))
    #- rain stats
    file_dist_pr = os.path.join(resultdir,'dist_pr_IL.pickle')
    pickle.dump(dist_pr_IL,open(file_dist_pr,'wb'))
    
    print('Success :)')
    
