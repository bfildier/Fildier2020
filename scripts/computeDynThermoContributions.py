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

## Add own library to path
# workdir = os.path.dirname(os.path.realpath(__file__))
workdir = '/Users/bfildier/Code/analyses/Fildier2020/scripts'
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
    data2D = xr.open_mfdataset(filepattern_2D,decode_cf=False,data_vars=varids)
    
    # Load 3D data
    NDmax = int(data2D.time.values[-1])
    stepmax = NDmax*24*60*60//15
    stepmin = stepmax - ndays*24*60*60//15
    offset = 1*60*60//15 # = 1 hr lag
    # offset = 0
    
    steprange = (stepmin-offset,stepmax-offset)
    files_in_steprange = get3DFilesBetweenSteps(archivedir,simname,steprange)
    
    data3D = xr.open_mfdataset(files_in_steprange,decode_cf=False)
    # data = []
    # for file in files_in_steprange:
    #     data.append(xr.open_dataset(file))
    # data3D = xr.concat(data,dim='time')

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
    pr = data2D.Prec.values[s_end,:,:]
    w = data3D.W.values
    tabs = data3D.TABS.values
    p_profile = np.array(np.mean(data3D.p,axis=0))
    z_coord = data3D.z.values
    
    
    # #- pr statistics
    # computePrStats = True
    # if computePrStats:
    #     dist_pr_IL = Distribution(name='pr',bintype='invlogQ',fill_last_decade=True,
    #                               distribution=None,overwrite=False)
    #     # Compute
    #     dist_pr_IL.computeDistribution(sample=pr)
    #     # Store locations of reference points in each bin
    #     dist_pr_IL.storeSamplePoints(sample=pr,verbose=True)
    #     # Compute inverse CDF on IL bins (fraction of rain mass above percentile)
    #     dist_pr_IL.computeInvCDF(sample=pr)
    #     # Compute bootstrapping to estimate uncertainty on percentiles
    #     nd_resample = 10*24 # time slices for 10 days
    #     dist_pr_IL.bootstrapPercentiles(sample=pr,nd_resample=nd_resample)
    #     # Compute individual percentiles for error bar on the mean
    #     dist_pr_IL.computeIndividualPercentiles(sample=pr,ranks=[5,25,50,75,95])
    
    
    # define
    cdist_w_on_pr_IL = ConditionalDistribution(name='W',is3D=True,isTime=True,on=dist_pr_IL)
    cdist_tabs_on_pr_IL = ConditionalDistribution(name='T',is3D=True,isTime=True,on=dist_pr_IL)
    # compute conditional distributions
    cdist_w_on_pr_IL.computeConditionalMeanAndVariance(w)
    cdist_tabs_on_pr_IL.computeConditionalMeanAndVariance(tabs)

    ##-- Compute OGS09 scaling approximation
    pr_OGS09 = np.nan*np.zeros(len(dist_pr_IL.ranks))
    for iQ in range(len(dist_pr_IL.ranks)):
        w_prof = cdist_w_on_pr_IL.cond_mean[:,iQ]
        tabs_prof = cdist_tabs_on_pr_IL.cond_mean[:,iQ]
        pr_OGS09[iQ] = scalingOGS09_zcoord(w_prof,tabs_prof,p_profile,z_coord,levdim=0)
    
    
    
    
    
    # TEST bin locations for pr
    coords_pr999 = dist_pr_IL.bin_locations[39]
    pr999_vals = []
    w999_vals = [[] for iz in range(len(p_profile))]
    for ind in coords_pr999:
        pr999_vals.append(np.take(pr.flatten(),ind))
        for iz in range(len(p_profile)):
            w999_vals[iz].append(np.take(w[:,iz,:,:].flatten(),ind))
    pr999_vals = np.array(pr999_vals)
    w999_vals = np.array(w999_vals)
    print(pr999_vals.mean(),pr999_vals.std())
    print(pr999_vals)
    w999_mean = np.mean(w999_vals,axis=1)
    plt.plot(w999_mean,p_profile)
    plt.plot(cdist_w_on_pr_IL.cond_mean[:,39],p_profile)
    
    
    # DEBUG
    sample = w.copy()
    sshape = sample.shape
    # Initialize default output
    sample_out = sample
    # Get dimensions and adjust output shape
    if self.is3D:
        # collapse dimensions other than z
        if len(sshape) > 2: # reshape
            # if time dimension (in 2nd dimension), reorder to have z in first dim
            if self.isTime:
                # nlev = sample_out.shape[0]
                sample_out = np.swapaxes(sample_out,0,1)
                sshape = sample_out.shape
            sample_out = np.reshape(sample_out,(sshape[0],np.prod(sshape[1:])))
        Nz,Npoints = sample_out.shape
        # if self.isTime:
        #     Npoints = Npoints/nlev
        #     print(Npoints)
    else:
        if len(sshape) > 1: # reshape
            sample_out = np.reshape(sample,np.prod(sshape))
        Npoints, = sample_out.shape
        Nz = None
    
    
    
    
    
    
    showfigure = True
    if showfigure:
        
        fig,ax = plt.subplots()
        subplotRanksILog(ax,dist_pr_IL.ranks,dist_pr_IL.percentiles,transformX=False)
        subplotRanksILog(ax,dist_pr_IL.ranks,pr_OGS09*100000,transformX=False)
        x = np.flipud(1./(1-dist_pr_IL.ranks/100.))
        transformXaxisIL(ax,x)
        plt.show()

        iQ = 42
        fig,ax = plt.subplots()
        ax.plot(cdist_w_on_pr_IL.cond_mean[:,iQ],p_profile)
        ax.fill_betweenx(p_profile,
                         cdist_w_on_pr_IL.cond_mean[:,iQ]-cdist_w_on_pr_IL.cond_std[:,iQ],
                         cdist_w_on_pr_IL.cond_mean[:,iQ]+cdist_w_on_pr_IL.cond_std[:,iQ],
                         alpha=0.2)
        ax.invert_yaxis()
        plt.show()
    
    
    