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
    data2D = xr.open_mfdataset(filepattern_2D,decode_cf=False,data_vars=varids,combine='by_coords')
    
    # Load 3D data
    NDmax = int(data2D.time.values[-1])
    stepmax = NDmax*24*60*60//15
    stepmin = stepmax - ndays*24*60*60//15
    offset = 1*60*60//15 # = 1 hr lag
    
    steprange = (stepmin-offset,stepmax-offset)
    files_in_steprange = get3DFilesBetweenSteps(archivedir,simname,steprange)
    data3D = xr.open_mfdataset(files_in_steprange,decode_cf=False,combine='by_coords')

    ##-- Loading rain statistics
    
    print('- loading rain statistics')
    print()
    
    # #- mean
    # file_mean_pr = os.path.join(resultdir,'mean_pr.pickle')
    # mean_pr = pickle.load(open(file_mean_pr,'rb'))
    # #- rain stats
    # file_dist_pr = os.path.join(resultdir,'dist_pr_IL.pickle')
    # dist_pr_IL = pickle.load(open(file_dist_pr,'rb'))
    

    ##-- Compute conditional statistics

    # slice to analyze
    NTmax = floor(len(data2D.time)/10/24)*10*24
    s_end = slice(NTmax-24*ndays,NTmax)    
    # s_end = slice(NTmax-1,NTmax)
    
    # data
    pr = data2D.Prec.values[s_end,:,:]
    w = data3D.W.values
    tabs = data3D.TABS.values
    qp = data3D.QP.values
    p_profile = np.array(np.mean(data3D.p,axis=0))
    # p_profile = np.array((data3D.p))
    z_coord = data3D.z.values
    
    
    #- pr statistics
    dist_pr_IL = Distribution(name='pr',bintype='invlogQ',fill_last_decade=True,
                              distribution=None,overwrite=False)
    # Compute
    dist_pr_IL.computeDistribution(sample=pr)
    # Store locations of reference points in each bin
    dist_pr_IL.storeSamplePoints(sample=pr,verbose=True,sizemax=500)
    # Compute inverse CDF on IL bins (fraction of rain mass above percentile)
    dist_pr_IL.computeInvCDF(sample=pr)
    # Compute bootstrapping to estimate uncertainty on percentiles
    nd_resample = 10*24 # time slices for 10 days
    dist_pr_IL.bootstrapPercentiles(sample=pr,nd_resample=nd_resample)
    # Compute individual percentiles for error bar on the mean
    dist_pr_IL.computeIndividualPercentiles(sample=pr,ranks=[5,25,50,75,95])
    
    #- Compute mean rainfall
    mean_pr = np.mean(pr)
    
    #- conditional statistics
    # define
    cdist_w_on_pr_IL = ConditionalDistribution(name='W',is3D=True,isTime=True,on=dist_pr_IL)
    cdist_tabs_on_pr_IL = ConditionalDistribution(name='T',is3D=True,isTime=True,on=dist_pr_IL)
    # cdist_qp_on_pr_IL = ConditionalDistribution(name='QP',is3D=True,isTime=True,on=dist_pr_IL)
    # compute
    cdist_w_on_pr_IL.computeConditionalMeanAndVariance(w)
    cdist_tabs_on_pr_IL.computeConditionalMeanAndVariance(tabs)
    # cdist_qp_on_pr_IL.computeConditionalMeanAndVariance(qp)

    ##-- Compute OGS09 scaling approximation and thermodynamic/dynamic components
    pr_OGS09 = np.nan*np.zeros(len(dist_pr_IL.ranks))
    M = np.nan*np.zeros(len(dist_pr_IL.ranks))
    Gamma = np.nan*np.zeros(len(dist_pr_IL.ranks))
    
    for iQ in range(len(dist_pr_IL.ranks)):
        
        # scaling
        w_prof = cdist_w_on_pr_IL.cond_mean[:,iQ]
        tabs_prof = cdist_tabs_on_pr_IL.cond_mean[:,iQ]
        pr_OGS09[iQ] = -scalingOGS09_zcoord(w_prof,tabs_prof,p_profile,z_coord,levdim=0)
        
        # Crop input profiles at the tropopause
        p_sc, tabs_sc, w_sc, z_sc = cropProfiles(p_profile*100, # Pa
                                                 tabs_prof,
                                                 [w_prof,z_coord.data],
                                                 levdim=0)
        
        # Average mass flux
        M[iQ] = verticalPressureIntegral(p_sc,w_sc) / verticalPressureIntegral(p_sc)
        
        # Thermodynamic component
        qvstar_sc = saturationSpecificHumidity(tabs_sc,p_sc)
        dqvs_dz = np.diff(qvstar_sc)/np.diff(z_sc)
        Gamma[iQ] = verticalPressureIntegral(p_sc,values=w_sc,dvdp=dqvs_dz)/M[iQ]
        

    #- precipitation efficiency
    def computePE(perc,scaling,sQ):

        x = scaling[sQ]
        y = perc[sQ]/86400
        f = lambda x,alpha: alpha*x
    
        eps, ecov = curve_fit(f,x,y,p0=1)
        return eps[0], ecov[0][0]
    
    sQ9999_99999 = slice(39,49)
    eps, ecov = computePE(dist_pr_IL.percentiles,pr_OGS09,sQ9999_99999)

    # combine in a single object
    scaling = {'OGS09': eps*pr_OGS09,
               'eps': eps,
               'ecov': ecov,
               'M': M,
               'Gamma': Gamma}

    ##-- Save results
    print('- saving')
    print()
    
    #- mean
    file_mean_pr = os.path.join(resultdir,'mean_pr.pickle')
    pickle.dump(mean_pr,open(file_mean_pr,'wb'))
    #- rain stats
    file_dist_pr = os.path.join(resultdir,'dist_pr_IL.pickle')
    pickle.dump(dist_pr_IL,open(file_dist_pr,'wb'))
    #- conditional statistics
    file_cdist_w = os.path.join(resultdir,'cdist_w_on_pr_IL.pickle')
    pickle.dump(cdist_w_on_pr_IL,open(file_cdist_w,'wb'))
    file_cdist_tabs = os.path.join(resultdir,'cdist_tabs_on_pr_IL.pickle')
    pickle.dump(cdist_tabs_on_pr_IL,open(file_cdist_tabs,'wb'))
    #- OGorman scaling approximation
    # file_eps= os.path.join(resultdir,'eps.pickle')
    # pickle.dump(eps,open(file_eps,'wb'))
    # file_ecov= os.path.join(resultdir,'ecov.pickle')
    # pickle.dump(ecov,open(file_ecov,'wb'))
    # file_pr_OGS09= os.path.join(resultdir,'pr_OGS09.pickle')
    # pickle.dump(pr_OGS09,open(file_pr_OGS09,'wb'))
    file_scaling = os.path.join(resultdir,'scaling.pickle')
    pickle.dump(scaling,open(file_scaling,'wb'))
    


    print('Success :)')

    sys.exit(0)
