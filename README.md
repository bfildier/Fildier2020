
This archive gathers scripts produced in the course of the work sent for publication 
as Fildier et al. (2020), Changing Skewness of the Rainfall Distribution with Warming, 
with and without Self-Aggregation


The structure is organized as follows:

- SAM_scripts/ contains the scripts used to launch the numerical simulations with 
the model SAM (System for Atmospheric Modeling), as used on the NERSC supercomputer 
Cori:
  
  . load_dirnames.sh 

  . build_and_run/ contains scripts to automate the compilation and launch of SAM
  with a variety of boundary conditions and parameters, for both the startup runs 
  and the branch runs
  
  . bash_specs/ contains subroutines for the launch scripts
  
  . process_outputs/ contains post processing scripts to convert SAM outputs to netcdf,
  move, archive the data and load it from the long-term storage HPSS


All other folders contain scripts or results related to the analysis:

- functions/ contain user-defined python modules for tasks such as computing statistical 
  distributions, thermodynamic and physical formulas, data manipulation and plotting

- scripts/ contains scripts used to calculate statistics and notebooks used to draw figures as well as scripts used to automate the analysis in scripts/batch_scripts/

- results/RCE_MPDATAxTKExCAMxSAM1MOM_4000x4000x15_256x256x64/ contains intermediate
results used for final display, classified per simulation type. Simulation names 
are chosen using the following nomenclature: 
  
  . SSTXXX specifies a SST forcing of XXX K
  . rX specifies that a run realization numbered X
  . no additional suffix means that the run was made using interactive radiative
   fluxes and interactive surface fluxes
  . radhomo means that radiative heating rates are homogenized in the horizontal
  dimension at each time step
  . radagg means that the vertical profile of radiative heating is prescibed uniformly
  from the reference aggregated simulation at the same SST (Oref)
  . sfcagg means that surface fluxes are prescibed uniformly from the reference
  aggregated simulation at the same SST (Oref)
  . sfcdisagg means surface fluxes are prescibed uniformly from the reference
  disaggregated simulation at the same SST (Dref)
  . bXXX specifies that the run has been branched after day XXX from the simulation
  whose name is used as root 

- figures/ contains intermediate and final figures used in the publication

** Steps performed in the analysis

- Run scripts/computeRainStatistics.py to only derive rainfall statistics for a given simulation
- Run scripts/computeRainStatisticsConditionalsAndScaling.py to derive rain statistics, conditional statistics (composites on extreme percentiles) and the decomposition using O'Gorman (2009)'s scaling approximation
- Calculations of enhancement factors are performed directly in the notebooks when figures are generated
