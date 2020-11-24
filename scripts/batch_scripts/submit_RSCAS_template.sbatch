#!/bin/bash
#
#SBATCH --mail-user=benjamin.fildier@lmd.ens.fr
#SBATCH --mail-type=ALL
#SBATCH --job-name=compute_stats
#SBATCH --output=logs/compute_stats_SIMNAME_%jobid.out
#SBATCH --error=logs/compute_stats_SIMNAME_%jobid.err
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
##SBATCH --mem-per-cpu=100
#SBATCH --constraint=knl,quad,cache

scriptdir=/global/cscratch1/sd/bfildier/dataAnalysis/Fildier2020/scripts
simname=SIMNAME

module load python/3.6-anaconda-4.4
conda activate pyPhD

srun python ${scriptdir}/computeRainStatisticsConditionalsAndScaling.py -s $simname
