#!/bin/bash
#
#SBATCH --mail-user=benjamin.fildier@lmd.ens.fr
#SBATCH --mail-type=ALL
#SBATCH --job-name=compute_stats
#SBATCH --output=compute_stats%jodid.out
#SBATCH --error=compute_stats%jodid.err
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=10:00
##SBATCH --mem-per-cpu=100
#SBATCH --constraint=knl,quad,cache

simname=SIMNAME

module load python/3.6-anaconda-4.4
conda activate pyPhD

srun python ../../computeRainStatistics.py -s $simname
