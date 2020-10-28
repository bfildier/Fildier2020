#!/bin/bash

conda activate pyLMD

for expname in `cat expnames`; do

    simname = RCE_MPDATAxTKExCAMxSAM1MOM_4000x4000x15_256x256x64_$expname
    
    echo calculate rain statistics for $simname
    python computeRainStatistics.py $simname

done
