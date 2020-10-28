#!/bin/bash


for expname in `cat ../expnames_test`; do

    runscript=submit_${expname}.sbatch
    cp submit_template.sbatch $runscript
    
    sed -i'' "s/SIMNAME/RCE_MPDATAxTKExCAMxSAM1MOM_4000x4000x15_256x256x64_$expname/" $runscript
    
   # sbatch $runscript

done
