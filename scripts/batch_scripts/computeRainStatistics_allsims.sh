#!/bin/bash

here=$(dirname $0)
scriptsdir=${here}/..

for expname in `cat ${scriptsdir}/expnames`; do

    runscript=${here}/temp/submit_RS_${expname}.sbatch
    cp ${here}/submit_RS_template.sbatch $runscript
    
    sed -i'' "s/SIMNAME/RCE_MPDATAxTKExCAMxSAM1MOM_4000x4000x15_256x256x64_$expname/" $runscript
    
    sbatch $runscript

done
