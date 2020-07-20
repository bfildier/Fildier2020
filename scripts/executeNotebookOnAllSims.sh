#!/bin/bash

# WARNING: the notebook used must have specific features:
# - It must include the following lines:
# $import matplotlib
# $# matplotlib.use("PDF")
# It must include an option 'fromsimname' that defines options from the simulation name
# - All figure instances, if displayed, must include a plt.show() statement, and otherwise a plt.close() statement


runmode="regular"
#runmode="debug"

##-- Set working directory --##
#SCRIPTDIR=`python -c "import os; print(os.getcwd())"`
SCRIPTDIR=/global/cscratch1/sd/bfildier/dataAnalysis/aggregation-and-extreme-rain/scripts
if [ "$runmode" == "regular" ]; then
    SCRIPTDIR_AUTO=${SCRIPTDIR}/automated_analysis
elif [ "$runmode" == "debug" ]; then
    SCRIPTDIR_AUTO=${SCRIPTDIR}/automated_analysis/debug
fi
[ -d ${SCRIPTDIR_AUTO} ] || mkdir ${SCRIPTDIR_AUTO}
[ -d ${SCRIPTDIR_AUTO}/scripts ] || mkdir ${SCRIPTDIR_AUTO}/scripts
[ -d ${SCRIPTDIR_AUTO}/logs ] || mkdir ${SCRIPTDIR_AUTO}/logs

# Simulations to consider

SST=308
#: '
simnames="RCE_MPDATAxTKExCAMxSAM1MOM_4000x4000x15_256x256x64_TKE-SST${SST}-r1-b150-sfcagg
RCE_MPDATAxTKExCAMxSAM1MOM_4000x4000x15_256x256x64_TKE-SST${SST}-r1-b150-sfcdisagg
RCE_MPDATAxTKExCAMxSAM1MOM_4000x4000x15_256x256x64_TKE-SST${SST}-radhomo-r1
RCE_MPDATAxTKExCAMxSAM1MOM_4000x4000x15_256x256x64_TKE-SST${SST}-radhomo-r1-b100-radagg
"
#'
: '
simnames="RCE_MPDATAxTKExCAMxSAM1MOM_4000x4000x15_256x256x64_TKE-SST${SST}-r1
RCE_MPDATAxTKExCAMxSAM1MOM_4000x4000x15_256x256x64_TKE-SST${SST}-r1-b150-sfcagg
RCE_MPDATAxTKExCAMxSAM1MOM_4000x4000x15_256x256x64_TKE-SST${SST}-r1-b150-sfcdisagg
RCE_MPDATAxTKExCAMxSAM1MOM_4000x4000x15_256x256x64_TKE-SST${SST}-radhomo-r1
RCE_MPDATAxTKExCAMxSAM1MOM_4000x4000x15_256x256x64_TKE-SST${SST}-radhomo-r1-b100-radagg
"
#'
: '
simnames="RCE_MPDATAxTKExCAMxSAM1MOM_4000x4000x15_256x256x64_TKE-SST${SST}-r1
RCE_MPDATAxTKExCAMxSAM1MOM_4000x4000x15_256x256x64_TKE-SST${SST}-r1-b150-sfcagg
RCE_MPDATAxTKExCAMxSAM1MOM_4000x4000x15_256x256x64_TKE-SST${SST}-r1-b150-sfcdisagg
RCE_MPDATAxTKExCAMxSAM1MOM_4000x4000x15_256x256x64_TKE-SST${SST}-radhomo-r1
RCE_MPDATAxTKExCAMxSAM1MOM_4000x4000x15_256x256x64_TKE-SST${SST}-radhomo-r1-b100-radagg
"
#'
: '
simnames="RCE_MPDATAxTKExCAMxSAM1MOM_4000x4000x15_256x256x64_TKE-SST302-r1-b150-sfcagg
RCE_MPDATAxTKExCAMxSAM1MOM_4000x4000x15_256x256x64_TKE-SST302-r1-b150-sfcdisagg
RCE_MPDATAxTKExCAMxSAM1MOM_4000x4000x15_256x256x64_TKE-SST302-radhomo-r1
RCE_MPDATAxTKExCAMxSAM1MOM_4000x4000x15_256x256x64_TKE-SST302-radhomo-r1-b100-radagg
"
#'
: '
simnames="RCE_MPDATAxTKExCAMxSAM1MOM_4000x4000x15_256x256x64_TKE-SST304-radhomo-sfchomo-r1
RCE_MPDATAxTKExCAMxSAM1MOM_4000x4000x15_256x256x64_TKE-SST302-radhomo-sfchomo-r1
RCE_MPDATAxTKExCAMxSAM1MOM_4000x4000x15_256x256x64_TKE-SST300-radhomo-sfchomo-r1
RCE_MPDATAxTKExCAMxSAM1MOM_4000x4000x15_256x256x64_TKE-SST300-r1-b150-sfchomo
RCE_MPDATAxTKExCAMxSAM1MOM_4000x4000x15_256x256x64_TKE-SST302-r1-b150-sfchomo
RCE_MPDATAxTKExCAMxSAM1MOM_4000x4000x15_256x256x64_TKE-SST304-r1-b150-sfchomo
"
#'
: '
simnames="RCE_MPDATAxTKExCAMxSAM1MOM_4000x4000x15_256x256x64_TKE-SST${SST}-r1
RCE_MPDATAxTKExCAMxSAM1MOM_4000x4000x15_256x256x64_TKE-SST${SST}-radhomo-r1"
#'
#simnames="RCE_MPDATAxTKExCAMxSAM1MOM_4000x4000x15_256x256x64_TKE-SST306-r1-b150-sfcagg
#RCE_MPDATAxTKExCAMxSAM1MOM_4000x4000x15_256x256x64_TKE-SST306-r1-b150-sfcdisagg"
#simnames="RCE_MPDATAxTKExCAMxSAM1MOM_4000x4000x15_256x256x64_TKE-SST302-radhomo-r1-b100-radagg"
#simnames="RCE_MPDATAxTKExCAMxSAM1MOM_4000x4000x15_256x256x64_TKE-SST304-radhomo-r1-b100-radagg"
#simnames="RCE_MPDATAxTKExCAMxSAM1MOM_4000x4000x15_256x256x64_TKE-SST306-r1"

##-- Batch script options --##

notebook=${1%.ipynb} # Remove .ipynb extension

##-- Run options--##
todaysdate=`date +"%Y%m%d-%H%M"`
template_nameroot=${notebook}_$todaysdate
template_analysis_script=${template_nameroot}.py
template_batch_script=sbatch_${template_nameroot}.sbatch

##-- Main --##

# Convert notebook to template python script 
./ipynb2py.sh ${notebook} ${SCRIPTDIR_AUTO}/scripts/${template_nameroot}
# Allow script to define simulation attributes by specifying simulation name
sed -i'' "s/^fromsimname.*/fromsimname = True/" ${SCRIPTDIR_AUTO}/scripts/${template_analysis_script}
# Change simulation name
sed -i'' "s/^simname.*/simname = ${simname}/" ${SCRIPTDIR_AUTO}/scripts/${template_analysis_script}

# Copy template batch script
cp sbatch_executeNotebookOnAllSims_template.sbatch ${SCRIPTDIR_AUTO}/scripts/${template_batch_script}

cd ${SCRIPTDIR_AUTO}/scripts/

# Change run option in batch script
if [ "$runmode" == "debug" ];
then
    sed -i'' 's/^#SBATCH --partition=.*/#SBATCH --partition=debug/' ${template_batch_script}
    sed -i'' 's/^#SBATCH --qos=/##SBATCH --qos=/' ${template_batch_script}
    sed -i'' 's/^#SBATCH --time=.*/#SBATCH --time=00:30:00/' ${template_batch_script}
    sed -i'' 's/^#SBATCH --mail-type=.*/#SBATCH --mail-type=ALL/' ${template_batch_script}
elif [ "$runmode" == "regular" ];
then
    sed -i'' 's/^#SBATCH --partition=.*/#SBATCH --partition=regular/' ${template_batch_script}
    sed -i'' 's/.*#SBATCH --qos=.*/#SBATCH --qos=premium/' ${template_batch_script}
    sed -i'' 's/^#SBATCH --time=.*/#SBATCH --time=05:00:00/' ${template_batch_script}
    sed -i'' 's/^#SBATCH --mail-type=.*/#SBATCH --mail-type=FAIL,TIME_LIMIT_90/' ${template_batch_script}
fi

# Set up and run batch script for each simulation name
for simname in `echo $simnames`; do

##-- Create scripts --##

    nameroot=${template_nameroot}_${simname}
    analysisscript=${nameroot}.py
    batchscript=sbatch_${nameroot}.sbatch

    # Duplicate analysis script
    cp ${template_analysis_script} ${analysisscript}
    # Duplicate batch script
    cp ${template_batch_script} $batchscript

##-- Analysis script options --#

    sed -i'' "s/simname =.*/simname = \"${simname}\"/" ${analysisscript}
    sed -i'' "s/fromsimname =.*/fromsimname = True/" ${analysisscript}

##-- Batch script options --#

    sed -i'' "s|SCRIPTDIR|${SCRIPTDIR_AUTO}/scripts|g" $batchscript
    sed -i'' "s/SCRIPTNAME/${nameroot}/g" $batchscript

##-- Launch batch script--##

    echo "$simname : "
    sbatch $batchscript

done

exit 0
