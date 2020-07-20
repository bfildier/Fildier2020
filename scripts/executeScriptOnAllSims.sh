#!/bin/bash

machine=coriknl

# Define UTILDIR and OUTPUTDIR
. /global/homes/b/bfildier/code/SAM_scripts/load_dirnames.sh ${machine}

#scriptname=drawFirstPlots
scriptname=computeAndDrawSpatialStats

for simname in `ls ${ARCHIVEDIR}/${machine}`; do
    simroot=${simname%_*}
    expname=${simname##*_}
    echo "------------------------------------------------------------------------"
    echo "    $simroot $expname"
    echo "------------------------------------------------------------------------"
    for keyword in simroot expname; do
        sed -i "s/${keyword}=.*/${keyword}=${!keyword}/" ${scriptname}.sh
    done
    source ${scriptname}.sh
done
