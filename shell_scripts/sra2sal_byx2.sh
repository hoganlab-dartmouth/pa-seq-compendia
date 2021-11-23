#!/bin/bash

########################################################################
#Script Name: sra2sal_byx.sh
#Description: linker script to run experiment based on job array index
#Args:        job array index
#Author:      Georgia Doing
#Date:        12-17-20
#Email:       Georgia.Doing.GR@Dartmouth.edu
########################################################################

cd $SLURM_SUBMIT_DIR
readarray SRXS < <(find $3 -mindepth 2 -maxdepth 2 -type d)

echo $3
echo $1
echo ${#SRXS[@]}

if  (( $1 >  ${#SRXS[@]} )) ; then
    echo 'no accession'
    exit 203
fi

../shell_scripts/salPa_byx.sh ${SRXS[$1]} $2

exit

