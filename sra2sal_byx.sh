#!/bin/bash

########################################################################
#Script Name: sra2sal_byx.sh
#Description: linker script to run experiment based on job array index
#Args:        job array index
#Author:      Georgia Doing
#Date:        12-17-20
#Email:       Georgia.Doing.GR@Dartmouth.edu
########################################################################


readarray SRXS < <(find sra_comp -mindepth 2 -maxdepth 2 -type d)
./salPa_byx.sh ${SRXS[$1]}
