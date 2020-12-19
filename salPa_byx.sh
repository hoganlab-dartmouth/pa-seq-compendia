#!/bin/bash -l

####################################################################
#Script Name: salPa_byx.sh
#Description: Runs sra experiments through salmon
#Args:        dir containing txt file with all run numbers
#Author:      Georgia Doing
#Date:        12-17-20
#Email:       Georgia.Doing.GR@Dartmouth.edu  
####################################################################

# set local scratch
echo "/repository/user/main/public/root = \"/scratch/sra_out\"" >> $HOME/.ncbi/user-settings.mkfg
echo "/repository/user/default-path = \"/scratch/sra_out\"" >> $HOME/.ncbi/user-settings.mkfg

# load modules
module load sratoolkit
module load salmon

# create temp dirs on local scratch
mkdir /scratch/sra_out
mkdir /scratch/salmon_out
mkdir -p /scratch/sra_out/$1
echo $1 # record sra experiment accession

# read in sra run accessions
readarray SRFILES < <(find $1 -type f -regex ".*txt")
echo ${SRFILES}
readarray SRRS < ${SRFILES}
echo ${SRRS}

# download each sra run, one at a time
for f in $SRRS
do
    fasterq-dump $f -O /scratch/sra_out/$1 -t /scratch/sra_out/tmp -f
done

# read in fastq files and run together
FQFILES=/scratch/sra_out/$1/*
salmon quant -i  /dartfs-hpc/scratch/f002bx6/pao1_index\
    -l A \
    -r $FQFILES \
    -o /scratch/salmon_out/$1.salmon 


# move from local scratch to dartfs scratch
mkdir -p /dartfs-hpc/scratch/f002bx6/salmon/$1
cp /scratch/salmon_out/$1.salmon/* /dartfs-hpc/scratch/f002bx6/salmon/$1/

# cleanup
rm -r /scratch/sra_out/$1*
rm -r /scratch/salmon_out/$1.salmon/*
