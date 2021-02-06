#!/bin/bash -l
####################################################################
#Script Name: salPa_byx.sh
#Description: Runs sra experiments through salmon
#Args:        dir containing txt file with all run numbers
#Author:      Georgia Doing
#Date:        12-17-20
#Email:       Georgia.Doing.GR@Dartmouth.edu  
####################################################################
#cd $PBS_O_WORKDIR

# set local scratch
#echo "/repository/user/main/public/root = \"/scratch/sra_out\"" >> $HOME/.ncbi/user-settings.mkfg
#echo "/repository/user/default-path = \"/scratch/sra_out\"" >> $HOME/.ncbi/user-settings.mkfg

# load modules
module load sratoolkit
module load salmon

# create temp dirs on local scratch
mkdir -p /scratch/$PBS_JOBID/salmon_out/$2
mkdir -p /scratch/$PBS_JOBID/sra_out/$1

# read in sra run accessions
readarray SRFILES < <(find $1 -type f -regex ".*txt")
echo SRFILES:
echo ${SRFILES}
readarray SRRS < ${SRFILES}
echo SRRS:
echo ${SRRS}

# download each sra run, one at a time
for f in $SRRS
do
    fasterq-dump $f -O /scratch/$PBS_JOBID/sra_out/$1 -t /scratch/$PBS_JOBID/sra_out/tmp -f
done

# read in fastq files and run together
#FQFILES=(/scratch/$PBS_JOBID/sra_out/$1/*)
readarray FQFILES < <(find /scratch/$PBS_JOBID/sra_out/$1 -type f -regex ".*fastq")
echo FQFILES:
echo ${FQFILES[*]}

for fq in $FQFILES
do
    wc -l $fq
done

salmon quant -i  ../t_indxs/$2 \
    --validateMappings \
    --writeUnmappedNames \
    --writeMappings=/scratch/$PBS_JOBID/salmon_out/$2/$1.salmon/mappings.sam \
    -l A \
    -r ${FQFILES[*]} \
    -o /scratch/$PBS_JOBID/salmon_out/$2/$1.salmon

# move from local scratch to dartfs scratch
mkdir -p /dartfs-hpc/scratch/f002bx6/salmon/$2/$1
cp -r /scratch/$PBS_JOBID/salmon_out/$2/$1.salmon /dartfs-hpc/scratch/f002bx6/salmon/$2/$1/

# cleanup
rm -r /scratch/$PBS_JOBID


exit

