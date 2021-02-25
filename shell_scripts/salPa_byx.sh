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
conda activate salmon



# check if this accession-index pair has already been processed
if du -hac -d 1 /dartfs-hpc/scratch/f002bx6/salmon/$2/$1/*.salmon ; then
    echo "/salmon/$2/$1 already processed"
    DT=$(date)
    echo "Attempt $PBS_JOBID $DT" >> /dartfs-hpc/scratch/f002bx6/salmon/$2/$1/processes.log
    rm -r /scratch/f002bx6/$PBS_JOBID
    exit 200
fi


mkdir -p /dartfs-hpc/scratch/f002bx6/salmon/$2/$1
DT=$(date)
echo "Processed $PBS_JOBID $DT" > /dartfs-hpc/scratch/f002bx6/salmon/$2/$1/processes.log

#echo ls /dartfs-hpc/scratch/f002bx6/salmon/$2/$1/
# create temp dirs on local scratch
mkdir -p /scratch/f002bx6/$PBS_JOBID/salmon_out/$2
mkdir -p /scratch/f002bx6/$PBS_JOBID/sra_out/$1

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
    if ! fasterq-dump $f -O /scratch/f002bx6/$PBS_JOBID/sra_out/$1 -t /scratch/f002bx6/$PBS_JOBID/sra_out/tmp -f ; then
        echo "fasterq-dump failed on $f"
        exit 201
    fi
done

# read in fastq files and run together
#FQFILES=(/scratch/$PBS_JOBID/sra_out/$1/*)
readarray FQFILES < <(find /scratch/f002bx6/$PBS_JOBID/sra_out/$1 -type f -regex ".*fastq")
echo FQFILES:
echo ${FQFILES[*]}

for fq in $FQFILES
do
    wc -l $fq
done

if ! \
salmon quant -i  ../t_indxs/$2 \
    --softclip \
    --softclipOverhangs \
    --minScoreFraction 0.65 \
    --fldMean 51 \
    --seqBias \
    -l A \
    -r ${FQFILES[*]} \
    -o /scratch/f002bx6/$PBS_JOBID/salmon_out/$2/$1.salmon
then
    echo 'salmon failed'
    exit 202
fi
# move from local scratch to dartfs scratch
#mkdir -p /dartfs-hpc/scratch/f002bx6/salmon/$2/$1
cp -r /scratch/f002bx6/$PBS_JOBID/salmon_out/$2/$1.salmon /dartfs-hpc/scratch/f002bx6/salmon/$2/$1/

# cleanup
du -h --max-depth 0 /scratch/f002bx6/$PBS_JOBID
rm -r /scratch/f002bx6/$PBS_JOBID
[ ! -e /scratch/f002bx6/$PBS_JOBID ] && echo "/scratch/f002bx6/$PBS_JOBID deleted"


exit

