#!/bin/bash

########################################################################
#Script Name: collect_jobIDs.sh
#Description: collects jobID of failed sra2sal runs for rerun
#Args:        t_inx for folder
#Author:      Georgia Doing
#Date:        2021-02-14
#Email:       Georgia.Doing.GR@Dartmouth.edu
########################################################################

cd $PBS_O_WORKDIR
readarray -t SRXS < <(find sra_comp -mindepth 2 -maxdepth 2 -type d)

IDS="0"

for ((IDX = 1; IDX <= ${#SRXS[@]} ; IDX++))
do
    echo /dartfs-hpc/scratch/f002bx6/salmon/$1/${SRXS[$IDX]}/*.salmon
    if ! ls /dartfs-hpc/scratch/f002bx6/salmon/$1/${SRXS[$IDX]}/*.salmon ; then
	IDS="$IDS,$IDX"
    fi
done

if $IDS = "" ; then
    echo 'all salmon experiments exists'
    exit 200
fi

IDS="$IDS%50"
echo $IDS

cat > collected_IDS_jobs_$1.pbs << EOF1
#!/bin/bash
#######################################################################
# Script Name:
# Description:
#Args:
#Author: Georgia Doing
#Date:
#Email: Georgia.Doing.GR@Dartmouth.edu
#######################################################################

#PBS -q largeq
#PBS -A NCCC
#PBS -N s2s14
#PBS -l walltime=02:00:00
#PBS -l nodes=1:ppn=6
#PBS -m bea
#PBS -M Georgia.Doing.GR@Dartmouth.edu
#PBS -t $IDS
#PBS -j oe
#PBS -e /dartfs-hpc/scratch/f002bx6/salmon/sra_comp/oe
#PBS -o /dartfs-hpc/scratch/f002bx6/salmon/sra_comp/oe
#PBS -l feature='cellk'

cd \$PBS_O_WORKDIR

less \$PBS_NODEFILE
mkdir -p /scratch/f002bx6
mkdir /scratch/f002bx6/\$PBS_JOBID

if ../shell_scripts/sra2sal_byx.sh \${PBS_ARRAYID} $1 ; then
    echo \$PBS_JOBID quant succeeded
    du -h --max-depth 1 /scratch/f002bx6
else
    rc=$?
    echo \$PBS_JOBID quant failed
    du -h --max-depth 1 /scratch/f002bx6
    rm -r /scratch/f002bx6/\$PBS_JOBID
    du -h --max-depth 1 /scratch/f002bx6
    exit $rc
    
fi

exit
EOF1

exit

