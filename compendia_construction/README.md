# Construction of compendia

This folder contains the bash, python, and R scripts used to construct the compendia.
To re-run these scripts, an [SRA run table](https://osf.io/pt7am/) is required.

Before using the PBS script, the headers will need to be altered to run through a specific user
account on Dartmouth's discovery cluster and/or dirs in bash scripts will need to be
restructured.

Relative paths are relative to this directory. 

## Installation

If these scripts are run locally they require the SRAtoolkit and salmon.
These software packages can be installed using miniconda.

All software necessary has been installed on the Dartmouth cluster discovery by admin.

## Usage

1. Upload a run table from the SRA to the working dir (e.g. use `scp` or `rsync` to upload the file if working on a remote computer).
2. Use the R script `run_table/run_table_dirs.R` to make the required directories for compendia processing. This script requires a run table from SRA.
3. Create a *P. aeruginosa* transcriptome index.
```
salmon index -t ./t_indxs/pao1_cnda.fa.gz -i ./t_indxs/pao1_cdna_k15 -k 15
```

4. Process samples in the dir against the transcriptome index.
```
mksub sra_to_salmon.pbs
```
or
```
./shell_scripts/sra2sal_byx.sh 1 # will fun first experiment found in sra_comp/
```
Iterate between the following scripts until all samples have sucessfully run
```
./pbs_scripts/collect_jobIDs.pbs
./pbs_scripts/collected_IDs_jobs.pbs
```
5. Collect salmon output and log data using the following scripts:
```
shape_comp/gene_names.R
shape_comp/quant_collect.R
shape_comp/logs_collect.py
```

6. Download the processed data (e.g. using `scp` or `rsync`).
