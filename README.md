# pa-seq-compendia

To build compendia of publicly available RNAseq data for *P. aeruginosa*, sepcifically.

This short pipeline (collection of bash scripts) is written to download *P. aeruginosa*
RNAseq data from the NBCI SRA database and process them with salmon. It works off-of
a pre-downloaded SRA run table and pre-determined salmon parameters.

Before using the PBS script, the headers will need to be altered to run through a specific user
account on Dartmouth's disocvery cluster and/or dirs in bash scripts will need to be
restructured.

## Installation

If these scripts are run locally they require the SRAtoolkit and salmon.

All software necesary has been installed on dscovery by admin.

## Usage

1. Upload a run table from sra to the working dir.

eg. scp

2. Assemble a directory structure.

/run_table/run_table_dirs.R

3. Create a *P. aeruginosa* transcriptome index.

salmon index -t ../t_indxs/pao1_cnda.fa.gz -i ../t_indxs/pao1_cdna_k15 -k 15


4. Process samples in the dir against the transcriptome index.

mksub sra_to_salmon.pbs

or

./sra2sal_byx.sh 1 # will fun first experiment found in sra_comp/

Iterate between the following scripts until all samples have sucessfully run

/collect_jobIDs.pbs
/collected_IDs_jobs.pbs

Check that all samples have been processed
find comp/ -mexdepth 2 -mindepth 2 | cat | wc -l
find /dartfs-hpc/scratch/f002bx6/salmon/pa14_cdna/comp -name "*.sf" | wc -l


5. Collect salmon output and log data

/shape_comp/gene_names.R
/shape_comp/quant_collect.R
/shape_comp/logs_collect.py


6. Download the processed data
eg. scp


## Support

Georgia.Doing.GR@Dartmouth.edu

## Roadmap

This is a proof-of-principle project that could help set the stage for
a versatile, scaleable pipeline.

## Contributing

## Authors

Georgia Doing (@georgiadoing)
## License

tba

## Status

development
