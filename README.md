# Computationally efficient assembly of a *Pseudomonas aeruginosa* gene expression compendium

This repository contains scripts used to build compendia of publicly available RNAseq data for *Pseudomonas aeruginosa*.

The construction and composition of these compendia is documented in [Doing et al. 2022](https://doi.org/10.1101/2022.01.24.477642).

Required input files and the main output files like [raw](https://osf.io/mn7tj/) and [normalized compendia](https://osf.io/vz42h/) are stored on OSF: https://osf.io/s9gyu/

The compendia were originally constructed using a collection of bash, python, and R scripts which are documented in `compendia_construction`. 
This pipeline required a pre-downloaded [SRA run table](https://osf.io/pt7am/) which was used to download FASTQ files of RNAseq data using `fasterq-dump`.
These files were then mapped against *P. aeruginosa* transcriptomes using salmon, and counts were post-processed in R and python.

To facilitate the addition of new RNAseq samples to the compendia, we have provided an automated pipeline that will process new SRA experiment accessions, add the new counts to the original raw compendia, and re-filter and re-normalize the compendia. 
This pipeline is provided in `auto_add_to_compendia`.
This pipeline automates the bash, R, and python scripts provided in `compendia_construction` with minor modifications to the original code. 

## Support

Please submit an issue to the issue tracker. 

## Roadmap

This is a proof-of-principle project that could help set the stage for a versatile, scaleable pipeline.

## Authors

Georgia Doing (@georgiadoing), Taylor Reiter (@taylorreiter)

## License

See the LICENSE file above.

## Status

Complete.

## Citation

Please cite

Computationally efficient assembly of a Pseudomonas aeruginosa gene expression compendium
Georgia Doing, Alexandra J. Lee, Samuel L. Neff, Jacob D. Holt, Bruce A. Stanton, Casey S. Greene, Deborah A. Hogan
bioRxiv 2022.01.24.477642; doi: https://doi.org/10.1101/2022.01.24.477642 
