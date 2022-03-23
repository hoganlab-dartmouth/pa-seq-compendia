## ---------------------------
##
## Script name: quant_collect.R
##
## Description: Gathers all quant.fs files into a csv
##
## Args: dir containing quant.fs files
##       genome fasta
##       output file name 
##
## Author: Georgia Doing
##
## Date Created: 2020-12-19
##
## Email: Georgia.Doing.GR@Dartmouth.edu
##
## ---------------------------
##
## Notes:
##
##
##  User must have write permissions
##
## ---------------------------
#library(Biostrings)

# to recieve commandline input
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)<3) {
  stop("Must supply dir, genome and out file", call.=FALSE)
}

## decend into dir to find all quant.fs files
sf_files <- list.files(args[1], pattern= '.sf', 
                       recursive = TRUE, full.names = TRUE)
#print(sf_files)
## parse file paths to get experiment accession
sf_names <- sapply(list.files(args[1], pattern= '.sf', 
                              recursive = TRUE, full.names = FALSE), function(x)
                                substring(x, regexpr('/[E,S]RX',x)+1, regexpr('/quant',x)-1))
# Read in the data from all quant files
sf_datasets <- lapply(sf_files, function(x) read.csv(x, sep = '\t',
                                                     stringsAsFactors = FALSE))
# combine into a dataframes of read numbers and TPM, merging by name
sf_all_reads <- Reduce(function(df1, df2) merge(df1, df2[,c(1,5)], by = 'Name'), sf_datasets)
rownames(sf_all_reads) <- sf_all_reads$Name
sf_all_reads <- sf_all_reads[,-c(2:4)]

sf_all_tpm <- Reduce(function(df1, df2) merge(df1, df2[,c(1,4)], by = 'Name'), sf_datasets)
rownames(sf_all_tpm) <- sf_all_tpm$Name
sf_all_tpm <- sf_all_tpm[,-c(2,3,5)]

# name columns as experiment accession
colnames(sf_all_reads) <- c('Name',sf_names)
colnames(sf_all_tpm) <- c('Name',sf_names)

# load in genome to convert t-index names to gene names
#genome_fasta <- readDNAStringSet(args[2])
anns_df <- read.csv(args[2], stringsAsFactors=F)

# convert t-index names to gene names
sf_all_reads[,1] <- as.character(sapply(rownames(sf_all_reads), function(x) anns_df$X2[anns_df$X1 == x]))
#colnames(sf_all_reads) <- c('X', colnames(sf_all_reads[-ncol(sf_all_reads)]))

sf_all_tpm[,1] <- as.character(sapply(rownames(sf_all_tpm), function(x) anns_df$X2[anns_df$X1 == x]))
#colnames(sf_all_tpm) <- c('X', colnames(sf_all_tpm[-ncol(sf_all_tpm)]))

# write out tables for num reads and TPM
write.csv(sf_all_reads, paste('num_reads_',args[3],sep=''))
write.csv(sf_all_tpm, paste('TPM_',args[3],sep=''))
