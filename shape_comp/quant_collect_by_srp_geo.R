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

system('export NCBI_API_KEY=ee22ceb4dead070e01f20a5e8cf948258908')
# to recieve commandline input
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)<3) {
  stop("Must supply dir, genome and out file", call.=FALSE)
}

## get directories for each experiment
srp_dirs <- list.files(args[1], full.names=TRUE)
#srp_dirs

srp_tables <- lapply(srp_dirs, function(x){
#print('x: ')
print(x)	   
## decend into dir to find all quant.fs files
	   sf_files <- list.files(x, pattern= '.sf', 
                       recursive = TRUE, full.names = TRUE)
print(length(sf_files))
## parse file paths to get experiment accession
   	   sf_names <- sapply(list.files(x, pattern= '.sf', 
                              recursive = TRUE, full.names = FALSE), function(x)
                                substring(x, regexpr('/[E,S]RX',x)+1, regexpr('/quant',x)-1))

gsm_names <- sapply(sf_names, function(x){
print('x: ')
print(x)
	  sra_to_gsm <- paste0(paste0('esearch -db gds -query "',sub('.salmon','',x)),'[ACCN" | efetch -format docsum | xtract -pattern DocumentSummary -element Accession')
	  system(sra_to_gsm, intern=TRUE)
})
#print(length(sf_names))
print(gsm_names)
sf_names <- paste(gsm_names, sf_names, sep='_')
#print(length(sf_names))
print(sf_names)
# Read in the data from all quant files
       	  sf_datasets <- lapply(sf_files, function(x) read.csv(x, sep = '\t',
                                                     stringsAsFactors = FALSE))
#print(head(sf_datasets[[1]]))
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
  	write.csv(sf_all_reads, paste0('',paste(x,paste('num_reads_',args[3],sep=''),sep='_')))
	write.csv(sf_all_tpm,paste0('',  paste(x,paste('TPM_',args[3],sep=''),sep='_')))
}
)
