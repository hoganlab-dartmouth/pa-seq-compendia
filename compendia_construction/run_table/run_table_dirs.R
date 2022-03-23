## ---------------------------
##
## Script name: run_table_dirs
##
## Description: Uses a run table from SRA to make dirs for compendia processing
##
## Args: SraRunTable.csv
##
## Author: Georgia Doing
##
## Date Created: 2020-12-18
##
## Email: Georgia.Doing.GR@Dartmouth.edu
##
## ---------------------------
##
## Notes:
##  SraRunTable.csv must have columns:
##      SRA_study - accession numbers SRP... or ERP...
##      Experiment - accession numbers SRX... or ERX...
##      Run - accession numbers SRR... or ERR...
##
##  User must have write permissions
##
## ---------------------------
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied SraRunTable.csv", call.=FALSE)
}
# load in SraRunTable from the first argument
run_table <- read.csv(args[1], stringsAsFactors = F)

print(colnames(run_table))
# Using the run table make dirs nested within sra_comp/
# of projects and experiments
# where leaf dirs contain txt files of run accessions
lapply(unique(run_table$SRA.Study), function(study){
  print(study)
  sdir <- paste('./sra_comp/',study, sep='')
  dir.create(sdir)
  exps <- unique(run_table$Experiment[run_table$SRA.Study == study])
  lapply(exps, function(exp){
    edir <- paste(sdir, exp, sep='/')
    dir.create(edir)
    rns <- unique(run_table$Run[run_table$Experiment == exp])
    erns <- file(paste(paste(edir,exp,sep='/'), '.txt',sep=''))
    writeLines(c(rns,''), erns)
    close(erns)
  })
})
