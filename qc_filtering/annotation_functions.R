## ---------------------------
##
## Script name: annotation_functions.R
##
## Description:
##
## Args:
##
## Author: Georgia Doing
##
## Date Created: 2021-09-28
##
## Email: Georgia.Doing.GR@Dartmouth.edu
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

# load in annotation files 
# (all originally sources from www.pseudomonas.com)


PAO1_orth <- read.csv('~/Dropbox (Hogan Lab)/Resources/Annotations/Pseudomonas_aeruginosa_PAO1_107_orthologs.csv', stringsAsFactors = FALSE)
#PAO1_orth <- tryCatch(read.csv('Pseudomonas_aeruginosa_PAO1_107_orthologs.csv', stringsAsFactors = F),
#                       .e = function(e){
#                         download.file('https://pseudomonas.com/downloads/pseudomonas/pgd_r_20_2/Pseudomonas_aeruginosa_PAO1_107/Pseudomonas_aeruginosa_PAO1_107_orthologs.csv.gz',
#                                       'PAO1_orthologs.csv.gz')
#                         read.csv('PAO1_orthologs.csv.gz', stringsAsFactors = F)
#                       }, finally = print("PAO1_orth")
#)

PA14_orth <- read.csv('~/Dropbox (Hogan Lab)/Resources/Annotations/Pseudomonas_aeruginosa_UCBPP-PA14_109_orthologs.csv', stringsAsFactors = FALSE)
#PA14_orth <- tryCatch(read.csv('Pseudomonas_aeruginosa_UCBPP-PA14_109_orthologs.csv', stringsAsFactors = F),
#                       error = function(e){
#                         download.file('https://pseudomonas.com/downloads/pseudomonas/pgd_r_20_2/Pseudomonas_aeruginosa_UCBPP-PA14_109/Pseudomonas_aeruginosa_UCBPP-PA14_109_orthologs.csv.gz',
#                                       'PA14_orthologs.csv.gz')
#                         read.csv('PA14_orthologs.csv.gz', stringsAsFactors = F)
#                       },finally = print("PAO1_orth")
#)

PA14_ann <- read.csv('~/Dropbox (Hogan Lab)/Resources/Annotations/Pseudomonas_aeruginosa_UCBPP-PA14_109.csv', stringsAsFactors = FALSE, skip = 2)
#PA14_ann <- tryCatch(read.csv('Pseudomonas_aeruginosa_UCBPP-PA14_109.csv', stringsAsFactors = F),
#                      error = function(e){
#                        download.file('https://pseudomonas.com/downloads/pseudomonas/pgd_r_20_2/Pseudomonas_aeruginosa_UCBPP-PA14_109/Pseudomonas_aeruginosa_UCBPP-PA14_109.csv.gz',
#                                      'PA14_annotations.csv.gz')
#                        read.csv('PA14_annotations.csv.gz', stringsAsFactors = F)
#                      },finally = print("PA14_ann")
#)

PAO1_ann <- read.csv('~/Dropbox (Hogan Lab)/Resources/Annotations/Pseudomonas_aeruginosa_PAO1_107.csv', stringsAsFactors = FALSE, skip = 2)
#PAO1_ann <- tryCatch(read.csv('Pseudomonas_aeruginosa_PAO1_107.csv', stringsAsFactors = F),
#                     error = function(e){
#                       download.file('https://pseudomonas.com/downloads/pseudomonas/pgd_r_20_2/Pseudomonas_aeruginosa_PAO1_107/Pseudomonas_aeruginosa_PAO1_107.csv.gz',
#                                     'PAO1_annotations.csv.gz')
#                       read.csv('PAO1_annotations.csv.gz', stringsAsFactors = F)
#                     },finally = print("PAO1_ann")
#)

PAO1_cdna <- read.csv('~/Dropbox (Hogan Lab)/Resources/Annotations/pao1_asm676v1cdna_gene_names.csv', stringsAsFactors = F)
#PAO1_cdna <- tryCatch(read.csv('pao1_asm676v1cdna_gene_names.csv', stringsAsFactors = F),
#                      error = function(e){
#                        print("need cDNA to PA num mapping")
#                        #read.csv('pao1_asm676v1cdna_gene_names.csv', stringsAsFactors = F)
#                      },finally = print("PAO1_cdna")
#)

PA14_cdna <- read.csv('~/Dropbox (Hogan Lab)/Resources/Annotations/pa14_cdna_gene_names.csv', stringsAsFactors = F)
#PA14_cdna <- tryCatch(read.csv('pa14_cdna_gene_names.csv', stringsAsFactors = F),
#                      error = function(e){
#                        print("need cDNA to PA num mapping")
#                        #read.csv('pa14_cdna_gene_names.csv', stringsAsFactors = F)
#                      }, finally = print("PA14_cdna")
#)

# functions

#' Convert from 3 or 4 letter name to PAO1 number syntax
#' 
#' @param x A 3 or 4 letter gene name
#' @return The corresponding PAO1 number
#' @example 
#' name_to_PAO1('dnaA')
name_to_PAO1 <- function(x){
  if(x==''){
    return(x)
  }
  out_num <- PAO1_ann$Locus.Tag[PAO1_ann$Name == x]
  if(identical(out_num, character(0))){
    return(x)
  }
  return(substring(out_num[1],1,6))
}


PAO1_to_name <- function(x){
  if(x==''){
    return(x)
  }
  out_num <- PAO1_ann$Name[substring(PAO1_ann$Locus.Tag,1,6) == x]
  if(identical(out_num, character(0))){
    return(x)
  }
  if(out_num[1] == ''){
    return(x)
  }
  if(grepl(' ',out_num[1])){
    return(x)
  }
  return(out_num[1])
}

name_to_PA14 <- function(x){
  if(x==''){
    return(x)
  }
  out_num <- PA14_ann$Locus.Tag[PA14_ann$Name == x]
  if(identical(out_num, character(0))){
    return(x)
  }
  return(substring(out_num[1],1,10))
}

PA14_to_name <- function(x){
  if(x==''){
    return(x)
  }
  out_num <- PA14_ann$Name[substring(PA14_ann$Locus.Tag,1,10) == x]
  if(identical(out_num, character(0))){
    return(x)
  }
  if(out_num == ''){
    return(x)
  }
  return(out_num)
}

PA14_to_PAO1 <- function(x){
  if(x==''){
    return(x)
  }
  out_num <- PA14_orth$Locus.Tag..Hit.[PA14_orth$Locus.Tag..Query. == x & PA14_orth$Strain..Hit. == 'Pseudomonas aeruginosa PAO1 (Reference)']
  if(identical(out_num, character(0))){
    return(x)
  }
  return(out_num[1])
}

PAO1_to_PA14 <- function(x){
  if(x==''){
    return(x)
  }
  out_num <- PAO1_orth$Locus.Tag..Hit.[PAO1_orth$Locus.Tag..Query. == x & PAO1_orth$Strain..Hit. == 'Pseudomonas aeruginosa UCBPP-PA14']
  if(identical(out_num, character(0))){
    return(x)
  }
  return(out_num[1])
} 

PAO1_to_cDNA <- function(x){
  if(x==''){
    return(x)
  }
  out_num <- PAO1_cdna$X1[PAO1_cdna$X2 == x]
  if(identical(out_num, character(0))){
    return(x)
  }
  return(out_num[1])
}

PA14_to_cDNA <- function(x){
  if(x==''){
    return(x)
  }
  out_num <- PA14_cdna$X1[PA14_cdna$X2 == x]
  if(identical(out_num, character(0))){
    return(x)
  }
  return(out_num[1])
}

cDNA_to_PAO1 <- function(x){
  if(x==''){
    return(x)
  }
  out_num <- PAO1_cdna$X2[PAO1_cdna$X1 == x]
  if(identical(out_num, character(0))){
    return(x)
  }
  return(out_num[1])
}

cDNA_to_PA14 <- function(x){
  if(x==''){
    return(x)
  }
  out_num <- PA14_cdna$X2[PA14_cdna$X1 == x]
  if(identical(out_num, character(0))){
    return(x)
  }
  return(out_num[1])
}