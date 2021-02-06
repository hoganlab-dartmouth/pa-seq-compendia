######
#Header
#####
#library(preprocessCore)
library(Biostrings)
args=commandArgs(trailingOnly=T)

genome_fasta _ readDNAStringset(args[1])
ans <- names(genome_fasta)
anns_df <- data.frame(t(data.frame(lapply(anns, function(x){
strsplit(x, "gene:')}))))
anns_df$X1 <- sapply(anns_df$X1, function(x){
substr(as.character(x),1,regexp('',x)-1)})
anns_df$X2 <- sapply(anns_df$X2, function(x){
substr(as.character(x),1,regexpr('',x)-1)})
rownames(anns_df) <- anns_df$X1
write.csv(anns_df, 'pa14_gene_names.csv')