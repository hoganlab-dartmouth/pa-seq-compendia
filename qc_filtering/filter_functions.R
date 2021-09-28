## ---------------------------
##
## Script name: filter_functions.R
##
## Description:
##
## Args:
##
## Author: Georgia Doing
##
## Date Created: 2021-04-27
##
## Email: Georgia.Doing.GR@Dartmouth.edu
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------
library(RVenn)
source('~/Dropbox (Hogan Lab)/Resources/Scripts/fsqn.R')

hlog <- function(df, base=exp(1)){
  # df: dataframe of gene expression to be logged
  # base: base for log calculation
  # returns: log of df with happy 0s for 0 or neg values
  dfl <- data.frame(log(data.matrix(df),base))
  dfl[dfl == -Inf] <- 0
  dfl[dfl < 0] <- 0
  dfl
}

MRnorm <- function(data){
  # data: dataframe of gene expression counts to be normalized via DESeq2
  meta <- data.frame(experiment = colnames(data))
  dds <- DESeqDataSetFromMatrix(countData = 
                                       ceiling(data), 
                                     colData = meta, design = ~ experiment)
  dds <- estimateSizeFactors(dds)
  counts(dds, normalized=TRUE)
}

get_zeros <- function(data, lt= .1, ut = .9){
  # data:
  # lt:
  # ut:
  quantile(apply(
    data, 2, 
    FUN = function(x) sum( x  == 0)), probs = c(lt,ut))
}

get_hks <- function(data, hk_genes, lt=.2,ut=.98){
  # data:
  # hk_genes:
  # lt: 
  # ut:
  quantile(apply(data[hks_genes,], 
                 2, FUN = function(x) median(x)), probs = c(lt,ut))
}

filter_hks <- function(data, hk_genes=NA, hk_min=4, hk_max=6){
  # data: dataframe of gene expression with rownames as gene names
  # hk_genes: array of gene names
  # hk_min: float of minimun mean expression, default 0
  # hk_max: float of maximum mean expression, deafult infinity
  # returns: array of bools for filtering columns of data
  if(is.na(hk_genes)){
    hk_genes <- sapply(c('ppiD','rpoD','rpoS','proC','recA','rpsL','rho','oprL','tipA','nadB','ampC'), function(x) name_to_PAO1(x))
  }
  hk_clean <- hk_genes %in% rownames(data)
  if(sum(hk_clean) < 1){
    print('HK genes not in dataframe rownames')
    break
  } else{
    hk_genes <- hk_genes[hk_clean]
  }
  hk_means <- apply(data[hk_genes,], 2, FUN = function(x) median(x))
  hk_filt <- (hk_means >= hk_min) & (hk_means <= hk_max)
  return(hk_filt)
}

filter_sparsity <- function(data, max_zeros=1000, min_zeros = 200, min_peak=0){
  # data: dataframe of gene expression with rownames as gene names
  # max_zeros: int of maximum unexpressed gene per sample, default infinity
  # min_peak: float [0-1] of multiplier for peak cutoff, default 1 (0.15 good)
  # returns: array of bools for filtering columns of data
  if(max_zeros < 0 | min_peak < 0){
    print('numer of zeros and peak coefficinet must be positive float values')
    break
  }
  toosp_test <- colSums(data == 0) < max_zeros
  unsparse_test_z <- colSums(data == 0) > min_zeros
  unsparse_test_p <- sapply(colnames(data), function(x){
    h1 <- hist(as.numeric(data[,x]), plot = F, breaks=50)
    if(h1$density[1] <= min_peak){
      F
    } else {T}
  })
  sparse_filt <- toosp_test & unsparse_test_z #& unsparse_test_p 
  return(sparse_filt)
}



array_match <- function(a_comp, s_comp, a_ref=T, fsqn=T){
  dfs <- list(a_comp, s_comp)
  genes <- lapply(dfs, function(x) rownames(x))
  genes_venn <- Venn(genes)
  common_genes <- overlap(genes_venn)
  if(fsqn){
    if(a_ref){
      s_comp_fsqn <- quantileNormalizeByFeature(t(s_comp[common_genes,]),t(a_comp[common_genes,]))
      comp_bind <- cbind(a_comp[common_genes,],t(s_comp_fsqn))
    } else{
      a_comp_fsqn <- quantileNormalizeByFeature(t(a_comp[common_genes,]),t(s_comp[common_genes,]))
      comp_bind <- cbind(t(a_comp_fsqn),s_comp[common_genes,])
    }
  }else{
    comp_bind <- cbind(a_comp[common_genes,],s_comp[common_genes,])
  }
  
  return(comp_bind)
}
