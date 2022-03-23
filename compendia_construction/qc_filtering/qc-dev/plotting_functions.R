## ---------------------------
##
## Script name: plotting_functions.R
##
## Description:
##
## Args:
##
## Author: Georgia Doing
##
## Date Created: 2021-06-21
##
## Email: Georgia.Doing.GR@Dartmouth.edu
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

annotated_heatmap <- function(data, count_data, sparse, hks_b,  anns=T, name = ""){
  plot_mat <- t((t(data.matrix(data))))
  plot_mat[is.na(plot_mat)] <- 1
  genes_sets <- sapply(rownames(plot_mat), function(x){
    if(x %in% hks){
      "hks"
    } else{
      if (x %in% lgs){
        "lgs"
      } else{"other_genes"}
    }})
  
  ha_row <- HeatmapAnnotation(reg = genes_sets, which='row', 
                              col= list(reg = c ('hks'='red', 'lgs' ='blue', 'other_genes'='gray')),
                              annotation_legend_param = list(
                                reg = list( title = "Gene Set",
                                            at = c("hks","lgs", "other_genes"),
                                            labels = c("Houskeeping","Lowly Expressed", "Random"),
                                            nrow=1)
                                
                              ))
  
  
  
  rownames(plot_mat) <- sapply(rownames(plot_mat), function(x) PAO1_to_name(x))
  if(anns){
    
    ha_col <- HeatmapAnnotation(depth = anno_barplot(log(colSums(count_data)), ylim = c(10,19), baseline = 10),#baseline='min'
                                #hk_gene_mean = anno_barplot(apply(data,2,FUN = function(x) median(x)), baseline='min'),
                                #compendium = c(rep('array',sum(colnames(plot_mat) %in% colnames(array_comp))), 
                                #               rep('seq',sum(colnames(plot_mat) %in% colnames(rnaseq_pao1_tpm)))),
                                sparsity = sparse,
                                hk_genes = hks_b,
                                #medians = medians,
                                annotation_legend_param = list(
                                  #compendium = list( title = "compendium",
                                  #                   at = c("array","seq"),
                                  #                   labels = c("Array","RNAseq")),
                                  sparsity = list( title = "sparsity",
                                                   at = c("TRUE","FALSE"),
                                                   labels = c("pass", "fail")),
                                  hk_genes = list( title = "hk_genes",
                                                   at = c("TRUE","FALSE"),
                                                   labels = c("pass", "fail"))#,
                                  #medians = list( title = "medians",
                                  #                at = c("TRUE","FALSE"),
                                  #                labels = c("pass", "fail"))
                                  
                                ),
                                col= list(#compendium = c('array' = 'green', 'seq'='orange'),
                                  sparsity = c('TRUE' = 'grey', 'FALSE'='red'),
                                  #medians = c('TRUE' = 'grey', 'FALSE'='red'),
                                  hk_genes = c('TRUE' = 'grey', 'FALSE'='red')
                                ),
                                show_annotation_name = T,
                                annotation_name_side='left'
    )
    
    ht <- Heatmap(plot_mat, 
                  show_column_names = F, 
                  show_row_names=T, 
                  row_labels = sapply(rownames(plot_mat), function(x) PAO1_to_name(x)),
                  row_names_gp = gpar(fontsize=9),
                  left_annotation = ha_row,
                  top_annotation = ha_col,
                  heatmap_legend_param = list(
                    title = "logTPM",
                    direction="horizontal"#,
                    #legend_width= unit(4,"cm")
                  ),
                  #column_split = c(rep('array',sum(colnames(plot_mat) %in% colnames(array_comp))), 
                  #                 rep('seq',sum(colnames(plot_mat) %in% colnames(rnaseq_pao1_tpm)))),
                  row_split = factor(genes_sets,
                                     levels=c('hks','lgs','other_genes')),
                  column_title = name
    )
  } else{
    ha_col <- HeatmapAnnotation(depth = anno_barplot(log(colSums(count_data)), ylim=c(10,19),axis=F, baseline=10), #baseline='min'
                                #hk_gene_mean = anno_barplot(apply(data,2,FUN = function(x) median(x)), baseline='min'),
                                #compendium = c(rep('array',sum(colnames(plot_mat) %in% colnames(array_comp))), 
                                #               rep('seq',sum(colnames(plot_mat) %in% colnames(rnaseq_pao1_tpm)))),
                                sparsity = sparse,
                                hk_genes = hks_b,
                                #medians = medians,
                                annotation_legend_param = list(
                                  #compendium = list( title = "compendium",
                                  #                   at = c("array","seq"),
                                  #                   labels = c("Array","RNAseq")),
                                  sparsity = list( title = "sparsity",
                                                   at = c("TRUE","FALSE"),
                                                   labels = c("pass", "fail"))
                                  
                                ),
                                col= list(#compendium = c('array' = 'green', 'seq'='orange'),
                                          sparsity = c('TRUE' = 'grey', 'FALSE'='red'),
                                          #medians = c('TRUE' = 'grey', 'FALSE'='red'),
                                          hk_genes = c('TRUE' = 'grey', 'FALSE'='red')
                                ),
                                show_annotation_name = F
    )
    
    ht <- Heatmap(plot_mat, 
                  show_column_names = F, 
                  #show_column_dend = F
                  show_row_names=T, 
                  row_labels = sapply(rownames(plot_mat), function(x) PAO1_to_name(x)),
                  row_names_gp = gpar(fontsize=9),
                  #left_annotation = ha_row,
                  top_annotation = ha_col,
                  #heatmap_legend_param = list(
                  #  title = "logTPM"
                  #),
                  show_heatmap_legend = F,
                  #column_split = c(rep('array',sum(colnames(plot_mat) %in% colnames(array_comp))), 
                  #                 rep('seq',sum(colnames(plot_mat) %in% colnames(rnaseq_pao1_tpm)))),
                  row_split = factor(genes_sets,
                                     levels=c('hks','lgs','other_genes')),
                  column_title=name
    )
  }
  
  
}


###################################################################################
annotated_heatmap_filt <- function(data, count_data, sparse, hks_b, medians){
  h1 <- annotated_heatmap(data, count_data, sparse, hks_b, medians, name="Unfiltered",
                          anns=T)
  
  h2 <- annotated_heatmap(data[,sparse & hks_b & medians], 
                          count_data[sparse & hks_b & medians], 
                          sparse[sparse & hks_b & medians], 
                          hks_b[sparse & hks_b & medians], 
                          medians[sparse & hks_b & medians], 
                          name="Filtered",
                          anns=F)
  draw(h1+h2, heatmap_legend_side = "bottom")
}



comparison_heatmap <- function(data, count_data, sparse, hks_b, medians, anns=T, name = ""){
  plot_mat <- t((t(data.matrix(data))))
  plot_mat[is.na(plot_mat)] <- 1
  genes_sets <- sapply(rownames(plot_mat), function(x){
    if(x %in% hks){
      "hks"
    } else{
      if (x %in% lgs){
        "lgs"
      } else{"other_genes"}
    }})
  
  ha_row <- HeatmapAnnotation(reg = genes_sets, which='row', 
                              col= list(reg = c ('hks'='red', 'lgs' ='blue', 'other_genes'='gray')),
                              annotation_legend_param = list(
                                reg = list( title = "Gene Set",
                                            at = c("hks","lgs", "other_genes"),
                                            labels = c("Houskeeping","Lowly Expressed", "Random"),
                                            nrow=1)
                                
                              ))
  
  
  ha_col <- HeatmapAnnotation(depth = anno_barplot(log(colSums(count_data)), ylim = c(10,19), baseline = 10),#baseline='min'
                              #hk_gene_mean = anno_barplot(apply(data,2,FUN = function(x) median(x)), baseline='min'),
                              compendium = c(rep('array',sum(colnames(plot_mat) %in% colnames(array_comp))), 
                                             rep('seq',sum(colnames(plot_mat) %in% colnames(rnaseq_pao1_tpm)))),
                              sparsity = sparse,
                              hk_genes = hks_b,
                              medians = medians,
                              annotation_legend_param = list(
                                compendium = list( title = "compendium",
                                                   at = c("array","seq"),
                                                   labels = c("Array","RNAseq")),
                                sparsity = list( title = "sparsity",
                                                 at = c("TRUE","FALSE"),
                                                 labels = c("pass", "fail")),
                                hk_genes = list( title = "hk_genes",
                                                 at = c("TRUE","FALSE"),
                                                 labels = c("pass", "fail")),
                                medians = list( title = "medians",
                                                at = c("TRUE","FALSE"),
                                                labels = c("pass", "fail"))
                                
                              ),
                              col= list(compendium = c('array' = 'green', 'seq'='orange'),
                                        sparsity = c('TRUE' = 'grey', 'FALSE'='red'),
                                        medians = c('TRUE' = 'grey', 'FALSE'='red'),
                                        hk_genes = c('TRUE' = 'grey', 'FALSE'='red')
                              ),
                              show_annotation_name = T,
                              annotation_name_side='left'
  )
  rownames(plot_mat) <- sapply(rownames(plot_mat), function(x) PAO1_to_name(x))
  if(anns){
    ht <- Heatmap(plot_mat, 
                  show_column_names = F, 
                  show_row_names=T, 
                  row_labels = sapply(rownames(plot_mat), function(x) PAO1_to_name(x)),
                  row_names_gp = gpar(fontsize=9),
                  left_annotation = ha_row,
                  top_annotation = ha_col,
                  heatmap_legend_param = list(
                    title = "logTPM",
                    direction="horizontal",
                    legend_width= unit(4,"cm")
                  ),
                  column_split = c(rep('array',sum(colnames(plot_mat) %in% colnames(array_comp))), 
                                   rep('seq',sum(colnames(plot_mat) %in% colnames(rnaseq_pao1_tpm)))),
                  row_split = factor(genes_sets,
                                     levels=c('hks','lgs','other_genes')),
                  column_title = name
    )
  } else{
    ha_col <- HeatmapAnnotation(depth = anno_barplot(log(colSums(count_data)), ylim=c(10,19),axis=F, baseline=10), #baseline='min'
                                #hk_gene_mean = anno_barplot(apply(data,2,FUN = function(x) median(x)), baseline='min'),
                                compendium = c(rep('array',sum(colnames(plot_mat) %in% colnames(array_comp))), 
                                               rep('seq',sum(colnames(plot_mat) %in% colnames(rnaseq_pao1_tpm)))),
                                sparsity = sparse,
                                hk_genes = hks_b,
                                medians = medians,
                                annotation_legend_param = list(
                                  compendium = list( title = "compendium",
                                                     at = c("array","seq"),
                                                     labels = c("Array","RNAseq")),
                                  sparsity = list( title = "sparsity",
                                                   at = c("TRUE","FALSE"),
                                                   labels = c("pass", "fail"))
                                  
                                ),
                                col= list(compendium = c('array' = 'green', 'seq'='orange'),
                                          sparsity = c('TRUE' = 'grey', 'FALSE'='red'),
                                          medians = c('TRUE' = 'grey', 'FALSE'='red'),
                                          hk_genes = c('TRUE' = 'grey', 'FALSE'='red')
                                ),
                                show_annotation_name = F
    )
    
    ht <- Heatmap(plot_mat, 
                  show_column_names = F, 
                  #show_column_dend = F
                  show_row_names=T, 
                  row_labels = sapply(rownames(plot_mat), function(x) PAO1_to_name(x)),
                  row_names_gp = gpar(fontsize=9),
                  #left_annotation = ha_row,
                  top_annotation = ha_col,
                  #heatmap_legend_param = list(
                  #  title = "logTPM"
                  #),
                  show_heatmap_legend = F,
                  column_split = c(rep('array',sum(colnames(plot_mat) %in% colnames(array_comp))), 
                                   rep('seq',sum(colnames(plot_mat) %in% colnames(rnaseq_pao1_tpm)))),
                  row_split = factor(genes_sets,
                                     levels=c('hks','lgs','other_genes')),
                  column_title=name
    )
  }
  
  
}

comparison_heatmap_filt <- function(data, count_data, sparse, hks_b, medians){
  h1 <- annotated_heatmap(data, count_data, sparse, hks_b, medians, name="Unfiltered",
                          anns=T)
  
  h2 <- annotated_heatmap(data[,sparse & hks_b & medians], 
                          count_data[sparse & hks_b & medians], 
                          sparse[sparse & hks_b & medians], 
                          hks_b[sparse & hks_b & medians], 
                          medians[sparse & hks_b & medians], 
                          name="Filtered",
                          anns=F)
  draw(h1+h2, heatmap_legend_side = "bottom")
}
