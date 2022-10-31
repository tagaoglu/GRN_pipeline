#' @description Mean gene filtering : Filter-out low-abundance genes (previous to GRN-inference) based on mean-based filtering
#' 
#' @param exprMat Expression matrix
#' @param aveCounts_thr average count threshold; by default, low-abundance genes are defined as those with an average count below a filter threshold of 0.25

geneFiltering_mean <- function(expMat,
                               aveCounts_thr=0.25)
{
  
  # Calculate stats
  aveCounts <- rowMeans(expMat, na.rm = T)
  
  ## Show info
  #message("Number of genes in the expression matrix: ", nrow(expMat))
  #message("Average counts per gene:")
  #print(summary(aveCounts))
  
  ## Filter
  #message("\nNumber of genes left after applying the filter:")
  genesLeft_aveCounts <- names(aveCounts)[which(aveCounts >= aveCounts_thr)]
  #message("\t", length(genesLeft_aveCounts), "\tgenes with average counts per gene >= ", aveCounts_thr)
  
  genesKept <- genesLeft_aveCounts
  
  return(genesKept)
}