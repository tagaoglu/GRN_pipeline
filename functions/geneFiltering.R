#' @description Gene filtering : Filter-out genes (previous to GRN-inference) based on the counts and number of cells in which they are detected.
#' @param exprMat Expression matrix
#' @param minCountsPerGene Minimum counts per gene required
#' @param minSamples Minimum number of samples (cells) in which the gene should be detected

# Gene filter/selection
#Before running GENIE3, it is recommended to apply soft gene filter, to remove genes that are expressed either at very low levels or in too few cells. Here we apply a filtering based on the total number of counts of the gene, and the number of cells in which it is detected.

#First filter: Filter by the total number of reads per gene. This filter is meant to remove genes that are most likely noise. By default it keeps only the genes with at least 6 UMI counts across all samples/cells (e.g. the total number the gene would have, if it was expressed with a value of 3 in 10% of the cells). Adjust this value (minCountsPerGene) according to the dataset (it will depend on the dataset units, e.g. UMI, TPMs…).

#Second filter: Filter by the number of cells in which the gene is detected** (e.g. >0 UMI, or >1 log2(TPM)). By default (minSamples), genes that are detected in at least 1% of the cells are kept. This filtering is meant to remove genes whose reads come from one a few ‘noisy’ cells (genes that are only expressed in one, or very few cells, gain a lot of weight if they happen to coincide in a given cell). To avoid removing small (but potentially interesting) cell populations, we recommend to set a percentage lower than the smallest population of cells to be detected.


geneFiltering <- function(expMat,
                          minCountsPerGene=3*.1*ncol(expMat),
                          minSamples=ncol(expMat)*.01)
{
  
  # Check expression matrix (e.g. not factor)
  if(is.data.frame(expMat)) 
  {
    supportedClasses <- paste(gsub("AUCell_buildRankings,", "", methods("AUCell_buildRankings")), collapse=", ")
    supportedClasses <- gsub("-method", "", supportedClasses)
    
    stop("'expMat' should be one of the following classes: ", supportedClasses, 
         "(data.frames are not supported. Please, convert the expression matrix to one of these classes.)")
  }
  if(any(table(rownames(expMat))>1))
    stop("The rownames (gene id/name) in the expression matrix should be unique.")
  
  
  # Calculate stats
  nCountsPerGene <- rowSums(expMat, na.rm = T)
  nCellsPerGene <- rowSums(expMat>0, na.rm = T)
  
  
  ## Show info
  message("Maximum value in the expression matrix: ", max(expMat, na.rm=T))
  message("Ratio of detected vs non-detected: ", signif(sum(expMat>0, na.rm=T) / sum(expMat==0, na.rm=T), 2))
  message("Number of counts (in the dataset units) per gene:")
  print(summary(nCountsPerGene))
  message("Number of cells in which each gene is detected:")
  print(summary(nCellsPerGene))
  
  ## Filter
  message("\nNumber of genes left after applying the following filters (sequential):")
  # First filter
  # minCountsPerGene <- 3*.01*ncol(expMat)
  genesLeft_minReads <- names(nCountsPerGene)[which(nCountsPerGene > minCountsPerGene)]
  message("\t", length(genesLeft_minReads), "\tgenes with counts per gene > ", minCountsPerGene)
  
  # Second filter
  # minSamples <- ncol(expMat)*.01
  nCellsPerGene2 <- nCellsPerGene[genesLeft_minReads]
  genesLeft_minCells <- names(nCellsPerGene2)[which(nCellsPerGene2 > minSamples)]
  message("\t", length(genesLeft_minCells), "\tgenes detected in more than ",minSamples," cells")
  
  genesKept <- genesLeft_minCells  # or genesLeft_minCells_inDatabases
  
  return(genesKept)
}
