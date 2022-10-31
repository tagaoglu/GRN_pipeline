#Custom Rcistarget function to select high confidence links
Custom.Rcis=function(Networks,
                     chosenDb,
                     MinGenesetSize=0, 
                     directed=T){
  library(tools)
  library(RcisTarget)
  #https://bioconductor.riken.jp/packages/3.9/bioc/vignettes/RcisTarget/inst/doc/RcisTarget.html
  ## Aibar et al. (2017) SCENIC: single-cell regulatory network 
  ##   inference and clustering. Nature Methods. doi: 10.1038/nmeth.4463
  
  ## Select motif database to use (i.e. organism and distance around TSS)
  #https://resources.aertslab.org/cistarget/
  setwd("~/p728/RSTUDIO/analysis/tagaoglu/data/cisTarget_databases")
  motifRankings <- importRankings(chosenDb)
  data(motifAnnotations_hgnc)
  
  #Filtering and selecting only the genes that are available in RcisTarget databases, 
  #i.e. excluding genes missing from database
  #and also Formatting the targets from GRN algorithms into co-expression modules
  TFtoTargetNetworks= lapply(Networks, function(x){
    if(directed==F){
      LinkstoKeep=c()
      LinkstoKeep.GeneStatus= data.frame(gene1= c(), gene2=c())
      for (i in 1:nrow(x)){
        Genes= as.character(x[i,1:2])
        GeneStatus= c()
        for (j in 1:length(Genes)){
          if (Genes[j] %in% colnames(getRanking(motifRankings))){
            GeneStatus[j]= 'target'
          }else{
            if (Genes[j] %in% motifAnnotations_hgnc$TF){
              GeneStatus[j]= 'TF'
            }else{GeneStatus[j]= NA}
          }
        }
        if (any(is.na(GeneStatus)) | all(GeneStatus=='target')){LinkstoKeep= LinkstoKeep}
        else{
          LinkstoKeep= c(LinkstoKeep, i)
          LinkstoKeep.GeneStatus[match(i, LinkstoKeep),1:2]= GeneStatus
        }
      }
      DbFilteredNetwork= x[LinkstoKeep,]
      for(i in 1:nrow(DbFilteredNetwork)){
        if (LinkstoKeep.GeneStatus[i,1] == 'target'){
          DbFilteredNetwork[i, 1:2]= rev( DbFilteredNetwork[i, 1:2])
        }
      }
      return(DbFilteredNetwork)
    }else{return(x)}
  })
  
  FilteredTFtoTargetNetworks= lapply(TFtoTargetNetworks, function(x){
    if(directed){
      networkFilteredforDbgenes= x[which(x$gene2 %in% colnames(getRanking(motifRankings))),]
      networkFilteredforDbgenes= networkFilteredforDbgenes[which(networkFilteredforDbgenes$gene1 %in% motifAnnotations_hgnc$TF),]
      return(networkFilteredforDbgenes)
    }else{return(x)}
  })
  
  #Create gene-sets associated to each TF for next step, RcisTarget
  NetworkGenesets= lapply(FilteredTFtoTargetNetworks, function(x){
    Geneset= list()
    for (i in unique(x$gene1)){
      table= x[which(x$gene1==i),]
      if (nrow(table) >= MinGenesetSize){
        Geneset[[i]]= as.character(table$gene2)
      }
    }
    return(Geneset)
  })
  
  ## Calculate motif enrichment for each TF-module (Run RcisTarget)
  #1.Calculate enrichment
  AUC= lapply(NetworkGenesets, calcAUC, rankings= motifRankings)
  
  #2.Select significant motifs and/or annotate to TFs
  #Convert to table, filter by NES & add the TFs to which the motif is annotated
  motifEnrichment= lapply(AUC, function(aucOutput){
    # Extract the TF of the gene-set name (i.e. MITF_w001):
    tf <- sapply(setNames(strsplit(rownames(aucOutput), "_"), rownames(aucOutput)), function(x) x[[1]])
    # Calculate NES and add motif annotation (provide tf in 'highlightTFs'):
    addMotifAnnotation(aucOutput, 
                       nesThreshold=3, digits=3, 
                       motifAnnot=motifAnnotations_hgnc,
                       motifAnnot_highConfCat=c("directAnnotation", "inferredBy_Orthology"),
                       motifAnnot_lowConfCat=c("inferredBy_MotifSimilarity",
                                               "inferredBy_MotifSimilarity_n_Orthology"), 
                       highlightTFs=tf)
  })
  
  #Keep only the motifs annotated to the initial TF
  motifEnrichment_selfMotifs = lapply(motifEnrichment, function(x){
    selfMotif= x[which(x$TFinDB != ""),]
    return(selfMotif)
  })
  
  #3.Prune targets
  #Identify the genes with the best enrichment for each Motif
  EnrichmentTable= list()
  for(i in 1:length(motifEnrichment_selfMotifs)){
    EnrichmentTable[[i]]= addSignificantGenes(motifEnrichment_selfMotifs[[i]],
                                              geneSets = NetworkGenesets[[i]],
                                              rankings = motifRankings)
  }
  names(EnrichmentTable)= names(motifEnrichment_selfMotifs)
  
  
  EnrichmentTableHighConf= lapply(EnrichmentTable, function(x){
    if(nrow(x) !=0){
      HighConf= x[which(x$TFinDB=='**'),]
      return(HighConf)
    }else{return(x)}
  })
  
  #Format regulons: In order to build the regulons, merge the genes from any of the enriched motifs for the same TF.
  HighConfRcisTargetNetwork= lapply(EnrichmentTableHighConf, function(x){
    regulons= data.frame(TF= c(), target=c())
    if(nrow(x) != 0){
      for(i in unique(x$geneSet)){
        str= x[ which(x$geneSet==i),'enrichedGenes']
        str= unlist(strsplit(str$enrichedGenes, ';'))
        str= unique(str)
        table=data.frame(TF= rep(i, length(str)), target= str)
        regulons= rbind(regulons, table)
      }
    }
    return(regulons)
  })
  
  CustomRcisTargetResult= list()
  for(i in 1:length(HighConfRcisTargetNetwork)){
    CustomRcisTargetResult[[i]]= list(AnalyzedLinks= FilteredTFtoTargetNetworks[[i]],
                                      FilteredLinks= HighConfRcisTargetNetwork[[i]])
  }
  names(CustomRcisTargetResult)= names(HighConfRcisTargetNetwork)
  
  return(CustomRcisTargetResult)
}
