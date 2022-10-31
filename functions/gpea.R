#' @description gpea:Gene pair enrichment analysis (GPEA)
#' #' The enrichment analysis is based on a one-sided Fisher's exact test.
#' 
#' Modification:
#' When a network G contains n interactions, of which k interactions are among genes from the given gene set S,
#' then a p-value for the enrichment of gene pairs of this gene set S can be calculated based on a e.g., one-sided Fisher's exact test. 
#' For p genes there is a total of N=p(p-1)/2 different gene pairs, if G is undirected, whereas 
#' For p genes there is a total of N=p(p-1) different gene pairs, if G is directed (clique graph) 
#' with the assumption that all genes within a gene set are associated to each other. 
#' If there are pS genes for a particular gene set (S) then the total number of gene pairs for this gene set is mS=pS(pS-1)/2, in case of undirected graph. If directed, mS=pS(pS-1))
#' 
#' @param gnet igraph object of a given network where the gene identifiers [V(net)$names] correspond to the provided gene identifiers in the reference gene sets.
#' @param genesets A named list object of a collection of gene sets. The identifiers used for the candidate and reference genes need to match the identifier types used for the gene sets.
#' @param verbose The default value is . If this option is set the number and name of the gene sets during their processing is reported.
#' @param cmax All provided genesets with more than cmax genes will be excluded from the analysis (default cmax=1000).
#' @param cmin All provided genesets with less than cmin genes will be excluded from the analysis (default cmin>=3).
#' @param adj The default value is (False discovery rate using the Benjamini-Hochberg approach). Multiple testing correction based on the function stats::p.adjust() with available options for "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr" and "none"
#' 
#' @output The function returns a data.frame object with the columns"TermID" the given name of a gene set i from the named gene set collection list object S. "edges" the number of connected gene pairs present a given geneset "genes" the number of candidate genes present in the gene set i "all" the number of all genes present in the gene set i "pval" the nominal p-value from a one-sided fisher's exact test "padj" the adjusted p-value to consider for multiple testing
#
# This function is a modified version of gpea function from bc3net package
#
#REFERENCE: 
#https://www.rdocumentation.org/packages/bc3net/versions/1.0.4/topics/gpea
#https://github.com/cran/bc3net


gpea <- function(gnet, genesets, directed=T, verbose=TRUE, cmax=1000, cmin=3, adj="bonferroni"){
  
  
  genes = V(gnet)$name  
  genesets = lapply(genesets, function(x) x[x%in%genes])
  genesets = genesets[sapply(genesets, function(x) length(x)>=cmin & length(x)<cmax)]
  
  
  # considering only genes defined in the genesets
  allgenes = unique(unlist(genesets))
  
  ref.net = induced.subgraph(gnet, allgenes)
  v.ref = vcount(ref.net)
  e.ref = ecount(ref.net)
  
  #THIS PART IS MODIFIED by me from its original function
  if(directed){
    ref.clique = (vcount(ref.net)*(vcount(ref.net)-1))
  }else{
    ref.clique = (vcount(ref.net)*(vcount(ref.net)-1))/2
  }
  
  res=matrix(0, ncol=4, nrow=length(genesets))
  rownames(res)=names(genesets)
  colnames(res)=c("edges", "genes", "pval", "padj")
  
  
  for(i in 1:nrow(res)){
    
    
    if(verbose==TRUE){
      cat(i," of ",nrow(res),"\n")
    }
    
    term = genesets[[i]]
    sub.net = induced.subgraph(gnet, term)
    
    v.sub=vcount(sub.net)
    e.sub=ecount(sub.net)
    
    #THIS PART IS MODIFIED by me from its original function
    if(directed){
      e.clique=(v.sub*(v.sub-1))
    }else{
      e.clique=(v.sub*(v.sub-1))/2
    }
    
    
    # sub of genesets[[i]] in GRN | sub not in GRN (clique)  
    # not in genesets[[i]] in GRN | not in genesets[[i]] (clique)
    
    # Example: #######################################
    #     170 | 2,680
    #  59,933 | 54,561,831
    ############################## p=8.15e-227
    
    mat = matrix(c(e.sub, e.ref-e.sub,
                   e.clique-e.sub, ref.clique-e.clique-e.ref+e.sub),
                 ncol=2, nrow=2)
    
    pval = fisher.test(mat, alternative="greater")$p.value
    
    res[i,1:3]=c(e.sub, v.sub, pval)
  }
  
  # res=res[res[,2]<=cmax & res[,2]>cmin,]
  res=res[order(res[,3]),]
  res[,4]=p.adjust(res[,3], method=adj)  
  
  
  tab=data.frame(TermID=rownames(res),res)
  
  return(tab)
}