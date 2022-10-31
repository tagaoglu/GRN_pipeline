# Set *.libPaths*  
.libPaths(c("~/p728/RSTUDIO/R/library/4.1/", .libPaths()))


#Set working directory and seed:
setwd("~/p728/RSTUDIO/analysis/tagaoglu")
set.seed(123)


#Load functions
source("/home/tagaoglu/p728/RSTUDIO/analysis/tagaoglu/functions/save_png_pdf.R")
source("/home/tagaoglu/p728/RSTUDIO/analysis/tagaoglu/functions/fgsea_bubblePlot.R")


##################################################################################################
###Load data
##################################################################################################

### First run "Critical_genes_analysis.R", then run codes below


##################################################################################################
## Perform enrichment analysis
##################################################################################################

#########
###GENIE3


### fgsea
library(msigdbr) #https://cran.r-project.org/web/packages/msigdbr/vignettes/msigdbr-intro.html
library(fgsea)  #https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html#


## 1. Generate ranked list
# since we have many samples, use sapply to make the rank genes in a named list (simplify=F will make sapply return results in a list)

# GENIE_net_int_rank <- sapply(names(GENIE_net_int), function(nam){
#   sapply(names(GENIE_net_int[[nam]]), function(group){
#     sapply(names(GENIE_net_int[[nam]][[group]]), function(cell){
#       setNames(GENIE_net_int[[nam]][[group]][[cell]]$minuslog10[order(GENIE_net_int[[nam]][[group]][[cell]]$minuslog10, decreasing = T)],
#                rownames(GENIE_net_int[[nam]][[group]][[cell]])[order(GENIE_net_int[[nam]][[group]][[cell]]$minuslog10, decreasing = T)])
#     }, simplify=F)
#   }, simplify=F)
# }, simplify=F)
      


# Generate ranked list, just for target genes-regulated by JUN for Viral Endothelial
groups <- c("viral")
cells <- c("Endothelial")

#Set TF, to extract neighbours list (target genes) of it
TF <- "JUN"


GENIE_net_target_genes <- list()
GENIE_net_targets <- list()

# Rank based on the integrated centrality measure (-log10 ??? or betweenness centrality(can find the code in the old version if needed)???)

#Code for -log10

#Get the first degree neighbors of a specific TF
for(nam in names(GENIE_net_igraph)){
  for(group in groups){
    for(cell in cells){
      GENIE_net_target_genes[[nam]][[group]][[cell]] <- names(neighbors(GENIE_net_igraph[[nam]][[group]][[cell]],TF)) ##or names(GENIE_net_igraph[[nam]][[group]][[cell]][[TF]])
      GENIE_net_target_genes[[nam]][[group]][[cell]] <- append(GENIE_net_target_genes[[nam]][[group]][[cell]], TF) # add TF
      GENIE_net_targets[[nam]][[group]][[cell]] <- subset(GENIE_net_int[[nam]][[group]][[cell]], rownames(GENIE_net_int[[nam]][[group]][[cell]]) %in% GENIE_net_target_genes[[nam]][[group]][[cell]])
    }
  }
}
      

#Get the graph of a specific TF vertexes and all their first 3 degree neighbors
for(nam in names(GENIE_net_igraph)){
  for(group in groups){
    for(cell in cells){
      sub_net <- induced.subgraph(graph=GENIE_net_igraph[[nam]][[group]][[cell]],
                                  vids=unlist(neighborhood(graph=GENIE_net_igraph[[nam]][[group]][[cell]],order=3,nodes=TF)))
      GENIE_net_target_genes[[nam]][[group]][[cell]] <- names(V(sub_net))
      GENIE_net_target_genes[[nam]][[group]][[cell]] <- append(GENIE_net_target_genes[[nam]][[group]][[cell]], TF) # add TF
      GENIE_net_targets[[nam]][[group]][[cell]] <- subset(GENIE_net_int[[nam]][[group]][[cell]], rownames(GENIE_net_int[[nam]][[group]][[cell]]) %in% GENIE_net_target_genes[[nam]][[group]][[cell]])
    }
  }
}
      

GENIE_net_int_rank <- sapply(names(GENIE_net_int), function(nam){
  sapply(groups, function(group){
    sapply(cells, function(cell){
      setNames(GENIE_net_targets[[nam]][[group]][[cell]]$minuslog10[order(GENIE_net_targets[[nam]][[group]][[cell]]$minuslog10, decreasing = T)],
               rownames(GENIE_net_targets[[nam]][[group]][[cell]])[order(GENIE_net_targets[[nam]][[group]][[cell]]$minuslog10, decreasing = T)])
    }, simplify=F)
  }, simplify=F)
}, simplify=F)




## 2. Download gene sets 
# Download gene sets and convert to gmt --> takes several min. you can save the final object for future use!
msig <- c("C2", "C4", "C5", "C6", "C8", "H") # choose your pathways. better more than less, you can also include all of them and omit from results what is not interesting.

msigdb.gmt <- sapply(msig,function(x) NULL)

for (i in msig) {
  tmp <- as.data.frame(msigdbr(species = "Homo sapiens", category = i))
  msigdb.gmt[[i]] <- lapply(unique(tmp$gs_name), function(x) {
    tmp[tmp$gs_name == x, "gene_symbol"]
  })
  names(msigdb.gmt[[i]]) <- unique(paste0(tmp$gs_name))
}



## 3. Run fgsea
# "loop" over samples and "loop" over every pathway in msig. Resulting object is a list of lists.

GENIE_net_int_fgseaRes <- sapply(names(GENIE_net_int_rank), function(nam){
  sapply(names(GENIE_net_int_rank[[nam]]), function(group){
    sapply(names(GENIE_net_int_rank[[nam]][[group]]), function(cell){
      sapply(msig, function(msig) {
        fgsea(pathways = msigdb.gmt[[msig]],
              stats = GENIE_net_int_rank[[nam]][[group]][[cell]],
              eps = 0.0, minSize = 15, maxSize = 500)
      }, simplify=F)
    }, simplify=F)
  }, simplify=F)
}, simplify=F)



## 4. Visualization

# Summary table (top 10 & bottom 10)
GENIE_net_int_fgseaRes_topPathways <- sapply(names(GENIE_net_int_rank), function(nam){
  sapply(names(GENIE_net_int_rank[[nam]]), function(group){
    sapply(names(GENIE_net_int_rank[[nam]][[group]]), function(cell){
      sapply(msig, function(msig) {
        c(head(GENIE_net_int_fgseaRes[[nam]][[group]][[cell]][[msig]][which(GENIE_net_int_fgseaRes[[nam]][[group]][[cell]][[msig]]$NES > 0),][order(padj)][order(-NES), pathway], n = 10),
          tail(GENIE_net_int_fgseaRes[[nam]][[group]][[cell]][[msig]][which(GENIE_net_int_fgseaRes[[nam]][[group]][[cell]][[msig]]$NES < 0),][order(-padj)][order(-NES), pathway], n = 10))
      }, simplify=F)
    }, simplify=F)
  }, simplify=F)
}, simplify=F)
        

#or 
#c(GENIE_net_int_fgseaRes[[nam]][[group]][[cell]][[msig]][NES > 0][head(order(pval), n=10), pathway],
#  rev(GENIE_net_int_fgseaRes[[nam]][[group]][[cell]][[msig]][NES < 0][head(order(pval), n=10), pathway]))



# Save GseaTable plot

##Create files if it does not exist
for (nam in names(GENIE_net_int_fgseaRes)) {
  plotdir_GENIE=paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/GseaTable/GENIE3/", nam)
  if(!dir.exists(plotdir_GENIE)){dir.create(plotdir_GENIE,  recursive = T)}
}


# cannot save TableGrob to object, plot directly to device
# using my save_png_pdf wrapper fro plotting. 

for(nam in names(GENIE_net_int_fgseaRes)){
  for(group in names(GENIE_net_int_fgseaRes[[nam]])){
    for(cell in names(GENIE_net_int_fgseaRes[[nam]][[group]])){
      for(sig in msig){
        save_png_pdf(plotGseaTable(msigdb.gmt[[sig]][GENIE_net_int_fgseaRes_topPathways[[nam]][[group]][[cell]][[sig]]], 
                                   GENIE_net_int_rank[[nam]][[group]][[cell]], 
                                   GENIE_net_int_fgseaRes[[nam]][[group]][[cell]][[sig]], 
                                   gseaParam=0.5), 
                     paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/GseaTable/GENIE3/", nam, "/", group, "_", cell, ".fgsea.", sig), height = 8, width = 18, res=250)
      }
    }
  }
}




# Save Bubble plot

##Create files if it does not exist
for (nam in names(GENIE_net_int_fgseaRes)) {
  plotdir_GENIE=paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/GseaBubblePlot/GENIE3/", nam)
  if(!dir.exists(plotdir_GENIE)){dir.create(plotdir_GENIE,  recursive = T)}
}


#Bubble plot of gene sets on y-axis and adjusted p-value (padj) on x-axis. 
#Bubble size indicates the number of genes in each gene set, and bubble color indicates the normalized enrichment score (NES). 
#Blue is for negative NES (enrichment of higher targeted genes in males), and red is for positive NES (enrichment of higher targeted genes in females).

for(nam in names(GENIE_net_int_fgseaRes)){
  for(group in names(GENIE_net_int_fgseaRes[[nam]])){
    for(cell in names(GENIE_net_int_fgseaRes[[nam]][[group]])){
      for(sig in msig){
        save_png_pdf(fgsea_bubblePlot(GENIE_net_int_fgseaRes[[nam]][[group]][[cell]][[sig]], fdrcut = 0.1), 
                     paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/GseaBubblePlot/GENIE3/", nam, "/", group, "_", cell, ".fgsea.", sig), height = 8, width = 18, res=250)
      }
    }
  }
}




#####
#may be of use later
#GO over-representation analysis
library(org.Hs.eg.db)
library(clusterProfiler)

GENIE_net_ego <- list()

for(nam in names(GENIE_net_targets)){
  for(group in groups){
    for(cell in cells){
      GENIE_net_ego[[nam]][[group]][[cell]] <- enrichGO(gene         = GENIE_net_targets[[nam]][[group]][[cell]][["Name"]],
                                                        OrgDb         = org.Hs.eg.db,
                                                        keyType       = 'SYMBOL',
                                                        ont           = "BP",
                                                        pAdjustMethod = "BH",
                                                        pvalueCutoff  = 0.01,
                                                        qvalueCutoff  = 0.05)
    }
  }
}


nam <- 1

head(GENIE_net_ego[[nam]][[group]][[cell]], 3) 

df <- summary(GENIE_net_ego[[nam]][[group]][[cell]]) 
names(df)



#########
###PIDC

# ...



