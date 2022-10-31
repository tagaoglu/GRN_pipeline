# Set *.libPaths*  
.libPaths(c("~/p728/RSTUDIO/R/library/4.1/", .libPaths()))


#Set working directory and seed:
setwd("~/p728/RSTUDIO/analysis/tagaoglu")
set.seed(123)


#Load functions
source("/home/tagaoglu/p728/RSTUDIO/analysis/tagaoglu/functions/save_png_pdf.R")


#Load libraries
library(org.Hs.eg.db)
library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(ggplot2)


# Following analyses will be done just for centrality measure: "integrated_rank" 
centrality_interest <- "integrated_rank"


##################################################################################################
### COMPARISON : NORMAL vs VIRAL vs NONVIRAL
##################################################################################################

comp <- "comp_NORMALvsVIRALvsNONVIRAL"


#or

# ##################################################################################################
# ### COMPARISON : NORMAL vs TUMOR
# ##################################################################################################
# 
# comp <- "comp_NORMALvsTUMOR"


##################################################################################################
###Load data
##################################################################################################

#top_centr_all <- readRDS(paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Network_analysis/networks_hub_lists.rds"))
#top_centr_ub_all <- readRDS(paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Network_analysis/networks_hub_ubiquitous_", comp, ".rds"))
top_centr_sp_all <- readRDS(paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Network_analysis/networks_hub_specific_", comp, ".rds"))



##################################################################################################
## Perform Gene Ontology & KEGG pathway & MSigDb Enrichment Analysis
##################################################################################################

#As we saw before, among all networks and centralities there is a large proportion of hubs that 
#are specific to that group (group-specific hubs, multiplicity 1, purple), 
#as well as a large proportion of ubiquitous hubs (multiplicity 3+, orange). 

# By conducting a functional enrichment analysis using the Bioconductor package ClusterProfiler
# We will further examine group-specific hubs, with multiplicity 1 (purple)

#NOTE for GO analysis and the universe set of genes:
#To perform The Gene Ontology Enrichment Analysis (GOEA) we need to create a gene background called the universe and we will use all genes with a GO term. 
#Normally the universe (the background gene list) is made of the list of genes that were actually assayed in your transcriptome analysis, 
#in our case, this can be all genes present in the corresponding network


#Functions I will be using:
#enrichGO     : ORA using GO 
#enrichKEGG   : ORA using KEGG pathway 
#enricher     : ORA using MSigDb


###############################################################################
#enrichGO     : ORA using GO 

# Performs a Gene Ontology enrichment analysis for a given target gene set and universe.
# Args:  target: character vector with the gene symbols corresponding to the target set.
#        universe: integer/character??? vector with the gene symbols corresponding to the gene universe. 
# Returns:
#        A data.frame with the GO id, the p-value, the Odds score and the description of every enriched GO term.


## 1. Get target genelist

#Create dataframe to be able to visualize enrichment results across groups in one plot
df_target_sp <- list()
for (algo in names(top_centr_sp_all)){ 
  for (nam in names(top_centr_sp_all[[algo]])){                                                    
    for (cell in names(top_centr_sp_all[[algo]][[nam]])){ 
      for (centr in centrality_interest){    
        mydf <- data.frame()
        for (org in names(top_centr_sp_all[[algo]][[nam]][[cell]][[centr]])){
          l <- top_centr_sp_all[[algo]][[nam]][[cell]][[centr]][[org]]
          if(length(l) != 0){
            tmp <- data.frame(Gene=l)
            tmp$group <- org
            mydf <- rbind(mydf, tmp)
          }
        }
        df_target_sp[[algo]][[nam]][[cell]][[centr]] <- mydf
      }
    }
  }
}



## 2. Run enrichGO via compareCluster
enrichGO_sp <- list()
for (algo in names(top_centr_sp_all)){ 
  for (nam in names(top_centr_sp_all[[algo]])){                                                    
    for (cell in names(top_centr_sp_all[[algo]][[nam]])){ 
      for (centr in centrality_interest){                                                    
        #group-specific hubs
        enrichGO_sp[[algo]][[nam]][[cell]][[centr]] <- compareCluster(Gene~group, 
                                                                      data=df_target_sp[[algo]][[nam]][[cell]][[centr]], 
                                                                      fun="enrichGO",
                                                                      OrgDb         = org.Hs.eg.db,
                                                                      keyType       = 'SYMBOL',
                                                                      ont           = "BP",
                                                                      pAdjustMethod = "BH",
                                                                      #qvalueCutoff  = 0.05,
                                                                      pvalueCutoff  = 0.05
                                                                      )
                                                                         
      }
    }
  }
}

saveRDS(enrichGO_sp, file = paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Network_analysis/enrichGO_specific_", comp, ".rds")) 

#enrichGO_sp <- readRDS(paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Network_analysis/enrichGO_specific_", comp, ".rds"))


###Visualization :  dotplots

#Create dotplots, and then Save
for (algo in names(enrichGO_sp)){
  for (nam in names(enrichGO_sp[[algo]])){
    for (cell in names(enrichGO_sp[[algo]][[nam]])){
      for (centr in names(enrichGO_sp[[algo]][[nam]][[cell]])){
        
          tmp_sp <- enrichGO_sp[[algo]][[nam]][[cell]][[centr]] 
          
          if(dim(tmp_sp)[1]!=0){
            out_sp <- clusterProfiler::dotplot(tmp_sp, showCategory = 10, title = paste0("GO enrichment analysis (BP)"))
          }
          
          ggsave(
            paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Network_analysis/", algo, "/", nam, "/", cell, "/", "enrichGO_sp_", comp, ".dotplot.pdf"), 
            plot = out_sp, 
            device = "pdf", 
            height = 30, 
            width = 15, 
            units = "cm"
          )  
        
      }
    }
  }
}




###############################################################################
#enricher     : ORA using MSigDb

# Performs a MSigDb enrichment analysis using different gene sets for a given target gene set and universe.
# Args:  target: character vector with the gene symbols corresponding to the target set.
# Returns:
#        A data.frame with the MSigDb id, the p-value, the Odds score and the description of every enriched MSigDb term.

# MSigDB gene sets are divided into 9 collections:
# Hallmark gene sets (H) are coherently expressed signatures derived by aggregating many MSigDB gene sets to represent well-defined biological states or processes.
# Positional gene sets (C1) for each human chromosome and cytogenetic band.
# Curated gene sets (C2) from online pathway databases, publications in PubMed, and knowledge of domain experts.
# Regulatory target gene sets (C3) based on gene target predictions for microRNA seed sequences and predicted transcription factor binding sites.
# Computational gene sets (C4) defined by mining large collections of cancer-oriented microarray data.
# Ontology gene sets (C5) consist of genes annotated by the same ontology term.
# Oncogenic signature gene (C6) sets defined directly from microarray gene expression data from cancer gene perturbations.
# Immunologic signature gene sets (C7) defined directly from microarray gene expression data from immunologic studies.
# Cell type signature gene sets (C8) curated from cluster markers identified in single-cell sequencing studies of human tissue.




## 1. Get target genelist
#done in the previous steps


## 2. Download gene sets 
# Download gene sets and convert to gmt --> takes several min. you can save the final object for future use!
msig <- c("C2", "C4", "C5", "C6", "C8", "H") # choose your pathways. better more than less, you can also include all of them and omit from results what is not interesting.

msigdb.gmt <- sapply(msig,function(x) NULL)

for (i in msig) {
  msigdb.gmt[[i]] <- msigdbr(species = "Homo sapiens", category = i) %>% dplyr::select(gs_name, gene_symbol)
}

saveRDS(msigdb.gmt, file = paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Network_analysis/msigdb.gmt.Rds"))



## 3. Run enricher via compareCluster

enricher_sp <- list()
for (algo in names(top_centr_sp_all)){ 
  for (nam in names(top_centr_sp_all[[algo]])){                                                    
    for (cell in names(top_centr_sp_all[[algo]][[nam]])){ 
      for (centr in centrality_interest){                                                      #names(top_centr_sp_all[[algo]][[nam]][[cell]])
        for (msig in names(msigdb.gmt)){
          #group-specific hubs
          enricher_sp[[algo]][[nam]][[cell]][[centr]][[msig]] <- compareCluster(Gene~group, 
                                                                                data = df_target_sp[[algo]][[nam]][[cell]][[centr]], 
                                                                                fun = "enricher",
                                                                                TERM2GENE = msigdb.gmt[[msig]],
                                                                                pvalueCutoff  = 0.05
          )
          
        }
      }
    }
  }
}

saveRDS(enricher_sp, file = paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Network_analysis/enricher_specific_", comp, ".rds")) 

#enricher_sp <- readRDS(paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Network_analysis/enricher_specific_", comp, ".rds"))

###Visualization :  dotplots

#Create dotplots, and then Save
for (algo in names(enricher_sp)){
  for (nam in names(enricher_sp[[algo]])){
    for (cell in names(enricher_sp[[algo]][[nam]])){
      for (centr in names(enricher_sp[[algo]][[nam]][[cell]])){
        for (msig in names(enricher_sp[[algo]][[nam]][[cell]][[centr]])){
          
          tmp_sp <- enricher_sp[[algo]][[nam]][[cell]][[centr]][[msig]] 
          
          if(dim(tmp_sp)[1]!=0){
            out_sp <- clusterProfiler::dotplot(tmp_sp, showCategory = 10, title = paste0("MSigDb enrichment analysis performed using gene set ", msig))
          }
          
          ggsave(
            paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Network_analysis/", algo, "/", nam, "/", cell, "/", "enricher_sp_", msig, "_", comp, ".dotplot.pdf"),
            plot = out_sp, 
            device = "pdf", 
            height = 30, 
            width = 20, 
            units = "cm"
          ) 
          
        }
      }
    }
  }
}





###############################################################################
#enrichKEGG   : ORA using KEGG pathway 

# Performs a KEGG pathway enrichment analysis for a given target gene set and universe.
# Args:  target: character vector with the gene symbols corresponding to the target set.
# Returns:
#        A data.frame with the KEGG id, the p-value, the Odds score and the description of every enriched KEGG term.


## 1. Get target genelist

# For enrichKEGG, Input "target" data should a vector of entrez gene ids, so first do the conversion
tmp_target_sp <- list()
for (algo in names(top_centr_sp_all)){
  for (nam in names(top_centr_sp_all[[algo]])){                                                    
    for (cell in names(top_centr_sp_all[[algo]][[nam]])){
      for (centr in centrality_interest){                                                      #names(top_centr_sp_all[[algo]][[nam]][[cell]])
        for (org in names(top_centr_sp_all[[algo]][[nam]][[cell]][[centr]])){
          
          #group-specific hubs
          x <- top_centr_sp_all[[algo]][[nam]][[cell]][[centr]][[org]]
          
          tmp <- bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
          
          tmp_target_sp[[algo]][[nam]][[cell]][[centr]][[org]] <- tmp[,2]
          
        }
      }
    }
  }
}


# Then, Create dataframe to be able to visualize enrichment results across groups in one plot
df_target_sp <- list()
for (algo in names(tmp_target_sp)){ 
  for (nam in names(tmp_target_sp[[algo]])){                                                    
    for (cell in names(tmp_target_sp[[algo]][[nam]])){ 
      for (centr in centrality_interest){    
        mydf <- data.frame()
        for (org in names(tmp_target_sp[[algo]][[nam]][[cell]][[centr]])){
          l <- tmp_target_sp[[algo]][[nam]][[cell]][[centr]][[org]]
          if(length(l) != 0){
            tmp <- data.frame(Gene=l)
            tmp$group <- org
            mydf <- rbind(mydf, tmp)
          }
        }
        df_target_sp[[algo]][[nam]][[cell]][[centr]] <- mydf
      }
    }
  }
}



## 2. Run enrichKEGG
tmp_enrichKEGG_sp <- list()
for (algo in names(top_centr_sp_all)){
  for (nam in names(top_centr_sp_all[[algo]])){                                                    
    for (cell in names(top_centr_sp_all[[algo]][[nam]])){
      for (centr in centrality_interest){                                                      #names(top_centr_sp_all[[algo]][[nam]][[cell]])
          
        #group-specific hubs
        tmp_enrichKEGG_sp[[algo]][[nam]][[cell]][[centr]] <- compareCluster(Gene~group, 
                                                                            data=df_target_sp[[algo]][[nam]][[cell]][[centr]], 
                                                                            fun="enrichKEGG",
                                                                            organism     = 'hsa',
                                                                            keyType = "kegg",
                                                                            pvalueCutoff = 0.05
                                                                            )
          
      }
    }
  }
}


## The geneID column is ENTREZID, So Using setReadable, The geneID column is translated to symbol
enrichKEGG_sp <- list()
for (algo in names(tmp_enrichKEGG_sp)){
  for (nam in names(tmp_enrichKEGG_sp[[algo]])){                                                    
    for (cell in names(tmp_enrichKEGG_sp[[algo]][[nam]])){
      for (centr in centrality_interest){                                                      #names(top_centr_sp_all[[algo]][[nam]][[cell]])
          #group-specific hubs
          enrichKEGG_sp[[algo]][[nam]][[cell]][[centr]] <- setReadable(tmp_enrichKEGG_sp[[algo]][[nam]][[cell]][[centr]], 
                                                                              OrgDb = org.Hs.eg.db, keyType="ENTREZID")
      }
    }
  }
}


saveRDS(enrichKEGG_sp, file = paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Network_analysis/enrichKEGG_specific_", comp, ".rds"))  

#enrichKEGG_sp <- readRDS(paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Network_analysis/enrichKEGG_specific_", comp, ".rds"))

###Visualization :  dotplots

#Create dotplots, and then Save
for (algo in names(enrichKEGG_sp)){
  for (nam in names(enrichKEGG_sp[[algo]])){
    for (cell in names(enrichKEGG_sp[[algo]][[nam]])){
      for (centr in names(enrichKEGG_sp[[algo]][[nam]][[cell]])){
        
        tmp_sp <- enrichKEGG_sp[[algo]][[nam]][[cell]][[centr]] 
        
        if(dim(tmp_sp)[1]!=0){
          out_sp <- clusterProfiler::dotplot(tmp_sp, showCategory = 10, title = paste0("KEGG pathway enrichment analysis"))
        }
        
        ggsave(
          paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Network_analysis/", algo, "/", nam, "/", cell, "/", "enrichKEGG_sp_", comp, ".dotplot.pdf"), 
          plot = out_sp, 
          device = "pdf", 
          height = 30, 
          width = 15, 
          units = "cm"
        )  
        
      }
    }
  }
}








###









##################################################################################################
### COMPARISON : hep NORMAL vs SAMPLES
##################################################################################################

comp <- "comp_hep_NORMALvsSAMPLES"

#...


##################################################################################################
### COMPARISON : hep NORMAL vs EACH SAMPLE separately
##################################################################################################

samples <- c("P123s2", "P123s8", "P193", "P199", "P108", "P191", "P207", "P207Gel", "P215", "P220", "P275", "P282C", "P295")

# Compare normal vs each sample separately
for(sample in samples){
  comp <- paste0("comp_hep_NORMALvs", sample)
  
  
  ##################################################################################################
  ###Load data
  ##################################################################################################
  
  #top_centr_all <- readRDS(paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Network_analysis/networks_hub_lists.rds"))
  top_centr_ub_all <- readRDS(paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Network_analysis/networks_hub_ubiquitous_", comp, ".rds"))
  top_centr_sp_all <- readRDS(paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Network_analysis/networks_hub_specific_", comp, ".rds"))
  
  
  ##################################################################################################
  ## Perform Gene Ontology & KEGG pathway & MSigDb Enrichment Analysis
  ##################################################################################################
  
  
  # ...
  
  
  
}
  






