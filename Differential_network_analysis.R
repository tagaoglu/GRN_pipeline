# Set *.libPaths*  
.libPaths(c("~/p728/RSTUDIO/R/library/4.1/", .libPaths()))


#Set working directory and seed:
setwd("~/p728/RSTUDIO/analysis/tagaoglu")
set.seed(123)


##################################################################################################
###Set parameters
##################################################################################################

nam_dir <- c("expData_mean_1","expData_mean_2")

#Set accordingly:Directed = T  or  Directed = F
GENIE_edge_prediction <- T
PIDC_edge_prediction <- F


##################################################################################################
###Load data
##################################################################################################

GENIE <- readRDS(paste("~/p728/RSTUDIO/analysis/tagaoglu/data/networks/GENIE3/GENIE3_all_10k.Rds", sep = '' ))
PIDC <- readRDS(paste("~/p728/RSTUDIO/analysis/tagaoglu/data/networks/PIDC/PIDC_all_10k.Rds", sep = '' ))

#CoDiNA analysis does not work with GENIE_wRcis because of lack of weight of edges

#exclude random data
for(nam in names(GENIE)){
  GENIE <- GENIE[names(GENIE) %in% nam_dir]
  PIDC <- PIDC[names(PIDC) %in% nam_dir]
}

#Sort cell names of list objects to compare easily
for(nam in names(GENIE)){
  for(group in names(GENIE[[nam]])){
    GENIE[[nam]][[group]] <- GENIE[[nam]][[group]][order(names(GENIE[[nam]][[group]]))]
    PIDC[[nam]][[group]] <- PIDC[[nam]][[group]][order(names(PIDC[[nam]][[group]]))]
  }
}

##################################################################################################


##Create files if it does not exist
plotdir=paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/CoDiNA")
if(!dir.exists(plotdir)){dir.create(plotdir,  recursive = T)}


#####################################################
#Creating the Differential Network (CoDiNA)
#####################################################

#CoDiNA: Co-Expression Differential Network Analysis
#Categorize links and nodes from multiple networks in 3 categories: 
#Common links (alpha) specific links (gamma), and different links (beta). 
#Also categorizes the links into sub-categories and groups. The package includes a visualization tool for the networks.

#INPUT: Input data for the CoDiNA R package can be any networks, filtered for containing only significant links (according to the network construction method used). Edge list is a list containing all the links and their weights. The user can assign a weight of zero to links for which the p-value is not significant. 
#STEPS: The function MakeDiffNet() classifies the links into the Φ and  categories, calculates and normalises the scores. Its output is used as input for assigning the nodes into categories by the function ClusterNodes(). The plot() function can be used on the output from MakeDiffNet() and automatically calls the function ClusterNodes().


#remove.packages('CoDiNA', lib="/home/tagaoglu/p728/RSTUDIO/R/library/4.1")
#install.packages('CoDiNA', version = "1.1.1")
#packageVersion('CoDiNA')

library(CoDiNA)

#install.packages("visNetwork")
library(visNetwork)
library(igraph)


########
#GENIE

# Differential Network Analysis for 3, "normal", "viral", "nonviral"
GENIE_DiffNet_Comp3 <- list()
for(nam in names(GENIE)){
  for(cell in names(GENIE[[nam]][["normal"]])){
    GENIE_DiffNet_Comp3[[nam]][[cell]] = MakeDiffNet(Data = list(GENIE[[nam]][["normal"]][[cell]],
                                                                 GENIE[[nam]][["viral"]][[cell]], 
                                                                 GENIE[[nam]][["nonviral"]][[cell]]), 
                                                     Code = c("normal", "viral", "nonviral"))
  }
}

# Differential Network Analysis for 2, "normal", "tumor"
GENIE_DiffNet_Comp2 <- list()
for(nam in names(GENIE)){
  for(cell in names(GENIE[[nam]][["normal"]])){
    GENIE_DiffNet_Comp2[[nam]][[cell]] = MakeDiffNet(Data = list(GENIE[[nam]][["normal"]][[cell]],
                                                                 GENIE[[nam]][["tumor"]][[cell]]), 
                                                     Code = c("normal", "tumor"))
  }
}    

GENIE_DiffNet <-list()
GENIE_DiffNet[["Comp3"]] <- GENIE_DiffNet_Comp3
rm(GENIE_DiffNet_Comp3)
GENIE_DiffNet[["Comp2"]] <- GENIE_DiffNet_Comp2
rm(GENIE_DiffNet_Comp2)


########
#PIDC

# Differential Network Analysis for 3, "normal", "viral", "nonviral"
PIDC_DiffNet_Comp3 <- list()
for(nam in names(PIDC)){
  for(cell in names(PIDC[[nam]][["normal"]])){
    PIDC_DiffNet_Comp3[[nam]][[cell]] = MakeDiffNet(Data = list(PIDC[[nam]][["normal"]][[cell]],
                                                                PIDC[[nam]][["viral"]][[cell]], 
                                                                PIDC[[nam]][["nonviral"]][[cell]]), 
                                                    Code = c("normal", "viral", "nonviral"))
  }
}

# Differential Network Analysis for 2, "normal", "tumor"
PIDC_DiffNet_Comp2 <- list()
for(nam in names(PIDC)){
  for(cell in names(PIDC[[nam]][["normal"]])){
    PIDC_DiffNet_Comp2[[nam]][[cell]] = MakeDiffNet(Data = list(PIDC[[nam]][["normal"]][[cell]],
                                                                PIDC[[nam]][["tumor"]][[cell]]), 
                                                    Code = c("normal", "tumor"))
  }
} 


PIDC_DiffNet <-list()
PIDC_DiffNet[["Comp3"]] <- PIDC_DiffNet_Comp3
rm(PIDC_DiffNet_Comp3)
PIDC_DiffNet[["Comp2"]] <- PIDC_DiffNet_Comp2
rm(PIDC_DiffNet_Comp2)


#
#head(print(GENIE_DiffNet_Comp3[[nam]][[cell]]))
#summary(GENIE_DiffNet_Comp3[[nam]][[cell]])
#summary(GENIE_DiffNet_Comp3[[nam]][[cell]])$Phi_tilda


# Save the data for CoDiNA
saveRDS(GENIE_DiffNet, file = paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/CoDiNA/GENIE3_all_10k_DiffNet.Rds"))  
saveRDS(PIDC_DiffNet, file = paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/CoDiNA/PIDC_all_10k_DiffNet.Rds"))  


# Load the data
#GENIE_DiffNet <- readRDS(paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/CoDiNA/GENIE3_all_10k_DiffNet.Rds"))
#PIDC_DiffNet <- readRDS(paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/CoDiNA/PIDC_all_10k_DiffNet.Rds"))



#####################################################
#Clustering the nodes into Φ and Φ̃ categories
#####################################################

#The suggested values for the internal and external cutoffs are the median or the first and third quantiles of the internal and Φ̃ scores, depending on how conservative the network should be.

#cutoff.external : The cut-off between the clusters (delta from the center to the edge coordinates), the closer to 1, the better.
#cutoff.internal : The cut-off inside the clusters (delta from the theoretical cluster to the edge coordinates), the closer to zero, the better.

#GENIE
GENIE_int_C <- list()
GENIE_ext_C <- list()
GENIE_Nodes_Groups <- list()

for (comp in names(GENIE_DiffNet)){
  for (nam in names(GENIE_DiffNet[[comp]])){
    for (cell in names(GENIE_DiffNet[[comp]][[nam]])){
      #Using the median:
      GENIE_int_C[[comp]][[nam]][[cell]] = quantile(GENIE_DiffNet[[comp]][[nam]][[cell]]$Score_internal, 0.5)
      GENIE_ext_C[[comp]][[nam]][[cell]] = quantile(GENIE_DiffNet[[comp]][[nam]][[cell]]$Score_Phi, 0.5)
      #or
      # #Using the first and third quantile:
      # GENIE_int_C[[comp]][[nam]][[cell]] = quantile(GENIE_DiffNet[[comp]][[nam]][[cell]]$Score_internal, 0.25)
      # GENIE_ext_C[[comp]][[nam]][[cell]] = quantile(GENIE_DiffNet[[comp]][[nam]][[cell]]$Score_Phi, 0.75)
      
      GENIE_Nodes_Groups[[comp]][[nam]][[cell]] = ClusterNodes(DiffNet = GENIE_DiffNet[[comp]][[nam]][[cell]], 
                                                               cutoff.external = GENIE_ext_C[[comp]][[nam]][[cell]], 
                                                               cutoff.internal = GENIE_int_C[[comp]][[nam]][[cell]])
    }
  }
}
  
  
  


#PIDC
PIDC_int_C <- list()
PIDC_ext_C <- list()
PIDC_Nodes_Groups <- list()

for (comp in names(PIDC_DiffNet)){
  for (nam in names(PIDC_DiffNet[[comp]])){
    for (cell in names(PIDC_DiffNet[[comp]][[nam]])){
      #Using the median:
      PIDC_int_C[[comp]][[nam]][[cell]] = quantile(PIDC_DiffNet[[comp]][[nam]][[cell]]$Score_internal, 0.5)
      PIDC_ext_C[[comp]][[nam]][[cell]] = quantile(PIDC_DiffNet[[comp]][[nam]][[cell]]$Score_Phi, 0.5)
      #or
      # #Using the first and third quantile:
      # PIDC_int_C[[comp]][[nam]][[cell]] = quantile(PIDC_DiffNet[[comp]][[nam]][[cell]]$Score_internal, 0.25)
      # PIDC_ext_C[[comp]][[nam]][[cell]] = quantile(PIDC_DiffNet[[comp]][[nam]][[cell]]$Score_Phi, 0.75)
      
      PIDC_Nodes_Groups[[comp]][[nam]][[cell]] = ClusterNodes(DiffNet = PIDC_DiffNet[[comp]][[nam]][[cell]], 
                                                              cutoff.external = PIDC_ext_C[[comp]][[nam]][[cell]], 
                                                              cutoff.internal = PIDC_int_C[[comp]][[nam]][[cell]])
    }
  }
}


#
#table(GENIE_Nodes_Groups[[comp]][[nam]][[cell]]$Phi_tilde)

#subset(GENIE_Nodes_Groups[[comp]][[nam]][[cell]], Phi_tilde=="g.nonviral")
#subset(GENIE_Nodes_Groups[[comp]][[nam]][[cell]], Phi_tilde=="a")



# Save the data for CoDiNA
saveRDS(GENIE_Nodes_Groups, file = paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/CoDiNA/GENIE3_all_10k_Nodes_Groups.Rds"))  
saveRDS(PIDC_Nodes_Groups, file = paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/CoDiNA/PIDC_all_10k_Groups.Rds"))  


# Load the data per sample
#GENIE_Nodes_Groups <- readRDS(paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/CoDiNA/GENIE3_all_10k_Nodes_Groups.Rds"))
#PIDC_Nodes_Groups <- readRDS(paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/CoDiNA/PIDC_all_10k_Groups.Rds"))




#####################################################
#Plotting the network
#####################################################

#The visualization of the final network can be quickly done with plot. The layout of the network can be also determined from a variety that is implemented in igraph package, the Make_Cluster argument allows the nodes to be clusterized according to many clustering algorithms that are implemented in igraph can be used. The final graph can be exported as an HTML or as png. The argument path saves the network in the given path.

#The plot returns the nodes and its information.


##Create files if it does not exist
for (comp in names(GENIE_DiffNet)){
  for (nam in names(GENIE_DiffNet[[comp]])){
    plotdir_GENIE=paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/CoDiNA/plots/Networks_html/", "GENIE/", comp, "/", nam)
  if(!dir.exists(plotdir_GENIE)){dir.create(plotdir_GENIE,  recursive = T)}
  }
}

for (comp in names(PIDC_DiffNet)){
  for (nam in names(PIDC_DiffNet[[comp]])){
    plotdir_PIDC=paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/CoDiNA/plots/Networks_html/", "PIDC/", comp, "/", nam)
    if(!dir.exists(plotdir_PIDC)){dir.create(plotdir_PIDC,  recursive = T)}
  }
}


#######
#GENIE

#Plot (html)
GENIE_Graph <- list()
for (comp in names(GENIE_DiffNet)){
  for (nam in names(GENIE_DiffNet[[comp]])){
    for (cell in names(GENIE_DiffNet[[comp]][[nam]])){
      GENIE_Graph[[comp]][[nam]][[cell]] = plot(GENIE_DiffNet[[comp]][[nam]][[cell]], 
                                                cutoff.external = GENIE_ext_C[[comp]][[nam]][[cell]], 
                                                cutoff.internal = GENIE_int_C[[comp]][[nam]][[cell]], 
                                                layout = 'layout_components', 
                                                path = paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/CoDiNA/plots/Networks_html/", "GENIE/", comp, "/", nam, "/", cell, ".html"))
    }
  }
}


#Plot (visNetwork)
GENIE_vis <- list()
for (comp in names(GENIE_Graph)){
  for (nam in names(GENIE_Graph[[comp]])){
    for (cell in names(GENIE_Graph[[comp]][[nam]])){
      net_fig<- GENIE_Graph[[comp]][[nam]][[cell]]
      
      e = net_fig$Edges 
      names(e) = c("from", "to", "group", "Phi", "width")
      e$weight = e$width
      
      n = net_fig$Nodes
      n$Phi_tilde %>% as.factor()
      n$group = n$Phi_tilde
      # n$color = RColorBrewer::brewer.pal(nlevels(n$Phi_tilde), "Pastel2")[unclass(n$Phi_tilde)]
      # n$frame.color = RColorBrewer::brewer.pal(nlevels(n$Phi_tilde), "Set2")[unclass(n$Phi_tilde)]
      # n$label.color = RColorBrewer::brewer.pal(nlevels(n$Phi_tilde), "Dark2")[unclass(n$Phi_tilde)]
      
      gDis = graph_from_data_frame(e, vertices = n, directed = F)
      V(gDis)$size = strength(gDis) %>% 
        CoDiNA::normalize()
      V(gDis)$size = (V(gDis)$size + 0.1)*30
      V(gDis)$label = ifelse(V(gDis)$size  > 4, V(gDis)$name, NA )
      V(gDis)$value = V(gDis)$size
      E(gDis)$weight = E(gDis)$width * 100
      
      E(gDis)$length = E(gDis)$width
      x = visIgraph(gDis)
      
      GENIE_vis[[comp]][[nam]][[cell]] <- visNetwork(nodes = x$x$nodes, edges = x$x$edges) %>%
        visPhysics(enabled = F) %>% 
        visIgraphLayout(layout = "layout_with_fr") %>%
        visLegend() %>% 
        visExport(type = "png", name = paste0("GENIE"), 
                  float = "left", label = "Save network", background = "white", style= "") 
      
      visSave(GENIE_vis[[comp]][[nam]][[cell]], 
              file = paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/CoDiNA/plots/Networks_html/", "GENIE/", comp, "/", nam, "/", cell, "_visNetwork.html"))
    }
  }
}




#######
#PIDC

#Plot (html)
PIDC_Graph <- list()
for (comp in names(PIDC_DiffNet)){
  for (nam in names(PIDC_DiffNet[[comp]])){
    for (cell in names(PIDC_DiffNet[[comp]][[nam]])){
      PIDC_Graph[[comp]][[nam]][[cell]] = plot(PIDC_DiffNet[[comp]][[nam]][[cell]], 
                                               cutoff.external = PIDC_ext_C[[comp]][[nam]][[cell]], 
                                               cutoff.internal = PIDC_int_C[[comp]][[nam]][[cell]], 
                                               layout = 'layout_components', 
                                               path = paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/CoDiNA/plots/Networks_html/", "PIDC/", comp, "/", nam, "/", cell, ".html"))
    }
  }
}


#Plot (visNetwork)
PIDC_vis <- list()
for (comp in names(PIDC_Graph)){
  for (nam in names(PIDC_Graph[[comp]])){
    for (cell in names(PIDC_Graph[[comp]][[nam]])){
      net_fig<- PIDC_Graph[[comp]][[nam]][[cell]]
      
      e = net_fig$Edges 
      names(e) = c("from", "to", "group", "Phi", "width")
      e$weight = e$width
      
      n = net_fig$Nodes
      n$Phi_tilde %>% as.factor()
      n$group = n$Phi_tilde
      # n$color = RColorBrewer::brewer.pal(nlevels(n$Phi_tilde), "Pastel2")[unclass(n$Phi_tilde)]
      # n$frame.color = RColorBrewer::brewer.pal(nlevels(n$Phi_tilde), "Set2")[unclass(n$Phi_tilde)]
      # n$label.color = RColorBrewer::brewer.pal(nlevels(n$Phi_tilde), "Dark2")[unclass(n$Phi_tilde)]
      
      gDis = graph_from_data_frame(e, vertices = n, directed = F)
      V(gDis)$size = strength(gDis) %>% 
        CoDiNA::normalize()
      V(gDis)$size = (V(gDis)$size + 0.1)*30
      V(gDis)$label = ifelse(V(gDis)$size  > 4, V(gDis)$name, NA )
      V(gDis)$value = V(gDis)$size
      E(gDis)$weight = E(gDis)$width * 100
      
      E(gDis)$length = E(gDis)$width
      x = visIgraph(gDis)
      
      PIDC_vis[[comp]][[nam]][[cell]] <- visNetwork(nodes = x$x$nodes, edges = x$x$edges) %>%
        visPhysics(enabled = F) %>% 
        visIgraphLayout(layout = "layout_with_fr") %>%
        visLegend() %>% 
        visExport(type = "png", name = paste0("PIDC"), 
                  float = "left", label = "Save network", background = "white", style= "") 
      
      visSave(PIDC_vis[[comp]][[nam]][[cell]], 
              file = paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/CoDiNA/plots/Networks_html/", "PIDC/", comp, "/", nam, "/", cell, "_visNetwork.html"))
    }
  }
}






#####################################################
#Create table
#####################################################

library(gt)
library(webshot)
library(scales)
library(readr)

########
#GENIE

GENIE_df <- list()
for (comp in names(GENIE_Nodes_Groups)){
  for (nam in names(GENIE_Nodes_Groups[[comp]])){
    for (cell in names(GENIE_Nodes_Groups[[comp]][[nam]])){
      df <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("Group", "Count", "Genes"))
      for (i in 1:length(names(table(GENIE_Nodes_Groups[[comp]][[nam]][[cell]]$Phi_tilde)))) {
        Group <- names(table(GENIE_Nodes_Groups[[comp]][[nam]][[cell]]$Phi_tilde))[i]
        Count <- as.numeric(table(GENIE_Nodes_Groups[[comp]][[nam]][[cell]]$Phi_tilde))[i]
        Gene_names <- subset(GENIE_Nodes_Groups[[comp]][[nam]][[cell]], Phi_tilde==Group)$Node
        Genes <- paste(Gene_names, collapse=', ')
        res= c(Group, Count, Genes)
        df[i,]=res
      }
      GENIE_df[[comp]][[nam]][[cell]] <- df
    }
  }
}


for (comp in names(GENIE_df)){
  for (nam in names(GENIE_df[[comp]])){
    for (cell in names(GENIE_df[[comp]][[nam]])){
      GENIE_df[[comp]][[nam]][[cell]] %>% 
        gt() %>%
        tab_header(title = cell) %>%
        gtsave(paste0(cell, "_SummaryTable.html"), 
               path = paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/CoDiNA/plots/Networks_html/", "GENIE/", comp, "/", nam, "/"))
    }
  }
}



########
#PIDC

PIDC_df <- list()
for (comp in names(PIDC_Nodes_Groups)){
  for (nam in names(PIDC_Nodes_Groups[[comp]])){
    for (cell in names(PIDC_Nodes_Groups[[comp]][[nam]])){
      df <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("Group", "Count", "Genes"))
      for (i in 1:length(names(table(PIDC_Nodes_Groups[[comp]][[nam]][[cell]]$Phi_tilde)))) {
        Group <- names(table(PIDC_Nodes_Groups[[comp]][[nam]][[cell]]$Phi_tilde))[i]
        Count <- as.numeric(table(PIDC_Nodes_Groups[[comp]][[nam]][[cell]]$Phi_tilde))[i]
        Gene_names <- subset(PIDC_Nodes_Groups[[comp]][[nam]][[cell]], Phi_tilde==Group)$Node
        Genes <- paste(Gene_names, collapse=', ')
        res= c(Group, Count, Genes)
        df[i,]=res
      }
      PIDC_df[[comp]][[nam]][[cell]] <- df
    }
  }
}

for (comp in names(PIDC_df)){
  for (nam in names(PIDC_df[[comp]])){
    for (cell in names(PIDC_df[[comp]][[nam]])){
      PIDC_df[[comp]][[nam]][[cell]] %>% 
        gt() %>%
        tab_header(title = cell) %>%
        gtsave(paste0(cell, "_SummaryTable.html"), 
               path = paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/CoDiNA/plots/Networks_html/", "PIDC/", comp, "/", nam, "/"))
    }
  }
}






##################################################################################################
## Perform Gene Ontology & KEGG pathway & MSigDb Enrichment Analysis
##################################################################################################

##################################################################################################
## for classified NODES
##################################################################################################


##Create files if it does not exist
for (comp in names(GENIE_Nodes_Groups)){
  for (nam in names(GENIE_Nodes_Groups[[comp]])){
    plotdir_GENIE=paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/CoDiNA/plots/ORA_plots/", "GENIE/", comp, "/", nam)
    if(!dir.exists(plotdir_GENIE)){dir.create(plotdir_GENIE,  recursive = T)}
  }
}

for (comp in names(PIDC_Nodes_Groups)){
  for (nam in names(PIDC_Nodes_Groups[[comp]])){
    plotdir_PIDC=paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/CoDiNA/plots/ORA_plots/", "PIDC/", comp, "/", nam)
    if(!dir.exists(plotdir_PIDC)){dir.create(plotdir_PIDC,  recursive = T)}
  }
}


#Load libraries
library(org.Hs.eg.db)
library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(ggplot2)



#Put all in one list
algo_list <- c("GENIE", "PIDC")

Nodes_Groups <- sapply(algo_list,function(x) NULL)
Nodes_Groups[["GENIE"]] <- GENIE_Nodes_Groups
Nodes_Groups[["PIDC"]] <- PIDC_Nodes_Groups


## 1. Get target genelist

#First subset correspondingly
Nodes_Groups_sub <- list()
for (algo in names(Nodes_Groups)){
  for (comp in names(Nodes_Groups[[algo]])){
    for (nam in names(Nodes_Groups[[algo]][[comp]])){
      for (cell in names(Nodes_Groups[[algo]][[comp]][[nam]])){
        if(comp=="Comp3"){
          groups <- c("normal", "viral", "nonviral")
          for(group in groups){
            Nodes_Groups_sub[[algo]][[comp]][[nam]][[cell]][[group]] <- subset(Nodes_Groups[[algo]][[comp]][[nam]][[cell]], Phi_tilde==paste0("g.", group))$Node
          }
        }else if(comp=="Comp2") {
          groups <- c("normal", "tumor")
          for(group in groups){
            Nodes_Groups_sub[[algo]][[comp]][[nam]][[cell]][[group]] <- subset(Nodes_Groups[[algo]][[comp]][[nam]][[cell]], Phi_tilde==paste0("g.", group))$Node
          }
        }
      }
    }
  }
}
          
#Then, Create dataframe to be able to visualize enrichment results across groups in one plot
df_target_sp <- list()
for (algo in names(Nodes_Groups_sub)){
  for (comp in names(Nodes_Groups_sub[[algo]])){
    for (nam in names(Nodes_Groups_sub[[algo]][[comp]])){
      for (cell in names(Nodes_Groups_sub[[algo]][[comp]][[nam]])){
        mydf <- data.frame()
        for (org in names(Nodes_Groups_sub[[algo]][[comp]][[nam]][[cell]])){
          l <- Nodes_Groups_sub[[algo]][[comp]][[nam]][[cell]][[org]]
          if(length(l) != 0){
            tmp <- data.frame(Gene=l)
            tmp$group <- org
            mydf <- rbind(mydf, tmp)
          }
        }
        df_target_sp[[algo]][[comp]][[nam]][[cell]] <- mydf
      }
    }
  }
}


## 2. Run enrichGO via compareCluster
enrichGO_sp <- list()
for (algo in names(Nodes_Groups_sub)){
  for (comp in names(Nodes_Groups_sub[[algo]])){
    for (nam in names(Nodes_Groups_sub[[algo]][[comp]])){
      for (cell in names(Nodes_Groups_sub[[algo]][[comp]][[nam]])){
        #group-specific hubs
        tmp <- df_target_sp[[algo]][[comp]][[nam]][[cell]]
        
        if(length(tmp) != 0){
          enrichGO_sp[[algo]][[comp]][[nam]][[cell]] <- compareCluster(Gene~group, 
                                                                       data=tmp, 
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
}


saveRDS(enrichGO_sp, file = paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/CoDiNA/enrichGO_specific.rds")) 

#enrichGO_sp = readRDS("~/p728/RSTUDIO/analysis/tagaoglu/data/CoDiNA/enrichGO_specific.rds")

###Visualization :  dotplots

#Create dotplots, and then Save

for (algo in names(enrichGO_sp)){
  for (comp in names(enrichGO_sp[[algo]])){
    for (nam in names(enrichGO_sp[[algo]][[comp]])){
      for (cell in names(enrichGO_sp[[algo]][[comp]][[nam]])){
        tmp_sp <- enrichGO_sp[[algo]][[comp]][[nam]][[cell]] 
        
        if(dim(tmp_sp)[1]!=0){
          out_sp <- clusterProfiler::dotplot(tmp_sp, showCategory = 10, title = paste0("GO enrichment analysis (BP)"))
        }
        
        ggsave(
          paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/CoDiNA/plots/ORA_plots/", algo, "/", comp, "/", nam, "/", cell, "_enrichGO_sp.dotplot.pdf"), 
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
#done in the previous steps for enrichGO


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
for (algo in names(Nodes_Groups_sub)){
  for (comp in names(Nodes_Groups_sub[[algo]])){
    for (nam in names(Nodes_Groups_sub[[algo]][[comp]])){
      for (cell in names(Nodes_Groups_sub[[algo]][[comp]][[nam]])){
        for (msig in names(msigdb.gmt)){
          #group-specific hubs
          tmp <- df_target_sp[[algo]][[comp]][[nam]][[cell]]
          
          if(length(tmp) != 0){
            enricher_sp[[algo]][[comp]][[nam]][[cell]][[msig]] <- compareCluster(Gene~group, 
                                                                                 data=tmp, 
                                                                                 fun = "enricher",
                                                                                 TERM2GENE = msigdb.gmt[[msig]],
                                                                                 pvalueCutoff  = 0.05
            )
         }
        }
      }
    }
  }
}

saveRDS(enricher_sp, file = paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/CoDiNA/enricher_specific.rds")) 

#enricher_sp = readRDS("~/p728/RSTUDIO/analysis/tagaoglu/data/CoDiNA/enricher_specific.rds")

###Visualization :  dotplots

#Create dotplots, and then Save

for (algo in names(enricher_sp)){
  for (comp in names(enricher_sp[[algo]])){
    for (nam in names(enricher_sp[[algo]][[comp]])){
      for (cell in names(enricher_sp[[algo]][[comp]][[nam]])){
        for(msig in names(enricher_sp[[algo]][[comp]][[nam]][[cell]])){
          tmp_sp <- enricher_sp[[algo]][[comp]][[nam]][[cell]][[msig]]
          
          if(dim(tmp_sp)[1]!=0){
            out_sp <- clusterProfiler::dotplot(tmp_sp, showCategory = 10, title = paste0("MSigDb enrichment analysis performed using gene set ", msig))
          }
          
          ggsave(
            paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/CoDiNA/plots/ORA_plots/", algo, "/", comp, "/", nam, "/", cell, "_enricher_sp_", msig, ".dotplot.pdf"), 
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

# ...

## 1. Get target genelist

# For enrichKEGG, Input "target" data should a vector of entrez gene ids, so first do the conversion
tmp_target_sp <- list()
for (algo in names(Nodes_Groups_sub)){
  for (comp in names(Nodes_Groups_sub[[algo]])){
    for (nam in names(Nodes_Groups_sub[[algo]][[comp]])){
      for (cell in names(Nodes_Groups_sub[[algo]][[comp]][[nam]])){
        for (org in names(Nodes_Groups_sub[[algo]][[comp]][[nam]][[cell]])){
          #group-specific hubs
          x <- Nodes_Groups_sub[[algo]][[comp]][[nam]][[cell]][[org]]
          
          tmp <- bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
          
          tmp_target_sp[[algo]][[comp]][[nam]][[cell]][[org]] <- tmp[,2]
          
        }
      }
    }
  }
}


# Then, Create dataframe to be able to visualize enrichment results across groups in one plot
#Then, Create dataframe to be able to visualize enrichment results across groups in one plot
df_target_sp <- list()
for (algo in names(tmp_target_sp)){
  for (comp in names(tmp_target_sp[[algo]])){
    for (nam in names(tmp_target_sp[[algo]][[comp]])){
      for (cell in names(tmp_target_sp[[algo]][[comp]][[nam]])){
        mydf <- data.frame()
        for (org in names(tmp_target_sp[[algo]][[comp]][[nam]][[cell]])){
          l <- tmp_target_sp[[algo]][[comp]][[nam]][[cell]][[org]]
          if(length(l) != 0){
            tmp <- data.frame(Gene=l)
            tmp$group <- org
            mydf <- rbind(mydf, tmp)
          }
        }
        df_target_sp[[algo]][[comp]][[nam]][[cell]] <- mydf
      }
    }
  }
}



## 2. Run enrichKEGG via compareCluster
tmp_enrichKEGG_sp <- list()
for (algo in names(Nodes_Groups_sub)){
  for (comp in names(Nodes_Groups_sub[[algo]])){
    for (nam in names(Nodes_Groups_sub[[algo]][[comp]])){
      for (cell in names(Nodes_Groups_sub[[algo]][[comp]][[nam]])){
        #group-specific hubs
        tmp <- df_target_sp[[algo]][[comp]][[nam]][[cell]]
        
        if(length(tmp) != 0){
          tmp_enrichKEGG_sp[[algo]][[comp]][[nam]][[cell]] <- compareCluster(Gene~group, 
                                                                         data=tmp,
                                                                         fun="enrichKEGG",
                                                                         organism  = 'hsa',
                                                                         keyType = "kegg",
                                                                         pvalueCutoff = 0.05
          )
        }
      }
    }
  }
}


## The geneID column is ENTREZID, So Using setReadable, The geneID column is translated to symbol
enrichKEGG_sp <- list()
for (algo in names(tmp_enrichKEGG_sp)){
  for (comp in names(tmp_enrichKEGG_sp[[algo]])){                                                    
    for (nam in names(tmp_enrichKEGG_sp[[algo]][[comp]])){ 
      for (cell in names(tmp_enrichKEGG_sp[[algo]][[comp]][[nam]])){ 
        #group-specific hubs
        enrichKEGG_sp[[algo]][[comp]][[nam]][[cell]] <- setReadable(tmp_enrichKEGG_sp[[algo]][[comp]][[nam]][[cell]], 
                                                                    OrgDb = org.Hs.eg.db, keyType="ENTREZID")
      }
    }
  }
}

saveRDS(enrichKEGG_sp, file = paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/CoDiNA/enrichKEGG_specific.rds")) 

#enrichKEGG_sp = readRDS("~/p728/RSTUDIO/analysis/tagaoglu/data/CoDiNA/enrichKEGG_specific.rds")

###Visualization :  dotplots

#Create dotplots, and then Save

for (algo in names(enrichKEGG_sp)){
  for (comp in names(enrichKEGG_sp[[algo]])){
    for (nam in names(enrichKEGG_sp[[algo]][[comp]])){
      for (cell in names(enrichKEGG_sp[[algo]][[comp]][[nam]])){
        tmp_sp <- enrichKEGG_sp[[algo]][[comp]][[nam]][[cell]] 
        
        if(dim(tmp_sp)[1]!=0){
          out_sp <- clusterProfiler::dotplot(tmp_sp, showCategory = 10, title = paste0("KEGG pathway enrichment analysis"))
        }
        
        ggsave(
          paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/CoDiNA/plots/ORA_plots/", algo, "/", comp, "/", nam, "/", cell, "_enrichKEGG_sp.dotplot.pdf"), 
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







##################################################################################################
## Perform Gene Ontology & KEGG pathway & MSigDb Enrichment Analysis
##################################################################################################

##################################################################################################
## for classified INTERACTIONS
##################################################################################################

##Create files if it does not exist
for (comp in names(GENIE_Nodes_Groups)){
  for (nam in names(GENIE_Nodes_Groups[[comp]])){
    plotdir_GENIE=paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/CoDiNA/plots/ORA_plots_forINTERACTIONS/", "GENIE/", comp, "/", nam)
    if(!dir.exists(plotdir_GENIE)){dir.create(plotdir_GENIE,  recursive = T)}
  }
}

for (comp in names(PIDC_Nodes_Groups)){
  for (nam in names(PIDC_Nodes_Groups[[comp]])){
    plotdir_PIDC=paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/CoDiNA/plots/ORA_plots_forINTERACTIONS/", "PIDC/", comp, "/", nam)
    if(!dir.exists(plotdir_PIDC)){dir.create(plotdir_PIDC,  recursive = T)}
  }
}


#Load libraries
library(org.Hs.eg.db)
library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(ggplot2)



#Put all in one list
algo_list <- c("GENIE", "PIDC")

DiffNet <- sapply(algo_list,function(x) NULL)
DiffNet[["GENIE"]] <- GENIE_DiffNet
DiffNet[["PIDC"]] <- PIDC_DiffNet


## 1. Get target genelist

#First subset correspondingly
DiffNet_sub <- list()
for (algo in names(DiffNet)){
  for (comp in names(DiffNet[[algo]])){
    for (nam in names(DiffNet[[algo]][[comp]])){
      for (cell in names(DiffNet[[algo]][[comp]][[nam]])){
        if(comp=="Comp3"){
          groups <- c("normal", "viral", "nonviral")
          for(group in groups){
            Node1_tmp <- subset(DiffNet[[algo]][[comp]][[nam]][[cell]], Phi_tilde==paste0("g.", group))$Node.1
            Node2_tmp <- subset(DiffNet[[algo]][[comp]][[nam]][[cell]], Phi_tilde==paste0("g.", group))$Node.2
            DiffNet_sub[[algo]][[comp]][[nam]][[cell]][[group]] <- unique(c(Node1_tmp,Node2_tmp))
          }
        }else if(comp=="Comp2") {
          groups <- c("normal", "tumor")
          for(group in groups){
            Node1_tmp <- subset(DiffNet[[algo]][[comp]][[nam]][[cell]], Phi_tilde==paste0("g.", group))$Node.1
            Node2_tmp <- subset(DiffNet[[algo]][[comp]][[nam]][[cell]], Phi_tilde==paste0("g.", group))$Node.2
            DiffNet_sub[[algo]][[comp]][[nam]][[cell]][[group]] <- unique(c(Node1_tmp,Node2_tmp))
          }
        }
      }
    }
  }
}

#Then, Create dataframe to be able to visualize enrichment results across groups in one plot
df_target_sp <- list()
for (algo in names(DiffNet_sub)){
  for (comp in names(DiffNet_sub[[algo]])){
    for (nam in names(DiffNet_sub[[algo]][[comp]])){
      for (cell in names(DiffNet_sub[[algo]][[comp]][[nam]])){
        mydf <- data.frame()
        for (org in names(DiffNet_sub[[algo]][[comp]][[nam]][[cell]])){
          l <- DiffNet_sub[[algo]][[comp]][[nam]][[cell]][[org]]
          if(length(l) != 0){
            tmp <- data.frame(Gene=l)
            tmp$group <- org
            mydf <- rbind(mydf, tmp)
          }
        }
        df_target_sp[[algo]][[comp]][[nam]][[cell]] <- mydf
      }
    }
  }
}


## 2. Run enrichGO via compareCluster
enrichGO_sp <- list()
for (algo in names(DiffNet_sub)){
  for (comp in names(DiffNet_sub[[algo]])){
    for (nam in names(DiffNet_sub[[algo]][[comp]])){
      for (cell in names(DiffNet_sub[[algo]][[comp]][[nam]])){
        #group-specific hubs
        tmp <- df_target_sp[[algo]][[comp]][[nam]][[cell]]
        
        if(length(tmp) != 0){
          enrichGO_sp[[algo]][[comp]][[nam]][[cell]] <- compareCluster(Gene~group, 
                                                                       data=tmp, 
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
}


saveRDS(enrichGO_sp, file = paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/CoDiNA/enrichGO_specific_forINTERACTIONS.rds")) 

#enrichGO_sp = readRDS("~/p728/RSTUDIO/analysis/tagaoglu/data/CoDiNA/enrichGO_specific_forINTERACTIONS.rds")

###Visualization :  dotplots

#Create dotplots, and then Save

for (algo in names(enrichGO_sp)){
  for (comp in names(enrichGO_sp[[algo]])){
    for (nam in names(enrichGO_sp[[algo]][[comp]])){
      for (cell in names(enrichGO_sp[[algo]][[comp]][[nam]])){
        tmp_sp <- enrichGO_sp[[algo]][[comp]][[nam]][[cell]] 
        
        if(dim(tmp_sp)[1]!=0){
          out_sp <- clusterProfiler::dotplot(tmp_sp, showCategory = 10, title = paste0("GO enrichment analysis (BP)"))
        }
        
        ggsave(
          paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/CoDiNA/plots/ORA_plots_forINTERACTIONS/", algo, "/", comp, "/", nam, "/", cell, "_enrichGO_sp.dotplot.pdf"), 
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
#done in the previous steps for enrichGO


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
for (algo in names(DiffNet_sub)){
  for (comp in names(DiffNet_sub[[algo]])){
    for (nam in names(DiffNet_sub[[algo]][[comp]])){
      for (cell in names(DiffNet_sub[[algo]][[comp]][[nam]])){
        for (msig in names(msigdb.gmt)){
          #group-specific hubs
          tmp <- df_target_sp[[algo]][[comp]][[nam]][[cell]]
          
          if(length(tmp) != 0){
            enricher_sp[[algo]][[comp]][[nam]][[cell]][[msig]] <- compareCluster(Gene~group, 
                                                                                 data=tmp, 
                                                                                 fun = "enricher",
                                                                                 TERM2GENE = msigdb.gmt[[msig]],
                                                                                 pvalueCutoff  = 0.05
            )
          }
        }
      }
    }
  }
}

saveRDS(enricher_sp, file = paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/CoDiNA/enricher_specific_forINTERACTIONS.rds")) 

#enricher_sp = readRDS("~/p728/RSTUDIO/analysis/tagaoglu/data/CoDiNA/enricher_specific_forINTERACTIONS.rds")

###Visualization :  dotplots

#Create dotplots, and then Save

for (algo in names(enricher_sp)){
  for (comp in names(enricher_sp[[algo]])){
    for (nam in names(enricher_sp[[algo]][[comp]])){
      for (cell in names(enricher_sp[[algo]][[comp]][[nam]])){
        for(msig in names(enricher_sp[[algo]][[comp]][[nam]][[cell]])){
          tmp_sp <- enricher_sp[[algo]][[comp]][[nam]][[cell]][[msig]]
          
          if(dim(tmp_sp)[1]!=0){
            out_sp <- clusterProfiler::dotplot(tmp_sp, showCategory = 10, title = paste0("MSigDb enrichment analysis performed using gene set ", msig))
          }
          
          ggsave(
            paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/CoDiNA/plots/ORA_plots_forINTERACTIONS/", algo, "/", comp, "/", nam, "/", cell, "_enricher_sp_", msig, ".dotplot.pdf"), 
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

# ...

## 1. Get target genelist

# For enrichKEGG, Input "target" data should a vector of entrez gene ids, so first do the conversion
tmp_target_sp <- list()
for (algo in names(DiffNet_sub)){
  for (comp in names(DiffNet_sub[[algo]])){
    for (nam in names(DiffNet_sub[[algo]][[comp]])){
      for (cell in names(DiffNet_sub[[algo]][[comp]][[nam]])){
        for (org in names(DiffNet_sub[[algo]][[comp]][[nam]][[cell]])){
          #group-specific hubs
          x <- DiffNet_sub[[algo]][[comp]][[nam]][[cell]][[org]]
          
          tmp <- bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
          
          tmp_target_sp[[algo]][[comp]][[nam]][[cell]][[org]] <- tmp[,2]
          
        }
      }
    }
  }
}


# Then, Create dataframe to be able to visualize enrichment results across groups in one plot
#Then, Create dataframe to be able to visualize enrichment results across groups in one plot
df_target_sp <- list()
for (algo in names(tmp_target_sp)){
  for (comp in names(tmp_target_sp[[algo]])){
    for (nam in names(tmp_target_sp[[algo]][[comp]])){
      for (cell in names(tmp_target_sp[[algo]][[comp]][[nam]])){
        mydf <- data.frame()
        for (org in names(tmp_target_sp[[algo]][[comp]][[nam]][[cell]])){
          l <- tmp_target_sp[[algo]][[comp]][[nam]][[cell]][[org]]
          if(length(l) != 0){
            tmp <- data.frame(Gene=l)
            tmp$group <- org
            mydf <- rbind(mydf, tmp)
          }
        }
        df_target_sp[[algo]][[comp]][[nam]][[cell]] <- mydf
      }
    }
  }
}



## 2. Run enrichKEGG via compareCluster
tmp_enrichKEGG_sp <- list()
for (algo in names(DiffNet_sub)){
  for (comp in names(DiffNet_sub[[algo]])){
    for (nam in names(DiffNet_sub[[algo]][[comp]])){
      for (cell in names(DiffNet_sub[[algo]][[comp]][[nam]])){
        #group-specific hubs
        tmp <- df_target_sp[[algo]][[comp]][[nam]][[cell]]
        
        if(length(tmp) != 0){
          tmp_enrichKEGG_sp[[algo]][[comp]][[nam]][[cell]] <- compareCluster(Gene~group, 
                                                                             data=tmp,
                                                                             fun="enrichKEGG",
                                                                             organism  = 'hsa',
                                                                             keyType = "kegg",
                                                                             pvalueCutoff = 0.05
          )
        }
      }
    }
  }
}


## The geneID column is ENTREZID, So Using setReadable, The geneID column is translated to symbol
enrichKEGG_sp <- list()
for (algo in names(tmp_enrichKEGG_sp)){
  for (comp in names(tmp_enrichKEGG_sp[[algo]])){                                                    
    for (nam in names(tmp_enrichKEGG_sp[[algo]][[comp]])){ 
      for (cell in names(tmp_enrichKEGG_sp[[algo]][[comp]][[nam]])){ 
        #group-specific hubs
        enrichKEGG_sp[[algo]][[comp]][[nam]][[cell]] <- setReadable(tmp_enrichKEGG_sp[[algo]][[comp]][[nam]][[cell]], 
                                                                    OrgDb = org.Hs.eg.db, keyType="ENTREZID")
      }
    }
  }
}

saveRDS(enrichKEGG_sp, file = paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/CoDiNA/enrichKEGG_specific_forINTERACTIONS.rds")) 

#enrichKEGG_sp = readRDS("~/p728/RSTUDIO/analysis/tagaoglu/data/CoDiNA/enrichKEGG_specific_forINTERACTIONS.rds")

###Visualization :  dotplots

#Create dotplots, and then Save

for (algo in names(enrichKEGG_sp)){
  for (comp in names(enrichKEGG_sp[[algo]])){
    for (nam in names(enrichKEGG_sp[[algo]][[comp]])){
      for (cell in names(enrichKEGG_sp[[algo]][[comp]][[nam]])){
        tmp_sp <- enrichKEGG_sp[[algo]][[comp]][[nam]][[cell]] 
        
        if(dim(tmp_sp)[1]!=0){
          out_sp <- clusterProfiler::dotplot(tmp_sp, showCategory = 10, title = paste0("KEGG pathway enrichment analysis"))
        }
        
        ggsave(
          paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/CoDiNA/plots/ORA_plots_forINTERACTIONS/", algo, "/", comp, "/", nam, "/", cell, "_enrichKEGG_sp.dotplot.pdf"), 
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













#############################################
### Run CoDiNA per sample
#############################################


#Put all in one list
algo_list <- c("GENIE", "PIDC")

all <- sapply(algo_list,function(x) NULL)
all[["GENIE"]] <- GENIE
all[["PIDC"]] <- PIDC

samples <- c("P123s2", "P123s8", "P193", "P199", "P108", "P191", "P207", "P207Gel", "P215", "P220", "P275", "P282C", "P295")

cell_interest <- "Hepatocytes"

#comp <- paste0("comp_hep_NORMALvs", sample)

#############################################


#####################################################
#Creating the Differential Network (CoDiNA)
#####################################################

#CoDiNA: Co-Expression Differential Network Analysis
#Categorize links and nodes from multiple networks in 3 categories: 
#Common links (alpha) specific links (gamma), and different links (beta). 
#Also categorizes the links into sub-categories and groups. The package includes a visualization tool for the networks.

#INPUT: Input data for the CoDiNA R package can be any networks, filtered for containing only significant links (according to the network construction method used). Edge list is a list containing all the links and their weights. The user can assign a weight of zero to links for which the p-value is not significant. 
#STEPS: The function MakeDiffNet() classifies the links into the Φ and  categories, calculates and normalises the scores. Its output is used as input for assigning the nodes into categories by the function ClusterNodes(). The plot() function can be used on the output from MakeDiffNet() and automatically calls the function ClusterNodes().


#perSample
# Differential Network Analysis
DiffNet <- list()
for(algo in names(all)){
  for(nam in names(all[[algo]])){
    for(cell in cell_interest){
      # Compare normal vs each sample separately
      for(sample in samples){
        DiffNet[[algo]][[nam]][[cell]][[sample]] = MakeDiffNet(Data = list(all[[algo]][[nam]][["normal"]][[cell]],
                                                                 all[[algo]][[nam]][[sample]][[cell]]), 
                                                     Code = c("normal", sample))
      }
    
    }
  }
}


# Save the data for CoDiNA
saveRDS(DiffNet, file = paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/CoDiNA/perSample_10k_DiffNet.Rds"))  




#####################################################
#Clustering the nodes into Φ and Φ̃ categories
#####################################################

#The suggested values for the internal and external cutoffs are the median or the first and third quantiles of the internal and Φ̃ scores, depending on how conservative the network should be.

#cutoff.external : The cut-off between the clusters (delta from the center to the edge coordinates), the closer to 1, the better.
#cutoff.internal : The cut-off inside the clusters (delta from the theoretical cluster to the edge coordinates), the closer to zero, the better.

int_C <- list()
ext_C <- list()
Nodes_Groups <- list()

for (algo in names(DiffNet)){
  for (nam in names(DiffNet[[algo]])){
    for (cell in names(DiffNet[[algo]][[nam]])){
      for (sample in names(DiffNet[[algo]][[nam]][[cell]])){
        #Using the median:
        int_C[[algo]][[nam]][[cell]][[sample]] = quantile(DiffNet[[algo]][[nam]][[cell]][[sample]]$Score_internal, 0.5)
        ext_C[[algo]][[nam]][[cell]][[sample]] = quantile(DiffNet[[algo]][[nam]][[cell]][[sample]]$Score_Phi, 0.5)
        #or
        # #Using the first and third quantile:
        # int_C[[algo]][[nam]][[cell]][[sample]] = quantile(DiffNet[[algo]][[nam]][[cell]][[sample]]$Score_internal, 0.25)
        # ext_C[[algo]][[nam]][[cell]][[sample]] = quantile(DiffNet[[algo]][[nam]][[cell]][[sample]]$Score_Phi, 0.75)
        
        Nodes_Groups[[algo]][[nam]][[cell]][[sample]] = ClusterNodes(DiffNet = DiffNet[[algo]][[nam]][[cell]][[sample]], 
                                                                     cutoff.external = ext_C[[algo]][[nam]][[cell]][[sample]], 
                                                                     cutoff.internal = int_C[[algo]][[nam]][[cell]][[sample]])
      }
    }
  }
}

# Save the data 
saveRDS(Nodes_Groups, file = paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/CoDiNA/perSample_10k_Nodes_Groups.Rds"))  


#####################################################
#Create table
#####################################################

##Create files if it does not exist
for (algo in names(Nodes_Groups)){
  for (nam in names(Nodes_Groups[[algo]])){
    plotdir=paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/CoDiNA/plots/Networks_html/perSample/", algo, "/", nam)
    if(!dir.exists(plotdir)){dir.create(plotdir,  recursive = T)}
  }
}


sample_df <- list()
for (algo in names(Nodes_Groups)){
  for (nam in names(Nodes_Groups[[algo]])){
    for (cell in names(Nodes_Groups[[algo]][[nam]])){
      for(sample in names(Nodes_Groups[[algo]][[nam]][[cell]])){
        df <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("Group", "Count", "Genes"))
        for (i in 1:length(names(table(Nodes_Groups[[algo]][[nam]][[cell]][[sample]]$Phi_tilde)))) {
          Group <- names(table(Nodes_Groups[[algo]][[nam]][[cell]][[sample]]$Phi_tilde))[i]
          Count <- as.numeric(table(Nodes_Groups[[algo]][[nam]][[cell]][[sample]]$Phi_tilde))[i]
          Gene_names <- subset(Nodes_Groups[[algo]][[nam]][[cell]][[sample]], Phi_tilde==Group)$Node
          Genes <- paste(Gene_names, collapse=', ')
          res= c(Group, Count, Genes)
          df[i,]=res
        }
        sample_df[[algo]][[nam]][[cell]][[sample]] <- df
      }
    }
  }
}


for (algo in names(sample_df)){
  for (nam in names(sample_df[[algo]])){
    for (cell in names(sample_df[[algo]][[nam]])){
      for (sample in names(sample_df[[algo]][[nam]][[cell]])){
        sample_df[[algo]][[nam]][[cell]][[sample]] %>% 
          gt() %>%
          tab_header(title = sample) %>%
          gtsave(paste0(sample, "_SummaryTable.html"), 
                 path = paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/CoDiNA/plots/Networks_html/perSample/", algo, "/", nam, "/"))
      }
    }
  }
}



