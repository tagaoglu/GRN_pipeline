# Set *.libPaths*  
.libPaths(c("~/p728/RSTUDIO/R/library/4.1/", .libPaths()))


#Set working directory and seed:
setwd("~/p728/RSTUDIO/analysis/tagaoglu")
set.seed(123)


#Load functions
source("/home/tagaoglu/p728/RSTUDIO/analysis/tagaoglu/functions/centralityMetrics.R")
source("/home/tagaoglu/p728/RSTUDIO/analysis/tagaoglu/functions/save_png_pdf.R")


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
GENIE_wRcis <- readRDS(paste("~/p728/RSTUDIO/analysis/tagaoglu/data/networks/GENIE3/GENIE3_all_10k_OnlyActivatedLinks_Rcis_5kb.Rds", sep = '' ))
PIDC <- readRDS(paste("~/p728/RSTUDIO/analysis/tagaoglu/data/networks/PIDC/PIDC_all_10k.Rds", sep = '' ))


#exclude random data
for(nam in names(GENIE)){
  GENIE <- GENIE[names(GENIE) %in% nam_dir]
  GENIE_wRcis <- GENIE_wRcis[names(GENIE_wRcis) %in% nam_dir]
  PIDC <- PIDC[names(PIDC) %in% nam_dir]
}

#Sort cell names of list objects to compare easily
for(nam in names(GENIE)){
  for(group in names(GENIE[[nam]])){
    GENIE[[nam]][[group]] <- GENIE[[nam]][[group]][order(names(GENIE[[nam]][[group]]))]
    GENIE_wRcis[[nam]][[group]] <- GENIE_wRcis[[nam]][[group]][order(names(GENIE_wRcis[[nam]][[group]]))]
    PIDC[[nam]][[group]] <- PIDC[[nam]][[group]][order(names(PIDC[[nam]][[group]]))]
  }
}


##################################################################################################
###Critical Genes Identification
##################################################################################################

##################################################################################################
###Network Statistics
#Analyze networks and extract certain statistics about the network that tell us about structural properties of networks
##################################################################################################

## Load package
library(igraph)


###Create igraph object for each edgelist
GENIE_igraph <- sapply(names(GENIE), function(nam){
  sapply(names(GENIE[[nam]]), function(group){
    sapply(names(GENIE[[nam]][[group]]), function(cell){
      graph_from_data_frame(GENIE[[nam]][[group]][[cell]], directed=GENIE_edge_prediction) ## Reconstructing the graph
    }, simplify=F)
  }, simplify=F)
}, simplify=F)

#Check next time
# Following stuff caused problem since the data frame should contain at least two columns
GENIE_wRcis[["expData_mean_1"]][["P275"]][["Hepatocytes"]] <- NULL
GENIE_wRcis[["expData_mean_2"]][["normal"]][["Hepatocytes"]] <- NULL
GENIE_wRcis[["expData_mean_2"]][["P275"]][["Hepatocytes"]] <- NULL
GENIE_wRcis[["expData_mean_2"]][["nonviral"]][["pDCs"]] <- NULL


GENIE_wRcis_igraph <- sapply(names(GENIE_wRcis), function(nam){
  sapply(names(GENIE_wRcis[[nam]]), function(group){
    sapply(names(GENIE_wRcis[[nam]][[group]]), function(cell){
      graph_from_data_frame(GENIE_wRcis[[nam]][[group]][[cell]][["FilteredLinks"]], directed=GENIE_edge_prediction) ## Reconstructing the graph
    }, simplify=F)
  }, simplify=F)
}, simplify=F)


PIDC_igraph <- sapply(names(PIDC), function(nam){
  sapply(names(PIDC[[nam]]), function(group){
    sapply(names(PIDC[[nam]][[group]]), function(cell){
      graph_from_data_frame(PIDC[[nam]][[group]][[cell]], directed=PIDC_edge_prediction) ## Reconstructing the graph
    }, simplify=F)
  }, simplify=F)
}, simplify=F)


# Calculate "5 centrality metrics" and then "integrated centrality metrics"
#degree, closeness, betweenness, eigenvector and PageRank
#node centrality metrics which are employed to measure the importance of nodes in the constructed networks.

#Note: I may need to exclude centrality metric: "eigenvalue", 
#because for GENIE_wRcis_centralityMetrics calculation I may encounter the following problem:
#ERROR message was: At core/centrality/centrality_other.c:328 : graph is directed and acyclic; eigenvector centralities will be zeros.
#PROBLEM explained here: http://users.dimi.uniud.it/~massimo.franceschet/teaching/datascience/network/katz.html

GENIE_centralityMetrics <- sapply(names(GENIE), function(nam){
  sapply(names(GENIE[[nam]]), function(group){
    sapply(names(GENIE[[nam]][[group]]), function(cell){
      centralityMetrics(GENIE[[nam]][[group]][[cell]], edge_prediction=GENIE_edge_prediction) #FilteredLinks with RcisTarget
    }, simplify=F)
  }, simplify=F)
}, simplify=F)


GENIE_wRcis_centralityMetrics <- sapply(names(GENIE_wRcis), function(nam){
  sapply(names(GENIE_wRcis[[nam]]), function(group){
    sapply(names(GENIE_wRcis[[nam]][[group]]), function(cell){
      centralityMetrics(GENIE_wRcis[[nam]][[group]][[cell]][["FilteredLinks"]], edge_prediction=GENIE_edge_prediction) #FilteredLinks with RcisTarget
    }, simplify=F)
  }, simplify=F)
}, simplify=F)


PIDC_centralityMetrics <-  sapply(names(PIDC), function(nam){
  sapply(names(PIDC[[nam]]), function(group){
    sapply(names(PIDC[[nam]][[group]]), function(cell){
      centralityMetrics(PIDC[[nam]][[group]][[cell]], edge_prediction=PIDC_edge_prediction) 
    }, simplify=F)
  }, simplify=F)
}, simplify=F)


###Extract all centrality metrics
GENIE_metrics <- sapply(names(GENIE), function(nam){
  sapply(names(GENIE[[nam]]), function(group){
    sapply(names(GENIE[[nam]][[group]]), function(cell){
      GENIE_centralityMetrics[[nam]][[group]][[cell]][["allMetrics"]]
    }, simplify=F)
  }, simplify=F)
}, simplify=F)


GENIE_wRcis_metrics <- sapply(names(GENIE_wRcis), function(nam){
  sapply(names(GENIE_wRcis[[nam]]), function(group){
    sapply(names(GENIE_wRcis[[nam]][[group]]), function(cell){
      GENIE_wRcis_centralityMetrics[[nam]][[group]][[cell]][["allMetrics"]]
    }, simplify=F)
  }, simplify=F)
}, simplify=F)


PIDC_metrics <- sapply(names(PIDC), function(nam){
  sapply(names(PIDC[[nam]]), function(group){
    sapply(names(PIDC[[nam]][[group]]), function(cell){
      PIDC_centralityMetrics[[nam]][[group]][[cell]][["allMetrics"]]
    }, simplify=F)
  }, simplify=F)
}, simplify=F)


###Extract integrated centrality metrics of the nodes for each GRN
GENIE_int <- sapply(names(GENIE), function(nam){
  sapply(names(GENIE[[nam]]), function(group){
    sapply(names(GENIE[[nam]][[group]]), function(cell){
      GENIE_centralityMetrics[[nam]][[group]][[cell]][["integratedMetrics"]]
    }, simplify=F)
  }, simplify=F)
}, simplify=F)


GENIE_wRcis_int <- sapply(names(GENIE_wRcis), function(nam){
  sapply(names(GENIE_wRcis[[nam]]), function(group){
    sapply(names(GENIE_wRcis[[nam]][[group]]), function(cell){
      GENIE_wRcis_centralityMetrics[[nam]][[group]][[cell]][["integratedMetrics"]]
    }, simplify=F)
  }, simplify=F)
}, simplify=F)


PIDC_int <- sapply(names(PIDC), function(nam){
  sapply(names(PIDC[[nam]]), function(group){
    sapply(names(PIDC[[nam]][[group]]), function(cell){
      PIDC_centralityMetrics[[nam]][[group]][[cell]][["integratedMetrics"]]
    }, simplify=F)
  }, simplify=F)
}, simplify=F)




##################################################################################################
###Visualization
##################################################################################################

#Put all in one list
algo_list <- c("GENIE", "GENIE_wRcis", "PIDC")


#Put all in one list

int_metrics <- sapply(algo_list,function(x) NULL)
int_metrics[["GENIE"]] <- GENIE_int
int_metrics[["GENIE_wRcis"]] <- GENIE_wRcis_int
int_metrics[["PIDC"]] <- PIDC_int


net_igraph <- sapply(algo_list,function(x) NULL)
net_igraph[["GENIE"]] <- GENIE_igraph
net_igraph[["GENIE_wRcis"]] <- GENIE_wRcis_igraph
net_igraph[["PIDC"]] <- PIDC_igraph


##################################################################################################
## For Top10 nodes, Visualization & Extract Common&OverlappingSpecificHubs
##################################################################################################

# library
library(ggplot2)
library(ggrepel)
library(ggpubr)
theme_set(theme_pubr())

library(gt)
library(webshot)
library(scales)
library(readr)

comparisons <- c("comp_NORMALvsVIRALvsNONVIRAL", "comp_NORMALvsTUMOR", "comp_hep_NORMALvsSAMPLES")


##Create files if it does not exist
for (algo in names(int_metrics)) {
  for(nam in names(int_metrics[[algo]])){
    for(comp in comparisons){
      plotdir=paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Top10_Nodes/", algo, "/", nam, "/", comp)
      if(!dir.exists(plotdir)){dir.create(plotdir,  recursive = T)}
    }
  }
}


####################
###Plot Ranks for Nodes based on integrated centrality measure with Top10 critical Nodes highlighted

#Due to error "(too many overlaps). Consider increasing max.overlaps", set the following:
options(ggrepel.max.overlaps = Inf)

#Plot Top10 Nodes

# Create a list to hold the plot objects.
pList <- list()
for (algo in names(int_metrics)) {
  for(nam in names(int_metrics[[algo]])){
    for(group in names(int_metrics[[algo]][[nam]])){
      for(cell in names(int_metrics[[algo]][[nam]][[group]])){
        data <- int_metrics[[algo]][[nam]][[group]][[cell]]
        data$Top10 <- ifelse(data$rank <= 10, "Top10", "others")
        p <- ggplot(data, aes(x=rank, y=minuslog10)) +
          geom_point(aes(color = Top10))  +
          scale_color_manual(values = c("grey", "red")) +
          #lims(x=c(0,100), y=c(0,100)) +
          labs(x = "Rank", y="-log10(pvalue)", title = paste0(group)) + 
          guides(col= guide_legend(title= "Nodes"), title.position = c(0.8, 0.8)) +
          geom_text_repel(
            data = subset(data, rank <= 10),
            aes(x=rank, y=minuslog10, label=Name),
            size = 5,
            box.padding = unit(0.35, "lines"),
            point.padding = unit(0.3, "lines")) + 
          theme_classic()
        pList[[algo]][[nam]][[group]][[cell]] <- p
      }
    }
  }
}


figure <- list()
for (algo in names(int_metrics)) {
  for(nam in names(int_metrics[[algo]])){
    for(group in names(int_metrics[[algo]][[nam]])){
      for(cell in names(int_metrics[[algo]][[nam]][[group]])){
        figure[[algo]][[nam]][[cell]][[group]] <- pList[[algo]][[nam]][[group]][[cell]]
      }
    }
  }
}



####################

#Save plots and Create SummaryTable for Common&OverlappingSpecificNodes

####################


##Extract Top10 Nodes
Top10_Nodes <- list()
for (algo in names(int_metrics)) {
  for(nam in names(int_metrics[[algo]])){
    for(group in names(int_metrics[[algo]][[nam]])){
      for(cell in names(int_metrics[[algo]][[nam]][[group]])){
        data <- int_metrics[[algo]][[nam]][[group]][[cell]]
        Top10_Nodes[[algo]][[nam]][[cell]][[group]] <- as.character(rownames(data[data$rank <= 10,]))
      }
    }
  }
}


# Save Plots for all groups
for (algo in names(figure)) {
  for(nam in names(figure[[algo]])){
    for(cell in names(figure[[algo]][[nam]])){
      tmp <- ggarrange(figure[[algo]][[nam]][[cell]][["nonviral"]],
                       figure[[algo]][[nam]][[cell]][["viral"]], 
                       figure[[algo]][[nam]][[cell]][["normal"]],
                       figure[[algo]][[nam]][[cell]][["tumor"]],
                       nrow = 1, ncol = 4)
      f <- annotate_figure(tmp, top = text_grob(paste0(algo, " || ",nam, " || ", cell), color = "blue", face = "bold", size = 14))
      save_png_pdf(f, paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Top10_Nodes/", algo, "/", nam, "/", cell, "_Top10"), 
                   height = 8, width = 18, res=250)
      
    }
  }
}




##################################################################################################
### COMPARISON : NORMAL vs VIRAL vs NONVIRAL
##################################################################################################

comp <- "comp_NORMALvsVIRALvsNONVIRAL"
groups <- c("normal", "viral", "nonviral")

#####################################################

# # Save Plots for corresponding groups
# for (algo in names(figure)) {
#   for(nam in names(figure[[algo]])){
#     for(cell in names(figure[[algo]][[nam]])){
#       tmp <- ggarrange(figure[[algo]][[nam]][[cell]][["normal"]],
#                        figure[[algo]][[nam]][[cell]][["viral"]], 
#                        figure[[algo]][[nam]][[cell]][["nonviral"]],
#                        nrow = 1, ncol = 3)
#       f <- annotate_figure(tmp, top = text_grob(paste0(algo, " || ",nam, " || ", cell), color = "blue", face = "bold", size = 14))
#       save_png_pdf(f, paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Top10_Nodes/", algo, "/", nam, "/", comp, "/", cell, "_Top10"), height = 10, width = 18, res=250)
#       
#     }
#   }
# }


# Extract Common&OverlappingSpecificNodes and Create Summary Table

#Keep only corresponding groups
Top10_Nodes_comp <- list()
for (algo in names(Top10_Nodes)) {
  for(nam in names(Top10_Nodes[[algo]])){
    for(cell in names(Top10_Nodes[[algo]][[nam]])){
      Top10_Nodes_comp[[algo]][[nam]][[cell]] <- Top10_Nodes[[algo]][[nam]][[cell]][names(Top10_Nodes[[algo]][[nam]][[cell]]) %in% groups]
    }
  }
}


#Get specific hubs for each group (which are identified in Network_analysis.R )
top_centr_sp_all <- readRDS(paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Network_analysis/networks_hub_specific_", comp, ".rds"))


#Extract common Top10 nodes among corresponding groups
common_Nodes <- list()
for (algo in names(Top10_Nodes_comp)) {
  for(nam in names(Top10_Nodes_comp[[algo]])){
    for(cell in names(Top10_Nodes_comp[[algo]][[nam]])){
      common_Nodes[[algo]][[nam]][[cell]] <- Reduce(intersect, Top10_Nodes_comp[[algo]][[nam]][[cell]])
    }
  }
}


#Extract overlapping nodes between Top10 nodes and specific hubs for each group (which are identified in Network_analysis.R))
Overlap_specific_Nodes <- list()
for (algo in names(Top10_Nodes_comp)) {
  for(nam in names(Top10_Nodes_comp[[algo]])){
    for(cell in names(Top10_Nodes_comp[[algo]][[nam]])){
      for(group in names(Top10_Nodes_comp[[algo]][[nam]][[cell]])){
        Overlap_specific_Nodes[[algo]][[nam]][[cell]][[group]] <- intersect(Top10_Nodes_comp[[algo]][[nam]][[cell]][[group]], 
                                                                            top_centr_sp_all[[algo]][[nam]][[cell]][["integrated_rank"]][[group]])
      }
    }
  }
}



#Create table
for (algo in names(Top10_Nodes_comp)) {
  for(nam in names(Top10_Nodes_comp[[algo]])){
    for(cell in names(Top10_Nodes_comp[[algo]][[nam]])){
      Group <- c("common", "normal", "viral", "nonviral")
      Count <- c(length(common_Nodes[[algo]][[nam]][[cell]]), 
                 length(Overlap_specific_Nodes[[algo]][[nam]][[cell]][["normal"]]), 
                 length(Overlap_specific_Nodes[[algo]][[nam]][[cell]][["viral"]]), 
                 length(Overlap_specific_Nodes[[algo]][[nam]][[cell]][["nonviral"]]))
      Genes <- c(paste(common_Nodes[[algo]][[nam]][[cell]], collapse=', '), 
                 paste(Overlap_specific_Nodes[[algo]][[nam]][[cell]][["normal"]], collapse=', '), 
                 paste(Overlap_specific_Nodes[[algo]][[nam]][[cell]][["viral"]], collapse=', '), 
                 paste(Overlap_specific_Nodes[[algo]][[nam]][[cell]][["nonviral"]], collapse=', '))
      df <- data.frame(Group, Count, Genes)
      
      df %>% 
        gt() %>%
        tab_header(title = paste0(algo, " | ", nam, " | ", cell)) %>%
        gtsave(paste0(cell, "_Common&OverlappingSpecificNodes_SummaryTable.html"), 
               path = paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Top10_Nodes/", algo, "/", nam, "/", comp,"/"))
    }
  }
}



##################################################################################################
### COMPARISON : NORMAL vs TUMOR
##################################################################################################

comp <- "comp_NORMALvsTUMOR"
groups <- c("normal", "tumor")

#####################################################

# # Save Plots for corresponding groups
# for (algo in names(figure)) {
#   for(nam in names(figure[[algo]])){
#     for(cell in names(figure[[algo]][[nam]])){
#       tmp <- ggarrange(figure[[algo]][[nam]][[cell]][["normal"]],
#                        figure[[algo]][[nam]][[cell]][["tumor"]],
#                        nrow = 1, ncol = 2)
#       f <- annotate_figure(tmp, top = text_grob(paste0(algo, " || ",nam, " || ", cell), color = "blue", face = "bold", size = 14))
#       save_png_pdf(f, paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Top10_Nodes/", algo, "/", nam, "/", comp, "/", cell, "_Top10"), 
#                    height = 10, width = 12, res=250)
#       
#     }
#   }
# }


# Extract Common&OverlappingSpecificNodes and Create Summary Table

#Keep only corresponding groups
Top10_Nodes_comp <- list()
for (algo in names(Top10_Nodes)) {
  for(nam in names(Top10_Nodes[[algo]])){
    for(cell in names(Top10_Nodes[[algo]][[nam]])){
      Top10_Nodes_comp[[algo]][[nam]][[cell]] <- Top10_Nodes[[algo]][[nam]][[cell]][names(Top10_Nodes[[algo]][[nam]][[cell]]) %in% groups]
    }
  }
}


#Get specific hubs for each group (which are identified in Network_analysis.R )
top_centr_sp_all <- readRDS(paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Network_analysis/networks_hub_specific_", comp, ".rds"))


#Extract common Top10 nodes among corresponding groups
common_Nodes <- list()
for (algo in names(Top10_Nodes_comp)) {
  for(nam in names(Top10_Nodes_comp[[algo]])){
    for(cell in names(Top10_Nodes_comp[[algo]][[nam]])){
      common_Nodes[[algo]][[nam]][[cell]] <- Reduce(intersect, Top10_Nodes_comp[[algo]][[nam]][[cell]])
    }
  }
}


#Extract overlapping nodes between Top10 nodes and specific hubs for each group (which are identified in Network_analysis.R))
Overlap_specific_Nodes <- list()
for (algo in names(Top10_Nodes_comp)) {
  for(nam in names(Top10_Nodes_comp[[algo]])){
    for(cell in names(Top10_Nodes_comp[[algo]][[nam]])){
      for(group in names(Top10_Nodes_comp[[algo]][[nam]][[cell]])){
        Overlap_specific_Nodes[[algo]][[nam]][[cell]][[group]] <- intersect(Top10_Nodes_comp[[algo]][[nam]][[cell]][[group]], 
                                                                            top_centr_sp_all[[algo]][[nam]][[cell]][["integrated_rank"]][[group]])
      }
    }
  }
}



#Create table
for (algo in names(Top10_Nodes_comp)) {
  for(nam in names(Top10_Nodes_comp[[algo]])){
    for(cell in names(Top10_Nodes_comp[[algo]][[nam]])){
      Group <- c("common", "normal", "tumor")
      Count <- c(length(common_Nodes[[algo]][[nam]][[cell]]), 
                 length(Overlap_specific_Nodes[[algo]][[nam]][[cell]][["normal"]]), 
                 length(Overlap_specific_Nodes[[algo]][[nam]][[cell]][["tumor"]]))
      Genes <- c(paste(common_Nodes[[algo]][[nam]][[cell]], collapse=', '), 
                 paste(Overlap_specific_Nodes[[algo]][[nam]][[cell]][["normal"]], collapse=', '), 
                 paste(Overlap_specific_Nodes[[algo]][[nam]][[cell]][["tumor"]], collapse=', '))
      df <- data.frame(Group, Count, Genes)
      
      df %>% 
        gt() %>%
        tab_header(title = paste0(algo, " | ", nam, " | ", cell)) %>%
        gtsave(paste0(cell, "_Common&OverlappingSpecificNodes_SummaryTable.html"), 
               path = paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Top10_Nodes/", algo, "/", nam, "/", comp,"/"))
    }
  }
}




##################################################################################################
### COMPARISON : hep NORMAL vs SAMPLES
##################################################################################################

comp <- "comp_hep_NORMALvsSAMPLES"
groups <- c("normal", "P123s2", "P123s8", "P193", "P199", "P108", "P191", "P207", "P207Gel", "P215", "P220", "P275", "P282C", "P295")

cell <- "Hepatocytes"


#####################################################

# Save Plots for corresponding groups

# Plot hep_perSample
for (algo in names(figure)) {
  for(nam in names(figure[[algo]])){
    tmp <- ggarrange(figure[[algo]][[nam]][[cell]][["normal"]], figure[[algo]][[nam]][[cell]][["P123s2"]], 
                     figure[[algo]][[nam]][[cell]][["P123s8"]], figure[[algo]][[nam]][[cell]][["P193"]], 
                     figure[[algo]][[nam]][[cell]][["P199"]], figure[[algo]][[nam]][[cell]][["P108"]], 
                     figure[[algo]][[nam]][[cell]][["P191"]], figure[[algo]][[nam]][[cell]][["P207"]], 
                     figure[[algo]][[nam]][[cell]][["P207Gel"]], figure[[algo]][[nam]][[cell]][["P215"]], 
                     figure[[algo]][[nam]][[cell]][["P220"]], figure[[algo]][[nam]][[cell]][["P275"]], 
                     figure[[algo]][[nam]][[cell]][["P282C"]], figure[[algo]][[nam]][[cell]][["P295"]], 
                     nrow = 2, ncol = 7)
    f <- annotate_figure(tmp, top = text_grob(paste0(algo, " || ",nam, " || ", cell), color = "blue", face = "bold", size = 14))
    save_png_pdf(f, paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Top10_Nodes/", algo, "/", nam, "/", comp, "/", cell, "_Top10_perSample"), 
                 height = 16, width = 30, res=250)
  }
}


# Extract Common&OverlappingSpecificNodes and Create Summary Table

#Keep only corresponding groups
Top10_Nodes_comp <- list()
for (algo in names(Top10_Nodes)) {
  for(nam in names(Top10_Nodes[[algo]])){
    for(cell in cell){
      Top10_Nodes_comp[[algo]][[nam]][[cell]] <- Top10_Nodes[[algo]][[nam]][[cell]][names(Top10_Nodes[[algo]][[nam]][[cell]]) %in% groups]
    }
  }
}


#Get specific hubs for each group (which are identified in Network_analysis.R )
top_centr_sp_all <- readRDS(paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Network_analysis/networks_hub_specific_", comp, ".rds"))


#Extract common Top10 nodes among corresponding groups
common_Nodes <- list()
for (algo in names(Top10_Nodes_comp)) {
  for(nam in names(Top10_Nodes_comp[[algo]])){
    for(cell in names(Top10_Nodes_comp[[algo]][[nam]])){
      common_Nodes[[algo]][[nam]][[cell]] <- Reduce(intersect, Top10_Nodes_comp[[algo]][[nam]][[cell]])
    }
  }
}

#Extract common Top10 nodes among only samples
samples <- c("P123s2", "P123s8", "P193", "P199", "P108", "P191", "P207", "P207Gel", "P215", "P220", "P275", "P282C", "P295")

#Keep only samples
Top10_Nodes_comp_samples <- list()
for (algo in names(Top10_Nodes)) {
  for(nam in names(Top10_Nodes[[algo]])){
    for(cell in cell){
      Top10_Nodes_comp_samples[[algo]][[nam]][[cell]] <- Top10_Nodes[[algo]][[nam]][[cell]][names(Top10_Nodes[[algo]][[nam]][[cell]]) %in% samples]
    }
  }
}

common_Nodes_onlyforSamples <- list()
for (algo in names(Top10_Nodes_comp_samples)) {
  for(nam in names(Top10_Nodes_comp_samples[[algo]])){
    for(cell in names(Top10_Nodes_comp_samples[[algo]][[nam]])){
      common_Nodes_onlyforSamples[[algo]][[nam]][[cell]] <- Reduce(intersect,Top10_Nodes_comp_samples[[algo]][[nam]][[cell]])
    }
  }
}


#Extract overlapping nodes between Top10 nodes and specific hubs for each group (which are identified in Network_analysis.R))
Overlap_specific_Nodes <- list()
for (algo in names(Top10_Nodes_comp)) {
  for(nam in names(Top10_Nodes_comp[[algo]])){
    for(cell in names(Top10_Nodes_comp[[algo]][[nam]])){
      for(group in names(Top10_Nodes_comp[[algo]][[nam]][[cell]])){
        Overlap_specific_Nodes[[algo]][[nam]][[cell]][[group]] <- intersect(Top10_Nodes_comp[[algo]][[nam]][[cell]][[group]], 
                                                                            top_centr_sp_all[[algo]][[nam]][[cell]][["integrated_rank"]][[group]])
      }
    }
  }
}


#Create table
for (algo in names(Top10_Nodes_comp)) {
  for(nam in names(Top10_Nodes_comp[[algo]])){
    for(cell in names(Top10_Nodes_comp[[algo]][[nam]])){
      Group <- c("common", "common (amongSamples)", "normal", "P123s2", "P123s8", "P193", "P199", "P108", "P191", "P207", "P207Gel", "P215", "P220", "P275", "P282C", "P295")
      Count <- c(length(common_Nodes[[algo]][[nam]][[cell]]), 
                 length(common_Nodes_onlyforSamples[[algo]][[nam]][[cell]]),
                 length(Overlap_specific_Nodes[[algo]][[nam]][[cell]][["normal"]]), 
                 length(Overlap_specific_Nodes[[algo]][[nam]][[cell]][["P123s2"]]),
                 length(Overlap_specific_Nodes[[algo]][[nam]][[cell]][["P123s8"]]),
                 length(Overlap_specific_Nodes[[algo]][[nam]][[cell]][["P193"]]),
                 length(Overlap_specific_Nodes[[algo]][[nam]][[cell]][["P199"]]),
                 length(Overlap_specific_Nodes[[algo]][[nam]][[cell]][["P108"]]),
                 length(Overlap_specific_Nodes[[algo]][[nam]][[cell]][["P191"]]),
                 length(Overlap_specific_Nodes[[algo]][[nam]][[cell]][["P207"]]),
                 length(Overlap_specific_Nodes[[algo]][[nam]][[cell]][["P207Gel"]]),
                 length(Overlap_specific_Nodes[[algo]][[nam]][[cell]][["P215"]]),
                 length(Overlap_specific_Nodes[[algo]][[nam]][[cell]][["P220"]]),
                 length(Overlap_specific_Nodes[[algo]][[nam]][[cell]][["P275"]]),
                 length(Overlap_specific_Nodes[[algo]][[nam]][[cell]][["P282C"]]),
                 length(Overlap_specific_Nodes[[algo]][[nam]][[cell]][["P295"]])
                 )
      Genes <- c(paste(common_Nodes[[algo]][[nam]][[cell]], collapse=', '), 
                 paste(common_Nodes_onlyforSamples[[algo]][[nam]][[cell]], collapse=', '),
                 paste(Overlap_specific_Nodes[[algo]][[nam]][[cell]][["normal"]], collapse=', '), 
                 paste(Overlap_specific_Nodes[[algo]][[nam]][[cell]][["P123s2"]], collapse=', '),
                 paste(Overlap_specific_Nodes[[algo]][[nam]][[cell]][["P123s8"]], collapse=', '),
                 paste(Overlap_specific_Nodes[[algo]][[nam]][[cell]][["P193"]], collapse=', '),
                 paste(Overlap_specific_Nodes[[algo]][[nam]][[cell]][["P199"]], collapse=', '),
                 paste(Overlap_specific_Nodes[[algo]][[nam]][[cell]][["P108"]], collapse=', '),
                 paste(Overlap_specific_Nodes[[algo]][[nam]][[cell]][["P191"]], collapse=', '),
                 paste(Overlap_specific_Nodes[[algo]][[nam]][[cell]][["P207"]], collapse=', '),
                 paste(Overlap_specific_Nodes[[algo]][[nam]][[cell]][["P207Gel"]], collapse=', '),
                 paste(Overlap_specific_Nodes[[algo]][[nam]][[cell]][["P215"]], collapse=', '),
                 paste(Overlap_specific_Nodes[[algo]][[nam]][[cell]][["P220"]], collapse=', '),
                 paste(Overlap_specific_Nodes[[algo]][[nam]][[cell]][["P275"]], collapse=', '),
                 paste(Overlap_specific_Nodes[[algo]][[nam]][[cell]][["P282C"]], collapse=', '),
                 paste(Overlap_specific_Nodes[[algo]][[nam]][[cell]][["P295"]], collapse=', ')
                 )
      df <- data.frame(Group, Count, Genes)
      
      df %>% 
        gt() %>%
        tab_header(title = paste0(algo, " | ", nam, " | ", cell)) %>%
        gtsave(paste0(cell, "_Common&OverlappingSpecificNodes_SummaryTable.html"), 
               path = paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Top10_Nodes/", algo, "/", nam, "/", comp,"/"))
    }
  }
}



##################################################################################################
### COMPARISON : hep NORMAL vs EACH SAMPLE separately
##################################################################################################

cell <- "Hepatocytes"

samples <- c("P123s2", "P123s8", "P193", "P199", "P108", "P191", "P207", "P207Gel", "P215", "P220", "P275", "P282C", "P295")

# Compare normal vs each sample separately
for(samp in samples){
  sample <- samp
  comp <- paste0("comp_hep_NORMALvs", sample)
  groups <- c("normal", sample)
  
  
  #####################################################
  
  # Save Plots for corresponding groups

  # # Plot hep_perSample
  # for (algo in names(figure)) {
  #   for(nam in names(figure[[algo]])){
  #     tmp <- ggarrange(figure[[algo]][[nam]][[cell]][["normal"]], figure[[algo]][[nam]][[cell]][[sample]], 
  #                      nrow = 1, ncol = 2)
  #     f <- annotate_figure(tmp, top = text_grob(paste0(algo, " || ",nam, " || ", cell), color = "blue", face = "bold", size = 14))
  #     save_png_pdf(f, paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Top10_Nodes/", algo, "/", nam, "/", comp, "_Top10"), 
  #                  height = 8, width = 9, res=250)
  #   }
  # }
  
  
  # Extract Common&OverlappingSpecificNodes and Create Summary Table
  
  #Keep only corresponding groups
  Top10_Nodes_comp <- list()
  for (algo in names(Top10_Nodes)) {
    for(nam in names(Top10_Nodes[[algo]])){
      for(cell in cell){
        Top10_Nodes_comp[[algo]][[nam]][[cell]] <- Top10_Nodes[[algo]][[nam]][[cell]][names(Top10_Nodes[[algo]][[nam]][[cell]]) %in% groups]
      }
    }
  }
  
  
  #Get specific hubs for each group (which are identified in Network_analysis.R )
  top_centr_sp_all <- readRDS(paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Network_analysis/networks_hub_specific_", comp, ".rds"))
  
  
  #Extract common Top10 nodes among corresponding groups
  common_Nodes <- list()
  for (algo in names(Top10_Nodes_comp)) {
    for(nam in names(Top10_Nodes_comp[[algo]])){
      for(cell in names(Top10_Nodes_comp[[algo]][[nam]])){
        common_Nodes[[algo]][[nam]][[cell]] <- Reduce(intersect, Top10_Nodes_comp[[algo]][[nam]][[cell]])
      }
    }
  }
  
  
  
  #Extract overlapping nodes between Top10 nodes and specific hubs for each group (which are identified in Network_analysis.R))
  Overlap_specific_Nodes <- list()
  for (algo in names(Top10_Nodes_comp)) {
    for(nam in names(Top10_Nodes_comp[[algo]])){
      for(cell in names(Top10_Nodes_comp[[algo]][[nam]])){
        for(group in names(Top10_Nodes_comp[[algo]][[nam]][[cell]])){
          Overlap_specific_Nodes[[algo]][[nam]][[cell]][[group]] <- intersect(Top10_Nodes_comp[[algo]][[nam]][[cell]][[group]], 
                                                                              top_centr_sp_all[[algo]][[nam]][[cell]][["integrated_rank"]][[group]])
        }
      }
    }
  }
  
  
  #Create table
  for (algo in names(Top10_Nodes_comp)) {
    for(nam in names(Top10_Nodes_comp[[algo]])){
      for(cell in names(Top10_Nodes_comp[[algo]][[nam]])){
        Group <- c("common", "normal", sample)
        Count <- c(length(common_Nodes[[algo]][[nam]][[cell]]), 
                   length(Overlap_specific_Nodes[[algo]][[nam]][[cell]][["normal"]]), 
                   length(Overlap_specific_Nodes[[algo]][[nam]][[cell]][[sample]])
        )
        Genes <- c(paste(common_Nodes[[algo]][[nam]][[cell]], collapse=', '),
                   paste(Overlap_specific_Nodes[[algo]][[nam]][[cell]][["normal"]], collapse=', '), 
                   paste(Overlap_specific_Nodes[[algo]][[nam]][[cell]][[sample]], collapse=', ')
        )
        df <- data.frame(Group, Count, Genes)
        
        df %>% 
          gt() %>%
          tab_header(title = paste0(algo, " | ", nam, " | ", cell)) %>%
          gtsave(paste0(comp, "_Common&OverlappingSpecificNodes_SummaryTable.html"), 
                 path = paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Top10_Nodes/", algo, "/", nam, "/"))
      }
    }
  }
  
  
}








#####################################################################################
#######
##################################################################################################


#####################################################################################
### Venn diagram to show the overlap of Nodes, and also to show the overlap of edges/regulations (TF-target pairs) across sample-specific GRNs
##################################################################################################

#Load library
library("ggvenn")
library(ggpubr)
theme_set(theme_pubr())


##Create files if it does not exist
for (algo in names(int_metrics)) {
  for(nam in names(int_metrics[[algo]])){
    plotdir=paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/VennDiagram/", algo, "/", nam)
    if(!dir.exists(plotdir)){dir.create(plotdir,  recursive = T)}
  }
}


####################
###Plot Venn diagram

Nodes <- list()
Edges <- list()
for (algo in names(net_igraph)) {
  for(nam in names(net_igraph[[algo]])){
    for(group in names(net_igraph[[algo]][[nam]])){
      for(cell in names(net_igraph[[algo]][[nam]][[group]])){
        net <- net_igraph[[algo]][[nam]][[group]][[cell]]
        Nodes[[algo]][[nam]][[cell]][[group]] <- names(V(net)) # Extracting the Nodes
        Edges[[algo]][[nam]][[cell]][[group]] <- unique(attr(E(net), "vnames")) # Extracting the edges, or #as_ids(E(net))
      }
    }
  }
}



#Exclude samples, and keep only corresponding groups
groups <- c("normal", "viral", "nonviral", "tumor") 

Nodes_perGroup <- list()
Edges_perGroup <- list()
for (algo in names(Nodes)) {
  for(nam in names(Nodes[[algo]])){
    for(cell in names(Nodes[[algo]][[nam]])){
      Nodes_perGroup[[algo]][[nam]][[cell]] <- Nodes[[algo]][[nam]][[cell]][names(Nodes[[algo]][[nam]][[cell]]) %in% groups]
      Edges_perGroup[[algo]][[nam]][[cell]] <- Edges[[algo]][[nam]][[cell]][names(Edges[[algo]][[nam]][[cell]]) %in% groups]
    }
  }
}


#Plot venn diagram
#perGroup
for (algo in names(Nodes_perGroup)) {
  for(nam in names(Nodes_perGroup[[algo]])){
    for(cell in names(Nodes_perGroup[[algo]][[nam]])){
      venn_Nodes <- ggvenn(Nodes_perGroup[[algo]][[nam]][[cell]], show_percentage = TRUE) + 
        labs(title = paste0("Nodes")) +
        theme(plot.caption=element_text(hjust = 0.5, size=8))
      
      venn_edges <- ggvenn(Edges_perGroup[[algo]][[nam]][[cell]], show_percentage = TRUE) + 
        labs(title = paste0("Edges")) + 
        theme(plot.caption=element_text(hjust = 0.5, size=8))
      
      tmp <- ggarrange(venn_Nodes, venn_edges, labels = c("A", "B"), nrow = 1, ncol = 2)
      
      f <- annotate_figure(tmp, top = text_grob(paste0(algo, " || ", nam, " || ", cell), color = "blue", face = "bold", size = 14))
      
      save_png_pdf(f, paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/VennDiagram/", algo, "/", nam, "/", cell, "_VennDiagram"), 
                   height = 5, width = 14, res=250)
    }
  }
}
      
      
# #Extract intersection information in a txt file
# #perGroup
# for (algo in names(Nodes_perGroup)) {
#   for(nam in names(Nodes_perGroup[[algo]])){
#     for(cell in names(Nodes_perGroup[[algo]][[nam]])){
#       fileConn<-file(paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/VennDiagram/", algo, "/", nam, "/", cell, "_Nodes.txt"))
#       overlap_Nodes <- Reduce(intersect, Nodes_perGroup[[algo]][[nam]][[cell]])
#       #overlap_Nodes_txt <- cat(paste0("Shared Nodes: \n", paste0(overlap_Nodes, collapse = "\n")))
#       writeLines(overlap_Nodes, fileConn)
#       close(fileConn)
#       
#       fileConn<-file(paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/VennDiagram/", algo, "/", nam, "/", cell, "_edges.txt"))
#       overlap_edges <- Reduce(intersect, Edges_perGroup[[algo]][[nam]][[cell]])
#       #overlap_edges_txt <- cat(paste0("Shared Edges: \n", paste0(overlap_edges, collapse = "\n")))
#       writeLines(overlap_edges, fileConn)
#       close(fileConn)
#     }
#   }
# }



# #Extract overlapping and unique Nodes, and write to a table
# #perGroup
# for (algo in names(Nodes_perGroup)) {
#   for(nam in names(Nodes_perGroup[[algo]])){
#     for(cell in names(Nodes_perGroup[[algo]][[nam]])){
#       # Extract intersection for all
#       overlap_Nodes <- Reduce(intersect,Nodes_perGroup[[algo]][[nam]][[cell]])
#       Genes_overlap_Nodes <- paste(overlap_Nodes, collapse=', ')
#       length_overlap_Nodes <- length(overlap_Nodes)
#       
#       # Extract union for all
#       union_Nodes <- Reduce(union,Nodes_perGroup[[algo]][[nam]][[cell]])
#       
#       #Extract Nodes that exist only in normal Nodes 
#       union_excl_normal_Nodes <- union(union(Nodes_perGroup[[algo]][[nam]][[cell]][["viral"]], 
#                                              Nodes_perGroup[[algo]][[nam]][[cell]][["nonviral"]]), 
#                                        Nodes_perGroup[[algo]][[nam]][[cell]][["tumor"]])
#       only_normal_Nodes <- union_Nodes[which(!union_Nodes %in% union_excl_normal_Nodes)]
#       Genes_only_normal_Nodes <- paste(only_normal_Nodes, collapse=', ')
#       length_only_normal_Nodes <- length(only_normal_Nodes)
#       
#       #Extract Nodes that exist only in viral Nodes 
#       union_excl_viral_Nodes <- union(union(Nodes_perGroup[[algo]][[nam]][[cell]][["normal"]], 
#                                             Nodes_perGroup[[algo]][[nam]][[cell]][["nonviral"]]), 
#                                       Nodes_perGroup[[algo]][[nam]][[cell]][["tumor"]])
#       only_viral_Nodes <- union_Nodes[which(!union_Nodes %in% union_excl_viral_Nodes)]
#       Genes_only_viral_Nodes <- paste(only_viral_Nodes, collapse=', ')
#       length_only_viral_Nodes <- length(only_viral_Nodes)
#       
#       #Extract Nodes that exist only in nonviral Nodes 
#       union_excl_nonviral_Nodes <- union(union(Nodes_perGroup[[algo]][[nam]][[cell]][["normal"]], 
#                                                Nodes_perGroup[[algo]][[nam]][[cell]][["viral"]]), 
#                                          Nodes_perGroup[[algo]][[nam]][[cell]][["tumor"]])
#       only_nonviral_Nodes <- union_Nodes[which(!union_Nodes %in% union_excl_nonviral_Nodes)]
#       Genes_only_nonviral_Nodes <- paste(only_nonviral_Nodes, collapse=', ')
#       length_only_nonviral_Nodes <- length(only_nonviral_Nodes)
#       
#       #Extract Nodes that exist only in tumor Nodes 
#       union_excl_tumor_Nodes <- union(union(Nodes_perGroup[[algo]][[nam]][[cell]][["viral"]], 
#                                             Nodes_perGroup[[algo]][[nam]][[cell]][["nonviral"]]), 
#                                       Nodes_perGroup[[algo]][[nam]][[cell]][["normal"]])
#       only_tumor_Nodes <- union_Nodes[which(!union_Nodes %in% union_excl_tumor_Nodes)]
#       Genes_only_tumor_Nodes <- paste(only_tumor_Nodes, collapse=', ')
#       length_only_tumor_Nodes <- length(only_tumor_Nodes)
#       
#       
#       #Create table
#       Group <- c("intersection", "normal", "viral", "nonviral", "tumor")
#       Count <- c(length_overlap_Nodes, length_only_normal_Nodes, length_only_viral_Nodes, length_only_nonviral_Nodes, length_only_tumor_Nodes)
#       Genes <- c(Genes_overlap_Nodes, Genes_only_normal_Nodes, Genes_only_viral_Nodes, Genes_only_nonviral_Nodes, Genes_only_tumor_Nodes)
#       df <- data.frame(Group, Count, Genes)
#       
#       df %>% 
#         gt() %>%
#         tab_header(title = paste0(algo, " | ", cell, " | ", nam)) %>%
#         gtsave(paste0(cell, "_overlap&unique_Nodes_SummaryTable.html"), 
#                path = paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/VennDiagram/", algo, "/", nam, "/"))
#     }
#   }
# }



##################################################################################################
### COMPARISON : NORMAL vs VIRAL vs NONVIRAL
##################################################################################################

comp <- "comp_NORMALvsVIRALvsNONVIRAL"
groups <- c("normal", "viral", "nonviral") #Exclude samples, and keep only corresponding groups

Nodes_perGroup <- list()
Edges_perGroup <- list()
for (algo in names(Nodes)) {
  for(nam in names(Nodes[[algo]])){
    for(cell in names(Nodes[[algo]][[nam]])){
      Nodes_perGroup[[algo]][[nam]][[cell]] <- Nodes[[algo]][[nam]][[cell]][names(Nodes[[algo]][[nam]][[cell]]) %in% groups]
      Edges_perGroup[[algo]][[nam]][[cell]] <- Edges[[algo]][[nam]][[cell]][names(Edges[[algo]][[nam]][[cell]]) %in% groups]
    }
  }
}


#Plot venn diagram
#perGroup
for (algo in names(Nodes_perGroup)) {
  for(nam in names(Nodes_perGroup[[algo]])){
    for(cell in names(Nodes_perGroup[[algo]][[nam]])){
      venn_Nodes <- ggvenn(Nodes_perGroup[[algo]][[nam]][[cell]], show_percentage = TRUE) + 
        labs(title = paste0("Nodes")) +
        theme(plot.caption=element_text(hjust = 0.5, size=8))
      
      venn_edges <- ggvenn(Edges_perGroup[[algo]][[nam]][[cell]], show_percentage = TRUE) + 
        labs(title = paste0("Edges")) + 
        theme(plot.caption=element_text(hjust = 0.5, size=8))
      
      tmp <- ggarrange(venn_Nodes, venn_edges, labels = c("A", "B"), nrow = 1, ncol = 2)
      
      f <- annotate_figure(tmp, top = text_grob(paste0(algo, " || ", nam, " || ", cell), color = "blue", face = "bold", size = 14))
      
      save_png_pdf(f, paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/VennDiagram/", algo, "/", nam, "/", cell, "_", comp, "_VennDiagram"), 
                   height = 5, width = 10, res=250)
    }
  }
}

##################################################################################################
### COMPARISON : NORMAL vs TUMOR
##################################################################################################

comp <- "comp_NORMALvsTUMOR"
groups <- c("normal", "tumor") #Exclude samples, and keep only corresponding groups

Nodes_perGroup <- list()
Edges_perGroup <- list()
for (algo in names(Nodes)) {
  for(nam in names(Nodes[[algo]])){
    for(cell in names(Nodes[[algo]][[nam]])){
      Nodes_perGroup[[algo]][[nam]][[cell]] <- Nodes[[algo]][[nam]][[cell]][names(Nodes[[algo]][[nam]][[cell]]) %in% groups]
      Edges_perGroup[[algo]][[nam]][[cell]] <- Edges[[algo]][[nam]][[cell]][names(Edges[[algo]][[nam]][[cell]]) %in% groups]
    }
  }
}

#Plot venn diagram
#perGroup
for (algo in names(Nodes_perGroup)) {
  for(nam in names(Nodes_perGroup[[algo]])){
    for(cell in names(Nodes_perGroup[[algo]][[nam]])){
      venn_Nodes <- ggvenn(Nodes_perGroup[[algo]][[nam]][[cell]], show_percentage = TRUE) + 
        labs(title = paste0("Nodes")) +
        theme(plot.caption=element_text(hjust = 0.5, size=8))
      
      venn_edges <- ggvenn(Edges_perGroup[[algo]][[nam]][[cell]], show_percentage = TRUE) + 
        labs(title = paste0("Edges")) + 
        theme(plot.caption=element_text(hjust = 0.5, size=8))
      
      tmp <- ggarrange(venn_Nodes, venn_edges, labels = c("A", "B"), nrow = 1, ncol = 2)
      
      f <- annotate_figure(tmp, top = text_grob(paste0(algo, " || ", nam, " || ", cell), color = "blue", face = "bold", size = 14))
      
      save_png_pdf(f, paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/VennDiagram/", algo, "/", nam, "/", cell, "_", comp, "_VennDiagram"), 
                   height = 5, width = 10, res=250)
    }
  }
}


##################################################################################################
###Network Visualization with Top10 critical Nodes highlighted
##################################################################################################

#Load library
library(igraph)
library(bc3net)

##Create files if it does not exist
for (algo in names(int_metrics)) {
  for(nam in names(int_metrics[[algo]])){
    plotdir=paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Networks_Top10_Nodes_highlighted/", algo, "/", nam)
    if(!dir.exists(plotdir)){dir.create(plotdir,  recursive = T)}
  }
}


##################################################################################################
#Plot
for (algo in names(int_metrics)) {
  for(nam in names(int_metrics[[algo]])){
    for(group in names(int_metrics[[algo]][[nam]])){
      for(cell in names(int_metrics[[algo]][[nam]][[group]])){
        net <- net_igraph[[algo]][[nam]][[group]][[cell]]
        net_connected = getgcc(net) ##getgcc() extracts the giant CONNECTED COMPONENT of the network and returns the igraph subnetwork. 
        
        data <- int_metrics[[algo]][[nam]][[group]][[cell]] #Get integrated centrality measures
        
        ToHighlight <- rownames(subset(data, rank <= 10)) #to color top10 nodes in red,and others in gray
        
        # ========================================================================================
        #All
        #The layout function layout.fruchterman.reingold() performs a force-based graph layout for a better visual distinction of the connectivity between nodes
        l = layout.fruchterman.reingold(net)
        #l = layout.fruchterman.reingold(net, niter=5000, area=30*vcount(net)^100000)
        #l <- layout.fruchterman.reingold(net, niter=5000, area=8*(vcount(net)^2),repulserad=(vcount(net)^10))
        #l = layout.norm(l)
        #other layouts to try: layout_with_lgl(net_connected, maxiter=20) , layout_in_circle(net_connected), layout_on_sphere(net_connected)
        
        png(paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Networks_Top10_Nodes_highlighted/", algo, "/", nam, "/", group, "_", cell, "_network.png"), 700, 700)
        #For the visualization, we consider the local point densities of the respective gene placements using the color functions colorRampPalette() and densCols().
        #color the nodes that are in the degree over 10 (magenta) and those that are not (green)
        par(mar=c(0,0,1,0))
        plot(net, 
             edge.width=1, edge.arrow.size=0.7, edge.curved = 0,
             vertex.label.cex=1, vertex.label.dist=1, vertex.label.font=2, vertex.label.color="black", 
             vertex.size=ifelse(V(net)$name%in%ToHighlight, 3, 1), 
             vertex.label=ifelse(V(net)$name%in%ToHighlight,V(net)$name,NA), 
             vertex.color=ifelse(V(net)$name%in%ToHighlight,"Red","Gray"), 
             layout=l)
        title(paste0(algo, " || ", nam, " || ", group, " || ", cell, " (Top10 nodes in red)"),cex.main=1,col.main="gray")
        dev.off()
        
        
        # ========================================================================================
        #CONNECTED COMPONENT of the network
        #The layout function layout.fruchterman.reingold() performs a force-based graph layout for a better visual distinction of the connectivity between nodes
        l = layout.fruchterman.reingold(net_connected)
        
        png(paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Networks_Top10_Nodes_highlighted/", algo, "/", nam, "/", group, "_", cell, "_connected_network.png"), 700, 700)
        #For the visualization, we consider the local point densities of the respective gene placements using the color functions colorRampPalette() and densCols().
        par(mar=c(0,0,1,0))
        plot(net_connected, 
             edge.width=1, edge.arrow.size=0.7, edge.curved = 0,
             vertex.label.cex=1, vertex.label.dist=1, vertex.label.font=2, vertex.label.color="black", 
             vertex.size=ifelse(V(net_connected)$name%in%ToHighlight, 3, 1), 
             vertex.label=ifelse(V(net_connected)$name%in%ToHighlight,V(net_connected)$name,NA), 
             vertex.color=ifelse(V(net_connected)$name%in%ToHighlight,"Red","Gray"), 
             layout=l)
        title(paste0(algo, " || ", nam, " || ", group, " || ", cell, "_connected", "(Top10 nodes in red)"),cex.main=1,col.main="gray")
        dev.off()
      }
    }
  }
}







