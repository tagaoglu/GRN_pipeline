# Set *.libPaths*  
.libPaths(c("~/p728/RSTUDIO/R/library/4.1/", .libPaths()))


#Set working directory and seed:
setwd("~/p728/RSTUDIO/analysis/tagaoglu")
set.seed(123)


## Package loading
library(tidyverse)
library(igraph)
library(ggpubr)
library(ggdendro)
library(ggalt)
library(gplots)
library(VennDiagram)
library(limma)
library(zoo)
library(DT)
library(grDevices)
library(extrafont)
library(purrr)
library(fmsb)


#Load functions
source("/home/tagaoglu/p728/RSTUDIO/analysis/tagaoglu/functions/centralityMetrics.R")

##################################################################################################
###Set parameters
##################################################################################################

nam_dir <- c("expData_mean_1","expData_mean_2")

#Set accordingly:Directed = T  or  Directed = F
GENIE_edge_prediction <- T
PIDC_edge_prediction <- F


##################################################################################################
###Load the networks (saved as edgelists)
##################################################################################################

#Put all in one list
algo_list <- c("GENIE", "GENIE_wRcis", "PIDC")

all <- sapply(algo_list,function(x) NULL)
all[["GENIE"]] <- readRDS(paste("~/p728/RSTUDIO/analysis/tagaoglu/data/networks/GENIE3/GENIE3_all_10k.Rds", sep = '' ))
all[["GENIE_wRcis"]] <- readRDS(paste("~/p728/RSTUDIO/analysis/tagaoglu/data/networks/GENIE3/GENIE3_all_10k_OnlyActivatedLinks_Rcis_5kb.Rds", sep = '' ))
all[["PIDC"]] <- readRDS(paste("~/p728/RSTUDIO/analysis/tagaoglu/data/networks/PIDC/PIDC_all_10k.Rds", sep = '' ))


#exclude random data
for(algo in names(all)){
  for(nam in names(all[[algo]])){
    all[[algo]] <- all[[algo]][names(all[[algo]]) %in% nam_dir]
  }
}


#Sort cell names of list objects to compare easily
for(algo in names(all)){
  for(nam in names(all[[algo]])){
    for(group in names(all[[algo]][[nam]])){
      all[[algo]][[nam]][[group]] <- all[[algo]][[nam]][[group]][order(names(all[[algo]][[nam]][[group]]))]
    }
  }
}


##################################################################################################

net_list_all <- list()
for(algo in names(all)){
  for(nam in names(all[[algo]])){
    for(group in names(all[[algo]][[nam]])){
      for(cell in names(all[[algo]][[nam]][[group]])){
        if(algo == "GENIE_wRcis"){
          net_list_all[[algo]][[nam]][[cell]][[group]] <- all[[algo]][[nam]][[group]][[cell]][["FilteredLinks"]]
        }else{
          net_list_all[[algo]][[nam]][[cell]][[group]] <- all[[algo]][[nam]][[group]][[cell]]
        }
      }
    }
  }
}


##################################################################################################
## Conversion to igraph objects
##################################################################################################

##Check next time
#net_list_all[["GENIE_wRcis"]][["expData_mean_1"]][["Hepatocytes"]][["P275"]] <- NULL
#net_list_all[["GENIE_wRcis"]][["expData_mean_2"]][["Hepatocytes"]][["normal"]] <- NULL
#net_list_all[["GENIE_wRcis"]][["expData_mean_2"]][["Hepatocytes"]][["P275"]] <- NULL
#net_list_all[["GENIE_wRcis"]][["expData_mean_2"]][["pDCs"]][["nonviral"]] <- NULL


###Create igraph object for each edgelist
net_graphs_all <- list()
for(algo in names(net_list_all)){
  for(nam in names(net_list_all[[algo]])){
    for(cell in names(net_list_all[[algo]][[nam]])){
        if(algo == "PIDC"){
          net_graphs_all[[algo]][[nam]][[cell]] <- map(net_list_all[[algo]][[nam]][[cell]], 
                                               ~ graph_from_edgelist(as.matrix(.[, 1:2]), directed = PIDC_edge_prediction)) # <= edge_prediction!
        }else{
          net_graphs_all[[algo]][[nam]][[cell]] <- map(net_list_all[[algo]][[nam]][[cell]], 
                                               ~ graph_from_edgelist(as.matrix(.[, 1:2]), directed = GENIE_edge_prediction)) # <= edge_prediction!
        }
    }
  }
}


##################################################################################################

##Create files if it does not exist
for(algo in names(net_graphs_all)){
  for(nam in names(net_graphs_all[[algo]])){
    for(cell in names(net_graphs_all[[algo]][[nam]])){
      plotdir=paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Network_analysis/", algo, "/", nam, "/", cell)
      if(!dir.exists(plotdir)){dir.create(plotdir,  recursive = T)}
    }
  }
}

##################################################################################################


##################################################################################################
###Network analysis
#Network centrality measures (CM) are node-centric, meaning that one can calculate each of the CM for each of the nodes. 
#These CM allow one to rank nodes from most to less important according to a given criterion.
## Centrality Measures:
#Degree
#Betweenness
#Closeness
#Eigenvector
#Pagerank
##################################################################################################


########################################################
#Compute the centralities for each node in each network:
net_k_all <- list()
net_bc_all <- list()
net_close_all <- list()
net_eigen_all <- list()
net_page_all <- list()
for(algo in names(net_graphs_all)){
  for(nam in names(net_graphs_all[[algo]])){
    for(cell in names(net_graphs_all[[algo]][[nam]])){
      # Degree centrality
      net_k_all[[algo]][[nam]][[cell]] <- map(net_graphs_all[[algo]][[nam]][[cell]], igraph::degree)
      # Betweenness centrality
      net_bc_all[[algo]][[nam]][[cell]] <- map(net_graphs_all[[algo]][[nam]][[cell]], betweenness)
      # Closeness centrality
      net_close_all[[algo]][[nam]][[cell]] <- map(net_graphs_all[[algo]][[nam]][[cell]], closeness)
      
      if(algo == "PIDC"){
        # Eigen centrality
        net_eigen_all[[algo]][[nam]][[cell]] <- map(net_graphs_all[[algo]][[nam]][[cell]], ~ eigen_centrality(.)$vector, directed = PIDC_edge_prediction) # <= edge_prediction
        # Pagerank centrality
        net_page_all[[algo]][[nam]][[cell]] <- map(net_graphs_all[[algo]][[nam]][[cell]], ~ page_rank(.)$vector, directed = PIDC_edge_prediction) # <= edge_prediction
      }else{
        # Eigen centrality
        net_eigen_all[[algo]][[nam]][[cell]] <- map(net_graphs_all[[algo]][[nam]][[cell]], ~ eigen_centrality(.)$vector, directed = GENIE_edge_prediction) # <= edge_prediction
        # Pagerank centrality
        net_page_all[[algo]][[nam]][[cell]] <- map(net_graphs_all[[algo]][[nam]][[cell]], ~ page_rank(.)$vector, directed = GENIE_edge_prediction) # <= edge_prediction
      }
    }
  }
}

#Compute integrated centrality metrics for each node in each network:
#Gives rank, the smaller the rank is, the better, because the more critical the gene is
net_int_all <- list()
for(algo in names(net_list_all)){
  for(nam in names(net_list_all[[algo]])){
    for(cell in names(net_list_all[[algo]][[nam]])){
      for(org in names(net_list_all[[algo]][[nam]][[cell]])){
        if(algo == "PIDC"){
          # Integrated centrality
          tmp <- centralityMetrics(net_list_all[[algo]][[nam]][[cell]][[org]], edge_prediction = PIDC_edge_prediction) # <= edge_prediction
          #Convert two columns of a data frame to a named vector
          net_int_all[[algo]][[nam]][[cell]][[org]] <- setNames(as.numeric(tmp[["integratedMetrics"]]$rank), tmp[["integratedMetrics"]]$Name)
          
        }else{
          # Integrated centrality
          tmp <- centralityMetrics(net_list_all[[algo]][[nam]][[cell]][[org]], edge_prediction = GENIE_edge_prediction) # <= edge_prediction
          #Convert two columns of a data frame to a named vector
          net_int_all[[algo]][[nam]][[cell]][[org]] <- setNames(as.numeric(tmp[["integratedMetrics"]]$rank), tmp[["integratedMetrics"]]$Name)
        }
      }
    }
  }
}




# Gather all metrics 
net_centr_l_all <- list()
for(algo in names(net_graphs_all)){
  for(nam in names(net_graphs_all[[algo]])){
    for(cell in names(net_graphs_all[[algo]][[nam]])){
      net_centr_l_all[[algo]][[nam]][[cell]] <- list(net_k_all[[algo]][[nam]][[cell]], net_bc_all[[algo]][[nam]][[cell]], 
                                                     net_close_all[[algo]][[nam]][[cell]], net_eigen_all[[algo]][[nam]][[cell]], 
                                                     net_page_all[[algo]][[nam]][[cell]], net_int_all[[algo]][[nam]][[cell]])
      names(net_centr_l_all[[algo]][[nam]][[cell]]) <- c("degree", "betweenness", 
                                                         "closeness", "eigen", 
                                                         "pagerank", "integrated_rank")
    }
  }
}




########################################################
### Correlation between CM
#We will use 5 different CM, each of which provides a different criterion for what a central node is. 
#As described in [this paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1184047/), there exists a significant correlation between CM 
#in biological networks. Specifically, it has been reported that nodes with a very high degree tend to have high betweenness and closeness. 
#However, nodes with a low degree have a large variance in betweenness and closeness, meaning that a gene can have a low degree 
#but a very high betweenness. For instance, these nodes are referred to as bottlenecks, and they are crucial to the flow of information 
#between networks modules. Thus, we will assess if we can depict similar correlations in our networks:

# arranged_gg_all: list of ggplot objects. Each object contains 4 scatterplots,
# corresponding to the correlation between degree and the other 4 CM for a 
# subset of 4 groups (viral, nonviral, normal, and tumor).

centr_interest <- c("closeness", "pagerank", "eigen", "betweenness")

arranged_gg_all <- list()

for(algo in names(net_centr_l_all)){
  for(nam in names(net_centr_l_all[[algo]])){
    for(cell in names(net_centr_l_all[[algo]][[nam]])){
      
        groups_interest <- names(net_graphs_all[[algo]][[nam]][[cell]])
        
        for (org in groups_interest) {
          org_gg_all <- list()
          for (centr in centr_interest) {
            org_gg_all[[algo]][[nam]][[cell]][[centr]] <- net_centr_l_all[[algo]][[nam]][[cell]] %>% 
              map(org) %>% 
              bind_cols() %>% 
              ggplot(aes_string("degree", centr)) +
              geom_point(shape = 1, color = "#0071bd", size = 1.5, alpha = 0.5) +
              scale_x_log10("DEGREE") +
              scale_y_log10(str_to_upper(centr)) +
              theme_classic() +
              theme(axis.title = element_text(size = 10))
          }
          arranged_gg_all[[algo]][[nam]][[cell]][[org]] <- ggarrange(plotlist = org_gg_all[[algo]][[nam]][[cell]], nrow = 2, ncol = 2)
          arranged_gg_all[[algo]][[nam]][[cell]][[org]] <- annotate_figure(arranged_gg_all[[algo]][[nam]][[cell]][[org]], top = text_grob(org, size = 18))
        }
      
    }
  }
}


# Save each plot object in arranged_gg
for(algo in names(arranged_gg_all)){
  for(nam in names(arranged_gg_all[[algo]])){
    for(cell in names(arranged_gg_all[[algo]][[nam]])){
      for(gg in names(arranged_gg_all[[algo]][[nam]][[cell]])){
        plot = arranged_gg_all[[algo]][[nam]][[cell]][[gg]]
        ggsave(paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Network_analysis/", algo, "/", nam, "/", cell, "/correlation_CM_", gg, ".pdf"), 
               plot = plot, 
               device = "pdf", 
               height = 10, 
               width = 18, 
               units = "cm")
      }
    }  
  }
}



### RESULTS: Indeed, wee see how the pattern we expected: correlated CM, with low-degree nodes having a larger variability of other centralities.



########################################################
## Definition of central genes
#By ranking nodes by decreasing centrality measure we can get an idea of which genes are the most important for that particular system. 
#We will define as "central" genes those within the top 20% of the ranking, as performed [here](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.0030059): 

#As for integrated metrics, it is already showing rank from the most central to the least central nodes (so in increasing order)

# top_centr_all: list of lists containing gene symbols of the top 20%  
# central nodes for each group (inner list) and centrality (outer list)

top_centr_all <- list()
for(algo in names(net_centr_l_all)){
  for(nam in names(net_centr_l_all[[algo]])){
    for(cell in names(net_centr_l_all[[algo]][[nam]])){
      for (centr in names(net_centr_l_all[[algo]][[nam]][[cell]])) { 
        if(centr == "integrated_rank"){
          top_centr_all[[algo]][[nam]][[cell]][[centr]] <- net_centr_l_all[[algo]][[nam]][[cell]][[centr]] %>% 
            map(sort, decreasing = FALSE) %>%
            map(~ names(.[1:(length(.) * 0.2)]))
        }else{
          top_centr_all[[algo]][[nam]][[cell]][[centr]] <- net_centr_l_all[[algo]][[nam]][[cell]][[centr]] %>% 
            map(sort, decreasing = TRUE) %>%
            map(~ names(.[1:(length(.) * 0.2)]))
        }
      }
    }
  }
}

saveRDS(top_centr_all, paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Network_analysis/networks_hub_lists.rds"))


### RESULTS: We now have sets of hubs


#Create the same without integrated metrics for some of the further steps
centralities <- c("degree", "betweenness", "closeness", "eigen", "pagerank")

top_centr_all_no_int <- list()
for(algo in names(net_centr_l_all)){
  for(nam in names(net_centr_l_all[[algo]])){
    for(cell in names(net_centr_l_all[[algo]][[nam]])){
      for (centr in centralities) { 
          top_centr_all_no_int[[algo]][[nam]][[cell]][[centr]] <- net_centr_l_all[[algo]][[nam]][[cell]][[centr]] %>% 
            map(sort, decreasing = TRUE) %>%
            map(~ names(.[1:(length(.) * 0.2)]))
      }
    }
  }
}




########################################################
## Different CM identify different sets of central genes
#In agreement with the previous scatter plots, we expect that there are genes that are ranked as central by more than one CM; 
#as well as others that are CM-specific. We can visualize this using Venn diagrams:

# Create venn diagrams (one per group), depicting the overlapping between
# each set of hubs (one per centrality). Save them in pdf format in the plots/Network_analysis folder

# for(algo in names(net_graphs_all)){
#   for(nam in names(net_graphs_all[[algo]])){
#     for(cell in names(net_graphs_all[[algo]][[nam]])){
#       groups_interest <- names(net_graphs_all[[algo]][[nam]][[cell]])
#     }
#   }
# }

loadfonts()

hubs_all <- list()
venn_centr_all <- list()
for(algo in names(net_graphs_all)){
  for(nam in names(net_graphs_all[[algo]])){
    for(cell in names(net_graphs_all[[algo]][[nam]])){
      for(org in names(net_graphs_all[[algo]][[nam]][[cell]])){
        hubs_all[[algo]][[nam]][[cell]][[org]] <- map(top_centr_all_no_int[[algo]][[nam]][[cell]], org)
        venn_centr_all[[algo]][[nam]][[cell]][[org]] <- venn.diagram(hubs_all[[algo]][[nam]][[cell]][[org]], 
                                                              fill = 2:6, alpha = 0.3, filename = NULL,
                                                              cat.just = list(c(0.6,1), c(0,0), c(0,0), c(1,1), c(1,0)),
                                                              main = org)
        venn_file <- str_c(paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Network_analysis/", algo, "/", nam, "/", cell, "/venn_cm_", org, ".pdf"))
        pdf(file = venn_file, height = 5, width = 5)
        grid.draw(venn_centr_all[[algo]][[nam]][[cell]][[org]])
        dev.off()
      }
      
      
    }
  }
}


### RESULTS: Indeed, we see how there are some genes that are in the top 20% nodes for all CM. Moreover, some genes are specific for each CM.



########################################################
## Single-cell derived regulatory networks are scale-free

#As reviewed by [Barabási AL, et al. (2004)](https://www.nature.com/articles/nrg1272), 
#biological networks are scale-free, as their degree distribution follows a power-law. 
#In other words, there is both a high proportion of nodes with a very low degree (few connections) 
#and a low proportion of nodes with high degree (referred to as hubs). 
#This makes biological networks robust against random perturbation, 
#as the whole system will fail only if hubs are compromised, which is unlikely.

#Thus, something we aim to test to validate our networks is whether they are scale-free. 
#To that end, we will use a Kolmogorov-Smirnov test to compare the degree distribution of our networks 
#to a theoretical power-law distribuiton, with the null hypothesis that they are equal. 
#Hence, p-values > 0.05 will be indicative of scale-free networks. Moreover, we will also compute the ⍺, that is, 
#the exponent of the power-law. It has been reported that biological networks have exponents between 2 and 3. 



## Visualizing scale-free distribution: histograms

#Visualize the degree distribution with a histogram:

# net_k_df_all: data.frame that contains the degree for all nodes in all groups.
net_k_df_all <- list()
for(algo in names(net_k_all)){
  for(nam in names(net_k_all[[algo]])){
    for(cell in names(net_k_all[[algo]][[nam]])){
      net_k_df_all[[algo]][[nam]][[cell]] <- net_k_all[[algo]][[nam]][[cell]] %>% 
        map(as.data.frame) %>% 
        bind_rows(.id = "group") %>% 
        set_names(c("group", "degree"))
    }
  }
}


# k_hists_all: list with histograms that show the degree distribution in every network.
k_hists_all <- list()
for(algo in names(net_k_df_all)){
  for(nam in names(net_k_df_all[[algo]])){
    for(cell in names(net_k_df_all[[algo]][[nam]])){
      
        groups_interest <- names(net_graphs_all[[algo]][[nam]][[cell]]) 
        
        k_hists_all[[algo]][[nam]][[cell]] <- groups_interest %>% map(~ filter(net_k_df_all[[algo]][[nam]][[cell]], group == .) %>% 
                                                                        ggplot(aes(x = degree)) +
                                                                        geom_histogram(binwidth = 15, color = "black") +
                                                                        ggtitle(str_to_title(.)) +
                                                                        scale_x_continuous("Degree", expand = c(0.025, 0.025)) +
                                                                        scale_y_continuous("Nodes", expand = c(0.025, 0.025)) +
                                                                        theme_classic() +
                                                                        theme(axis.title.x = element_text(size = 11),
                                                                              plot.title = element_text(hjust = 0.5, vjust = 0, size = 14, face = "bold")) 
        )
        names(k_hists_all[[algo]][[nam]][[cell]]) <- groups_interest
      
    }
  }
}


# Save plots: k_hists as a single graph with all the plots distributed in 2 rows and 2 cols.
k_hists_out_all <- list()
for(algo in names(k_hists_all)){
  for(nam in names(k_hists_all[[algo]])){
    for(cell in names(k_hists_all[[algo]][[nam]])){
        k_hists_out_all[[algo]][[nam]][[cell]] <- ggarrange(plotlist = k_hists_all[[algo]][[nam]][[cell]], nrow = 5, ncol = 4, vjust = 0)
        #k_hists_out_all[[algo]][[nam]][[cell]] <- annotate_figure(k_hists_out_all[[algo]][[nam]][[cell]], top = text_grob(cell, size = 15))
        ggsave(paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Network_analysis/", algo, "/", nam, "/", cell, "/scale_free_histograms", ".pdf"), 
               plot = k_hists_out_all[[algo]][[nam]][[cell]], 
               device = "pdf", 
               height = 24, 
               width = 18, 
               units = "cm")
    }
  }
}



########################################################
## Power-law fit to the degree distibution: scatter plots

#With both axes in log-scale the power-law fit is visualized as a decreasing line. This time, we also include the p-value and alpha:
  
# net_k_df2_all: data.frame that contains the degree distribution for all groups
# degree distribution (P(k)): probability that a selected node has exactly k links
net_k_distr_all <- list()
len_dist_all <- list()
net_k_df2_all <- list()
for(algo in names(net_graphs_all)){
  for(nam in names(net_graphs_all[[algo]])){
    for(cell in names(net_graphs_all[[algo]][[nam]])){
      net_k_distr_all[[algo]][[nam]][[cell]] <- map(net_graphs_all[[algo]][[nam]][[cell]], degree_distribution)
      len_dist_all[[algo]][[nam]][[cell]] <- map_dbl(net_k_distr_all[[algo]][[nam]][[cell]], length)   
      net_k_df2_all[[algo]][[nam]][[cell]] <- list(len_dist_all[[algo]][[nam]][[cell]], net_k_distr_all[[algo]][[nam]][[cell]]) %>% 
        pmap(~ cbind(0:(..1 - 1), ..2)) %>% 
        map(as.data.frame) %>% 
        bind_rows(.id = "group") %>% 
        set_names(c("group", "degree", "freq")) %>% 
        filter(degree != 0)
    }
  }
}


# ks_test_all: list with the results of the power law fit to the all degree distributions. 
#KS.p contains the p-value of the Kolmogorov-Smirnov test, and alpha contains the exponent of the power law: P(k) = k^-alpha
ks_test_all <- list()
for(algo in names(net_k_all)){
  for(nam in names(net_k_all[[algo]])){
    for(cell in names(net_k_all[[algo]][[nam]])){
      ks_test_all[[algo]][[nam]][[cell]] <- map(net_k_all[[algo]][[nam]][[cell]], fit_power_law)
    }
  }
}


# k_logs_gg_all: list with scatter plots, showing the degree distributions in log-log scale.
k_logs_gg_all <- list()
for(algo in names(net_k_df2_all)){
  for(nam in names(net_k_df2_all[[algo]])){
    for(cell in names(net_k_df2_all[[algo]][[nam]])){
      
      groups_interest <- names(net_graphs_all[[algo]][[nam]][[cell]])
      
      get_ks_stat <- function(group, stat){
        # Returns the desired statistic for the power law fit to the degree 
       # distribution of the network of a certain group
       # 
       # Args:
       #   group: string specifies the group of interest
       #   stat: string of the statistic of interest (p-value or alpha)
        stat1 <- ifelse(stat == "p-value", "KS.p", "alpha")
        stat2 <- formatC(ks_test_all[[algo]][[nam]][[cell]][[group]][[stat1]], digits = 2)
        ifelse(stat == "p-value", "KS.p", "\u03B1") %>% 
          str_c(" = ", stat2, collapse = TRUE)
      }
      
      x_labels <- expression("10", "10"^2)
      y_labels <- expression("10"^-3, "10"^-2, "10"^-1) 
      
      
      k_logs_gg_all[[algo]][[nam]][[cell]] <- groups_interest %>% 
        map(~ filter(net_k_df2_all[[algo]][[nam]][[cell]], group == .) %>% 
              ggplot(aes(degree, freq)) +
              geom_point(color = "#003CFF", size = 1) +
              geom_smooth(method = "lm", se = FALSE, color = "#003CFF") +
              annotate("text", 
                       x = 2, #11
                       y = 0.3, 
                       label = str_c(get_ks_stat(., "p-value"), 
                                     "\n",
                                     get_ks_stat(., "alpha")),
                       size = 14, #4
                       fontface =2,
                       hjust = 0) +
              scale_x_log10("Degree", expand = c(0,0), breaks = c(10, 100), 
                            labels = x_labels, limits = c(1, 200)) +
              scale_y_log10("Density", breaks = c(0.001, 0.01, 0.1), 
                            labels = y_labels, limits = c(0.001, 1)) +
              expand_limits(x = 1, y = 0) +
              ggtitle(.) +
              theme_classic() + 
              theme(panel.grid.major = element_line(color = "lightgrey"),
                    panel.border = element_rect(colour = "black", fill = NA),
                    axis.title = element_text(size = 50)) #13
            )
      
      names(k_logs_gg_all[[algo]][[nam]][[cell]]) <- groups_interest
    }
  }
}


# Save each scatterplot as pdf file
for(algo in names(k_logs_gg_all)){
  for(nam in names(k_logs_gg_all[[algo]])){
    for(cell in names(k_logs_gg_all[[algo]][[nam]])){
      
      groups_interest <- names(net_graphs_all[[algo]][[nam]][[cell]])
      
      walk2(k_logs_gg_all[[algo]][[nam]][[cell]], groups_interest, function(gg, org) {
        ggsave(str_c(paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Network_analysis/", algo, "/", nam, "/", cell, "/scale_free_scater_", org, ".png")),
               plot = gg, 
               device = "png", 
               width = 15, 
               height = 15, 
               units = "cm"
               )
        })
      
    }
  }
}



##################################################################################################
### COMPARISON : NORMAL vs VIRAL vs NONVIRAL
##################################################################################################

comp <- "comp_NORMALvsVIRALvsNONVIRAL"
groups_interest <- c("normal", "viral", "nonviral")


#Extract only groups_interest for correct comparison
top_centr_all_tmp <- list()
for(algo in names(top_centr_all)){
  for(nam in names(top_centr_all[[algo]])){
    for(cell in names(top_centr_all[[algo]][[nam]])){
      for (centr in names(top_centr_all[[algo]][[nam]][[cell]])) {
        #Extract only groups_interest for correct comparison
        for(org in groups_interest){
          top_centr_all_tmp[[algo]][[nam]][[cell]][[centr]][[org]] <- top_centr_all[[algo]][[nam]][[cell]][[centr]][[org]]
        }
      }
    }
  }
}


########################################################
# Multiplicity of the top central nodes

#Another essential question we aim to answer is whether these hubs are specific to a particular group (group-specific hubs) 
#or whether they are hubs in several groups (ubiquitous hubs). To assess this question, let us compute and plot the multiplicity of each hub 
#(i.e. # groups a given gene acts as a hub):


##  Inspect multiplicity

# hub_multi_df_all: df that contains, for each group and CM, the percentage of hubs that 
# have a multiplicity of 1, 2 or 3+ (i.e. are hubs in 3 or more groups).

centralities <- c("degree", "betweenness", "closeness", "eigen", "pagerank", "integrated_rank")

hub_multi_df_all <- list()
for(algo in names(top_centr_all_tmp)){
  for(nam in names(top_centr_all_tmp[[algo]])){
    for(cell in names(top_centr_all_tmp[[algo]][[nam]])){
      hub_multi_df_all[[algo]][[nam]][[cell]] <- data.frame(
        centrality = c(), 
        group = c(), 
        multiplicity = c(), 
        percentage = c()
        )
      
      for (centr in centralities) {
        centr_counts <- table(unlist(top_centr_all_tmp[[algo]][[nam]][[cell]][[centr]]))  
        
        for (org in groups_interest) {
          curr_hubs <- top_centr_all_tmp[[algo]][[nam]][[cell]][[centr]][[org]]
          multipl <- table(centr_counts[curr_hubs]) 
          mult_perc <- c(
            "1" = multipl[1] / sum(multipl) * 100,
            "2" = multipl[2] / sum(multipl) * 100,
            "3" = multipl[3] / sum(multipl) * 100   #"3+" = sum(multipl[3:length(multipl)]) / sum(multipl) * 100
          )
          curr_mult_df <- data.frame(
            centrality = rep(centr, 3), 
            group = rep(org, 3), 
            multiplicity = c("1", "2", "3"), 
            percentage = mult_perc
          )
          hub_multi_df_all[[algo]][[nam]][[cell]] <- bind_rows(hub_multi_df_all[[algo]][[nam]][[cell]], curr_mult_df)
        }
        
      }
      
    }
  }
}


# hub_multi_gg_all: stacked bar plot showing the multiplicity distribition per group and centrality.

hub_multi_gg_all <- list()
for(algo in names(hub_multi_df_all)){
  for(nam in names(hub_multi_df_all[[algo]])){
    for(cell in names(hub_multi_df_all[[algo]][[nam]])){
      bar_colors <- c("#ff8533", "#599ad3", "#9e66ab")
      hub_multi_df_all[[algo]][[nam]][[cell]]$multiplicity <- factor(hub_multi_df_all[[algo]][[nam]][[cell]]$multiplicity,
                                                                     levels = c("3", "2", "1")) 
      hub_multi_gg_all[[algo]][[nam]][[cell]] <- hub_multi_df_all[[algo]][[nam]][[cell]] %>% 
        ggplot(aes(group, percentage, 
                   fill = multiplicity)) +
        geom_bar(stat = "identity", colour = "black") +
        scale_x_discrete(name = "") +
        scale_y_continuous(name = "Percentage", expand = c(0,0)) +
        scale_fill_manual(name = "multiplicity", values = bar_colors) +
        facet_grid(. ~ centrality) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                         vjust = 0.5, lineheight = 1,
                                         size = 8),
              axis.title.y = element_text(size = 10),
              legend.title = element_text(size = 10),
              legend.text = element_text(size = 8))
      
    }
  }
}
      

# Save the hub_multiplicity_gg plot as a pdf file and the hub_multi_df as csv
for(algo in names(hub_multi_gg_all)){
  for(nam in names(hub_multi_gg_all[[algo]])){
    for(cell in names(hub_multi_gg_all[[algo]][[nam]])){
      ggsave(
        paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Network_analysis/", algo, "/", nam, "/", cell, "/hub_multiplicity_", comp, ".pdf"), 
        plot = hub_multi_gg_all[[algo]][[nam]][[cell]], 
        device = "pdf", 
        width = 7, 
        height = 3
      )
    }
  }
}



########################################################
## Define group-specific and shared hubs.
#Based on the previous stacked bar plot, we can see that, indeed, there are hubs that are specific for each group, 
#while others are ubiquitous. Let us separate them for downstream analysis:

# Separate ubiquitious (multiplicity == 3+; top_centr_ub) from group-specific hubs (multiplicity = 1; top_centr_sp). 
# Both top_centr_ub and top_centr_sp are lists of lists (outer list = CM, inner list = groups).

top_all_all <- list()
top_centr_ub_all <- list()
top_centr_sp_all <- list()
for(algo in names(top_centr_all_tmp)){
  for(nam in names(top_centr_all_tmp[[algo]])){
    for(cell in names(top_centr_all_tmp[[algo]][[nam]])){
      for (centr in names(top_centr_all_tmp[[algo]][[nam]][[cell]])) {
        
        top_all_all[[algo]][[nam]][[cell]] <- table(unlist(top_centr_all_tmp[[algo]][[nam]][[cell]][[centr]]))
        top_centr_ub_all[[algo]][[nam]][[cell]][[centr]] <- map(top_centr_all_tmp[[algo]][[nam]][[cell]][[centr]], ~ .[top_all_all[[algo]][[nam]][[cell]][.] >= 3])
        top_centr_sp_all[[algo]][[nam]][[cell]][[centr]] <- map(top_centr_all_tmp[[algo]][[nam]][[cell]][[centr]], ~ .[top_all_all[[algo]][[nam]][[cell]][.] == 1])
      
      }
    }
  }
}

saveRDS(top_centr_ub_all, paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Network_analysis/networks_hub_ubiquitous_", comp, ".rds"))
saveRDS(top_centr_sp_all, paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Network_analysis/networks_hub_specific_", comp, ".rds"))






##################################################################################################
### COMPARISON : NORMAL vs TUMOR
##################################################################################################

comp <- "comp_NORMALvsTUMOR"
groups_interest <- c("normal", "tumor")


#Extract only groups_interest for correct comparison
top_centr_all_tmp <- list()
for(algo in names(top_centr_all)){
  for(nam in names(top_centr_all[[algo]])){
    for(cell in names(top_centr_all[[algo]][[nam]])){
      for (centr in names(top_centr_all[[algo]][[nam]][[cell]])) {
        #Extract only groups_interest for correct comparison
        for(org in groups_interest){
          top_centr_all_tmp[[algo]][[nam]][[cell]][[centr]][[org]] <- top_centr_all[[algo]][[nam]][[cell]][[centr]][[org]]
        }
      }
    }
  }
}


########################################################
# Multiplicity of the top central nodes

#Another essential question we aim to answer is whether these hubs are specific to a particular group (group-specific hubs) 
#or whether they are hubs in several groups (ubiquitous hubs). To assess this question, let us compute and plot the multiplicity of each hub 
#(i.e. # groups a given gene acts as a hub):


##  Inspect multiplicity

# hub_multi_df_all: df that contains, for each group and CM, the percentage of hubs that 
# have a multiplicity of 1, 2 or 3+ (i.e. are hubs in 3 or more groups).

centralities <- c("degree", "betweenness", "closeness", "eigen", "pagerank", "integrated_rank")

hub_multi_df_all <- list()
for(algo in names(top_centr_all_tmp)){
  for(nam in names(top_centr_all_tmp[[algo]])){
    for(cell in names(top_centr_all_tmp[[algo]][[nam]])){
      hub_multi_df_all[[algo]][[nam]][[cell]] <- data.frame(
        centrality = c(), 
        group = c(), 
        multiplicity = c(), 
        percentage = c()
      )
      
      for (centr in centralities) {
        centr_counts <- table(unlist(top_centr_all_tmp[[algo]][[nam]][[cell]][[centr]]))  
        
        for (org in groups_interest) {
          curr_hubs <- top_centr_all_tmp[[algo]][[nam]][[cell]][[centr]][[org]]
          multipl <- table(centr_counts[curr_hubs]) 
          mult_perc <- c(
            "1" = multipl[1] / sum(multipl) * 100,
            "2" = multipl[2] / sum(multipl) * 100
          )
          curr_mult_df <- data.frame(
            centrality = rep(centr, 2), 
            group = rep(org, 2), 
            multiplicity = c("1", "2"), 
            percentage = mult_perc
          )
          hub_multi_df_all[[algo]][[nam]][[cell]] <- bind_rows(hub_multi_df_all[[algo]][[nam]][[cell]], curr_mult_df)
        }
        
      }
      
    }
  }
}


# hub_multi_gg_all: stacked bar plot showing the multiplicity distribition per group and centrality.

hub_multi_gg_all <- list()
for(algo in names(hub_multi_df_all)){
  for(nam in names(hub_multi_df_all[[algo]])){
    for(cell in names(hub_multi_df_all[[algo]][[nam]])){
      bar_colors <- c("#599ad3", "#9e66ab") #"#ff8533"
      hub_multi_df_all[[algo]][[nam]][[cell]]$multiplicity <- factor(hub_multi_df_all[[algo]][[nam]][[cell]]$multiplicity,
                                                                     levels = c("2", "1")) 
      hub_multi_gg_all[[algo]][[nam]][[cell]] <- hub_multi_df_all[[algo]][[nam]][[cell]] %>% 
        ggplot(aes(group, percentage, 
                   fill = multiplicity)) +
        geom_bar(stat = "identity", colour = "black") +
        scale_x_discrete(name = "") +
        scale_y_continuous(name = "Percentage", expand = c(0,0)) +
        scale_fill_manual(name = "multiplicity", values = bar_colors) +
        facet_grid(. ~ centrality) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                         vjust = 0.5, lineheight = 1,
                                         size = 8),
              axis.title.y = element_text(size = 10),
              legend.title = element_text(size = 10),
              legend.text = element_text(size = 8))
      
    }
  }
}


# Save the hub_multiplicity_gg plot as a pdf file and the hub_multi_df as csv
for(algo in names(hub_multi_gg_all)){
  for(nam in names(hub_multi_gg_all[[algo]])){
    for(cell in names(hub_multi_gg_all[[algo]][[nam]])){
      ggsave(
        paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Network_analysis/", algo, "/", nam, "/", cell, "/hub_multiplicity_", comp, ".pdf"), 
        plot = hub_multi_gg_all[[algo]][[nam]][[cell]], 
        device = "pdf", 
        width = 7, 
        height = 3
      )
    }
  }
}



########################################################
## Define group-specific and shared hubs.
#Based on the previous stacked bar plot, we can see that, indeed, there are hubs that are specific for each group, 
#while others are ubiquitous. Let us separate them for downstream analysis:

# Separate ubiquitious (multiplicity == 3+; top_centr_ub) from group-specific hubs (multiplicity = 1; top_centr_sp). 
# Both top_centr_ub and top_centr_sp are lists of lists (outer list = CM, inner list = groups).

top_all_all <- list()
top_centr_ub_all <- list()
top_centr_sp_all <- list()
for(algo in names(top_centr_all_tmp)){
  for(nam in names(top_centr_all_tmp[[algo]])){
    for(cell in names(top_centr_all_tmp[[algo]][[nam]])){
      for (centr in names(top_centr_all_tmp[[algo]][[nam]][[cell]])) {
        
        top_all_all[[algo]][[nam]][[cell]] <- table(unlist(top_centr_all_tmp[[algo]][[nam]][[cell]][[centr]]))
        top_centr_ub_all[[algo]][[nam]][[cell]][[centr]] <- map(top_centr_all_tmp[[algo]][[nam]][[cell]][[centr]], ~ .[top_all_all[[algo]][[nam]][[cell]][.] >= 2])
        top_centr_sp_all[[algo]][[nam]][[cell]][[centr]] <- map(top_centr_all_tmp[[algo]][[nam]][[cell]][[centr]], ~ .[top_all_all[[algo]][[nam]][[cell]][.] == 1])
        
      }
    }
  }
}

saveRDS(top_centr_ub_all, paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Network_analysis/networks_hub_ubiquitous_", comp, ".rds"))
saveRDS(top_centr_sp_all, paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Network_analysis/networks_hub_specific_", comp, ".rds"))





##################################################################################################
### COMPARISON : hep NORMAL vs SAMPLES
##################################################################################################


comp <- "comp_hep_NORMALvsSAMPLES"
groups_interest <- c("normal", "P123s2", "P123s8", "P193", "P199", "P108", "P191", "P207", "P207Gel", "P215", "P220", "P275", "P282C", "P295")

cells <- c("Hepatocytes")

#Extract only groups_interest for correct comparison
top_centr_all_tmp <- list()
for(algo in names(top_centr_all)){
  for(nam in names(top_centr_all[[algo]])){
    for(cell in cells){
      for (centr in names(top_centr_all[[algo]][[nam]][[cell]])) {
        #Extract only groups_interest for correct comparison
        for(org in groups_interest){
          top_centr_all_tmp[[algo]][[nam]][[cell]][[centr]][[org]] <- top_centr_all[[algo]][[nam]][[cell]][[centr]][[org]]
        }
      }
    }
  }
}


########################################################
# Multiplicity of the top central nodes

#Another essential question we aim to answer is whether these hubs are specific to a particular group (group-specific hubs) 
#or whether they are hubs in several groups (ubiquitous hubs). To assess this question, let us compute and plot the multiplicity of each hub 
#(i.e. # groups a given gene acts as a hub):


##  Inspect multiplicity

# hub_multi_df_all: df that contains, for each group and CM, the percentage of hubs that 
# have a multiplicity of 1, 2 or 3+ (i.e. are hubs in 3 or more groups).

centralities <- c("degree", "betweenness", "closeness", "eigen", "pagerank", "integrated_rank")

hub_multi_df_all <- list()
for(algo in names(top_centr_all_tmp)){
  for(nam in names(top_centr_all_tmp[[algo]])){
    for(cell in names(top_centr_all_tmp[[algo]][[nam]])){
      hub_multi_df_all[[algo]][[nam]][[cell]] <- data.frame(
        centrality = c(), 
        group = c(), 
        multiplicity = c(), 
        percentage = c()
      )
      
      for (centr in centralities) {
        centr_counts <- table(unlist(top_centr_all_tmp[[algo]][[nam]][[cell]][[centr]]))  
        
        for (org in groups_interest) {
          curr_hubs <- top_centr_all_tmp[[algo]][[nam]][[cell]][[centr]][[org]]
          multipl <- table(centr_counts[curr_hubs]) 
          mult_perc <- c(
            "1" = multipl[1] / sum(multipl) * 100,
            "2" = multipl[2] / sum(multipl) * 100,
            "3+" = sum(multipl[3:length(multipl)]) / sum(multipl) * 100   
          )
          curr_mult_df <- data.frame(
            centrality = rep(centr, 3), 
            group = rep(org, 3), 
            multiplicity = c("1", "2", "3+"), 
            percentage = mult_perc
          )
          hub_multi_df_all[[algo]][[nam]][[cell]] <- bind_rows(hub_multi_df_all[[algo]][[nam]][[cell]], curr_mult_df)
        }
        
      }
      
    }
  }
}


# hub_multi_gg_all: stacked bar plot showing the multiplicity distribition per group and centrality.

hub_multi_gg_all <- list()
for(algo in names(hub_multi_df_all)){
  for(nam in names(hub_multi_df_all[[algo]])){
    for(cell in names(hub_multi_df_all[[algo]][[nam]])){
      bar_colors <- c("#ff8533", "#599ad3", "#9e66ab")
      hub_multi_df_all[[algo]][[nam]][[cell]]$multiplicity <- factor(hub_multi_df_all[[algo]][[nam]][[cell]]$multiplicity,
                                                                     levels = c("3+", "2", "1")) 
      hub_multi_gg_all[[algo]][[nam]][[cell]] <- hub_multi_df_all[[algo]][[nam]][[cell]] %>% 
        ggplot(aes(group, percentage, 
                   fill = multiplicity)) +
        geom_bar(stat = "identity", colour = "black") +
        scale_x_discrete(name = "") +
        scale_y_continuous(name = "Percentage", expand = c(0,0)) +
        scale_fill_manual(name = "multiplicity", values = bar_colors) +
        facet_grid(. ~ centrality) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                         vjust = 0.5, lineheight = 1,
                                         size = 8),
              axis.title.y = element_text(size = 10),
              legend.title = element_text(size = 10),
              legend.text = element_text(size = 8))
      
    }
  }
}


# Save the hub_multiplicity_gg plot as a pdf file and the hub_multi_df as csv
for(algo in names(hub_multi_gg_all)){
  for(nam in names(hub_multi_gg_all[[algo]])){
    for(cell in names(hub_multi_gg_all[[algo]][[nam]])){
      ggsave(
        paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Network_analysis/", algo, "/", nam, "/", cell, "/hub_multiplicity_", comp, ".pdf"), 
        plot = hub_multi_gg_all[[algo]][[nam]][[cell]], 
        device = "pdf", 
        width = 14, 
        height = 5
      )
    }
  }
}



########################################################
## Define group-specific and shared hubs.
#Based on the previous stacked bar plot, we can see that, indeed, there are hubs that are specific for each group, 
#while others are ubiquitous. Let us separate them for downstream analysis:

# Separate ubiquitious (multiplicity == 3+; top_centr_ub) from group-specific hubs (multiplicity = 1; top_centr_sp). 
# Both top_centr_ub and top_centr_sp are lists of lists (outer list = CM, inner list = groups).

top_all_all <- list()
top_centr_ub_all <- list()
top_centr_sp_all <- list()
for(algo in names(top_centr_all_tmp)){
  for(nam in names(top_centr_all_tmp[[algo]])){
    for(cell in names(top_centr_all_tmp[[algo]][[nam]])){
      for (centr in names(top_centr_all_tmp[[algo]][[nam]][[cell]])) {
        
        top_all_all[[algo]][[nam]][[cell]] <- table(unlist(top_centr_all_tmp[[algo]][[nam]][[cell]][[centr]]))
        top_centr_ub_all[[algo]][[nam]][[cell]][[centr]] <- map(top_centr_all_tmp[[algo]][[nam]][[cell]][[centr]], ~ .[top_all_all[[algo]][[nam]][[cell]][.] >= 3])
        top_centr_sp_all[[algo]][[nam]][[cell]][[centr]] <- map(top_centr_all_tmp[[algo]][[nam]][[cell]][[centr]], ~ .[top_all_all[[algo]][[nam]][[cell]][.] == 1])
        
      }
    }
  }
}

saveRDS(top_centr_ub_all, paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Network_analysis/networks_hub_ubiquitous_", comp, ".rds"))
saveRDS(top_centr_sp_all, paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Network_analysis/networks_hub_specific_", comp, ".rds"))




##################################################################################################
### COMPARISON : hep NORMAL vs EACH SAMPLE separately
##################################################################################################

cells <- c("Hepatocytes")

samples <- c("P123s2", "P123s8", "P193", "P199", "P108", "P191", "P207", "P207Gel", "P215", "P220", "P275", "P282C", "P295")

# Compare normal vs each sample separately
for(samp in samples){
  sample <- samp
  comp <- paste0("comp_hep_NORMALvs", sample)
  groups_interest <- c("normal", sample)
  
  
  
  #Extract only groups_interest for correct comparison
  top_centr_all_tmp <- list()
  for(algo in names(top_centr_all)){
    for(nam in names(top_centr_all[[algo]])){
      for(cell in cells){
        for (centr in names(top_centr_all[[algo]][[nam]][[cell]])) {
          #Extract only groups_interest for correct comparison
          for(org in groups_interest){
            top_centr_all_tmp[[algo]][[nam]][[cell]][[centr]][[org]] <- top_centr_all[[algo]][[nam]][[cell]][[centr]][[org]]
          }
        }
      }
    }
  }
  
  
  ########################################################
  # Multiplicity of the top central nodes
  
  #Another essential question we aim to answer is whether these hubs are specific to a particular group (group-specific hubs) 
  #or whether they are hubs in several groups (ubiquitous hubs). To assess this question, let us compute and plot the multiplicity of each hub 
  #(i.e. # groups a given gene acts as a hub):
  
  
  ##  Inspect multiplicity
  
  # hub_multi_df_all: df that contains, for each group and CM, the percentage of hubs that 
  # have a multiplicity of 1, 2 or 3+ (i.e. are hubs in 3 or more groups).
  
  centralities <- c("degree", "betweenness", "closeness", "eigen", "pagerank", "integrated_rank")
  
  hub_multi_df_all <- list()
  for(algo in names(top_centr_all_tmp)){
    for(nam in names(top_centr_all_tmp[[algo]])){
      for(cell in names(top_centr_all_tmp[[algo]][[nam]])){
        hub_multi_df_all[[algo]][[nam]][[cell]] <- data.frame(
          centrality = c(), 
          group = c(), 
          multiplicity = c(), 
          percentage = c()
        )
        
        for (centr in centralities) {
          centr_counts <- table(unlist(top_centr_all_tmp[[algo]][[nam]][[cell]][[centr]]))  
          
          for (org in groups_interest) {
            curr_hubs <- top_centr_all_tmp[[algo]][[nam]][[cell]][[centr]][[org]]
            multipl <- table(centr_counts[curr_hubs]) 
            mult_perc <- c(
              "1" = multipl[1] / sum(multipl) * 100,
              "2" = multipl[2] / sum(multipl) * 100
            )
            curr_mult_df <- data.frame(
              centrality = rep(centr, 2), 
              group = rep(org, 2), 
              multiplicity = c("1", "2"), 
              percentage = mult_perc
            )
            hub_multi_df_all[[algo]][[nam]][[cell]] <- bind_rows(hub_multi_df_all[[algo]][[nam]][[cell]], curr_mult_df)
          }
          
        }
        
      }
    }
  }
  
  
  # hub_multi_gg_all: stacked bar plot showing the multiplicity distribition per group and centrality.
  
  hub_multi_gg_all <- list()
  for(algo in names(hub_multi_df_all)){
    for(nam in names(hub_multi_df_all[[algo]])){
      for(cell in names(hub_multi_df_all[[algo]][[nam]])){
        bar_colors <- c("#599ad3", "#9e66ab")
        hub_multi_df_all[[algo]][[nam]][[cell]]$multiplicity <- factor(hub_multi_df_all[[algo]][[nam]][[cell]]$multiplicity,
                                                                       levels = c("2", "1")) 
        hub_multi_gg_all[[algo]][[nam]][[cell]] <- hub_multi_df_all[[algo]][[nam]][[cell]] %>% 
          ggplot(aes(group, percentage, 
                     fill = multiplicity)) +
          geom_bar(stat = "identity", colour = "black") +
          scale_x_discrete(name = "") +
          scale_y_continuous(name = "Percentage", expand = c(0,0)) +
          scale_fill_manual(name = "multiplicity", values = bar_colors) +
          facet_grid(. ~ centrality) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                           vjust = 0.5, lineheight = 1,
                                           size = 8),
                axis.title.y = element_text(size = 10),
                legend.title = element_text(size = 10),
                legend.text = element_text(size = 8))
        
      }
    }
  }
  
  
  # Save the hub_multiplicity_gg plot as a pdf file and the hub_multi_df as csv
  for(algo in names(hub_multi_gg_all)){
    for(nam in names(hub_multi_gg_all[[algo]])){
      for(cell in names(hub_multi_gg_all[[algo]][[nam]])){
        ggsave(
          paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Network_analysis/", algo, "/", nam, "/", cell, "/hub_multiplicity_", comp, ".pdf"), 
          plot = hub_multi_gg_all[[algo]][[nam]][[cell]], 
          device = "pdf", 
          width = 7, 
          height = 3
        )
      }
    }
  }
  
  
  
  ########################################################
  ## Define group-specific and shared hubs.
  #Based on the previous stacked bar plot, we can see that, indeed, there are hubs that are specific for each group, 
  #while others are ubiquitous. Let us separate them for downstream analysis:
  
  # Separate ubiquitious (multiplicity == 3+; top_centr_ub) from group-specific hubs (multiplicity = 1; top_centr_sp). 
  # Both top_centr_ub and top_centr_sp are lists of lists (outer list = CM, inner list = groups).
  
  top_all_all <- list()
  top_centr_ub_all <- list()
  top_centr_sp_all <- list()
  for(algo in names(top_centr_all_tmp)){
    for(nam in names(top_centr_all_tmp[[algo]])){
      for(cell in names(top_centr_all_tmp[[algo]][[nam]])){
        for (centr in names(top_centr_all_tmp[[algo]][[nam]][[cell]])) {
          
          top_all_all[[algo]][[nam]][[cell]] <- table(unlist(top_centr_all_tmp[[algo]][[nam]][[cell]][[centr]]))
          top_centr_ub_all[[algo]][[nam]][[cell]][[centr]] <- map(top_centr_all_tmp[[algo]][[nam]][[cell]][[centr]], ~ .[top_all_all[[algo]][[nam]][[cell]][.] >= 2])
          top_centr_sp_all[[algo]][[nam]][[cell]][[centr]] <- map(top_centr_all_tmp[[algo]][[nam]][[cell]][[centr]], ~ .[top_all_all[[algo]][[nam]][[cell]][.] == 1])
          
        }
      }
    }
  }
  
  saveRDS(top_centr_ub_all, paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Network_analysis/networks_hub_ubiquitous_", comp, ".rds"))
  saveRDS(top_centr_sp_all, paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Network_analysis/networks_hub_specific_", comp, ".rds"))
  
  
}













