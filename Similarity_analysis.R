#Evaluation of similarity between networks, by computing the overlapping index, and Jaccard index and the WJS between these two networks

# Set *.libPaths*  
.libPaths(c("~/p728/RSTUDIO/R/library/4.1/", .libPaths()))


#Set working directory and seed:
setwd("~/p728/RSTUDIO/analysis/tagaoglu")
set.seed(123)

#Load functions
source("/home/tagaoglu/p728/RSTUDIO/analysis/tagaoglu/functions/Similarity_stats.R")
source("/home/tagaoglu/p728/RSTUDIO/analysis/tagaoglu/functions/save_png_pdf.R")


#Set accordingly:Directed = T  or  Directed = F
GENIE_edge_prediction <- T

#Set accordingly:Directed = T  or  Directed = F
PIDC_edge_prediction <- F


##################################################################################################
###Load data 
##################################################################################################

##################################################################################################
#GENIE3

GENIE_filt_type <-  c("all", "all_10k", "all_5k", "all_2k", "all_1k", "all_0.01_weightThr", "all_0.001_weightThr")

# Load the data
GENIE_net <- sapply(GENIE_filt_type,function(x) NULL)
for (thr in GENIE_filt_type){
  GENIE_net[[thr]] <- readRDS(paste("~/p728/RSTUDIO/analysis/tagaoglu/data/networks/GENIE3/GENIE3_", thr, ".Rds", sep = '' ))
}

##################################################################################################
#PIDC

PIDC_filt_type <-  c("all", "all_10k", "all_5k", "all_2k", "all_1k", "all_1_weightThr", "all_0.1_weightThr")

# Load the data
PIDC_net <- sapply(PIDC_filt_type,function(x) NULL)
for (thr in PIDC_filt_type){
  PIDC_net[[thr]] <- readRDS(paste("~/p728/RSTUDIO/analysis/tagaoglu/data/networks/PIDC/PIDC_", thr, ".Rds", sep = '' ))
}

##################################################################################################

#Sort cell names of list objects to compare easily
#GENIE
for(filt in names(GENIE_net)){
  for(nam in names(GENIE_net[[filt]])){
    for(group in names(GENIE_net[[filt]][[nam]])){
      GENIE_net[[filt]][[nam]][[group]] = GENIE_net[[filt]][[nam]][[group]][order(names(GENIE_net[[filt]][[nam]][[group]]))] # order cell names
    }
  }
}

#PIDC
for(filt in names(PIDC_net)){
  for(nam in names(PIDC_net[[filt]])){
    for(group in names(PIDC_net[[filt]][[nam]])){
      PIDC_net[[filt]][[nam]][[group]] = PIDC_net[[filt]][[nam]][[group]][order(names(PIDC_net[[filt]][[nam]][[group]]))] # order cell names
    }
  }
}




##################################################################################################
#Similarity Comparison
##################################################################################################

##################################
#All cells vs Random cells (25%)
##################################

#Calculate similarity index between networks inferred from the data with all cells and the data with sampling 25% of cells by the same algorithm
#Side note: I took just 25% of random cells and then try to construct network, to see whether I can get comparable answer or I get nothing. To see how reliable our results are.)

comp <- "AllvsRandom"

########
#GENIE3

#Rearrange lists to easily compare
nam_dir <- c("expData_mean_1", "expData_mean_2")
GENIE <- list()
for (filt in names(GENIE_net)){
  for (nam in nam_dir){
  GENIE[[filt]][[nam]] <- GENIE_net[[filt]][[nam]]
  }
}

nam_dir_RANDOM <- c("expData_mean_1_RANDOM", "expData_mean_2_RANDOM")
GENIE_RANDOM <- list()
for (filt in names(GENIE_net)){
  for (nam in nam_dir_RANDOM){
    GENIE_RANDOM[[filt]][[nam]] <- GENIE_net[[filt]][[nam]]
  }
}

#Rename to compare easily in the next steps
for (filt in names(GENIE_RANDOM)){
  names(GENIE_RANDOM[[filt]]) <- c("expData_mean_1", "expData_mean_2")
}

######
#PIDC

#Rearrange lists to easily compare
nam_dir <- c("expData_mean_1", "expData_mean_2")
PIDC <- list()
for (filt in names(PIDC_net)){
  for (nam in nam_dir){
    PIDC[[filt]][[nam]] <- PIDC_net[[filt]][[nam]]
  }
}

nam_dir_RANDOM <- c("expData_mean_1_RANDOM","expData_mean_2_RANDOM")
PIDC_RANDOM <- list()
for (filt in names(PIDC_net)){
  for (nam in nam_dir_RANDOM){
    PIDC_RANDOM[[filt]][[nam]] <- PIDC_net[[filt]][[nam]]
  }
}

#Rename to compare easily in the next steps
for (filt in names(PIDC_RANDOM)){
  names(PIDC_RANDOM[[filt]]) <- c("expData_mean_1", "expData_mean_2")
}



#####################################################################################
###Similarity Stat results
#GENIE3
#####################################################################################

## Load package
library(igraph)

GENIE_similarity <- list()
for (filt in names(GENIE_RANDOM)){
  for(nam in names(GENIE_RANDOM[[filt]])){
    for(group in names(GENIE_RANDOM[[filt]][[nam]])){
      res.df= setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("overlap.index", "jaccard.index", "cell_type", "group", "expData_thr", "links_thr"))
      for(i in 1:length(GENIE_RANDOM[[filt]][[nam]][[group]])){
        if(names(GENIE[[filt]][[nam]][[group]])[i]==names(GENIE_RANDOM[[filt]][[nam]][[group]])[i]){
          overlap.index= Intersection.index(GENIE[[filt]][[nam]][[group]][[i]], GENIE_RANDOM[[filt]][[nam]][[group]][[i]], 
                                            directed1=GENIE_edge_prediction, directed2=GENIE_edge_prediction)[1] 
          jaccard.index= Intersection.index(GENIE[[filt]][[nam]][[group]][[i]], GENIE_RANDOM[[filt]][[nam]][[group]][[i]], 
                                            directed1=GENIE_edge_prediction, directed2=GENIE_edge_prediction)[2] 
          #WJS= wjs(GENIE[[filt]][[nam]][[group]][[i]], GENIE_RANDOM[[filt]][[nam]][[group]][[i]], directed1=GENIE_edge_prediction, directed2=GENIE_edge_prediction)
          #jes= jaccard_edgeset_similarity(GENIE[[filt]][[nam]][[group]][[i]], GENIE_RANDOM[[filt]][[nam]][[group]][[i]], directed1=GENIE_edge_prediction, directed2=GENIE_edge_prediction)
          res= c(overlap.index, jaccard.index, names(GENIE[[filt]][[nam]][[group]])[i], group, nam, filt)
          res.df[i,]=res
        }
      }
      
      #Format numbers as percentages with one decimal place
      res.df$overlap.index <- as.numeric(res.df$overlap.index) * 100
      res.df$jaccard.index <- as.numeric(res.df$jaccard.index) * 100
      #Change the variable class to numeric
      #res.df$WJS <- as.numeric(res.df$WJS)
      
      GENIE_similarity[[filt]][[nam]][[group]] <- res.df
    }
  }
}


#Following steps for boxplot for all JI results;

#Combine a list of data frames into one data frame by row
GENIE_JI_merged <- data.frame()
for (filt in names(GENIE_similarity)){
  for(nam in names(GENIE_similarity[[filt]])){
    for(group in names(GENIE_similarity[[filt]][[nam]])){
      GENIE_JI_merged <- rbind(GENIE_JI_merged, GENIE_similarity[[filt]][[nam]][[group]])
    }
  }
}

GENIE_JI_merged$overlap.index <- as.numeric(GENIE_JI_merged$overlap.index) / 100
GENIE_JI_merged$jaccard.index <- as.numeric(GENIE_JI_merged$jaccard.index) / 100
GENIE_JI_merged$jaccard.index <- round(as.numeric(GENIE_JI_merged$jaccard.index),digits = 3)

#
str(GENIE_JI_merged)
#
quantile(GENIE_JI_merged$jaccard.index, probs=c(0, 0.25,0.5, 0.75, 1))



#####################################################################################
###Plotting results

##Create files if it does not exist
for (filt in names(GENIE_similarity)){
  for(nam in names(GENIE_similarity[[filt]])){
    plotdir_GENIE=paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Similarity/GENIE3/", filt, "/", nam)
    if(!dir.exists(plotdir_GENIE)){dir.create(plotdir_GENIE,  recursive = T)}
  }
}

plotdir_GENIE=paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Similarity/GENIE3/")

#####################################################################################
#Boxplot
#Keep in mind that data distribution is hidden behind each boxin box plot
#Since, the exact distribution of each group is hidden behind boxes, I added individual observations on top of boxes

# Libraries
library(ggplot2)
library(tidyverse)
library(viridis) 

# Simple boxplot
#boxplot(jaccard.index ~ links_thr, data = GENIE_JI_merged, xlab = "Links threshold", ylab = "Jaccard index", main = comp)

# Plot boxplot
#https://r-graph-gallery.com/89-box-and-scatter-plot-with-ggplot2.html
myplot <- GENIE_JI_merged %>%
  ggplot( aes(x=links_thr, y=jaccard.index, fill=links_thr)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) + 
  theme_linedraw() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle(comp) +
  xlab("Links threshold") + 
  ylab("Jaccard index") 


png(paste0(plotdir_GENIE, comp, ".png"), 600, 500)
print(myplot)
dev.off()

#ggsave(paste0(plotdir_GENIE, comp, "_2.png"), plot = myplot)



#####################################################################################
#Heat table
library(gt)
library(scales)
library(readr)
library(webshot)


# Visualize and save results
for (filt in names(GENIE_similarity)){
  for(nam in names(GENIE_similarity[[filt]])){
    for(group in names(GENIE_similarity[[filt]][[nam]])){
      res.df <- GENIE_similarity[[filt]][[nam]][[group]]
      #Make a heatmap table
      p <- res.df %>%
        gt() %>%
        data_color(columns = 1:2, 
                   colors = col_numeric(palette = c("red","yellow","green"),
                                        domain = c(0,100))) %>%
        fmt_percent(columns = 1:2, scale_values = FALSE, decimals = 1) %>%
        tab_header(title = md("all vs random"),
                   subtitle = md(paste0(filt, "_", nam, "_", group))
        )
      #%>%
      #data_color(columns = 3, 
      #           colors = col_numeric(palette = c("white","red"),
      #                                domain = c(0,1))) %>%
      #fmt_number(columns = 3, decimals = 3)
      gtsave(p, paste0(group, ".html"), path = paste0(plotdir_GENIE, filt, "/", nam)) 
      
    }
  }
}



#####################################################################################
###Similarity Stat results
#PIDC
#####################################################################################

## Load package
library(igraph)

PIDC_similarity <- list()
for (filt in names(PIDC_RANDOM)){
  for(nam in names(PIDC_RANDOM[[filt]])){
    for(group in names(PIDC_RANDOM[[filt]][[nam]])){
      res.df= setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("overlap.index", "jaccard.index", "cell_type", "group", "expData_thr", "links_thr"))
      for(i in 1:length(PIDC_RANDOM[[filt]][[nam]][[group]])){
        if(names(PIDC[[filt]][[nam]][[group]])[i]==names(PIDC_RANDOM[[filt]][[nam]][[group]])[i]){
          overlap.index= Intersection.index(PIDC[[filt]][[nam]][[group]][[i]], PIDC_RANDOM[[filt]][[nam]][[group]][[i]], 
                                            directed1=PIDC_edge_prediction, directed2=PIDC_edge_prediction)[1] 
          jaccard.index= Intersection.index(PIDC[[filt]][[nam]][[group]][[i]], PIDC_RANDOM[[filt]][[nam]][[group]][[i]], 
                                            directed1=PIDC_edge_prediction, directed2=PIDC_edge_prediction)[2] 
          #WJS= wjs(PIDC[[filt]][[nam]][[group]][[i]], PIDC_RANDOM[[filt]][[nam]][[group]][[i]], directed1=PIDC_edge_prediction, directed2=PIDC_edge_prediction)
          #jes= jaccard_edgeset_similarity(PIDC[[filt]][[nam]][[group]][[i]], PIDC_RANDOM[[filt]][[nam]][[group]][[i]], directed1=PIDC_edge_prediction, directed2=PIDC_edge_prediction)
          res= c(overlap.index, jaccard.index, names(PIDC[[filt]][[nam]][[group]])[i], group, nam, filt)
          res.df[i,]=res
        }
      }
      
      #Format numbers as percentages with one decimal place
      res.df$overlap.index <- as.numeric(res.df$overlap.index) * 100
      res.df$jaccard.index <- as.numeric(res.df$jaccard.index) * 100
      #Change the variable class to numeric
      #res.df$WJS <- as.numeric(res.df$WJS)
      
      PIDC_similarity[[filt]][[nam]][[group]] <- res.df
    }
  }
}


#Following steps for boxplot for all JI results;

#Combine a list of data frames into one data frame by row
PIDC_JI_merged <- data.frame()
for (filt in names(PIDC_similarity)){
  for(nam in names(PIDC_similarity[[filt]])){
    for(group in names(PIDC_similarity[[filt]][[nam]])){
      PIDC_JI_merged <- rbind(PIDC_JI_merged, PIDC_similarity[[filt]][[nam]][[group]])
    }
  }
}

PIDC_JI_merged$overlap.index <- as.numeric(PIDC_JI_merged$overlap.index) / 100
PIDC_JI_merged$jaccard.index <- as.numeric(PIDC_JI_merged$jaccard.index) / 100
PIDC_JI_merged$jaccard.index <- round(as.numeric(PIDC_JI_merged$jaccard.index),digits = 3)

#
str(PIDC_JI_merged)
#
quantile(PIDC_JI_merged$jaccard.index, probs=c(0, 0.25,0.5, 0.75, 1))



#####################################################################################
###Plotting results

##Create files if it does not exist
for (filt in names(PIDC_similarity)){
  for(nam in names(PIDC_similarity[[filt]])){
    plotdir_PIDC=paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Similarity/PIDC/", filt, "/", nam)
    if(!dir.exists(plotdir_PIDC)){dir.create(plotdir_PIDC,  recursive = T)}
  }
}

plotdir_PIDC=paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Similarity/PIDC/")

#####################################################################################
#Boxplot
#Keep in mind that data distribution is hidden behind each boxin box plot
#Since, the exact distribution of each group is hidden behind boxes, I added individual observations on top of boxes

# Libraries
library(ggplot2)
library(tidyverse)
library(viridis) 

# Simple boxplot
#boxplot(jaccard.index ~ links_thr, data = PIDC_JI_merged, xlab = "Links threshold", ylab = "Jaccard index", main = comp)

# Plot boxplot
#https://r-graph-gallery.com/89-box-and-scatter-plot-with-ggplot2.html
myplot <- PIDC_JI_merged %>%
  ggplot( aes(x=links_thr, y=jaccard.index, fill=links_thr)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) + 
  theme_linedraw() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle(comp) +
  xlab("Links threshold") + 
  ylab("Jaccard index") 


png(paste0(plotdir_PIDC, comp, ".png"), 600, 500)
print(myplot)
dev.off()

#ggsave(paste0(plotdir_PIDC, comp, "_2.png"), plot = myplot)



#####################################################################################
#Heat table
library(gt)
library(scales)
library(readr)
library(webshot)


# Visualize and save results
for (filt in names(PIDC_similarity)){
  for(nam in names(PIDC_similarity[[filt]])){
    for(group in names(PIDC_similarity[[filt]][[nam]])){
      res.df <- PIDC_similarity[[filt]][[nam]][[group]]
      #Make a heatmap table
      p <- res.df %>%
        gt() %>%
        data_color(columns = 1:2, 
                   colors = col_numeric(palette = c("red","yellow","green"),
                                        domain = c(0,100))) %>%
        fmt_percent(columns = 1:2, scale_values = FALSE, decimals = 1) %>%
        tab_header(title = md("all vs random"),
                   subtitle = md(paste0(filt, "_", nam, "_", group))
        )
      #%>%
      #data_color(columns = 3, 
      #           colors = col_numeric(palette = c("white","red"),
      #                                domain = c(0,1))) %>%
      #fmt_number(columns = 3, decimals = 3)
      gtsave(p, paste0(group, ".html"), path = paste0(plotdir_PIDC, filt, "/", nam)) 
      
    }
  }
}




##################################################################################################
#Similarity Comparison
##################################################################################################

##################################
#GENIE3 vs PIDC
##################################

comp <- "GENIE_vs_PIDC"

#####################################################################################
###Similarity Stat results
#####################################################################################

## Load package
library(igraph)


#Rearrange lists to easily compare
nam_dir <- c("expData_mean_1", "expData_mean_2")

GENIE <- list()
for (filt in names(GENIE_net)){
  for (nam in nam_dir){
    GENIE[[filt]][[nam]] <- GENIE_net[[filt]][[nam]]
  }
}

PIDC <- list()
for (filt in names(PIDC_net)){
  for (nam in nam_dir){
    PIDC[[filt]][[nam]] <- PIDC_net[[filt]][[nam]]
  }
}

# ###
# #In case of need, this is to modify code above
# #Define to exclude perSample data
# groups <- c("normal", "viral", "nonviral", "tumor")
# GENIE[[filt]][[nam]] <- GENIE_net[[filt]][[nam]][names(GENIE_net[[filt]][[nam]]) %in% groups]
# PIDC[[filt]][[nam]] <- PIDC_net[[filt]][[nam]][names(PIDC_net[[filt]][[nam]]) %in% groups]
# ###


#Second, Extract network that will be compared, and calculated similarity index, 
#and also Rename list items to easily compare
net_filt_type <-  c("all", 
                    "10k_10k", "5k_5k", 
                    "2k_2k", "1k_1k",
                    "0.01w_1w", "0.001w_0.1w")

net1 <- sapply(net_filt_type,function(x) NULL)
net1[["all"]] <- GENIE[["all"]]
net1[["10k_10k"]] <- GENIE[["all_10k"]]
net1[["5k_5k"]] <- GENIE[["all_5k"]]
net1[["2k_2k"]] <- GENIE[["all_2k"]]
net1[["1k_1k"]] <- GENIE[["all_1k"]]
net1[["0.01w_1w"]] <- GENIE[["all_0.01_weightThr"]]
net1[["0.001w_0.1w"]] <- GENIE[["all_0.001_weightThr"]]


net2 <- sapply(net_filt_type,function(x) NULL)
net2[["all"]] <- PIDC[["all"]]
net2[["10k_10k"]] <- PIDC[["all_10k"]]
net2[["5k_5k"]] <- PIDC[["all_5k"]]
net2[["2k_2k"]] <- PIDC[["all_2k"]]
net2[["1k_1k"]] <- PIDC[["all_1k"]]
net2[["0.01w_1w"]] <- PIDC[["all_1_weightThr"]]
net2[["0.001w_0.1w"]] <- PIDC[["all_0.1_weightThr"]]


GENIE_vs_PIDC <- list()
for (filt in names(net1)){
  for(nam in names(net1[[filt]])){
    for(group in names(net1[[filt]][[nam]])){
      res.df= setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("overlap.index", "jaccard.index", "cell_type", "group", "expData_thr", "links_thr"))
      for(i in 1:length(net1[[filt]][[nam]][[group]])){
        if(names(net1[[filt]][[nam]][[group]])[i]==names(net2[[filt]][[nam]][[group]])[i]){
          overlap.index= Intersection.index(net1[[filt]][[nam]][[group]][[i]], net2[[filt]][[nam]][[group]][[i]], 
                                            directed1=GENIE_edge_prediction, directed2=PIDC_edge_prediction)[1] 
          jaccard.index= Intersection.index(net1[[filt]][[nam]][[group]][[i]], net2[[filt]][[nam]][[group]][[i]], 
                                            directed1=GENIE_edge_prediction, directed2=PIDC_edge_prediction)[2] 
          #WJS= wjs(net1[[filt]][[nam]][[group]][[i]], net2[[filt]][[nam]][[group]][[i]], directed1=GENIE_edge_prediction, directed2=PIDC_edge_prediction)
          #jes= jaccard_edgeset_similarity(net1[[filt]][[nam]][[group]][[i]], net2[[filt]][[nam]][[group]][[i]], directed1=GENIE_edge_prediction, directed2=PIDC_edge_prediction)
          res= c(overlap.index, jaccard.index, names(net1[[filt]][[nam]][[group]])[i], group, nam, filt)
          res.df[i,]=res
        }
      }
      
      #Format numbers as percentages with one decimal place
      res.df$overlap.index <- as.numeric(res.df$overlap.index) * 100
      res.df$jaccard.index <- as.numeric(res.df$jaccard.index) * 100
      #Change the variable class to numeric
      #res.df$WJS <- as.numeric(res.df$WJS)
      
      GENIE_vs_PIDC[[filt]][[nam]][[group]] <- res.df
    }
  }
}


#Following steps for boxplot for all JI results;

#Combine a list of data frames into one data frame by row
GENIE_vs_PIDC_merged <- data.frame()
for (filt in names(GENIE_vs_PIDC)){
  for(nam in names(GENIE_vs_PIDC[[filt]])){
    for(group in names(GENIE_vs_PIDC[[filt]][[nam]])){
      GENIE_vs_PIDC_merged <- rbind(GENIE_vs_PIDC_merged, GENIE_vs_PIDC[[filt]][[nam]][[group]])
    }
  }
}

GENIE_vs_PIDC_merged$overlap.index <- as.numeric(GENIE_vs_PIDC_merged$overlap.index) / 100
GENIE_vs_PIDC_merged$jaccard.index <- as.numeric(GENIE_vs_PIDC_merged$jaccard.index) / 100
GENIE_vs_PIDC_merged$jaccard.index <- round(as.numeric(GENIE_vs_PIDC_merged$jaccard.index),digits = 3)

#
str(GENIE_vs_PIDC_merged)
#
quantile(GENIE_vs_PIDC_merged$jaccard.index, probs=c(0, 0.25, 0.5, 0.75, 1))



#####################################################################################
#Boxplot
#Keep in mind that data distribution is hidden behind each boxin box plot
#Since, the exact distribution of each group is hidden behind boxes, I added individual observations on top of boxes

plotdir=paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Similarity/")

# Libraries
library(ggplot2)
library(tidyverse)
library(viridis) 

# Simple boxplot
#boxplot(jaccard.index ~ links_thr, data = GENIE_vs_PIDC_merged, xlab = "Links threshold", ylab = "Jaccard index", main = comp)

# Plot boxplot
#https://r-graph-gallery.com/89-box-and-scatter-plot-with-ggplot2.html
myplot <- GENIE_vs_PIDC_merged %>%
  ggplot( aes(x=links_thr, y=jaccard.index, fill=links_thr)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) + 
  theme_linedraw() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle(comp) +
  xlab("Links threshold") + 
  ylab("Jaccard index") 


png(paste0(plotdir, comp, ".png"), 600, 500)
print(myplot)
dev.off()

#ggsave(paste0(plotdir, comp, "_2.png"), plot = myplot)









#################
#Additional analysis
#BACKUP

#####################################################################################
#Just to check quickly if they were following the same trend for different expData thresholds
#which they did, that's why I do not save them in cluster

str(GENIE_JI_merged)

GENIE_JI_merged_2 <- subset(GENIE_JI_merged, expData_thr=="expData_mean_1")

GENIE_JI_merged_3 <- subset(GENIE_JI_merged, expData_thr=="expData_mean_2")


# Simple boxplot
boxplot(jaccard.index ~ links_thr, data = GENIE_JI_merged_2, xlab = "Links threshold", ylab = "Jaccard index", main = "All vs Random(25%) || expData_mean_1")
boxplot(jaccard.index ~ links_thr, data = GENIE_JI_merged_3, xlab = "Links threshold", ylab = "Jaccard index", main = "All vs Random(25%) || expData_mean_2")

# Plot boxplot
#https://r-graph-gallery.com/89-box-and-scatter-plot-with-ggplot2.html
GENIE_JI_merged_2 %>%
  ggplot( aes(x=links_thr, y=jaccard.index, fill=links_thr)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) + 
  theme_linedraw() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("All vs Random(25%) || expData_mean_1") +
  xlab("Links threshold") + 
  ylab("Jaccard index") 


GENIE_JI_merged_3 %>%
  ggplot( aes(x=links_thr, y=jaccard.index, fill=links_thr)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) + 
  theme_linedraw() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("All vs Random(25%) || expData_mean_2") +
  xlab("Links threshold") + 
  ylab("Jaccard index") 











