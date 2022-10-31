# Set *.libPaths*  
.libPaths(c("~/p728/RSTUDIO/R/library/4.1/", .libPaths()))


#Set working directory and seed:
setwd("~/p728/RSTUDIO/analysis/tagaoglu")
set.seed(123)


#Load functions
source("/home/tagaoglu/p728/RSTUDIO/analysis/tagaoglu/functions/Custom_Rcis.R")


##################################################################################################
###Set parameters first
##################################################################################################

#Set accordingly:Directed = T  or  Directed = F
GENIE_edge_prediction <- T

#Set accordingly:Directed = T  or  Directed = F
PIDC_edge_prediction <- F

##################################################################################################

#Determine positivecor, default 0, calculating the correlation between TF and target, if the correlation more than positivecor, we consider the edge is up regulated
positivecor <- 0

#Set distance
#https://resources.aertslab.org/cistarget/
#"hg19-tss-centered-5kb-7species.mc9nr.feather"
#TSS+/-5kb: 5kb around the TSS (total: 10kb)
distance <- "5kb"

##################################################################################################

##Either
#Determine k, in order to extract maximum number of most likely links to be analyzed
k <- 10000
tmp <- k/1000
net_path_to_file <- paste("~/p728/RSTUDIO/analysis/tagaoglu/data/networks/GENIE3/GENIE3_all_", tmp, "k.Rds", sep = '' )
act_path_to_file <- paste("~/p728/RSTUDIO/analysis/tagaoglu/data/networks/GENIE3/GENIE3_all_", tmp, "k_OnlyActivatedLinks.Rds", sep = '' )
Rcis_path_to_file <- paste("~/p728/RSTUDIO/analysis/tagaoglu/data/networks/GENIE3/GENIE3_all_", tmp, "k_OnlyActivatedLinks_Rcis_", distance, ".Rds", sep = '' )

# ##or
# #Determine global filter. Weight is greater than the threshold (weightThr)
# weightThr <- 0.005
# net_path_to_file <- paste("~/p728/RSTUDIO/analysis/tagaoglu/data/networks/GENIE3/GENIE3_all_", weightThr, "_weightThr.Rds", sep = '' )
# act_path_to_file <- paste("~/p728/RSTUDIO/analysis/tagaoglu/data/networks/GENIE3/GENIE3_all_", weightThr, "_weightThr_OnlyActivatedLinks.Rds", sep = '' )
# Rcis_path_to_file <- paste("~/p728/RSTUDIO/analysis/tagaoglu/data/networks/GENIE3/GENIE3_all_", weightThr, "_weightThr_OnlyActivatedLinks_Rcis_", distance, ".Rds", sep = '' )


##################################################################################################


##################################################################################################
###Load data (LinkList and ExpData) 
##################################################################################################

nam_dir <- c("expData_mean_1","expData_mean_2")

#perGroup
#Load expData_filtered_log 
expData_filtered_log <- sapply(nam_dir,function(x) NULL)
# Load the data a list
for (nam in nam_dir) {
  inputs_dir= paste("~/p728/RSTUDIO/analysis/tagaoglu/data/GRN_inputs/", nam, "/", sep = "")
  expData_filtered_log[[nam]] <- readRDS(paste0(inputs_dir, "expData_filtered_log.Rds"))
}

#hep_perSample
#Load expData_filtered_log 
expData_filtered_log_perSample <- sapply(nam_dir,function(x) NULL)
# Load the data a list
for (nam in nam_dir) {
  inputs_dir= paste("~/p728/RSTUDIO/analysis/tagaoglu/data/GRN_inputs/hep_perSample/", nam, "/", sep = "")
  expData_filtered_log_perSample[[nam]] <- readRDS(paste0(inputs_dir, "expData_filtered_log.Rds"))
  expData_filtered_log_perSample[[nam]][["normal"]] <- NULL
}


# Merge all in one
for (nam in names(expData_filtered_log_perSample)) {
  expData_filtered_log[[nam]] <- append(expData_filtered_log[[nam]], expData_filtered_log_perSample[[nam]])
}


# Load edgelist
GENIE_net_filtered <- readRDS(net_path_to_file)
#exclude random data
GENIE_net_filtered <- GENIE_net_filtered[names(GENIE_net_filtered) %in% nam_dir]


##################################################################################################
###Correlation
###GENIE
##################################################################################################
#GENIE3 can detect both positive and negative associations. In order to distinguish potential activation from repression, 
#we will split the targets into positive- and negative-correlated targets (i.e. Spearman correlation between the TF and the potential target).
#(This step can be run either before/after or simultaneously to GENIE3)


#! Taking long time (~30min)
#Step 1: extract the up- and down- regulate information 
GENIE_net_cor <- list()
for(nam in names(GENIE_net_filtered)){
  for(group in names(GENIE_net_filtered[[nam]])){
    for (cell in names(GENIE_net_filtered[[nam]][[group]])) {
      Exp <- expData_filtered_log[[nam]][[group]][[cell]]
      spearmanCor <- apply(GENIE_net_filtered[[nam]][[group]][[cell]], 1, function(y){
        cor(Exp[y[1], ], Exp[y[2], ], method = "spearman")
      })
      direct <- rep(0, length(spearmanCor))
      direct[spearmanCor > positivecor] <- 1
      direct[spearmanCor < (0-positivecor)] <- -1
      GENIE_net_cor[[nam]][[group]][[cell]] <- data.frame(GENIE_net_filtered[[nam]][[group]][[cell]], spearmanCor, direct)
    }
  }
}


#Step 2: extract active edges
GENIE_net_active <- list()
for(nam in names(GENIE_net_cor)){
  for(group in names(GENIE_net_cor[[nam]])){
    for(cell in names(GENIE_net_cor[[nam]][[group]])){
      x <- GENIE_net_cor[[nam]][[group]][[cell]]
      rownames(x) <- paste0(x[,1], "_", x[,2])
      GENIE_net_active[[nam]][[group]][[cell]] <- x[x[,"spearmanCor"]>0, ] 
    }
  }
}



# Save the data for later use
saveRDS(GENIE_net_active, file = act_path_to_file)

# Load the data
#GENIE_net_active <- readRDS(act_path_to_file)



##################################################################################################
#Rcis Target 
###GENIE
##################################################################################################

#RcisTarget has been customized to select only High Confidence links, that are supported by strong motif evidence

#https://resources.aertslab.org/cistarget/
#"hg19-tss-centered-5kb-7species.mc9nr.feather"
#TSS+/-5kb: 5kb around the TSS (total: 10kb)


#Taking long time! (~4h)
GENIE_net_wRcis <- list()
for(nam in names(GENIE_net_active)){
  for(group in names(GENIE_net_active[[nam]])){
  GENIE_net_wRcis[[nam]][[group]] <- Custom.Rcis(Networks = GENIE_net_active[[nam]][[group]],
                                          chosenDb = "hg19-tss-centered-5kb-7species.mc9nr.feather",
                                          MinGenesetSize = 0,
                                          directed = GENIE_edge_prediction)
  }
}


# Save the data for later use
saveRDS(GENIE_net_wRcis, file = Rcis_path_to_file)

# Load the data
#GENIE_net_wRcis <- readRDS(Rcis_path_to_file)





#########
###PIDC

##################################################################################################
#Correlation and Rcis Target :
#does not work for PIDC
#Because PIDC estimate undirected networks
##################################################################################################

