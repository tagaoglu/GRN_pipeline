#RUN IN CLUSTER

#Set working directory and seed:
setwd("/data/projects/p728_scRNA_HCC/RSTUDIO/analysis/tagaoglu")
set.seed(123)


library(GENIE3)

#########################################################################

# Set accordingly
group = c("viral") #group = c("normal", "viral", "nonviral", "tumor")
sample = c("P123s2", "P123s8", "P193", "P199")
cell = c("Hepatocytes") #cell = c("Hepatocytes", "Endothelial", "Fibroblasts", "T_NK_cells", "B_cells", "Macrophages", "Mast_cells", "pDCs")
algo = c("GENIE3") #algo = c("GENIE3", "PIDC")

#########################################################################

# Set accordingly

thr <- 2

# Load data
nam <- paste("expData_mean_", thr, sep = "")
# or 
#nam <- paste("expData_mean_", thr, "_RANDOM", sep = "")
inputs_dir= paste("data/GRN_inputs/hep_perSample/", nam, "/", sep = "")
expData_filtered_log <- readRDS(paste0(inputs_dir, "expData_filtered_log.Rds"))

# or 

# minC <- 3*.5
# minS <- .01
# 
# # Load data
# nam <- paste("expData_minC_", minC, "_minS_", minS, sep = "") 
# # or 
# #nam <- paste("expData_minC_", minC, "_minS_", minS, "_RANDOM", sep = "")
# inputs_dir= paste("data/GRN_inputs/", nam, "/", sep = "")
# expData_filtered_log <- readRDS(paste0(inputs_dir, "expData_filtered_log.Rds"))


################################################################################


##create networks file if it does not exist
networks_dir= paste("data/networks/", algo, "/hep_perSample/", nam, "/", sep = "") 
if(!dir.exists(networks_dir)){dir.create(networks_dir,  recursive = T)}


# Run GENIE3 with the default parameters
LinkList <- list()
weightMat <- list()

for (s in sample) {
  for (i in cell) {
    weightMat[[s]][[i]] <- GENIE3(expData_filtered_log[[s]][[i]])
    LinkList[[s]][[i]] <- getLinkList(weightMat[[s]][[i]])
  }
}

# save RDS for later use
saveRDS(weightMat, file=paste0(networks_dir, group, "_hep_weightMat.Rds"))
saveRDS(LinkList, file=paste0(networks_dir, group, "_hep_LinkList.Rds"))  




