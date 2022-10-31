#RUN IN CLUSTER

#Set working directory and seed:
setwd("/data/projects/p728_scRNA_HCC/RSTUDIO/analysis/tagaoglu")
set.seed(123)


library(GENIE3)

#########################################################################

# Set accordingly

group = c("normal") #group = c("normal", "viral", "nonviral", "tumor")
cell = c("Endothelial", "Fibroblasts", "T_NK_cells", "B_cells", "Macrophages", "Mast_cells", "pDCs", "Plasma_cells") #cell = c("Hepatocytes", "Endothelial", "Fibroblasts", "T_NK_cells", "B_cells", "Macrophages", "Mast_cells", "pDCs", "Plasma_cells", "Cholangiocytes")
algo = c("GENIE3") #algo = c("GENIE3", "PIDC")

#########################################################################

# Set accordingly

thr <- 2

# Load data
nam <- paste("expData_mean_", thr, sep = "")
# or 
#nam <- paste("expData_mean_", thr, "_RANDOM", sep = "")
inputs_dir= paste("data/GRN_inputs/", nam, "/", sep = "")
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
networks_dir= paste("data/networks/", algo, "/", nam, "/", sep = "") 
if(!dir.exists(networks_dir)){dir.create(networks_dir,  recursive = T)}


# Run GENIE3 with the default parameters
LinkList <- list()
weightMat <- list()

for (i in cell) {
  weightMat[[i]] <- GENIE3(expData_filtered_log[[group]][[i]])
  LinkList[[i]] <- getLinkList(weightMat[[i]])
  #write.table(weightMat[[i]], file=paste0(networks_dir, group, "_", i, "_weightMat.tsv"), quote=FALSE, sep='\t', col.names = TRUE)
  #write.table(LinkList[[i]], file=paste0(networks_dir, group, "_", i, "_LinkList.tsv"), quote=FALSE, sep='\t', col.names = TRUE)
}

# save RDS for later use
saveRDS(weightMat, file=paste0(networks_dir, group, "_weightMat.Rds"))
saveRDS(LinkList, file=paste0(networks_dir, group, "_LinkList.Rds"))  




