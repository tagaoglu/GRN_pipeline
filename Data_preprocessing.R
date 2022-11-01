# Set *.libPaths*  
.libPaths(c("~/p728/RSTUDIO/R/library/4.1/", .libPaths()))


#Set working directory and seed:
setwd("~/p728/RSTUDIO/analysis/tagaoglu")
set.seed(123)


library(Seurat)


##Load Data
integrated.seurat <- readRDS("data/merged.seurat.rds")


# Get cell identity classes
Idents(object = integrated.seurat)
levels(x = integrated.seurat)
#[1] "T/NK-cells"     "Hepatocytes"    "Macrophages"    "Endothelial"    "Fibroblasts"    "Kupffer_cells"  "B-cells"   "Cholangiocytes"  "erythrocytes"  
#[10] "pDCs"           "Plasma_cells"   "Mast_cells"



#To use later on for sanity-check new objects
table(integrated.seurat@meta.data$sample)
# P67    P108   P123s2  P123s8   P191   P193    P199    P207    P207Gel     P210   P215    P220    P269    P275    P282C    P295    P297
# 14312  4310   1779    1364     2646   1816    1004    2356    1499        1003   2907    4017    938     5170    5759     4140    4381
# Andrews_2022_Donor_4         Guilliams_2022_H02         Guilliams_2022_H06            Guilliams_2022_H16            Payen_2021_158   
#                 7892                      17369                      11633                         11782                     21655    
# Ramachandran_2019_Healthy1 Ramachandran_2019_Healthy2 Ramachandran_2019_Healthy3      Ramachandran_2019_Healthy4     Ramachandran_2019_Healthy5 
#                       2837                       9372                       8982                            8017                           5342
# Wang_2022_Con_1     Wang_2022_Con_2       Wang_2022_Con_3       Wang_2022_Con_4
#           10120                5735                  6575                  6628


load("data/full_info_samples.RData")

#Check sample names for viral HCC, non-viral HCC and normal liver
#consider HCV, steatosis and normal livers and omit those mixed and idiopathic

#viral HCC
full_info$Sample[which(full_info$HCV_vs_steatosis == "HCV")]
# P123s2  P123s8  P193  P199

#non-viral HCC
full_info$Sample[which(full_info$HCV_vs_steatosis == "Steatosis")]
# P108    P140    P142    P191    P207    P207Gel
# P215    P220    P275    P282C   P295    P307

#normal liver
full_info$Sample[which(full_info$Sample_type == "Normal_liver")]
# P210 P269


# IMPORTANT NOTE:
# ###viral HCC 
# # P123s2  P123s8  P193  P199 
#
# ###non-viral HCC
# # P108  P191  P207 P207Gel  P215    P220    P275   P282C    P295
# 
# ###normal liver
# #P210 P269
# 
# # Side note  
# #Others
# #P67  P297 =  HCV+steatosis !!!!!!
# #Charlotte: "only consider HCV, steatosis and normal livers and omit those mixed and idiopathic." 
# 
# #P140  P142  P307 : excluded !!!!!


############################################################################
### Data extraction per group
############################################################################

# Create an empty list in order to load all in one list , 
groups <- c("tumor", "nonviral", "viral", "normal")


# Since viral HCC include no cell for "Kupffer_cells", "erythrocytes", "Cholangiocytes", and nonviral HCC include only 1 cell for "erythrocytes", "Cholangiocytes" and no cell for "Kupffer_cells"
# I will exclude these cell types for the following steps 

levels <- c("Hepatocytes", "Endothelial", "Fibroblasts", "T/NK-cells", "B-cells", "Macrophages", "Mast_cells", "pDCs", "Plasma_cells")



############################################################################
### 1.Subset data for each group: viral HCC, non-viral HCC, normal liver

expData_tmp <- sapply(groups,function(x) NULL)

#viral HCC 
expData_tmp[["viral"]] <- subset(x = integrated.seurat, subset = (sample == "P123s2" | sample == "P123s8"| sample == "P193"| sample == "P199"))
#sanity-check new object
table(expData_tmp[["viral"]]@meta.data$sample)

#non-viral HCC
expData_tmp[["nonviral"]] <- subset(x = integrated.seurat, subset = (sample == "P108" | sample == "P191"| sample == "P207"| sample == "P207Gel"| sample == "P215"| sample == "P220"| sample == "P275"| sample == "P282C"| sample == "P295"))
#sanity-check new object
table(expData_tmp[["nonviral"]]@meta.data$sample)


#tumor: viral HCC + nonviral HCC
expData_tmp[["tumor"]] <- subset(x = integrated.seurat, subset = (sample == "P123s2" | sample == "P123s8"| sample == "P193"| sample == "P199" | sample == "P108" | sample == "P191"| sample == "P207"| sample == "P207Gel"| sample == "P215"| sample == "P220"| sample == "P275"| sample == "P282C"| sample == "P295"))
#sanity-check new object
table(expData_tmp[["tumor"]]@meta.data$sample)


#normal liver
expData_tmp[["normal"]] <- subset(x = integrated.seurat, subset = (sample == "P210" | sample == "P269" | sample == "Andrews_2022_Donor_4"  | sample == "Guilliams_2022_H02"  |
                                                                     sample == "Guilliams_2022_H06"  | sample == "Guilliams_2022_H16"  | sample == "Payen_2021_158"  |
                                                                     sample == "Ramachandran_2019_Healthy1"  | sample == "Ramachandran_2019_Healthy2"  | sample == "Ramachandran_2019_Healthy3"  |
                                                                     sample == "Ramachandran_2019_Healthy4"  | sample == "Ramachandran_2019_Healthy5"  | sample == "Wang_2022_Con_1"  |
                                                                     sample == "Wang_2022_Con_2"  | sample == "Wang_2022_Con_3"  | sample == "Wang_2022_Con_4"))
#sanity-check new object
table(expData_tmp[["normal"]]@meta.data$sample)

############################################################################
### 2. Extract raw expression data for each cell cluster and subcluster

expData <- list() # to put expression data in a list for easier manipulation
for(g in groups){
  for(l in levels){
    expData[[g]][[l]] <- as.matrix(GetAssayData(object = expData_tmp[[g]], assay = 'RNA', slot = "counts")[, WhichCells(expData_tmp[[g]], ident = l)])
  }
  #Rename list items, B-cells and T/NK-cells
  names(expData[[g]])[4] <- "T_NK_cells"
  names(expData[[g]])[5] <- "B_cells"
}


############################################################################
### Filtering out genes and cells
############################################################################

############################################################################
### 1. Filter out the mitochondrial (MT- ) genes

#The list of top central genes is dominated by mitochondrial (MT-) genes.
#Therefore, these genes will not be retained for subsequent analysis steps.

expData_noMT <- list()
for(group in names(expData)){
  for(cell in names(expData[[group]])){
    mt_genes <- grepl("^MT-",  rownames(expData[[group]][[cell]]))
    mt_genes_list <- rownames(expData[[group]][[cell]])[mt_genes]
    expData_noMT[[group]][[cell]] <- as.matrix(expData[[group]][[cell]][!(rownames(expData[[group]][[cell]]) %in% mt_genes_list), ])
  }
}


# #in case I will need to exclude ribosomal genes
# rps_genes <- grepl("^RPS",  rownames(expData[[group]][[cell]]))
# rpl_genes <- grepl("^RPL",  rownames(expData[[group]][[cell]]))
# rps_genes_list <- rownames(expData[[group]][[cell]])[rps_genes]
# rpl_genes_list <- rownames(expData[[group]][[cell]])[rpl_genes]


############################################################################
### 2. Gene Filtering and Normalization

######################
# Gene filtering based on mean-based filtering

source("functions/geneFiltering_mean.R")
# Mean gene filtering : Filter-out low-abundance genes (previous to GRN-inference) based on mean-based filtering
# low-abundance genes are defined as those with an average count below a filter threshold of, by default, 0.25
# You can adjust threshold according to your dataset
# Default value: aveCounts_thr=0.25

thr <- 1

##create an input file if it does not exist
nam <- paste("expData_mean_", thr, sep = "")
inputs_dir= paste("~/p728/RSTUDIO/analysis/tagaoglu/data/GRN_inputs/", nam, "/", sep = "")
if(!dir.exists(inputs_dir)){dir.create(inputs_dir,  recursive = T)}


expData_filtered_log <- list()
for(group in names(expData_noMT)){
  for(cell in names(expData_noMT[[group]])){
    genesKept <- geneFiltering_mean(expData_noMT[[group]][[cell]], aveCounts_thr=thr)
    expMat_filtered <- expData_noMT[[group]][[cell]][genesKept, ]
    expMat_filtered_log <- log2(expMat_filtered+1)
    expData_filtered_log[[group]][[cell]] <- as.matrix(expMat_filtered_log)
    message("\nNumber of genes left after applying the filter for ", cell, ":")
    message("\t", length(genesKept), "\tgenes with average counts per gene >= ", thr)
    nam <- paste(group, "_", cell, "_expMat", sep = "")
    write.table(expData_filtered_log[[group]][[cell]], file=paste0(inputs_dir, nam, ".tsv"), quote=FALSE, sep='\t', col.names = TRUE) # required for PIDC
  }
}


# Save RDS for later use, required for GENIE3
saveRDS(expData_filtered_log, file = paste0(inputs_dir, "expData_filtered_log.Rds"))  







# ######################
# # Gene filter/selection
# #Before running GENIE3/GRNBoost, it is recommended to apply soft gene filter, to remove genes that are 
# #expressed either at very low levels or in too few cells. Here we apply a filtering based on 
# #the total number of counts of the gene, and the number of cells in which it is detected.
# #First filter: Filter by the total number of reads per gene. This filter is meant to remove genes that are most likely noise. By default it keeps only the genes with at least 6 UMI counts across all samples/cells (e.g. the total number the gene would have, if it was expressed with a value of 3 in 1% of the cells). Adjust this value (minCountsPerGene) according to the dataset (it will depend on the dataset units, e.g. UMI, TPMs…).
# #Second filter: Filter by the number of cells in which the gene is detected** (e.g. >0 UMI, or >1 log2(TPM)). By default (minSamples), genes that are detected in at least 1% of the cells are kept. This filtering is meant to remove genes whose reads come from one a few ‘noisy’ cells (genes that are only expressed in one, or very few cells, gain a lot of weight if they happen to coincide in a given cell). To avoid removing small (but potentially interesting) cell populations, we recommend to set a percentage lower than the smallest population of cells to be detected.
# 
# 
# source("~/p728/RSTUDIO/analysis/tagaoglu/functions/geneFiltering.R")
# #Gene filtering : Filter-out genes (previous to GRN-inference) based on the counts and number of cells in which they are detected.
# #minCountsPerGene : Minimum counts per gene required
# #minSamples : Minimum number of samples (cells) in which the gene should be detected
# 
# # You can adjust threshold according to your dataset
# # Default values: 
# #minCountsPerGene=3*.1*ncol(expMat)
# #minSamples=ncol(expMat)*.01
# 
# minC <- 3*1
# minS <- .01
# 
# ##create an input file if it does not exist
# nam <- paste("expData_minC_", minC, "_minS_", minS, sep = "")
# inputs_dir= paste("~/p728/RSTUDIO/analysis/tagaoglu/data/GRN_inputs/", nam, "/", sep = "")
# if(!dir.exists(inputs_dir)){dir.create(inputs_dir,  recursive = T)}
# 
# 
# expData_filtered_log <- list()
# for(group in names(expData_noMT)){
#   for(cell in names(expData_noMT[[group]])){
#     genesKept <- geneFiltering(expData_noMT[[group]][[cell]], 
#                                minCountsPerGene = minC * ncol(expData_noMT[[group]][[cell]]),
#                                minSamples = ncol(expData_noMT[[group]][[cell]]) * minS)
#     expMat_filtered <- expData_noMT[[group]][[cell]][genesKept, ]
#     expMat_filtered_log <- log2(expMat_filtered+1)
#     expData_filtered_log[[group]][[cell]] <- expMat_filtered_log
#     nam <- paste(group, "_", cell, "_expMat", sep = "")
#     write.table(expData_filtered_log[[group]][[cell]], file=paste0(inputs_dir, nam, ".tsv"), quote=FALSE, sep='\t', col.names = TRUE) # required for PIDC
#   }
# }
# 
# # Save RDS for later use, required for GENIE3
# saveRDS(expData_filtered_log, file = paste0(inputs_dir, "expData_filtered_log.Rds"))  
# 










############################################################################
#BACK-UP

############################################################################
# Take just 25% of cells randomly for each 
# In order to see how reliable our results are, and to see if the results of GRN algorithms depends on the number of cells

thr <- 1

# Load data
nam <- paste("expData_mean_", thr, sep = "")
inputs_dir= paste("~/p728/RSTUDIO/analysis/tagaoglu/data/GRN_inputs/", nam, "/", sep = "")
expData_filtered_log <- readRDS(paste0(inputs_dir, "expData_filtered_log.Rds"))


##Create an input file if it does not exist
nam <- paste("expData_mean_", thr, "_RANDOM", sep = "")
inputs_dir= paste("~/p728/RSTUDIO/analysis/tagaoglu/data/GRN_inputs/", nam, "/", sep = "")
if(!dir.exists(inputs_dir)){dir.create(inputs_dir,  recursive = T)}


# #or

# minC <- 3*.5
# minS <- .01
# 
# # Load data
# nam <- paste("expData_minC_", minC, "_minS_", minS, sep = "")
# inputs_dir= paste("~/p728/RSTUDIO/analysis/tagaoglu/data/GRN_inputs/", nam, "/", sep = "")
# expData_filtered_log <- readRDS(paste0(inputs_dir, "expData_filtered_log.Rds"))
# 
# 
# ##Create an input file if it does not exist
# nam <- paste("expData_minC_", minC, "_minS_", minS, "_RANDOM", sep = "")
# inputs_dir= paste("~/p728/RSTUDIO/analysis/tagaoglu/data/GRN_inputs/", nam, "/", sep = "")
# if(!dir.exists(inputs_dir)){dir.create(inputs_dir,  recursive = T)}

expData_filtered_log_RANDOM <- list()
for(group in names(expData_filtered_log)){
  for(cell in names(expData_filtered_log[[group]])){
    tmp <- expData_filtered_log[[group]][[cell]]
    randomCells <- tmp[ ,sample(ncol(tmp), round(0.25*ncol(tmp)))] #Randomly sampling 25% of columns
    expData_filtered_log_RANDOM[[group]][[cell]] <- as.matrix(randomCells)
  }
}

expData_filtered_log_RANDOM[["viral"]][["Cholangiocytes"]] <- expData_filtered_log[["viral"]][["Cholangiocytes"]]
expData_filtered_log_RANDOM[["nonviral"]][["Cholangiocytes"]] <- expData_filtered_log[["nonviral"]][["Cholangiocytes"]]


for(group in names(expData_filtered_log_RANDOM)){
  for(cell in names(expData_filtered_log_RANDOM[[group]])){
    nam <- paste(group, "_", cell, "_expMat", sep = "")
    write.table(expData_filtered_log_RANDOM[[group]][[cell]], file=paste0(inputs_dir, nam, ".tsv"), quote=FALSE, sep='\t', col.names = TRUE) # required for PIDC
  }
}


# Save RDS for later use, required for GENIE3
saveRDS(expData_filtered_log_RANDOM, file = paste0(inputs_dir, "expData_filtered_log.Rds"))  


