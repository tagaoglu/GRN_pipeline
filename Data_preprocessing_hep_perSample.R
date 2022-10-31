# Set *.libPaths*  
.libPaths(c("~/p728/RSTUDIO/R/library/4.1/", .libPaths()))


#Set working directory and seed:
setwd("~/p728/RSTUDIO/analysis/tagaoglu")
set.seed(123)


library(Seurat)


##Load Data
#load("old_data/data_v20220928/merged.seurat.RData")
integrated.seurat <- readRDS("data/merged.seurat.rds")


# Get cell identity classes
Idents(object = integrated.seurat)
levels(x = integrated.seurat)
#OLD:
#[1] "Hepatocytes"    "T/NK-cells"     "Macrophages"    "Endothelial"    "Fibroblasts"    "Plasma_cells"   "B-cells"       
#[8] "Mast_cells"     "pDCs"           "Cholangiocytes"

#UP TO DATE:
#[1] "T/NK-cells"     "Hepatocytes"    "Macrophages"    "Endothelial"    "Fibroblasts"    "Kupffer_cells"  "B-cells"   "Cholangiocytes"  "erythrocytes"  
#[10] "pDCs"           "Plasma_cells"   "Mast_cells"


table(integrated.seurat@meta.data$cluster_main_ctype)
#OLD:
#Hepatocytes     T/NK-cells    Macrophages    Endothelial    Fibroblasts   Plasma_cells        B-cells     Mast_cells 
#19279          15025          11079           6021           2139            451            366            189 
#pDCs Cholangiocytes 
#115             87 

#UP TO DATE:
#T/NK-cells    Hepatocytes    Macrophages    Endothelial    Fibroblasts  Kupffer_cells        B-cells Cholangiocytes   erythrocytes           pDCs 
#     76291          45252          42096          14513           4488           3062           2790           1892           1118            947 
#Plasma_cells     Mast_cells 
#         652            239 

#To use later on for sanity-check new objects
table(integrated.seurat@meta.data$sample)
#OLD:
# P67    P108  P123s2  P123s8    P191    P193    P199    P207     P207Gel   P210    P215   P220    P269    P275    P282C    P295    P297 
#13823    3713    1771    1310    2627    1525     975    2240    1391      972    2298    3922    919     5048    5104     3934    3179 

#UP TO DATE:
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
# #P140  P142  P307 : excluded in the upstream analysis !!!!!


############################################################################
### Data extraction per group
############################################################################

# Create an empty list in order to load all in one list ,
samples <- c("P123s2", "P123s8", "P193", "P199", "P108", "P191", "P207", "P207Gel", "P215", "P220", "P275", "P282C", "P295")
expData_tmp <- sapply(samples,function(x) NULL)

############################################################################
### 1.Subset data for each group: viral HCC, non-viral HCC

#viral&nonviral HCC samples
for(i in samples){
  expData_tmp[[i]] <- subset(x = integrated.seurat, subset = (sample == i))
}


#sanity-check new object
table(expData_tmp[["P123s2"]]@meta.data$sample)
table(expData_tmp[["P207"]]@meta.data$sample)

############################################################################
### 2. Extract raw expression data for each cell cluster and subcluster

#OLD:
#levels <- c("Hepatocytes", "Endothelial", "Fibroblasts", "T/NK-cells", "B-cells", "Macrophages", "Mast_cells", "pDCs")
#"Undefined_1", "Undefined_2", "Undefined_3" : not included in the list above

#UP TO DATE:
levels <- c("Hepatocytes")

expData <- list() # to put expression data in a list for easier manipulation
for(g in names(expData_tmp)){
  for(l in levels){
    expData[[g]][[l]] <- as.matrix(GetAssayData(object = expData_tmp[[g]], assay = 'RNA', slot = "counts")[, WhichCells(expData_tmp[[g]], ident = l)])
  }
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

thr <- 2

##create an input file if it does not exist
nam <- paste("expData_mean_", thr, sep = "")
inputs_dir= paste("~/p728/RSTUDIO/analysis/tagaoglu/data/GRN_inputs/hep_perSample/", nam, "/", sep = "")
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
















