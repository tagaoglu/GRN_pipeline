# Set *.libPaths*  
.libPaths(c("~/p728/RSTUDIO/R/library/4.1/", .libPaths()))


#Set working directory and seed:
setwd("~/p728/RSTUDIO/analysis/tagaoglu")
set.seed(123)


library(Seurat)


##############################################################################

# levels <- c("Hepatocytes", "Endothelial", "Fibroblasts", "T_NK_cells", "B_cells", "Macrophages", "Mast_cells", "pDCs", "Plasma_cells")
# 
# ##Create files if it does not exist
# for(cell in levels){
#   plotdir=paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Expr_level/", cell)
#   if(!dir.exists(plotdir)){dir.create(plotdir,  recursive = T)}
# }

##############################################################################

##Load Data
#load("old_data/data_v20220928/merged.seurat.RData")
integrated.seurat <- readRDS("data/merged.seurat.rds")


# IMPORTANT NOTE:
# ###viral HCC 
# # P123s2  P123s8  P193  P199 
#
# ###non-viral HCC
# # P108  P191  P207 P207Gel  P215    P220    P275   P282C    P295
# 
# ###normal liver
# #P210 P269


levels(x = integrated.seurat)
#[1] "T/NK-cells"     "Hepatocytes"    "Macrophages"    "Endothelial"    "Fibroblasts"    "Kupffer_cells"  "B-cells"       
#[8] "Cholangiocytes" "erythrocytes"   "pDCs"           "Plasma_cells"   "Mast_cells"


xx <- subset(x = integrated.seurat, subset = (sample == "P123s2" | sample == "P123s8"| sample == "P193"| sample == "P199" | 
                                               sample == "P108" | sample == "P191"| sample == "P207"| sample == "P207Gel"| 
                                               sample == "P215"| sample == "P220"| sample == "P275"| sample == "P282C"| sample == "P295" |
                                               sample == "P210" | sample == "P269" | sample == "Andrews_2022_Donor_4"  | sample == "Guilliams_2022_H02"  |
                                               sample == "Guilliams_2022_H06"  | sample == "Guilliams_2022_H16"  | sample == "Payen_2021_158"  |
                                               sample == "Ramachandran_2019_Healthy1"  | sample == "Ramachandran_2019_Healthy2"  | sample == "Ramachandran_2019_Healthy3"  |
                                               sample == "Ramachandran_2019_Healthy4"  | sample == "Ramachandran_2019_Healthy5"  | sample == "Wang_2022_Con_1"  |
                                               sample == "Wang_2022_Con_2"  | sample == "Wang_2022_Con_3"  | sample == "Wang_2022_Con_4"))



###############################################################################

cell_interest <- c("Hepatocytes")
features <- c("MT2A", "HINT1", "SERF2", "ATF3", "JUND", "DUSP1", 
              "CPS1", "FOXP1", "DST", "LENG8", "MCRIP1", "MTRNR2L12", 
              "NFIA", "FKBP5", "ALDOB", "SULT2A1", "TOR1AIP2",
              "FOS", "JUN")

xxx <- subset(xx, idents = c(cell_interest)) #Create a Seurat object of just cell of interest

for(feature in features){
  v<- VlnPlot(xxx, features = feature, group.by = 'sample') # Violin plot - Visualize single cell expression distributions in each cluster
  print(v)
  #ggsave(paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Expr_level/", cell_interest, "/Expr_level_", feature, ".png"), plot = v, device = "png")
}


###############################################################################


cell_interest <- c("Endothelial")
features <- c("APP", 
              "HNRNPU", "KLF6", "BTG1", "TAGLN2",
              "ZFP36", "BHLHE40", "MEG3", "ETS1", 
              "CDKN1A", "NAMPT", "PPP1R15A", "GADD45B",
              "NFKBIZ", "DUSP1", "IER2", "MAFF", "NR4A2",
              "STC1", "HSPA1B", "DNAJB1", "DNAJB4", "CD74",
              "MACF1", "HSPG2",
              "FCN2", "CLEC1B", "SELENOP", "VIM") 

xxx <- subset(xx, idents = c(cell_interest)) #Create a Seurat object of just cell of interest

for(feature in features){
  v<- VlnPlot(xxx, features = feature, group.by = 'sample') # Violin plot - Visualize single cell expression distributions in each cluster
  print(v)
  #ggsave(paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Expr_level/", cell_interest, "/Expr_level_", feature, ".png"), plot = v, device = "png")
}


###############################################################################


cell_interest <- c("Macrophages")
features <- c("ITM2B", "RPL19", "PABPC1", "HLA-DRA", "PSAP", 
              "NEAT1", "JUN", "ATF3", "CEBPB", "CXCL10", "MAFF", "DUSP1", "IER2",
              "HNRNPU", "KLF6", "NR4A2", "KLF4", "MCL1",
              "S100A6", "FCGR3A", "HLA-B", "REL",
              "CTSS", "FOS", "NFKBIA") 

xxx <- subset(xx, idents = c(cell_interest)) #Create a Seurat object of just cell of interest

for(feature in features){
  v<- VlnPlot(xxx, features = feature, group.by = 'sample') # Violin plot - Visualize single cell expression distributions in each cluster
  print(v)
  #ggsave(paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Expr_level/", cell_interest, "/Expr_level_", feature, ".png"), plot = v, device = "png")
}


###############################################################################

cell_interest <- c("Fibroblasts")
features <- c("RPL15", "SPARC", "TMSB4X", 
              "BTG1", "ZFP36L1", "JUND", "DUSP1", "IER2", "PPP1R15A", "CCNL1", "MYLK",
              "RPL10", "EFEMP1", "NEAT1", "FSTL1", "COL4A1") 

xxx <- subset(xx, idents = c(cell_interest)) #Create a Seurat object of just cell of interest

for(feature in features){
  v<- VlnPlot(xxx, features = feature, group.by = 'sample') # Violin plot - Visualize single cell expression distributions in each cluster
  print(v)
  #ggsave(paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Expr_level/", cell_interest, "/Expr_level_", feature, ".png"), plot = v, device = "png")
}


###############################################################################

cell_interest <- c("T/NK-cells")
features <- c("RPL18A", "H3F3B", "RPL30", "RPL37A",
              "JUN", "IER2", "NR4A2", "YPEL5", "REL", "LMNA", "BTG1",
              "NR4A2", "YPEL5", "LMNA", "UBC",
              "JUND", "HSP90AB1", "HSP90AA1", "DNAJA1", "DDX5", "HSPH1", "H3F3B",
              "SFPQ", "FUS", "SRGN", "CNOT6L", "RPL37", "RPLP1") 

xxx <- subset(xx, idents = c(cell_interest)) #Create a Seurat object of just cell of interest

for(feature in features){
  v<- VlnPlot(xxx, features = feature, group.by = 'sample') # Violin plot - Visualize single cell expression distributions in each cluster
  print(v)
  #ggsave(paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Expr_level/", cell_interest, "/Expr_level_", feature, ".png"), plot = v, device = "png")
}


###############################################################################


cell_interest <- c("T/NK-cells")
features <- c("RPL18A", "H3F3B", "RPL30", "RPL37A",
              "JUN", "IER2", "NR4A2", "YPEL5", "REL", "LMNA", "BTG1",
              "NR4A2", "YPEL5", "LMNA", "UBC",
              "JUND", "HSP90AB1", "HSP90AA1", "DNAJA1", "DDX5", "HSPH1", "H3F3B",
              "SFPQ", "FUS", "SRGN", "CNOT6L", "RPL37", "RPLP1") 

xxx <- subset(xx, idents = c(cell_interest)) #Create a Seurat object of just cell of interest

for(feature in features){
  v<- VlnPlot(xxx, features = feature, group.by = 'sample') # Violin plot - Visualize single cell expression distributions in each cluster
  print(v)
  #ggsave(paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Expr_level/", cell_interest, "/Expr_level_", feature, ".png"), plot = v, device = "png")
}


###############################################################################

cell_interest <- c("B-cells")
features <- c("RPS10", "RPL30", "TMSB4X", "NR4A3", "RPL37A",
              "JUN", "DUSP1", "CREM", "NR4A2", "YPEL5", "CYCS", "CYTIP", "JUND",
              "NR3C1", "ETS1", "ANKRD11", "ATP6V0C", "REL", "CHD1", "PAX5", "BRD2",
              "MEF2C", "ZFP36", "JUNB", "NR4A1", "ZNF331", "H3F3B",
              "JUND", "DDX5", "HSP90AB1", "RPS19") 

xxx <- subset(xx, idents = c(cell_interest)) #Create a Seurat object of just cell of interest

for(feature in features){
  v<- VlnPlot(xxx, features = feature, group.by = 'sample') # Violin plot - Visualize single cell expression distributions in each cluster
  print(v)
  #ggsave(paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Expr_level/", cell_interest, "/Expr_level_", feature, ".png"), plot = v, device = "png")
}

###############################################################################


cell_interest <- c("Mast_cells")
features <- c("HPGDS", "RPS26", "RPS21", "VIM", "TMSB4X", "TPT1", "RPS29",
              "UBC", "KLF6", "CDC42", "CD9", "JUND", 
              "RUNX3", "XBP1", "PABPC4", "REL", "MITF", "CD82", 
              "IER2", "EMP3", "ZEB2", "CD69", "NFKB1", "JUNB",
              "RPL21", "TPSAB1", "RPL13A", "RPL7", "RPL27A", "RPS25",
              "TMSB10", "RPS23", "RPS21", "RPS3A", "RPS29") 

xxx <- subset(xx, idents = c(cell_interest)) #Create a Seurat object of just cell of interest

for(feature in features){
  v<- VlnPlot(xxx, features = feature, group.by = 'sample') # Violin plot - Visualize single cell expression distributions in each cluster
  print(v)
  #ggsave(paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Expr_level/", cell_interest, "/Expr_level_", feature, ".png"), plot = v, device = "png")
}

###############################################################################


cell_interest <- c("pDCs")
features <- c("RPS7","RPL21", "RPL11", "PABPC1", "HSP90AA1", "CLN8", "NABP1", "UBC",
              "FOSB", "RAB11FIP1", "UBC", "KLF6", "FOS",
              "IRF4", "DUSP1", "NR3C1", "SPIB", "COTL1", "TIPARP", "RUNX2", "ETV6", "REL",
              "IKZF1", "JUNB", "CREM", "GRASP", "BRD2", "HIVEP1",
              "RPL10", "RPL19") 

xxx <- subset(xx, idents = c(cell_interest)) #Create a Seurat object of just cell of interest

for(feature in features){
  v<- VlnPlot(xxx, features = feature, group.by = 'sample') # Violin plot - Visualize single cell expression distributions in each cluster
  print(v)
  #ggsave(paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Expr_level/", cell_interest, "/Expr_level_", feature, ".png"), plot = v, device = "png")
}

###############################################################################


cell_interest <- c("Plasma_cells")
features <- c("RPS15A", "RPS5", "RPL15", "RPL19", "PPIB", 
              "RGS1", "IGHGP", "SRRM2", "RPS11", "NEAT1", "RPL3",
              "YPEL5", "NEAT1", "IER2", "CHD2", "PDE4B", "AHNAK", "H3F3B", "ZFP36", "JUND", "FBXW7",
              "RPS23", "IGHG1", "RGS1", "FOS", "PTMA") 

xxx <- subset(xx, idents = c(cell_interest)) #Create a Seurat object of just cell of interest

for(feature in features){
  v<- VlnPlot(xxx, features = feature, group.by = 'sample') # Violin plot - Visualize single cell expression distributions in each cluster
  print(v)
  #ggsave(paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/plots/Expr_level/", cell_interest, "/Expr_level_", feature, ".png"), plot = v, device = "png")
}

###############################################################################

