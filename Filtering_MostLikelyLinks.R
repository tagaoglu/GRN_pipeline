# Set *.libPaths*  
.libPaths(c("~/p728/RSTUDIO/R/library/4.1/", .libPaths()))


#Set working directory and seed:
setwd("~/p728/RSTUDIO/analysis/tagaoglu")
set.seed(123)

##################################################################################################
###Set parameters first
##################################################################################################

#Set accordingly:Directed = T  or  Directed = F
GENIE_edge_prediction <- T

#Set accordingly:Directed = T  or  Directed = F
PIDC_edge_prediction <- F

#perGroup
nam_dir <- c("expData_mean_1", "expData_mean_2")
groups <- c("normal", "viral", "nonviral", "tumor")

#perSample
nam_dir_perSample <- c("expData_mean_1", "expData_mean_2")
groups_perSample_GENIE <- c("viral", "nonviral")
groups_perSample_PIDC <- c("P123s2", "P123s8", "P193", "P199", "P108", "P191", "P207", "P207Gel", "P215", "P220", "P275", "P282C", "P295")


##################################################################################################
###Load data (LinkList and ExpData) 
##################################################################################################

########################################
###GENIE3

#perGroup
# Load the data "edgelist" to a list
GENIE_net <- sapply(nam_dir,function(x) NULL)
for (nam in nam_dir) {
  for (group in groups) {
    GENIE_net[[nam]][[group]] <- readRDS(paste0("data/networks/GENIE3/", nam, "/", group, "_LinkList.Rds"))
    GENIE_net[[nam]][[group]] <- append(GENIE_net[[nam]][[group]], 
                                        readRDS(paste0("data/networks/GENIE3/", nam, "/", group, "_hep_LinkList.Rds")))
  }
}

#hep_perSample
GENIE_net_perSample <- sapply(nam_dir_perSample,function(x) NULL)
# Load the data "edgelist" to a list
for (nam in nam_dir_perSample) {
  for (group in groups_perSample_GENIE) {
    GENIE_net_perSample[[nam]][[group]] <- readRDS(paste0("data/networks/GENIE3/hep_perSample/", nam, "/", group, "_hep_LinkList.Rds"))
  }
}


for (nam in names(GENIE_net_perSample)) {
  GENIE_net_perSample[[nam]] <- append(GENIE_net_perSample[[nam]][["nonviral"]], GENIE_net_perSample[[nam]])
  GENIE_net_perSample[[nam]] <- append(GENIE_net_perSample[[nam]][["viral"]], GENIE_net_perSample[[nam]])
  GENIE_net_perSample[[nam]][["viral"]] <- NULL
  GENIE_net_perSample[[nam]][["nonviral"]] <- NULL
}


# Merge all in one
for (nam in names(GENIE_net_perSample)) {
  GENIE_net[[nam]] <- append(GENIE_net[[nam]], GENIE_net_perSample[[nam]])
}

path_to_file <- paste("~/p728/RSTUDIO/analysis/tagaoglu/data/networks/GENIE3/GENIE3_all.Rds", sep = '' )

# Save the data for later use
saveRDS(GENIE_net, file = path_to_file)

# Load the data
#GENIE_net <- readRDS(path_to_file)


########################################
###PIDC

#perGroup
# Load the data a list
PIDC_net <- sapply(nam_dir,function(x) NULL)
for (nam in nam_dir) {
  for (group in groups) {
    setwd(paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/networks/PIDC/", nam))
    sample_files <- list.files(pattern=glob2rx(paste0(group, "_", "*.tsv")))
    PIDC_net[[nam]][[group]] <- lapply(sample_files, read.table) #load
    names(PIDC_net[[nam]][[group]]) <- gsub(paste0(group, "_|_LinkList.tsv"), "", sample_files) #replaces the prefix and the suffix with the empty string in one go
  }
}

#hep_perSample
# Load the data a list
PIDC_net_perSample <- sapply(nam_dir_perSample,function(x) NULL)
for (nam in nam_dir_perSample) {
  for (group in groups_perSample_PIDC) {
    setwd(paste0("~/p728/RSTUDIO/analysis/tagaoglu/data/networks/PIDC/hep_perSample/", nam))
    sample_files <- list.files(pattern=glob2rx(paste0(group, "_", "*.tsv")))
    PIDC_net_perSample[[nam]][[group]] <- lapply(sample_files, read.table) #load
    names(PIDC_net_perSample[[nam]][[group]]) <- gsub(paste0(group, "_|_LinkList.tsv"), "", sample_files) #replaces the prefix and the suffix with the empty string in one go
  }
}


# Merge all in one
for (nam in names(PIDC_net_perSample)) {
  PIDC_net[[nam]] <- append(PIDC_net[[nam]], PIDC_net_perSample[[nam]])
}

path_to_file <- paste("~/p728/RSTUDIO/analysis/tagaoglu/data/networks/PIDC/PIDC_all.Rds", sep = '' )

# Save the data for later use
saveRDS(PIDC_net, file = path_to_file)

# Load the data
#PIDC_net <- readRDS(path_to_file)




##################################################################################################
#Filter the edges
##################################################################################################

########################################
###GENIE3

##################################################################################################
###Set parameters first

##############################################################
##Either
#Determine k, in order to extract maximum number of most likely links to be analyzed
k <- 10000
tmp <- k/1000
GENIE_path_to_file <- paste("~/p728/RSTUDIO/analysis/tagaoglu/data/networks/GENIE3/GENIE3_all_", tmp, "k.Rds", sep = '' )

GENIE_net_filtered <- list()
for(nam in names(GENIE_net)){
  for(group in names(GENIE_net[[nam]])){
    for(cell in names(GENIE_net[[nam]][[group]])){
      colnames(GENIE_net[[nam]][[group]][[cell]]) <- c('gene1', 'gene2', 'weight')
      GENIE_net[[nam]][[group]][[cell]] <- GENIE_net[[nam]][[group]][[cell]][, 1:3]
      #Either
      #Extract maximum number of most likely links to be analyzed, according to k
      if(dim(GENIE_net[[nam]][[group]][[cell]])[1] < k){
        thr <- dim(GENIE_net[[nam]][[group]][[cell]])[1]
        GENIE_net_filtered[[nam]][[group]][[cell]] <- GENIE_net[[nam]][[group]][[cell]][order(GENIE_net[[nam]][[group]][[cell]]$weight, decreasing = TRUE),][1:thr,]
      }else{
        thr <- k
        GENIE_net_filtered[[nam]][[group]][[cell]] <- GENIE_net[[nam]][[group]][[cell]][order(GENIE_net[[nam]][[group]][[cell]]$weight, decreasing = TRUE),][1:thr,]
      }
    }
  }
}


# Save the data for later use
saveRDS(GENIE_net_filtered, file = GENIE_path_to_file)

# Load the data
#GENIE_net_filtered <- readRDS(GENIE_path_to_file)


##############################################################
##or
#Determine global filter. Weight is greater than the threshold (weightThr)
weightThr <- 0.001
GENIE_path_to_file <- paste("~/p728/RSTUDIO/analysis/tagaoglu/data/networks/GENIE3/GENIE3_all_", weightThr, "_weightThr.Rds", sep = '' )

GENIE_net_filtered <- list()
for(nam in names(GENIE_net)){
  for(group in names(GENIE_net[[nam]])){
    for(cell in names(GENIE_net[[nam]][[group]])){
      colnames(GENIE_net[[nam]][[group]][[cell]]) <- c('gene1', 'gene2', 'weight')
      GENIE_net[[nam]][[group]][[cell]] <- GENIE_net[[nam]][[group]][[cell]][, 1:3]
      
      ##or
      ##Extract links according to global filter. Weight is greater than the threshold weightThr
      ## # Keep only genes with weight > threshold
      GENIE_net_filtered[[nam]][[group]][[cell]] <- GENIE_net[[nam]][[group]][[cell]][GENIE_net[[nam]][[group]][[cell]][,3]>weightThr,]
    }
  }
}


# Save the data for later use
saveRDS(GENIE_net_filtered, file = GENIE_path_to_file)

# Load the data
#GENIE_net_filtered <- readRDS(GENIE_path_to_file)




##################################################################################################
#Filter the edges
##################################################################################################

########################################
###PIDC

##################################################################################################
###Set parameters first

#############################################################
##Either
#Determine k, in order to extract maximum number of most likely links to be analyzed
k <- 10000
tmp <- k/1000
PIDC_path_to_file <- paste("~/p728/RSTUDIO/analysis/tagaoglu/data/networks/PIDC/PIDC_all_", tmp, "k.Rds", sep = '' )


PIDC_net_filtered <- list()
for(nam in names(PIDC_net)){
  for(group in names(PIDC_net[[nam]])){
    for(cell in names(PIDC_net[[nam]][[group]])){
      colnames(PIDC_net[[nam]][[group]][[cell]]) <- c('gene1', 'gene2', 'weight')
      PIDC_net[[nam]][[group]][[cell]] <- PIDC_net[[nam]][[group]][[cell]][, 1:3]
      
      #Either
      #Extract maximum number of most likely links to be analyzed, according to k
      if(dim(PIDC_net[[nam]][[group]][[cell]])[1] < k){
        thr <- dim(PIDC_net[[nam]][[group]][[cell]])[1]
        PIDC_net_filtered[[nam]][[group]][[cell]] <- PIDC_net[[nam]][[group]][[cell]][order(PIDC_net[[nam]][[group]][[cell]]$weight, decreasing = TRUE),][1:thr,]
      }else{
        thr <- k
        PIDC_net_filtered[[nam]][[group]][[cell]] <- PIDC_net[[nam]][[group]][[cell]][order(PIDC_net[[nam]][[group]][[cell]]$weight, decreasing = TRUE),][1:thr,]
      }
    }
  }
}

# Save the data for later use
saveRDS(PIDC_net_filtered, file = PIDC_path_to_file)

# Load the data
#PIDC_net_filtered <- readRDS(PIDC_path_to_file)



#############################################################
##or
#Determine global filter. Weight is greater than the threshold (weightThr)
weightThr <- 1 
PIDC_path_to_file <- paste("~/p728/RSTUDIO/analysis/tagaoglu/data/networks/PIDC/PIDC_all_", weightThr, "_weightThr.Rds", sep = '' )

PIDC_net_filtered <- list()
for(nam in names(PIDC_net)){
  for(group in names(PIDC_net[[nam]])){
    for(cell in names(PIDC_net[[nam]][[group]])){
      colnames(PIDC_net[[nam]][[group]][[cell]]) <- c('gene1', 'gene2', 'weight')
      PIDC_net[[nam]][[group]][[cell]] <- PIDC_net[[nam]][[group]][[cell]][, 1:3]
      
      #or
      #Extract links according to global filter. Weight is greater than the threshold weightThr
      # # Keep only genes with weight > threshold
      PIDC_net_filtered[[nam]][[group]][[cell]] <- PIDC_net[[nam]][[group]][[cell]][PIDC_net[[nam]][[group]][[cell]][,3]>weightThr,]
    }
  }
}


# Save the data for later use
saveRDS(PIDC_net_filtered, file = PIDC_path_to_file)

# Load the data
#PIDC_net_filtered <- readRDS(PIDC_path_to_file)


