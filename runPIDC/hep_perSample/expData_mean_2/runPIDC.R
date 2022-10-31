#RUN IN CLUSTER

#Set working directory and seed:
setwd("/data/projects/p728_scRNA_HCC/RSTUDIO/analysis/tagaoglu")
set.seed(123)


#install JuliaCall, Package JuliaCall is an R interface to Julia
#install.packages("JuliaCall")

#to use JuliaCall you must have a working installation of Julia
library(JuliaCall)
julia_setup(installJulia = TRUE)

# Install packages, if you don't already have the package installed
#julia_install_package_if_needed("NetworkInference")
#or to include version info:
#julia_command('Pkg.add(PackageSpec(name = "NetworkInference", version = "0.1.0"))')

# Include packages
julia_library("NetworkInference")

#########################################################################

# Set accordingly

groups <- c("P123s2", "P123s8", "P193", "P199", "P108", "P191", "P207", "P207Gel", "P215", "P220", "P275", "P282C", "P295")
julia_assign("groups", groups)

algo = c("PIDC") #algo = c("GENIE3", "PIDC")
julia_assign("algo", algo)

#########################################################################

# Set accordingly

thr = 2 # <-----------------------------------------------------------------------------------------
julia_assign("thr", thr)


# Load data
nam = paste("expData_mean_", thr, sep = "") # <-----------------------------------------------------
# or 
#nam = paste("expData_mean_", thr, "_RANDOM", sep = "") # <-------------------------------------------
julia_assign("nam", nam)
inputs_dir= paste("data/GRN_inputs/hep_perSample/", nam, "/", sep = "")
julia_assign("inputs_dir", inputs_dir)

# or 

# minC = 3*.5 # <------------------------------------------------------------------------------------
# minS = .01  # <------------------------------------------------------------------------------------
# julia_assign("minC", minC)
# julia_assign("minS", minS)
# 
#  
# # Load data
# nam = paste("expData_minC_", minC, "_minS_", minS, sep = "") # <------------------------------------
# # or 
# nam = paste("expData_minC_", minC, "_minS_", minS, "_RANDOM", sep = "") # <-------------------------
# julia_assign("nam", nam)
# inputs_dir= paste("data/GRN_inputs/", nam, "/", sep = "")
# julia_assign("inputs_dir", inputs_dir)


################################################################################


##create networks file if it does not exist
networks_dir= paste("data/networks/", algo, "/hep_perSample/", nam, "/", sep = "") 
julia_assign("networks_dir", networks_dir)
if(!dir.exists(networks_dir)){dir.create(networks_dir,  recursive = T)}


################################################################################


# Customize the dataset, and algorithm

# Choose an algorithm
julia_command('algorithm = PIDCNetworkInference()')


# Run PIDC 

julia_command('cell = "Hepatocytes"')

julia_command('
for group in groups
  expMat = string(inputs_dir, group, "_", cell, "_expMat.tsv")
  genes = get_nodes(expMat, discretizer = "uniform_width")
  # Infer the network
  network = InferredNetwork(algorithm, genes)
  network_name = string(networks_dir, group, "_", cell, "_LinkList.tsv")
  write_network_file(network_name, network)
end')

# Description of functions above
# get_nodes : Get the genes and discretize the expression levels
# InferredNetwork : Infer the network
# write_network_file : Write the network to file:

# OUTPUT: An InferredNetwork has an array of nodes and an array of edges 
#between all possible node pairs (sorted in descending order of edge weight, 
#i.e. confidence of the edge existing in the true network).

