#' @description This function calculates 5 centrality metrics, and also integrates the five centrality metrics of the nodes for each GRN
#'              degree, closeness, betweenness, eigenvalue, and PageRank :
#'              node centrality metrics which are employed to measure the importance of nodes in the constructed networks.
#' @param network_edge_list network in the form of "edge list"
#' @param edge_prediction directed: TRUE or FALSE
#'  
#' @output The function returns a data.frame object with integrated centrality metric with corresponding p-values for each node


centralityMetrics <- function(network_edge_list, edge_prediction)
{
  ## Load package
  library(igraph)
  
  ## Reconstructing the graph
  net <- graph_from_data_frame(network_edge_list, directed=edge_prediction)
  
  # Extracting the vertices
  graph_vertices <- V(net)
  
  # Calculate node centrality metrics which are employed to measure the importance of nodes in the constructed networks
  
  # Degree centrality (number of ties)
  #Degree is the number of adjacent nodes of the corresponding node. Node with high degree is usually considered a hub with essential functions.
  deg <- degree(net, mode="all",loops = F)
  deg_rank <- deg%>%sort(decreasing = TRUE)
  names(deg_rank)
  
  # Closeness centrality (centrality based on distance to others in the graph)
  #Closeness measures the average distance to all nodes for the corresponding node. Node with high closeness denotes that the node is located in the central location of network. 
  closeness <- closeness(net, mode="all", weights=NA)
  closeness_rank <- closeness%>%sort(decreasing = TRUE)
  names(closeness_rank)
  
  # Vertex Betweenness centrality (centrality based on a broker position connecting others)
  #Betweenness is calculated based on the number of shortest paths passed through the corresponding node. Node with high betweenness could be the bottleneck of the GRN.
  betweenness <- betweenness(net, directed=edge_prediction, weights=NA)
  betweenness_rank <- betweenness%>%sort(decreasing = TRUE)
  names(betweenness_rank)
  
  ###!!!
  #Note: I may need to exclude eigenvalue, because for GENIE_wRcis_centralityMetrics calculation I may encounter the following problem:
  #ERROR message was: At core/centrality/centrality_other.c:328 : graph is directed and acyclic; eigenvector centralities will be zeros.
  #PROBLEM explained here: http://users.dimi.uniud.it/~massimo.franceschet/teaching/datascience/network/katz.html
  
  ##!!!
  # Eigenvector centrality (centrality proportional to the sum of connection centralities)
  #The eigenvalue measures the node importance by taking the importance of neighbors into consideration.
  eigen <- eigen_centrality(net, directed=edge_prediction, weights=NA)$vector
  eigen_rank <- eigen%>%sort(decreasing = TRUE)
  names(eigen_rank)
  
  #PageRank
  #PageRank is the probability of the random walk of the corresponding node. PageRank is similar to eigenvalue, whereas PageRank introduces the damping factor (default 0.85)
  pageRank <- page_rank(net)$vector
  pageRank_rank <- pageRank%>%sort(decreasing = TRUE)
  names(pageRank_rank)
  
  ##!!!
  #Put it all together in one data frame.
  allMetrics <- data.frame(row.names   = names(graph_vertices),
                           degree      = deg,
                           closeness   = closeness,
                           betweenness = betweenness,
                           eigenvector = eigen,
                           pageRank    = pageRank)
  
  # ##!!!
  # #Put it all together in one data frame.
  # allMetrics <- data.frame(row.names   = names(graph_vertices),
  #                          degree      = deg,
  #                          closeness   = closeness,
  #                          betweenness = betweenness,
  #                          pageRank    = pageRank)
  
  
  #create rank list for the next step
  rank_list <- list()
  rank_list[[1]] <- names(deg_rank)
  names(rank_list)[1] <- "degree"
  
  rank_list[[2]] <- names(closeness_rank)
  names(rank_list)[2] <- "closeness"
  
  rank_list[[3]] <- names(betweenness_rank)
  names(rank_list)[3] <- "betweenness"
  
  ##!!!
  rank_list[[4]] <- names(eigen_rank)
  names(rank_list)[4] <- "eigenvalue"

  rank_list[[5]] <- names(pageRank_rank)
  names(rank_list)[5] <- "pageRank"
  
  # ##!!!
  # rank_list[[4]] <- names(pageRank_rank)
  # names(rank_list)[4] <- "pageRank"
  
  
  ### Integrate the five centrality metrics of the nodes for each GRN
  
  #Apply a rank aggregation method "RRA" to integrate all centrality metrics. 
  #Finally, the rank integration score for each TF is calculated as –log10(p value) to indicate the final centrality measure.
  #https://cran.r-project.org/web/packages/RobustRankAggreg/RobustRankAggreg.pdf
  
  library(RobustRankAggreg) 
  
  # Aggregate the inputs
  r = RobustRankAggreg::rankMatrix(rank_list, full = TRUE)
  integratedMetrics = RobustRankAggreg::aggregateRanks(rmat = r, method = "RRA")
  integratedMetrics$minuslog10 <- -log10(integratedMetrics$Score)
  
  integratedMetrics$rank <- 1:dim(integratedMetrics)[1]
  
  centralityMetrics <- list("allMetrics" = allMetrics, "integratedMetrics" = integratedMetrics)
  
  return(centralityMetrics)
}



