#Algorithm comparison/evaluation
#Metrics are calculated: Intersection index, Weighted Jaccard Similarity (WJS)

# These functions are modified versions from:
#
#REFERENCE: 
#https://github.com/ComputationalSystemsBiology/scNET



library(igraph)
library(tools)


#Percentage of intersection (perINT) is used to detect the presence of links shared between two compared networks
#i.e. to test the amount of common links between the two networks
#Overlap Coefficient : the size of the intersection of set A and set B over the size of the smaller set between A and B.
#Jaccard Index or Jaccard Similarity Coefficient : the size of the intersection of set A and set B (i.e. the number of common elements) over the size of the union of set A and set B (i.e. the number of unique elements).
Intersection.index= function(netA, netB, directed1=T, directed2=F){
  if(nrow(netA) > nrow(netB)){
    trans=netA
    netA=netB
    netB=trans
  }
  
  charA= paste(netA[,1], ':', netA[,2])
  charB= paste(netB[,1], ':', netB[,2])
  
  if (directed1==T & directed2==T){
    inter= intersect(charA, charB)
  }else if (directed1==F & directed2==T) {
    charA2= paste(netA[,2], ':', netA[,1])
    charA2= c(charA, charA2)
    inter= intersect(charA2, charB)
  }else if (directed1==T & directed2==F) {
    charB2= paste(netB[,2], ':', netB[,1])
    charB2= c(charB, charB2)
    inter= intersect(charA, charB2)
  }else if (directed1==F & directed2==F) {
    charA2= paste(netA[,2], ':', netA[,1])
    charA2= c(charA, charA2)
    charB2= paste(netB[,2], ':', netB[,1])
    charB2= c(charB, charB2)
    inter= intersect(charA2, charB2)
  }
  
  overlap.index = length(inter)/min(length(charA), length(charB))
  jaccard.index  = length(inter)/length(union(charA, charB))
  index = c(overlap.index, jaccard.index)
  return(index)
}


#Weighted Jaccard Similarity (WJS) takes into account the similarity of the weights 
#associated with the links shared between the compared networks 
wjs= function(netA, netB, directed1=T, directed2=F){
  colnames(netA)= c('gene1', 'gene2', 'weight')
  colnames(netB)= c('gene1', 'gene2', 'weight')
  #create weighted adjacency matrixes if input is an edgelist
  AG= graph.data.frame(netA, directed = directed1) 
  Aadj= as_adjacency_matrix(AG, attr = 'weight')
  BG= graph.data.frame(netB, directed = directed2) 
  Badj= as_adjacency_matrix(BG, attr = 'weight')
  
  #match rows and columns for both adjacency lists
  comm= intersect(row.names(Aadj), row.names(Badj))
  Badj=Badj[comm,comm]
  Aadj=Aadj[comm,comm]
  
  #ensure that everything matches
  if(all(rownames(Aadj)==rownames(Badj)) &
     all(colnames(Aadj)==colnames(Badj)) ){
    Aadj= abs(Aadj)
    Badj= abs(Badj)
    mn= pmin(Aadj, Badj)
    mx= pmax(Aadj, Badj)
    
    #calculate
    mn= sum(mn)
    mx= sum(mx)
    wjs= mn/mx
  }
  return(wjs)
}

#Function for jaccard similarity of edge sets
#to calculate the Jaccard similarity of the edge sets of two graphs. 
#Important to note that This function can be used to compare graphs where the vertex set is identical(where the vertices are the same set of genes) but where the edge set may be different(where the edge sets are different according to interactions)
#Note that the igraph package has special operators %s% for the intersection of two graphs (that is, all vertices and edges that both graphs have in common), 
#and %u% for the union of two graphs (that is, all unique vertices and edges in both graphs combined).
#Note that a Jaccard similarity of 1 means that both graphs have identical edge sets and so are identical in structure, and a similarity of 0 means that both graphs have no edges in common. 
#For instance, getting a jaccard_edgeset_similarity of 0.09727385 means that there is only about 10% similarity between 2 networks/edgelists
jaccard_edgeset_similarity <- function(netA, netB, directed1=T, directed2=F) {
  #Convert the data to an igraph object
  G1= graph_from_data_frame(netA, directed = directed1) 
  G2= graph_from_data_frame(netB, directed = directed2) 
  
  if(directed1==F & directed2==T){
    #To do a fair comparison, redefine our graph to be a directed graph where edges always go in both directions.
    G1 <- as.directed(G1) #create directed version of network
  }else if(directed1==T & directed2==F){
    G2 <- as.directed(G2) 
  }else if(directed1==F & directed2==F){
    G1 <- as.directed(G1) 
    G2 <- as.directed(G2) 
  }
  
  inter <- length(E(G1 %s% G2))
  un <- length(E(G1 %u% G2))
  
  if (un == 0) {
    0
  } else {
    inter/un
  }
}











