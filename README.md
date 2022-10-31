# Gene Regulatory Network Inference from Single-cell Transcriptomic Data in Hepatocellular Carcinoma

This repository contains the code for all the analysis relative to my master thesis « Gene Regulatory Network Inference from Single-cell Transcriptomic Data in Hepatocellular Carcinoma », in order to ensure the reproducibility of all the results and figures presented in this study.


**/functions**
- additional functions


***Data_preprocessing.R & Data_preprocessing_hep_perSample.R***
- Gene Filtering and Normalization
- Create expression matrix to use as input for GRN algorithms


***/runGENIE3 
- codes to submit jobs in cluster
- codes to run GENIE3


***/runPIDC 
- codes to submit jobs in cluster
- codes to run PIDC


***Filtering_MostLikelyLinks.R
- Extract and Save filtered versions of networks with the most likely links for later use


***Similarity_analysis.R 
- Evaluation of similarity between networks, with all cells and with random cells (25%)
- Evaluation of similarity between networks, constructed by GENIE3 and PIDC
- Plot boxplot for JI
- Plot heat table


***TF_Motif_Analysis.R
- Run correlation, and extract active links and save
- Run RcisTarget and save


***Network_analysis.R 
- Compute Centrality Measures (CM) for each node in each network : Degree, Betweenness, Closeness, Eigenvector, Pagerank
- then also compute integrated centrality measure (aggregated rank)
- Compute Correlation between CM
- Rank nodes by decreasing centrality measure, and then Define as "central" genes those within the top 20% of the ranking
- Compute the overlapping between each set of hubs, visualizing with Venn diagrams
- Test to validate the networks whether they are scale-free, by using a Kolmogorov-Smirnov test to compare the degree distribution of the networks to a theoretical power-law distribuiton
  - 1. Visualizing scale-free distribution: histograms
  - 2. Power-law fit to the degree distibution: scatter plots
- Compute and plot the multiplicity of each hub (i.e. # groups a given gene acts as a hub), in order to assess whether these hubs are specific to a       particular group (group-specific hubs) or whether they are hubs in several groups (ubiquitous hubs)
- Define group-specific and shared hubs


***Critical_genes_analysis.R
- Critical Genes Identification
- Calculate 5 centrality metrics and then integrated centrality metrics
- Plot Top10 Nodes 
- Create SummaryTable for Common&OverlappingSpecificNodes 
- Create Venn diagram to show the overlap of Nodes, and also to show the overlap of edges (gene pairs) across group-specific GRNs
- Network Visualization with Top10 critical TFs highlighted


***goea_analysis.R 
- Perform Gene Ontology (using enrichGO) & KEGG pathway (using enrichKEGG) & MSigDb (using enricher) Enrichment Analysis
- Visualize results using dotplots
    

***Differential_network_analysis.R
- Run CoDiNA algorithm to perform differential network analysis
- Plot the network
- Create table to summarize candidate key nodes
- Perform Gene Ontology (using enrichGO) & KEGG pathway (using enrichKEGG) & MSigDb (using enricher) Enrichment Analysis
- Visualize results using dotplots

