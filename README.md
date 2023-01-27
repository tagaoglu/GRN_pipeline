# Gene Regulatory Network Inference from Single-cell Transcriptomic Data in Hepatocellular Carcinoma

This repository contains the code for all the analysis relative to my master thesis « Gene Regulatory Network Inference from Single-cell Transcriptomic Data in Hepatocellular Carcinoma », in order to ensure the reproducibility of all the results and figures presented in this study.

<br>

**/functions**
- additional functions

<br>

**Data_preprocessing.R   &   Data_preprocessing_hep_perSample.R**
- Gene Filtering and Normalization
- Create expression matrix to use as input for GRN algorithms

<br>

**/runGENIE3**
- scripts to submit jobs in cluster
- codes to run GENIE3

<br>

**/runPIDC**
- scripts to submit jobs in cluster
- codes to run PIDC

<br>

**Filtering_MostLikelyLinks.R**
- Extract and Save filtered versions of networks with the most likely links for later use

<br>

**Similarity_analysis.R** 
- Evaluation of similarity between networks, with all cells and with random cells (25%)
- Evaluation of similarity between networks, constructed by GENIE3 and PIDC
- Plot boxplot for JI
- Plot heat table

<br>

**TF_Motif_Analysis.R**
- Run correlation, and extract active links and save
- Run RcisTarget and save

<br>

**Network_analysis.R** 
- Compute Centrality Measures (CM) for each node in each network : Degree, Betweenness, Closeness, Eigenvector, Pagerank
- then also compute integrated centrality measure (aggregated rank)
- Compute Correlation between CM
- Rank nodes by decreasing centrality measure, and then Define as "central" genes those within the top 20% of the ranking
- Compute the overlapping between each set of hubs, visualizing with Venn diagrams
- Test to validate the networks whether they are scale-free, by using a Kolmogorov-Smirnov test to compare the degree distribution of the networks to a theoretical power-law distribuiton
  - Visualizing scale-free distribution: histograms
  - Power-law fit to the degree distibution: scatter plots
- Compute and plot the multiplicity of each hub (i.e. # groups a given gene acts as a hub), in order to assess whether these hubs are specific to a       particular group (group-specific hubs) or whether they are hubs in several groups (ubiquitous hubs)
- Define group-specific and shared hubs

<br>

**Critical_genes_analysis.R**
- Critical Genes Identification
- Calculate 5 centrality metrics and then integrated centrality metrics
- Plot Top10 Nodes 
- Create SummaryTable for Common&OverlappingSpecificNodes 
- Create Venn diagram to show the overlap of Nodes, and also to show the overlap of edges (gene pairs) across group-specific GRNs
- Network Visualization with Top10 critical TFs highlighted

<br>

**Expr_level.R**
- Create violin plot gene expression

<br>

**goea_analysis.R** 
- Perform Gene Ontology (using enrichGO) & KEGG pathway (using enrichKEGG) & MSigDb (using enricher) Enrichment Analysis
- Visualize results using dotplots
    
<br>

**Differential_network_analysis.R**
- Run CoDiNA algorithm to perform differential network analysis
- Plot the network
- Create table to summarize candidate key nodes
- Perform Gene Ontology (using enrichGO) & KEGG pathway (using enrichKEGG) & MSigDb (using enricher) Enrichment Analysis
- Visualize results using dotplots

<br>

**/figures**
- figures published in this study

## List of Figures
- **Figure 1: Summary of the workflow followed in this study**
- **Figure 2: Schematic diagram of GRN inference followed in this study**
- **Figure 3: Schematic diagram of network centrality analysis followed in this study**
- **Figure 4: Schematic diagram of differential network analysis followed in this study**
- **Figure 5: Single-cell RNA-sequencing data of HCC and non-tumoral/normal livers with annotated cell clusters**
- **Figure 6: Boxplots showing Jaccard Index of edges in two compared networks, one inferred by GENIE3 and the other inferred by PIDC, for different link thresholds**
- **Figure 7: Single-cell gene regulatory networks inferred by GENIE3 are scale-free.**         
The degree distribution of the single-cell gene regulatory networks inferred by GENIE3, derived for four phenotypes, nonviral, viral, normal and tumor, and for each cell type, displayed in linear (histogram) and logarithmic scale (scatter plot). Each distribution was fitted to a theoretical power-law distribution, and the p-value of the Kolmogorov-Smirnov test (KS.p) and the degree exponent of the power-law (alpha) are labelled for each network. 
- **Figure 8: GRNs of each phenotype in Hepatocytes, inferred by GENIE3**       
Colored nodes imply the top 10 critical genes. (A-D) The GRNs of nonviral (A), viral (B), normal (C), tumor (D) 
- **Figure 9: Relationship between degree and other centralities for the GRNs inferred by GENIE3, for Hepatocytes**        
Scatter plots displaying relationship between degree and the other centralities (log-scale) for the GRNs inferred by GENIE3, derived for four phenotypes, nonviral, viral, normal and tumor, and for Hepatocytes. The p-value of the spearman correlation test (p) and the correlation coefficient (R) are labelled. 
- **Figure 10: The central genes of different metrics show overlap, but there are also some central genes which are specific for each metric**     
Venn diagrams showing the overlap of central genes (top 20% of ranked genes) of different centralities in the networks of four phenotypes, nonviral, viral, normal and tumor, for each cell type separately. 
- **Figure 11: Results for Hepatocytes, obtained from GRNs inferred by GENIE3**       
(A - B) Venn diagrams showing the overlap of network nodes (A) and network edges (B) across three phenotypes, normal, viral and nonviral. (C - D) Venn diagrams showing the overlap of network nodes (C) and network edges (D) across two phenotypes, normal and tumor. (E) The central genes (top 20%) of each of three phenotypes, normal, viral and nonviral, classified by their multiplicity for all tested centrality metrics. (F) The central genes (top 20%) of each of two phenotypes, normal and tumor, classified by their multiplicity for all tested centrality metrics. Multiplicity=1 means that they are central only in that particular phenotype, whereas multiplicity=2 and 3 means that they are central also in other phenotypes. (G - J) Ranks for central genes in the four phenotype-specific networks based on the integrated centrality measure. For each network, the top 10 TFs are labeled and colored in red. The blue arrow marks the critical genes present in both Top10 genes list and phenotype-specific hubs. 
- **Figure 12: Violin plots displaying expression levels of MT2A in each sample for Hepatocytes**     
- **Figure 13: Functional exploration of phenotype-specific hubs in Hepatocytes**       
(A - B) Dot plots displaying enriched GO terms (BP) (A) and enriched Hallmark terms from MSigDB (B) in viral HCC, nonviral HCC and normal liver specific hubs. (C - D) Dot plots showing enriched GO terms (BP) (C) and enriched Hallmark terms from MSigDB (D) in tumor and normal liver specific hubs. The size of nodes represents the ratio of enriched genes; the color of the nodes corresponds to the adjusted p values. 

