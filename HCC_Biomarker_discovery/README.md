
## Data-Driven (In Silico) Biomarker Discovery in Hepatocellular Carcinoma (HCC)

This repository contains scripts for the differential expression analysis and enrichment analysis of genes associated with Hepatocellular Carcinoma (HCC) using publicly available GEO data (GSE62232). The analysis is performed using R, STRING and two Cytoscape apps (MCODE and CytoHubba). The R packages used were `limma`, `clusterProfiler`, and `org.Hs.eg.db`.



### Introduction

Hepatocellular carcinoma (HCC) is a major form of liver cancer with high mortality rates. This analysis aims to identify differentially expressed genes (DEGs) between tumor and normal liver tissues. Highly ranked DEGs were subjected to PPI Network Analysis to obtain subcluster and hub genes as potential biomarkers implicated in HCC.  Enrichment analysis were executed to uncover biological processes, cellular components, molecular functions, and pathways these potential biomarkers were associated with.

## Requirements

The following R packages are required to run the scripts:

- `GEOquery`
- `limma`
- `clusterProfiler`
- `org.Hs.eg.db`

The following tools are required for Biomarker discovery:
- `STRING`
- `Cytoscape (MCODE & CytoHubba)` 

You can install these packages using the following commands in R:

```r
install.packages("BiocManager")
BiocManager::install(c("GEOquery", "limma", "clusterProfiler", "org.Hs.eg.db"))
```



## Workflow


![workflow (2)](https://github.com/user-attachments/assets/6fa3c5ff-b341-40f0-b751-869a006cadc4)



## Data Acquisition

The GEO dataset `GSE62232` is downloaded using the `GEOquery` package. The expression set is extracted for further analysis:

```r
geo_data <- getGEO("GSE62232", GSEMatrix = TRUE)
expression_set <- geo_data[[1]]
```

## Differential Expression Analysis

1. **Design Matrix Creation**: A design matrix is constructed based on the groupings in the phenotypic data.
2. **Contrasts Definition**: Contrasts are defined to compare tumor samples against normal liver tissues.
3. **Linear Model Fitting**: The `lmFit` function from the `limma` package is used to fit the linear model.
4. **Significance Testing**: eBayes is applied to estimate the differential expression, and DEGs are filtered based on adjusted p-value and log fold change.

```r
fit <- lmFit(exprs(expression_set), design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
significant_degs <- topTable(fit2, adjust="fdr", number=Inf)
```


## PPI Network Analysis

Protein-protein interaction (PPI) network analysis is carried out using STRING and Cytoscape. DEGs are mapped onto the PPI network to explore potential interactions and key hub proteins.


## Gene Ontology (GO) Enrichment

Gene Ontology enrichment analysis is conducted using the `clusterProfiler` package to identify the biological processes (BP), cellular components (CC), and molecular functions (MF) associated with the identified biomarkers.

```r
GO_results_BP <- enrichGO(gene = biomarkers$Entrez.id, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", readable = TRUE)
```

## KEGG Pathway Enrichment

KEGG pathway enrichment is performed to identify pathways enriched in the DEGs using the `clusterProfiler` package:

```r
kegg_enrich <- enrichKEGG(gene = biomarkers$Entrez.id, organism = 'hsa')
```



## Visualization

The PPI network is visualized using Cytoscape. Bar plots and dot plots are generated to visualize the enriched GO terms and KEGG pathways. 
```r
barplot(GO_results_BP, showCategory = 10)
dotplot(kegg_enrich, showCategory = 20, title = "KEGG Pathway Enriched Terms")
```

## Author

Nana Ansah Adomako
