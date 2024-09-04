```markdown
# Hepatocellular Carcinoma (HCC) Biomarker Discovery and Enrichment Analysis

This repository contains scripts for the differential expression analysis and enrichment analysis of genes associated with Hepatocellular Carcinoma (HCC) using publicly available GEO data (GSE62232). The analysis is performed using R and involves various bioinformatics packages such as `limma`, `clusterProfiler`, and `org.Hs.eg.db`.

## Table of Contents

- [Introduction](#introduction)
- [Requirements](#requirements)
- [Data Acquisition](#data-acquisition)
- [Differential Expression Analysis](#differential-expression-analysis)
- [Gene Ontology (GO) Enrichment](#gene-ontology-go-enrichment)
- [KEGG Pathway Enrichment](#kegg-pathway-enrichment)
- [PPI Network Analysis](#ppi-network-analysis)
- [Files](#files)

## Introduction

Hepatocellular carcinoma (HCC) is a major form of liver cancer with high mortality rates. This analysis aims to identify differentially expressed genes (DEGs) between tumor and normal liver tissues, and to perform enrichment analysis to uncover biological processes, cellular components, molecular functions, and pathways associated with these DEGs.

## Requirements

The following R packages are required to run the scripts:

- `GEOquery`
- `limma`
- `clusterProfiler`
- `org.Hs.eg.db`

You can install these packages using the following commands in R:

```r
install.packages("BiocManager")
BiocManager::install(c("GEOquery", "limma", "clusterProfiler", "org.Hs.eg.db"))
```

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

## Gene Ontology (GO) Enrichment

Gene Ontology enrichment analysis is conducted using the `clusterProfiler` package to identify the biological processes (BP), cellular components (CC), and molecular functions (MF) associated with the DEGs.

```r
GO_results_BP <- enrichGO(gene = biomarkers$Entrez.id, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", readable = TRUE)
```

## KEGG Pathway Enrichment

KEGG pathway enrichment is performed to identify pathways enriched in the DEGs using the `clusterProfiler` package:

```r
kegg_enrich <- enrichKEGG(gene = biomarkers$Entrez.id, organism = 'hsa')
```

## PPI Network Analysis

Protein-protein interaction (PPI) network analysis is carried out using STRING and Cytoscape. DEGs are mapped onto the PPI network to explore potential interactions and key hub proteins.

## Files

- **HCC_DEGs.csv**: List of differentially expressed genes.
- **upregulated.csv**: Upregulated genes in tumor samples.
- **downregulated.csv**: Downregulated genes in tumor samples.
- **hcc_biomarkers.csv**: Filtered list of biomarkers with unique Entrez IDs.
- **top_10_categories_BP.csv**: Top 10 biological processes enriched in the DEGs.
- **top_10_categories_CC.csv**: Top 10 cellular components enriched in the DEGs.
- **top_10_categories_MF.csv**: Top 10 molecular functions enriched in the DEGs.

## Visualization

Bar plots and dot plots are generated to visualize the enriched GO terms and KEGG pathways. The PPI network is visualized using Cytoscape.

```r
barplot(GO_results_BP, showCategory = 10)
dotplot(kegg_enrich, showCategory = 20, title = "KEGG Pathway Enriched Terms")
```

## Author

Nana Ansah Adomako

## License

This project is licensed under the MIT License - see the LICENSE file for details.
```

This `README.md` file provides an overview of the project, including the purpose, methods, and required files. It is structured to guide users through the analysis steps and offers details on how to replicate the study.
