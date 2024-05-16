# Osteoarthritis PBMCs Gene Expression Analysis

This project focuses on the analysis of gene expression profiles in peripheral blood mononuclear cells (PBMCs) from patients with osteoarthritis (OA). The objective is to gain insights into the feasibility of using gene expression profiling of PBMCs to detect the onset of osteoarthritis. The analysis involves several steps including data download, exploratory data analysis (EDA), preprocessing, differential expression analysis (DEA), and gene set enrichment analysis.

## Data Download and Preprocessing

The first step involves downloading the gene expression dataset from the GEO database with the accession number "GSE48556". The dataset is processed for background correction, normalization, and log2 transformation. Exploratory data analysis is performed to assign groups and set up the design matrix for subsequent analyses.

## Differential Expression Analysis (DEA)

DEA is conducted using the t-test and limma package to identify genes that are differentially expressed between osteoarthritis and control groups. Holm correction is applied to adjust p-values, and volcano plots are generated to visualize the results.

## Gene Set Enrichment Analysis

Gene set enrichment analysis is performed to identify enriched biological pathways associated with differentially expressed genes. The clusterProfiler package is utilized for enrichment analysis, and various plots such as bar plots, dot plots, enrichment map plots, and category net plots are generated to visualize the results.

## Results and Observations

### Differential Expression Analysis:

- The t-test and limma package are used to identify significantly differentially expressed genes between osteoarthritis and control groups.
- Volcano plots show genes with significant fold changes and adjusted p-values.

### Gene Set Enrichment Analysis:

- Enrichment analysis reveals biological pathways enriched with differentially expressed genes.
- Bar plots, dot plots, enrichment map plots, and category net plots are used to visualize enriched pathways and functional clusters.

## Interpretation of Results

- Significance cutoffs of p-value ≤ 0.05 and log fold change (logFC) ≥ 1 or ≤ -1 are chosen based on statistical significance and biological relevance.
- Enriched pathways provide insights into the molecular mechanisms underlying osteoarthritis pathogenesis and potential therapeutic targets.

## Conclusion

This analysis demonstrates the utility of gene expression profiling of PBMCs in identifying biomarkers and pathways associated with osteoarthritis. The identified genes and pathways may serve as potential targets for further research and drug development in the field of osteoarthritis.
