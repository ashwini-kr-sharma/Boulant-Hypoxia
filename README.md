# Ileum organoids grown under hypoxic and normal oxygen conditions

### Methodology
- This repository hosts the souce code used to analyze the data from the paper - [XXX]()
- The companion website for this paper is hosted [here](https://ashwini-kr-sharma.github.io/Boulant-Hypoxia/)

### Scripts

All script can be found in the `/src/` directory
0. `/src/00_submitJobs.sh` - This is the `bash` script to run analysis `1,2,3,4` in the cluster environment
1. `/src/01_fastQC.sh` - Script to perform `fastQC` analysis on the raw `fastq` files
2. `/src/02_multiQC.sh` - Script to perform `QC analysis` of raw reads and alignment statistics
3. `/src/03_downloadGenomes.sh` - Download the `hg38` reference genome from Ensembl
4. `/src/04_indexRsubread.R` - Indexing of the reference genome
5. `/src/05_alignmentRsubread.R` - Script to perform `read alignment` analysis on the raw `fastq` files using `Rsubread`
6. `/src/06_countsRsubread.R` - Script to `count reads` overlapping genes
7. `/src/07_postProcessing.R` - Script to preprocess and filter the raw count data
8. `/src/08_exploratoryAnalysis.R` - Script to perform various `Exploratory analysis`
9. `/src/09_DSeq2Analysis.R` - Script to perform `Differential gene expression analysis (DGE)` using DSeq2
10. `/src/09_DSeq2_Rmarkdown/` - Script to `visualize` the results of Dseq2 as interactive data tables
11. `/src/10_differentialGeneExpression.R` - Script to `analyze DGE results` focusing on specific genes of interest
12. `/src/11_enrichmentAnalysis.R` -  Script to perform `Enrichment analyses`
13. `/src/12_plotGeneModulesExpression.R` - Script to plot the expression of `Hypoxia`, `EMT` and `WNT` signalling genes
14. `/src/13_computeCellFractions.R` - Using CIBERSORTx to identify Ileum cell types using single-cell data
