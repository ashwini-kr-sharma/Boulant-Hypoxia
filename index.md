## Ileum organoids grown under hypoxic conditions

### Experimental design
Ileum derived organoids were grown in hypoxic and normal oxygen concentrations over 6hr, 12hr, 1 and 2 days. Each time point and oxygen concentration levels were replicated four times. RNA sequencing was performed on a total of 32 samples consisting of ileum organoids grown in hypoxia and normoxia conditions.

![experimental design](/Hypoxia_ExpDesign.png)

### Methodology

Raw RNA sequencing reads (fastq) were aligned to the ensembl human GRCh38 genome reference using [Rsubread](https://doi.org/10.1093/nar/gkz114) with default settings. Read summariztion was done using [featureCounts](https://doi.org/10.1093/bioinformatics/btt656). Various quality metrics of the the raw reads and alignment statistics were analysed using [MultiQC](https://doi.org/10.1093/bioinformatics/btw354). Differential gene expression analysis was performed uing DSeq2, rlog transformed data was used for multi dimensional scaling and clustering analyses. Signalling programs were quantified using [PROGENy](https://doi.org/10.1038/s41467-017-02391-6). Transcription factor activities were computed using [DoRothEA](https://doi.org/10.1101/gr.240663.118) and [VIPER](https://doi.org/10.1038/ng.3593). Enrichment analysis on the most differentially expressed genes (-1< logFC >+1 and adjusted p value < 0.05) was performed using [enrichR](https://doi.org/10.1093/nar/gkw377).

### Data Download
Raw and processsed data from this study can be downloaded below. The source code for this study is available [here](https://github.com/ashwini-kr-sharma/Ileum-Hypoxia)

- [Raw fastq files](https://www.ncbi.nlm.nih.gov/gds)
- Raw counts - [RDS](/data/hypoxia_filtered_counts.RDS), [csv](/data/hypoxia_filtered_counts.csv)
- rlog transformed data - [RDS](/data/rlogTransformation.RDS)
- [Quality control - MultiQC](/results/MultiQC/multiqc_report.html)
- Differentially expresed genes (All) - [RDS](/data/diffExpGenes.RDS), [xls](/data/DGEtables.xls)
- Differentially expressed genes -
  - [Hypoxia vs Normoxia (6hr)](/src/09_DSeq2_Rmarkdown/Hypoxia_vs_Normoxia_6hr.html)
  - [Hypoxia vs Normoxia (12hr)](/src/09_DSeq2_Rmarkdown/Hypoxia_vs_Normoxia_12hr.html)
  - [Hypoxia vs Normoxia (D1)](/src/09_DSeq2_Rmarkdown/Hypoxia_vs_Normoxia_D1.html)
  - [Hypoxia vs Normoxia (D2)](/src/09_DSeq2_Rmarkdown/Hypoxia_vs_Normoxia_D2.html)
- [log2 fold changes](/src/09_DSeq2_Rmarkdown/log2_foldchange.html) [(RDS)](/data/diffExpLogFCmatrix.RDS)
- [Signalling programs](/data/progeny_all_results.csv)
- [Transcription factor activities](data/tfactivity_all_results.csv)
- CIBERSORTx deconvolution - [ABSOLUTE](/results/Deconvolution/CIBERSORTx_Absolute_Results.txt), [RELATIVE](/results/Deconvolution/CIBERSORTx_Relative_Results.txt)

