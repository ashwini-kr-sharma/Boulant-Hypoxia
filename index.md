## Ileum organoids grown under hypoxic conditions

### Experimental design

![experimental design](/Hypoxia_ExpDesign.png)

### Methodology

Raw RNA sequencing reads (fastq) were aligned to the ensembl human GRCh38 genome reference using [Rsubread](https://doi.org/10.1093/nar/gkz114) with default settings. Read summariztion was done using [featureCounts](https://doi.org/10.1093/bioinformatics/btt656). Various quality metrics of the the raw reads and alignment statistics were analysed using [MultiQC](https://doi.org/10.1093/bioinformatics/btw354). Differential gene expression analysis was performed uing DSeq2, rlog transformed data was used for multi dimensional scaling and clustering analyses. Signalling programs were quantified using [PROGENy](https://doi.org/10.1038/s41467-017-02391-6). Transcription factor activities were computed using [DoRothEA](https://doi.org/10.1101/gr.240663.118) and [VIPER](https://doi.org/10.1038/ng.3593). Enrichment analysis on the most differentially expressed genes (-1< logFC >+1 and adjusted p value < 0.05) was performed using [enrichR](https://doi.org/10.1093/nar/gkw377).

### Data Download
Raw and processsed data from this study can be downloaded below. The source code for this study is available [here](https://github.com/ashwini-kr-sharma/Boulant-Hypoxia)

- [Raw fastq files](https://www.ncbi.nlm.nih.gov/gds)
- Raw counts - [RDS](/data/T84_IL22_INFL_filtered_counts.RDS), [csv](/data/T84_IL22_INFL_filtered_counts.csv)
- [Quality control - MultiQC](/data/multiqc_report.html)
- Differentially expressed genes -
  - [Hypoxia vs Normoxia (6hr)](/data/DGE/IL22_3hr_vs_Mock_3hr.html)
  - [Hypoxia vs Normoxia (12hr)](/data/DGE/IL22_6hr_vs_Mock_6hr.html)
  - [Hypoxia vs Normoxia (D1)](/data/DGE/IL22_12hr_vs_Mock_12hr.html)
  - [Hypoxia vs Normoxia (D2)](/data/DGE/IL22_24hr_vs_Mock_24hr.html)
- [log2 fold changes](/data/DGE/log2_fold_change.html)
- [Signalling programs](/data/progeny_all_results.csv)
- [Transcription factor activities](data/tfactivity_all_results.csv)

