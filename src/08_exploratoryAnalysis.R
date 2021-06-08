#-------------------------------------------------------------------------------
# Required packages, setting root directory and reading the count data
#-------------------------------------------------------------------------------

library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)
library(progeny)
library(reshape2)
library(GSVA)
library(GO.db)
library(org.Hs.eg.db)
library(ggvenn)
library(patchwork)
library(readxl)
library(ggpubr)
library(dorothea)

path = "/icgc/dkfzlsdf/analysis/B080/sharma/boulantLab/bulk_hypoxia/"
dat  = readRDS(paste0(path, "data/counts/hypoxia_filtered_counts.RDS"))

#-------------------------------------------------------------------------------
# Create directories and download two independent signatures for stemness
#-------------------------------------------------------------------------------

# Create output directory 
if( ! dir.exists(paste0(path, "analysis/Exploratory"))){
  dir.create(paste0(path, "analysis/Exploratory"),recursive = T)
}

# Create directory to store the downloaded signature
if( ! dir.exists(paste0(path, "data/extdata"))){
dir.create(paste0(path, "data/extdata"),recursive = T)
}

# Stemness signature https://doi.org/10.1016/j.cell.2018.03.034
if( ! file.exists(paste0(path, "data/extdata/Stemness_Signature_PMID_29625051.xlsx"))){
  download.file(url = "https://api.gdc.cancer.gov/data/401e9e48-12d2-4177-81a7-10043a32867c",
                destfile = paste0(path, "data/extdata/Stemness_Signature_PMID_29625051.xlsx"))
}

# Stemness signature https://doi.org/10.1073/pnas.1818210116
if( ! file.exists(paste0(path, "data/extdata/Stemness_Signature_PMID_30996127.xlsx"))){
download.file(url = "https://www.pnas.org/highwire/filestream/860118/field_highwire_adjunct_files/1/pnas.1818210116.sd01.xlsx",
              destfile = paste0(path, "data/extdata/Stemness_Signature_PMID_30996127.xlsx"))
}

#-------------------------------------------------------------------------------
# DESeq2 object
#-------------------------------------------------------------------------------

dds = DESeqDataSetFromMatrix(countData = dat$counts,
                             colData   = dat$sampanno,
                             rowData   = dat$geneanno,
                             design    = ~ Type + Time)

rm(dat)

#-------------------------------------------------------------------------------
# Data transformation - rlog
#-------------------------------------------------------------------------------

rld = rlog(dds, blind = TRUE)

#-------------------------------------------------------------------------------
# Sample similarity heatmap based on Euclidean distance
#-------------------------------------------------------------------------------

# Compute distance
sampleDists = dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) = colData(dds)$SampleID
colnames(sampleDistMatrix) = colData(dds)$SampleID

# Annotations for plotting
colors <- colorRampPalette( rev(brewer.pal(9, "Greens")) )(10)
col_anno = as.data.frame(colData(dds), row.names = 1 )[,1:2]
col_anno$Library_Size = colSums(dds@assays@data$counts)/10^6

# Sanity check !!
identical(colnames(sampleDistMatrix), rownames(col_anno))
#[1] TRUE

ann_colors = list(
  Type = c(Normoxic = "#cb181d", Hypoxic = "#0570b0"),
  Time = c(`6h` = "#cccccc", `12h` = "#969696", `D1` = "#636363", `D2` = "#252525")
)

# Heatmap
pheatmap(sampleDistMatrix, clustering_method = "ward.D2",
         clustering_distance_rows = sampleDists,
         cutree_rows = 6, cutree_cols = 8, fontsize = 6,
         clustering_distance_cols = sampleDists,
         annotation_col = col_anno,
         annotation_colors = ann_colors,
         border_color = NA, 
         col = colors,
         treeheight_row = 10, treeheight_col = 10,
         filename = paste0(path, "analysis/Exploratory/sample_similarity_heatmap_cluster.pdf"), width = 4.8, height = 3.8)

# Heatmap
sample_order = c(paste0("N_6h_",1:4), paste0("H_6h_",1:4),
                 paste0("N_12h_",1:3), paste0("H_12h_",1:2),
                 paste0("N_D1_",1:4), paste0("H_D1_",1:4),
                 paste0("N_D2_",1:4), paste0("H_D2_",1:4))

# Sanity check !!
identical(sort(sample_order), sort(colnames(sampleDistMatrix)))
#[1] TRUE

sampleDistMatrix = sampleDistMatrix[sample_order, sample_order]

# Sanity check !!
identical(sort(sample_order), sort(rownames(col_anno)))
#[1] TRUE

col_anno = col_anno[sample_order,]

# Sanity check !!
identical(colnames(sampleDistMatrix), rownames(col_anno))
#[1] TRUE

pheatmap(sampleDistMatrix, 
         cluster_rows = F, cluster_cols = F,
         fontsize = 6,
         annotation_col = col_anno,
         annotation_colors = ann_colors,
         border_color = "white", 
         col = colors,
         filename = paste0(path, "analysis/Exploratory/sample_similarity_heatmap_orderedBy_Name.pdf"), width = 6, height = 5)

rm(colors, ann_colors, sample_order, sampleDists)

#-------------------------------------------------------------------------------
# Multidimensional scaling (MDS) plot
#-------------------------------------------------------------------------------

mds <- cbind(col_anno, cmdscale(sampleDistMatrix))

ggplot(mds, aes(x = `1`, y = `2`, color = Time, shape = Type)) + 
theme_bw(base_size = 8) + labs(x = "MDS dimension 1", y = "MDS dimension 2") +
geom_point(size = 2) +  ggtitle("MDS - rlog transformed data") +
scale_color_manual(values = c(`6h` = "#cccccc", `12h` = "#969696", `D1` = "#636363", `D2` = "#252525")) +
theme(panel.grid = element_blank(),
      axis.text = element_text(colour="black"), 
      axis.line = element_blank(),
      legend.position = "right") +
ggsave(filename = paste0(path, "analysis/Exploratory/sample_similarity_MDS.pdf"), width = 3.2, height = 2.1)

rm(sampleDistMatrix, mds, col_anno)

#-------------------------------------------------------------------------------
# Normalized expression data to use for all the analysis below
#-------------------------------------------------------------------------------

expr = assay(rld)

# Sanity check !!
identical(colnames(expr), rownames(colData(rld)))
#[1] TRUE

colnames(expr) = colData(rld)$SampleID
saveRDS(object = expr, file = paste0(path, "analysis/Exploratory/rlogTransformation.RDS"))

#-------------------------------------------------------------------------------
# Remove duplicated genes
#-------------------------------------------------------------------------------

# Sanity check !!
if(identical(rownames(expr), rowData(rld)$GeneID))
{
  rowID = rowData(rld)[, c("GeneID", "gene_name")]
  expr  = expr[- which(duplicated(rowID$gene_name)),]
  rowID = rowID[- which(duplicated(rowID$gene_name)),]
  
  # Sanity check !!
  if(identical(rownames(expr), rowID$GeneID))
  {
    rownames(expr) = rowID$gene_name
  }else{
    expr = NULL
  }
  rm(rowID)
}

#-------------------------------------------------------------------------------
# Progeny analysis for pathway activities - Compute pathway activity scores
#-------------------------------------------------------------------------------

pr_res = progeny(expr = expr)

write.csv(data.frame(Samples = rownames(pr_res), pr_res, stringsAsFactors = F),
          paste0(path, "analysis/Exploratory/progeny_all_results.csv"), 
          quote = FALSE, row.names = F)

pr_res_summ = split(data.frame(pr_res), 
                    sapply(strsplit(rownames(pr_res), "_"), function(x) paste0(x[1:2], collapse ="_")))

pr_res_summ = t(sapply(pr_res_summ, function(x){
  apply(x, 2, median)
}))

pr_res_summ = pr_res_summ[c("N_6h", "N_12h", "N_D1", "N_D2",
                            "H_6h", "H_12h", "H_D1", "H_D2"),]

paletteLength = 100
cols = colorRampPalette(rev(brewer.pal(9, "PiYG")))(paletteLength)
brk = c(seq(min(pr_res_summ), 0, length.out=ceiling(paletteLength/2) + 1),
        seq(max(pr_res_summ)/paletteLength,  max(pr_res_summ), length.out=floor(paletteLength/2)))

col_anno <- data.frame(do.call("rbind", strsplit(rownames(pr_res_summ), "_")), row.names = rownames(pr_res_summ))
colnames(col_anno) = c("Treatment", "Time")
col_anno$Treatment = ifelse(col_anno$Treatment == "H", "Hypoxic", "Normoxic")

# Sanity check
identical(rownames(col_anno), rownames(pr_res_summ))
#[1] TRUE

ann_colors = list(
  Treatment = c(Normoxic = "#cb181d", Hypoxic = "#0570b0"),
  Time = c(`6h` = "#cccccc", `12h` = "#969696", `D1` = "#636363", `D2` = "#252525")
)

#-------------------------------------------------------------------------------
# Progeny analysis for pathway activities - Heatmap of activity scores
#-------------------------------------------------------------------------------

pheatmap(t(pr_res_summ), clustering_method = "ward.D2",
         cluster_cols  = F,
         cutree_rows = 2, 
         fontsize = 6,
         annotation_col = col_anno,
         annotation_colors = ann_colors,
         border_color = "white", 
         col = cols,
         breaks = brk,
         treeheight_row = 15, treeheight_col = 15,
         filename = paste0(path, "analysis/Exploratory/progeny_pathway_scores_heatmap.pdf"), width = 2.8, height = 2.5)

rm(cols, brk, paletteLength, col_anno, ann_colors, pr_res_summ)

#-------------------------------------------------------------------------------
# Progeny analysis for pathway activities - Violin plots of activity scores
#-------------------------------------------------------------------------------

if(identical(rownames(pr_res), colData(rld)$SampleID))
{
  pr_res_dat = cbind(pr_res, colData(rld)[,2:3])
  pr_res_dat = melt(as.data.frame(pr_res_dat))
  
  p1 = ggplot(pr_res_dat, aes(x = Time, y = value, fill = Type)) + 
    geom_boxplot(lwd = 0.1, outlier.size = 0.5) +
    labs(y = "Signalling score", x = "", fill = "") +
    scale_fill_manual(values = c("#0570b0", "#cb181d")) +
    #geom_boxplot(width=0.1, outlier.shape = NA) +
    theme_classic(base_size = 8) + facet_wrap(~variable , scales = "free", ncol = 4) + 
    theme(strip.background = element_rect(size = 0.1), 
          axis.ticks = element_line(size = 0.1), 
          axis.line = element_line(size = 0.1), 
          legend.position = "top")
  
  ggsave(filename = paste0(path, "analysis/Exploratory/progeny_pathway_scores_violin_byTreatment.pdf"), plot = p1, width = 7.5, height = 6)
  
  p2 = ggplot(pr_res_dat, aes(x = Type, y = value, fill = Time)) + 
    geom_boxplot(lwd = 0.1, outlier.size = 0.5) +
    labs(y = "Signalling score - PROGENy", x = "", fill = "") +
    scale_fill_manual(values = c("#cccccc", "#969696", "#636363", "#252525")) +
    #geom_boxplot(width=0.1, outlier.shape = NA) +
    theme_classic(base_size = 8) + facet_wrap(~variable , scales = "free", ncol = 4) + 
    theme(strip.background = element_rect(size = 0.1), 
          axis.ticks = element_line(size = 0.1), 
          axis.line = element_line(size = 0.1), 
          legend.position = "top")
  
  ggsave(filename = paste0(path, "analysis/Exploratory/progeny_pathway_scores_violin_byTime.pdf"), plot = p2, width = 7.5, height = 6)
}

rm(p1, p2, pr_res_dat)

#-------------------------------------------------------------------------------
# Secondary validation of Hypoxia scores using related GO terms
#-------------------------------------------------------------------------------

# Find Hypoxia related GO terms

# xx <- as.list(GOTERM)
# selGO = sapply(xx, function(x) ifelse(length(grep("hypoxia", x@Definition, ignore.case = T)) > 0, TRUE, FALSE))
# hypGO = data.frame(GOID = names(xx[selGO]), Term = Term(names(xx[selGO])))
# tmp = lapply(hypGO$GOID, function(x) unique(unlist(mget(get(x, org.Hs.egGO2ALLEGS), org.Hs.egSYMBOL))))
# names(tmp) = hypGO$Term
# tmp = tmp[sapply(tmp, length) > 10]
# venn(tmp)

# These analysis showed that "response to hypoxia | GO:0001666" was the only meaningful Hypoxia related GO term
hypGO = unique(unlist(mget(get("GO:0001666", org.Hs.egGO2ALLEGS), org.Hs.egSYMBOL)))

length(hypGO)
#[1] 359

# Overlap between Progeny Hypoxia genes and "response to hypoxia | GO:0001666" terms
hyp.overlap = progeny::model_human_full
hyp.overlap = as.character(hyp.overlap$gene[hyp.overlap$pathway %in% "Hypoxia" & hyp.overlap$p.value < 0.01])
hyp.overlap = list(`GO:0001666 | Response to Hypoxia` = hypGO, `Progeny | Hypoxia genes (p < 0.01)` = hyp.overlap)

p1 = ggvenn(hyp.overlap, fill_color = c("#0073C2FF", "#EFC000FF"), 
            stroke_size = 0.3, set_name_size = 2, text_size = 2.5, show_percentage = F)

# Pathway activities for Hypoxia computed using ssGSEA based on genes from GO:0001666
hypoGO_res_gse = gsva(expr = as.matrix(expr), gset.idx.list = list(hypGO), method = "gsva")[1,]

col_anno <- data.frame(do.call("rbind", strsplit(names(hypoGO_res_gse), "_")), row.names = names(hypoGO_res_gse))
colnames(col_anno) = c("Treatment", "Time", "Replicate")
col_anno$Treatment = ifelse(col_anno$Treatment == "H", "Hypoxic", "Normoxic")
col_anno$Time = factor(col_anno$Time, levels = c("6h", "12h", "D1", "D2"))

hypoGO_res_gse = data.frame(Hypoxia = hypoGO_res_gse, col_anno)

p2 = ggplot(hypoGO_res_gse, aes(x = Time, y = Hypoxia, fill = Treatment)) + 
  geom_boxplot(lwd = 0.1, outlier.size = 0.5) +
  labs(y = "Signalling score - GSVA", x = "", fill = "") +
  scale_fill_manual(values = c("#0570b0", "#cb181d")) +
  #geom_boxplot(width=0.1, outlier.shape = NA) +
  theme_classic(base_size = 8) + 
  theme(strip.background = element_rect(size = 0.1), 
        axis.ticks = element_line(size = 0.1), 
        axis.line = element_line(size = 0.1), 
        legend.position = "top")

p = p1 + p2 + plot_layout(widths = c(1, 0.6))
ggsave(filename = paste0(path, "analysis/Exploratory/secondary_validation_hypoxia_scores.pdf"), plot = p, width = 4, height = 2.5)

rm(hypGO, hyp.overlap, hypoGO_res_gse, col_anno, p1, p2, p)

#-------------------------------------------------------------------------------
# TF activity using DoroTHEA
#-------------------------------------------------------------------------------

# Extracting literature curatded TF regulons from Dorothea
data(dorothea_hs, package = "dorothea")

regulons = dorothea_hs %>%
  filter(confidence %in% c("A", "B", "C"))

# Computing TF activities using viver
tf_activities <- run_viper(expr, regulons, 
                           options =  list(method = "scale", minsize = 25, 
                                           eset.filter = FALSE, cores = 1, 
                                           verbose = FALSE, nes=T))
nrow(tf_activities)
#[1] 181

write.csv(data.frame(Samples = rownames(tf_activities), tf_activities, stringsAsFactors = F),
          paste0(path, "analysis/Exploratory/tfactivity_all_results.csv"), 
          quote = FALSE, row.names = F)

# Selecting the top 25% most variable TFs
v = apply(tf_activities, 1, mad)
tf_activities = t(tf_activities[v > quantile(v, 0.75),])

# Summarization
tf_activities.summ = split(data.frame(tf_activities), 
                     sapply(strsplit(rownames(tf_activities), "_"), function(x) paste0(x[1:2], collapse ="_")))

tf_activities.summ = t(sapply(tf_activities.summ, function(x){
  apply(x, 2, median)
}))

tf_activities.summ = tf_activities.summ[c("N_6h", "N_12h", "N_D1", "N_D2",
                                          "H_6h", "H_12h", "H_D1", "H_D2"),]

ncol(tf_activities.summ)
#[1] 45

paletteLength = 100
cols = colorRampPalette(rev(brewer.pal(9, "PiYG")))(paletteLength)
brk = c(seq(min(tf_activities.summ), 0, length.out=ceiling(paletteLength/2) + 1),
        seq(max(tf_activities.summ)/paletteLength,  max(tf_activities.summ), length.out=floor(paletteLength/2)))

col_anno <- data.frame(do.call("rbind", strsplit(rownames(tf_activities.summ), "_")), row.names = rownames(tf_activities.summ))
colnames(col_anno) = c("Treatment", "Time")
col_anno$Treatment = ifelse(col_anno$Treatment == "H", "Hypoxic", "Normoxic")

# Sanity check
identical(rownames(col_anno), rownames(tf_activities.summ))
#[1] TRUE

ann_colors = list(
  Treatment = c(Normoxic = "#cb181d", Hypoxic = "#0570b0"),
  Time = c(`6h` = "#cccccc", `12h` = "#969696", `D1` = "#636363", `D2` = "#252525")
)

pheatmap(tf_activities.summ, clustering_method = "ward.D2",
         cluster_rows  = F,
         cutree_cols = 5, 
         fontsize = 6,
         annotation_row = col_anno,
         annotation_colors = ann_colors,
         border_color = "white", 
         col = cols,
         breaks = brk,
         treeheight_row = 15, treeheight_col = 15,
         filename = paste0(path, "analysis/Exploratory/TFactivity_Dorothea_Viper.pdf"), width = 7.8, height = 2.5)

rm(cols, brk, paletteLength, col_anno, ann_colors, tf_activities.summ, v, dorothea_hs, regulons)

#-------------------------------------------------------------------------------
# Load signatures for stemness
#-------------------------------------------------------------------------------

stem1 = data.frame(read_excel(paste0(path, "data/extdata/Stemness_Signature_PMID_29625051.xlsx"), sheet = 2), stringsAsFactors = F)
stem1 = stem1[stem1$HUGO %in% rownames(expr),][,2:3]
stem1$Weight = as.numeric(stem1$Weight)

stem2 = data.frame(read_excel(paste0(path, "data/extdata/Stemness_Signature_PMID_30996127.xlsx"), sheet = 1), stringsAsFactors = F)
stem2 = as.character(na.omit(stem2[,3]))

#-------------------------------------------------------------------------------
# Compute stemness score from signature 1 and 2 - PMID 29625051, PMID 30996127
#-------------------------------------------------------------------------------

stem1_rho  = apply(expr[stem1$HUGO,], 2, function(x) cor(x, stem1$Weight, method = "spearman"))
stem2_gsva = gsva(expr = as.matrix(expr), gset.idx.list = list(stem2), method = "gsva")[1,]

stem = data.frame(Stemness_Malta_etal = stem1_rho, Stemness_Miranda_etal = stem2_gsva)

rm(stem1, stem2, stem1_rho, stem2_gsva)

#-------------------------------------------------------------------------------
# Correlation between Hypoxia and Stemness
#-------------------------------------------------------------------------------

if(identical(rownames(stem), rownames(pr_res)) & 
   identical(rownames(stem), colnames(expr))){
  
  hyp_stem = data.frame(stem, Hypoxia = pr_res[,"Hypoxia"])
  
  col_anno <- data.frame(do.call("rbind", strsplit(rownames(hyp_stem), "_")), row.names = rownames(hyp_stem))
  colnames(col_anno) = c("Treatment", "Time", "Replicate")
  col_anno$Treatment = ifelse(col_anno$Treatment == "H", "Hypoxic", "Normoxic")
  col_anno$Time = factor(col_anno$Time, levels = c("6h", "12h", "D1", "D2"))
  
  hyp_stem = cbind(hyp_stem, col_anno)
}

p1 = ggscatter(hyp_stem, x = "Hypoxia", y = "Stemness_Malta_etal",
                 color = "black", size = 1, # Points color, shape and size
                 add = "reg.line",  # Add regressin line
                 add.params = list(color = "firebrick"), # Customize reg. line
                 conf.int = TRUE, # Add confidence interval
                 cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                 cor.coeff.args = list(method = "spearman", label.x = 0.5, label.y = 0, label.sep = "\n"))
p1 = p1 + theme_classic(base_size = 9)

p2 = ggscatter(hyp_stem, x = "Hypoxia", y = "Stemness_Miranda_etal",
               color = "black", size = 1, # Points color, shape and size
               add = "reg.line",  # Add regressin line
               add.params = list(color = "firebrick"), # Customize reg. line
               conf.int = TRUE, # Add confidence interval
               cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
               cor.coeff.args = list(method = "spearman", label.x = 0.5, label.y = 0.4, label.sep = "\n"))
p2 = p2 + theme_classic(base_size = 9)
p = p1 + p2
ggsave(filename = paste0(path, "analysis/Exploratory/hypoxia_stemness_correlation.pdf"), plot = p, width = 6, height = 2.5)

rm(hyp_stem, stem, col_anno, p1, p2, p)
