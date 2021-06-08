library(Seurat)
library(tidyverse)
library(data.table)
library(pheatmap)
library(RColorBrewer)

path = "/icgc/dkfzlsdf/analysis/B080/sharma/boulantLab/bulk_hypoxia/"

#-------------------------------------------------------------------------------
# Create directories
#-------------------------------------------------------------------------------

# Create output directory 
if( ! dir.exists(paste0(path, "analysis/Deconvolution"))){
  dir.create(paste0(path, "analysis/Deconvolution"),recursive = T)
}

#-------------------------------------------------------------------------------
# Hypoxia - CPM
#-------------------------------------------------------------------------------

hyp.cpm = readRDS(paste0(path, "data/counts/hypoxia_filtered_counts.RDS"))
anno = hyp.cpm$sampanno
genes = hyp.cpm$geneanno[,c("GeneID", "gene_name")]

hyp.cpm = data.frame(GeneID = rownames(hyp.cpm$counts), hyp.cpm$counts, stringsAsFactors = F)
hyp.cpm = merge(genes, hyp.cpm)
hyp.cpm = hyp.cpm[! duplicated(hyp.cpm$gene_name),]
hyp.cpm = data.frame(hyp.cpm[,3:ncol(hyp.cpm)], row.names = hyp.cpm$gene_name)
hyp.cpm = RelativeCounts(data = hyp.cpm,
                         scale.factor = 1e6)
hyp.cpm = round(as.matrix(hyp.cpm), 3)
hyp.cpm = hyp.cpm[rowSums(hyp.cpm) > 5,] #remove lowly expressed genes
rm(genes)

#-------------------------------------------------------------------------------
# 2D Ileum single cell organoid model - get CPM values for CIBERSORTx
#-------------------------------------------------------------------------------

# Data downloaded from https://doi.org/10.6084/m9.figshare.13703752.v1

load("/icgc/dkfzlsdf/analysis/B080/sharma/boulantLab/sc_interferon/data/COVID19_July.rda")

il2d = subset(Illeum_H_T, idents = "Illeum_Mock")
tmp = CreateSeuratObject(il2d@assays$RNA@counts)
if(identical(rownames(tmp@meta.data), rownames(il2d@meta.data))){
  tmp@meta.data = cbind(tmp@meta.data, il2d@meta.data[,c("Smillie", "Ileum_Tissue", "Ileum_org", "CellTypes")])
  tmp = PercentageFeatureSet(tmp, pattern = "^MT-", col.name = "percent.mt")
  tmp@meta.data$organoid_model = "2D"
  il2d = tmp
  rm(tmp)
}

il.cpm = RelativeCounts(data = GetAssayData(il2d, slot = "counts"),
                        scale.factor = 1e6)

if(identical(colnames(il2d), rownames(il2d@meta.data))){
  colnames(il.cpm) = gsub(" ", "_", il2d$CellTypes)
  il.cpm = round(as.matrix(il.cpm), 3)
}

il.cpm = il.cpm[rowSums(il.cpm) > 5,] #remove lowly expressed genes

# co2d = subset(Colon_H_T, idents = "Colon_Mock")
# tmp = CreateSeuratObject(co2d@assays$RNA@counts)
# if(identical(rownames(tmp@meta.data), rownames(co2d@meta.data))){
#   tmp@meta.data = cbind(tmp@meta.data, co2d@meta.data[,c("Smillie", "Ileum_Tissue", "Ileum_org", "CellTypes")])
#   tmp = PercentageFeatureSet(tmp, pattern = "^MT-", col.name = "percent.mt")
#   tmp@meta.data$organoid_model = "2D"
#   co2d = tmp
#   rm(tmp)
# }
# 
rm(Colon_H_T, Illeum_H_T, il2d)

#-------------------------------------------------------------------------------
# 2D Ileum single cell organoid model - get CPM values for CIBERSORTx
#-------------------------------------------------------------------------------

write.table(data.frame(Gene = rownames(hyp.cpm), hyp.cpm),  
            paste0(path, "analysis/Deconvolution/hypoxia_mixture.txt"),
            sep= "\t", quote = F, row.names = F)

write.table(data.frame(GeneSymbol = rownames(il.cpm), il.cpm),  
            paste0(path, "analysis/Deconvolution/ileum_reference.txt"),
            sep= "\t", quote = F, row.names = F)

rm(hyp.cpm, il.cpm)
#-------------------------------------------------------------------------------
#
#
# RUN CIBERSORTx ANALYSIS at https://cibersortx.stanford.edu/
# COPY RESULTS AND READ THEM BACK FOR VISUALIZATION AND ANALYSIS
#
#
#-------------------------------------------------------------------------------

# # Deconvolution results from CIBERSORTx:ABSOLUTE
# decov_abs = fread(paste0(path, "analysis/Deconvolution/CIBERSORTx_Absolute_Results.txt"), header=T, stringsAsFactors = F)
# decov_abs = t(data.frame(decov_abs, row.names = 1)[,1:9])
# colSums(decov_abs)

# Deconvolution results from CIBERSORTx:RELATIVE
decov_rel = fread(paste0(path, "analysis/Deconvolution/CIBERSORTx_Relative_Results.txt"), header=T, stringsAsFactors = F)
decov_rel = t(data.frame(decov_rel, row.names = 1)[,1:9])
colSums(decov_rel)

# identical(colnames(decov_rel), colnames(decov_abs))

decov_rel = reshape2::melt(decov_rel)

anno = data.frame(do.call("rbind", strsplit(as.character(decov_rel$Var2), "_", fixed=T)))[,1:2]
colnames(anno) = c("State", "Time")
anno$State = ifelse(anno$State == "H", "Hypoxia", "Normoxia")

decov_rel = cbind(decov_rel, anno)

# Since D2 has the highest hypoxic effect, only focussing on these samples
decov_rel = decov_rel[decov_rel$Time %in% c("D2", "12h"),]

decov_rel$Var1 = factor(as.character(decov_rel$Var1) , levels =  c("Stem_Cells", "Cycling_TA", "Inmature_Enterocyte_1", "Goblet_Cells", 
  "Secretory_TA", "Enteroendocrine_cells", "Inmature_Enterocyte_2", "Enterocyte_1", "TA"))

p1 = ggplot(decov_rel, aes(x = Time, y = value, fill = State)) + 
  geom_boxplot(lwd = 0.1, outlier.size = 0.5) +
  labs(y = "Cell type fraction", x = "", fill = "") +
  scale_fill_manual(values = c("#0570b0", "#cb181d")) +
  theme_classic(base_size = 8) + facet_wrap(~Var1 , scales = "free", ncol = 3) + 
  theme(strip.background = element_rect(size = 0.1), 
        axis.ticks = element_line(size = 0.1), 
        axis.line = element_line(size = 0.1), 
        legend.position = "top")

ggsave(filename = paste0(path, "analysis/Deconvolution/CIBERSORTx_relative_cell_fraction_analysis.pdf"), 
       plot = p1, width = 4, height = 6)
