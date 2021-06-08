#-------------------------------------------------------------------------------
# Required packages, setting root directory and reading the count data
#-------------------------------------------------------------------------------

library(DESeq2)
library(WriteXLS)
library(ashr)

path = "/icgc/dkfzlsdf/analysis/B080/sharma/boulantLab/bulk_hypoxia/"
dat = readRDS(paste0(path, "data/counts/hypoxia_filtered_counts.RDS"))

#-------------------------------------------------------------------------------
# Create directories and download two independent signatures for stemness
#-------------------------------------------------------------------------------

# Create output directory 
if( ! dir.exists(paste0(path, "analysis/DiffExp"))){
  dir.create(paste0(path, "analysis/DiffExp"),recursive = T)
}

#-------------------------------------------------------------------------------
# DESeq2 object - All conditions logFC
#-------------------------------------------------------------------------------

ddsFC = DESeqDataSetFromMatrix(countData = dat$counts,
                               colData = dat$sampanno,
                               rowData = dat$geneanno,
                               design = ~ Condition)
ddsFC = DESeq(ddsFC)

# Function to compute all possible contrasts
diffexp = function(cond1, cond2){
  dat = results(ddsFC, contrast = c("Condition", cond1, cond2), alpha = 0.01)
  dat = lfcShrink(ddsFC, res = dat, contrast = c("Condition", cond1, cond2), type="ashr")
}

# Hypoxia vs Normoxia
res_6h  = diffexp(cond1 = "H_6h",  cond2 = "N_6h")
res_12h = diffexp(cond1 = "H_12h", cond2 = "N_12h")
res_D1  = diffexp(cond1 = "H_D1",  cond2 = "N_D1")
res_D2  = diffexp(cond1 = "H_D2",  cond2 = "N_D2")

# Sanity checks !!
all(identical(rownames(res_6h), rownames(res_12h)),
    identical(rownames(res_D1), rownames(res_D2)),
    identical(rownames(res_D1), rownames(res_6h)),
    identical(rownames(res_D1), rowData(ddsFC)$GeneID))

# Add gene name
res_6h$Gene  = rowData(ddsFC)$gene_name
res_12h$Gene = rowData(ddsFC)$gene_name
res_D1$Gene  = rowData(ddsFC)$gene_name
res_D2$Gene  = rowData(ddsFC)$gene_name

dge = list(T_6hr = res_6h, T_12hr = res_12h, T_D1 = res_D1, T_D2 = res_D2)

saveRDS(dge, paste0(path, "analysis/DiffExp/diffExpGenes.RDS"))

dge = lapply(dge, function(x){
  x = data.frame(x)
  x[,c(1:3)] = round(x[,c(1:3)], 2)
  x$pvalue = signif(x$pvalue, 2)
  x$padj = signif(x$padj, 2)
  return (x)
})

WriteXLS(lapply(dge, function(x){data.frame(Gene = rownames(x), x)}), 
         ExcelFileName = paste0(path, "analysis/DiffExp/DGEtables.xls"),
         AdjWidth = TRUE, BoldHeaderRow = TRUE, FreezeRow = 1, 
         SheetNames = gsub("res_", "", names(dge)))

rm(dge)

# Consolidated FC 
resFC = data.frame(T_6hr  = res_6h$log2FoldChange,
                   T_12hr = res_12h$log2FoldChange,
                   T_D1   = res_D1$log2FoldChange,
                   T_D2   = res_D2$log2FoldChange,
                   row.names = paste(rowData(ddsFC)$GeneID, rowData(ddsFC)$gene_name, sep="|"))

resFC = resFC[!duplicated(rowData(ddsFC)$gene_name),]
rownames(resFC) = sapply(strsplit(rownames(resFC), "|", fixed=T), function(x)x[2])

saveRDS(resFC, paste0(path, "analysis/DiffExp/diffExpLogFCmatrix.RDS"))

#-------------------------------------------------------------------------------
# DESeq2 object - all Time course
#-------------------------------------------------------------------------------

# ddsTC = DESeqDataSetFromMatrix(countData = dat$counts,
#                                colData   = dat$sampanno,
#                                rowData   = dat$geneanno,
#                                design    = ~ Type + Time + Type:Time)
# ddsTC = DESeq(ddsTC, test = "LRT", reduced = ~ Type + Time)
# 
# if(identical(rownames(resTC), rowData(ddsTC)$GeneID)){
#   resTC = results(ddsTC, alpha = 0.01)
#   resTC$Symbol = rowData(ddsTC)$gene_name
#   resTC = as.data.frame(resTC)
#   resTC = resTC[order(resTC$log2FoldChange, resTC$padj),]
# }
# 
# resTC[resTC$Symbol %in% hif, ]
# 
# sel = resTC[which(resTC$Symbol %in% hypGO &
#             resTC$log2FoldChange > 1 & 
#             resTC$padj < 0.01), ]
# 
# fiss <- plotCounts(ddsTC, gene = "ENSG00000112715", 
#                    intgroup = c("Time","Type"), returnData = TRUE)
# 
# ggplot(fiss,
#        aes(x = Time, y = count, color = Type, group = Type)) + 
#   geom_point() + stat_summary(fun=median, geom="line") +
#   scale_y_log10()
