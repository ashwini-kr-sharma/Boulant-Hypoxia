#-------------------------------------------------------------------------------
# Required packages, setting root directory and reading the count data
#-------------------------------------------------------------------------------

library(DESeq2)
library(ggvenn)
library(tidyverse)
library(patchwork)
library(reshape2)
library(enrichR)
library(pheatmap)
library(RColorBrewer)

path = "/icgc/dkfzlsdf/analysis/B080/sharma/boulantLab/bulk_hypoxia/"
dat = readRDS(paste0(path, "analysis/DiffExp/diffExpGenes.RDS"))
resFC = readRDS(paste0(path, "analysis/DiffExp/diffExpLogFCmatrix.RDS"))

#-------------------------------------------------------------------------------
# Create directories and download two independent signatures for stemness
#-------------------------------------------------------------------------------

# Create output directory 
if( ! dir.exists(paste0(path, "analysis/DiffExp"))){
  dir.create(paste0(path, "analysis/DiffExp"),recursive = T)
}

#-------------------------------------------------------------------------------
# Top up-regulated or down-regulated genes
#-------------------------------------------------------------------------------

topgenes = lapply(dat, function(x){
  up = x$Gene[x$log2FoldChange > 1 & x$padj < 0.01]
  up = up[!is.na(up)]
  
  dw = x$Gene[x$log2FoldChange < -1 & x$padj < 0.01]
  dw = dw[!is.na(dw)]
  
  return(list(up = up, down = dw))
})

up = lapply(topgenes, function(x) x$up)
dw = lapply(topgenes, function(x) x$down)

rm(topgenes)

#-------------------------------------------------------------------------------
# Venn diagrams
#-------------------------------------------------------------------------------

# UP
p_up = ggvenn(up, fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"), 
              stroke_size = 0.3, set_name_size = 2, text_size = 2, show_percentage = F) + 
       labs(subtitle = "Upregulated genes")

# DOWN
p_dw = ggvenn(dw, fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"), 
              stroke_size = 0.3, set_name_size = 2, text_size = 2, show_percentage = F) + 
       labs(subtitle = "Downregulated genes")

p1 = p_up + p_dw

#-------------------------------------------------------------------------------
# Enrichment of common core genes
#-------------------------------------------------------------------------------

dge_up = Reduce(intersect, up)
dge_dw = Reduce(intersect, dw)

setEnrichrSite("Enrichr")
dbs <- listEnrichrDbs()
dbs = c("TRANSFAC_and_JASPAR_PWMs",
        "TRRUST_Transcription_Factors_2019", 
        "MSigDB_Hallmark_2020")

enrichDGE = function(dat, dbs)
{
  res <- melt(lapply(enrichr(dat, dbs), function(x){
    x = x[,c(1,4,7)]
    x = x[x$Adjusted.P.value < 0.001,]
  }), id = c("Term", "Adjusted.P.value", "Odds.Ratio"))
  res = res[grep("mouse", res$Term, invert = T),]
  
  res$Adjusted.P.value = -log10(res$Adjusted.P.value)
  res$Odds.Ratio = log10(res$Odds.Ratio)
  res$Adjusted.P.value[res$Adjusted.P.value > 10] = 10
  
  p.enr <- ggplot(data = res, aes(x = Adjusted.P.value, y = reorder(Term, Adjusted.P.value))) + 
            theme_bw(base_size = 8) + labs(y = "") + 
            geom_bar(stat="identity") + xlim(0,10) +
            facet_wrap(~ L1, scales ="free", ncol = 1) + geom_vline(xintercept = -log10(0.01)) +
            theme(panel.grid = element_blank(), axis.line = element_blank(), strip.background = element_blank())
  
  return(p.enr)
  
}

p2 = enrichDGE(dat = dge_up, dbs = dbs)
p3 = enrichDGE(dat = dge_dw, dbs = dbs)

#-------------------------------------------------------------------------------
# Plot heatmap of the top 50 genes from each condition
#-------------------------------------------------------------------------------

topgenes = lapply(dat, function(x){

  x = as.data.frame(x)
  x = na.omit(x)

  up = x[x$log2FoldChange > 1 & x$padj < 0.01,]
  up = up[order(up$log2FoldChange, decreasing = T),]
  up = head(up, 50)
  up = up$Gene

  dw = x[x$log2FoldChange < -1 & x$padj < 0.01,]
  dw = dw[order(dw$log2FoldChange, decreasing = F),]
  dw = head(dw, 50)
  dw = dw$Gene
  
  selgenes = c(up, dw)

  return(selgenes)
})

tmpdat = resFC[rownames(resFC) %in% unique(unlist(topgenes)),]

dim(tmpdat)
#[1] 249   4

# Heatmap design 
myColor = colorRampPalette(c("red", "white","darkblue"))(100)
br = c(seq(min(tmpdat), 0, length.out=ceiling(100/2) + 1), 
       seq(max(tmpdat)/100, max(tmpdat), length.out=floor(100/2)))

p4 = pheatmap(tmpdat, cluster_cols = F, treeheight_row = 20, show_rownames = F, 
         breaks = br, color = myColor, cutree_rows = 1)

p = p1 / (p2 | (p3 / p4[[4]])) 
ggsave(filename = paste0(path, "analysis/DiffExp/differentiallyExpressedGenes.pdf"), plot = p, width = 7, height = 10)

#rm(up, dw, p_dw, p_up, p1, p2, p3, p4, dge_up, dge_dw, p, dbs, enrichDGE, myColor, br, tmpdat)
