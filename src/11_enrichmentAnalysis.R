#-------------------------------------------------------------------------------
# Required packages, setting root directory and reading the count data
#-------------------------------------------------------------------------------

library(tidyverse)
library(fgsea)
library(msigdbr)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(enrichR)

path = "/icgc/dkfzlsdf/analysis/B080/sharma/boulantLab/bulk_hypoxia/"
resFC = readRDS(paste0(path, "analysis/DiffExp/diffExpLogFCmatrix.RDS"))
dge = readRDS(paste0(path, "analysis/DiffExp/diffExpGenes.RDS"))
  
#-------------------------------------------------------------------------------
# Enrichment analysis - Hallmarks and GO-BP - GSEA
#-------------------------------------------------------------------------------

# Enrichment analysis function -------------------------------------------------

performEnrichment = function(gene_sets, type)
{
  gene_sets = split(x = gene_sets$gene_symbol, f = gene_sets$gs_name)
  
  if(type == "Hallmarks"){
    pval_cutoff = 0.01
    w = 2.6
    h = 3
  }else if(type == "GOBP"){
    pval_cutoff = 0.001
    w = 5.8
    h = 35
  }
  
  # Enrichment analysis
  fgseaRes <- lapply(resFC_list, function(x) {
    enr = fgsea(pathways = gene_sets, 
                stats = x,
                eps = 0,
                minSize = 15,
                maxSize = 500)
    enr = enr[enr$padj < pval_cutoff,]
    enr = enr[order(enr$padj),][,c(1,2,3,6,8)]
    return(enr)
  })
  
  # Enrichment heatmap
  gseaHM = melt(lapply(fgseaRes, function(x) x[,c("pathway", "NES")]))
  gseaHM = gseaHM[,c(1,4,3)]
  gseaHM = dcast(data = gseaHM, formula = pathway ~ L1)
  
  gseaHM = gseaHM[apply(gseaHM, 1, function(x) sum(is.na(x))) < 11,]
  
  gseaHM = data.frame(gseaHM, row.names = 1)
  gseaHM[is.na(gseaHM)] = 0
  
  # Arranging the columns and rows as we want
  gseaHM = gseaHM[,c("T_6hr", "T_12hr", "T_D1", "T_D2")]
  gseaHM = gseaHM[order(rowSums(gseaHM), decreasing = T),]
  
  # Heatmap design 
  myColor = colorRampPalette(c("Darkblue", "white","red"))(100)
  br = c(seq(min(gseaHM), 0, length.out=ceiling(100/2) + 1), 
         seq(max(gseaHM)/100, max(gseaHM), length.out=floor(100/2)))
  
  if(type == "Hallmarks"){
    rownames(gseaHM) = gsub("HALLMARK_", "", rownames(gseaHM))
  }
  
  if(type == "GOBP"){
    rownames(gseaHM) = gsub("GO_", "", rownames(gseaHM))
  }
  
  pheatmap(gseaHM, breaks = br,
           cluster_cols = F,
           cluster_rows = F,
           fontsize = 6,
           labels_row = gsub("_", " ", rownames(gseaHM), fixed = T),
           border_color = "white", 
           col = myColor,
           filename = paste0(path, "analysis/DiffExp/enriched_pathways_", type, ".pdf"), width = w, height = h)
  
  return(NULL)
}

# ------------------------------------------------------------------------------

# Ranked list for GSEA
resFC_list = apply(resFC, 2, function(x) {
  x = sort(setNames(x, rownames(resFC)))
  return(list(x))
})
resFC_list = lapply(resFC_list, unlist)

# Gene sets
hl_gene_sets = msigdbr(species = "Homo sapiens", category = "H")
go_gene_sets = msigdbr(species = "Homo sapiens", subcategory = "GO:BP")

performEnrichment(gene_sets = hl_gene_sets, type ="Hallmarks")
performEnrichment(gene_sets = go_gene_sets, type ="GOBP")
