library(GO.db)
library(org.Hs.eg.db)
library(msigdbr)
library(pheatmap)
library(RColorBrewer)

path = "/icgc/dkfzlsdf/analysis/B080/sharma/boulantLab/bulk_hypoxia/"
resFC = readRDS(paste0(path, "analysis/DiffExp/diffExpLogFCmatrix.RDS"))

#-------------------------------------------------------------------------------
# Selected gene modules to plot
#-------------------------------------------------------------------------------

# Response to hypoxia | GO:0001666
hypGO = unique(unlist(mget(get("GO:0001666", org.Hs.egGO2ALLEGS), org.Hs.egSYMBOL)))

genesets = msigdbr(species = "Homo sapiens", category = "H")
hl_hifsig = split(x = genesets$gene_symbol, f = genesets$gs_name)$HALLMARK_HYPOXIA

hif = union(hypGO, hl_hifsig); rm(hypGO, hl_hifsig)
emt = split(x = genesets$gene_symbol, f = genesets$gs_name)$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
wnt = split(x = genesets$gene_symbol, f = genesets$gs_name)$HALLMARK_WNT_BETA_CATENIN_SIGNALING

#-------------------------------------------------------------------------------
# Function to plot gene modules
#-------------------------------------------------------------------------------

plotGeneModule = function(dat, name, abs_logFC_cutoff, palleteSize, rowclustsize, plotwidth, plotheight)
{
  if(name == "hypoxia_gene"){
    dat = dat[rowSums(abs(dat) > abs_logFC_cutoff) > 0 | rownames(dat) == "HIF1A",]
  }else{
    dat = dat[rowSums(abs(dat) > abs_logFC_cutoff) > 0 ,]
  }
  
  # Coloring
  myColor = colorRampPalette(c("red", "white","darkblue"))(palleteSize)
  br = c(seq(min(dat), 0, length.out=ceiling(palleteSize/2) + 1), 
         seq(max(dat)/palleteSize, max(dat), length.out=floor(palleteSize/2)))
  
  pheatmap(dat, breaks = br,
           clustering_method = "ward.D2",
           cluster_cols = FALSE, 
           cutree_rows = rowclustsize,
           treeheight_row = 0,
           fontsize = 5,
           border_color = "white", 
           col = myColor,
           filename = paste0(path, "analysis/DiffExp/", name, ".pdf"), width = plotwidth, height = plotheight)
  
  return(NULL)
}

#-------------------------------------------------------------------------------
# Hypoxia genes : only those genes with -1 > logFC > 1 in atleast one condition
#-------------------------------------------------------------------------------

hifsig = resFC[rownames(resFC) %in% hif,]
plotGeneModule(dat = hifsig, name = "logFC_hypoxia_genes", abs_logFC_cutoff = 1, 
               palleteSize = 20, rowclustsize = 5, plotwidth = 1.3, plotheight = 25)

#-------------------------------------------------------------------------------
# WNT genes : only those genes with -1 > logFC > 1 in atleast one condition
#-------------------------------------------------------------------------------

wntsig = resFC[rownames(resFC) %in% wnt,]
plotGeneModule(dat = wntsig, name = "logFC_WNT_genes", abs_logFC_cutoff = 1,
               palleteSize = 20, rowclustsize = 2, plotwidth = 1.2, plotheight = 2.5)

#-------------------------------------------------------------------------------
# EMT genes : only those genes with -1 > logFC > 1 in atleast one condition
#-------------------------------------------------------------------------------

emtsig = resFC[rownames(resFC) %in% emt,]
plotGeneModule(dat = emtsig, name = "logFC_EMT_genes", abs_logFC_cutoff = 1,
               palleteSize = 20, rowclustsize = 3, plotwidth = 1.3, plotheight = 10)
