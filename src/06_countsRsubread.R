options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(TRUE)

bamdir  = as.character(args[1])
gtffile = as.character(args[2])
outdir  = as.character(args[3])

library(Rsubread)

counts = featureCounts(files = list.files(bamdir, full.names = T, pattern=".*bam$"),
                       annot.ext = gtffile,
                       isGTFAnnotationFile = T,
                       GTF.attrType.extra = c("gene_name", "gene_biotype"),
                       nthreads = 15)

saveRDS(counts, file = paste0(outdir, "hypoxia_counts.RDS"))

