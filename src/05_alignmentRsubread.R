options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(TRUE)

fq     = as.character(args[1])
index  = as.character(args[2])
outdir = as.character(args[3])

library(Rsubread)

align.stat = align(index = paste0(index,"GRCh38index"),
                   readfile1 = fq,
                   type="rna",
                   phredOffset = 33,
                   nthreads = 15,
                   output_file = paste0(outdir, gsub(".fastq.gz", "", basename(fq)), ".bam")
                  )

# Default values used
# args(align)
# function (index, readfile1, readfile2 = NULL, type = "rna", input_format = "gzFASTQ", 
#           output_format = "BAM", output_file = paste(readfile1, "subread", output_format, sep = "."), 
#           isBCLinput = FALSE, phredOffset = 33, 
#           nsubreads = 10, TH1 = 3, TH2 = 1, maxMismatches = 3, unique = FALSE, 
#           nBestLocations = 1, indels = 5, complexIndels = FALSE, nTrim5 = 0, 
#           nTrim3 = 0, minFragLength = 50, maxFragLength = 600, PE_orientation = "fr", 
#           nthreads = 1, readGroupID = NULL, readGroup = NULL, keepReadOrder = FALSE, 
#           sortReadsByCoordinates = FALSE, color2base = FALSE, DP_GapOpenPenalty = -1, 
#           DP_GapExtPenalty = 0, DP_MismatchPenalty = 0, DP_MatchScore = 2, 
#           detectSV = FALSE, useAnnotation = FALSE, annot.inbuilt = "mm10", 
#           annot.ext = NULL, isGTF = FALSE, GTF.featureType = "exon", 
#           GTF.attrType = "gene_id", chrAliases = NULL) 

#saveRDS(align.stat, file = paste0(outdir, "alignStats/RDS/", gsub(".fastq.gz", "", basename(fq)), ".alignstats.RDS"))
