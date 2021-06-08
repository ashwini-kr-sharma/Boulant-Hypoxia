options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(TRUE)

path = as.character(args[1])
print(path)

library(Rsubread)

print("Starting to build index ...")

buildindex(basename = paste0(path, "/index/GRCh38index"), reference = paste0(path, "/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"))

print("Index building done ...")
