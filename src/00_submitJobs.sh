#!/bin/bash

#-------------------------------------------------------------------------------

# Scripts for generating counts from the single-end bulk hypoxia data from Boulant Lab
# Analysis done by - Ashwini Kumar Sharma, PhD
# Last update - 9th Dec, 2020

# Starting folder structure requirement inside the root directory
# NOTE: the names should be exact as shown below and dont change any script names etc

# <ROOTDIR>/data/fq/   > containing all the fastq files (can be downloaded from GEO)
# <ROOTDIR>/data/anno/ > containing the annotation information (can be downloaded from GEO/Supplementary data)
# <ROOTDIR>/src/       > containing all the analysis scripts (can be downloaded from GitHub)

# We assume that the user has conda installed and has created a conda environment -
# conda create --name hypoxia -c bioconda r fastqc multiqc bioconductor-rsubread r-base
# A detailed environment .yaml file is also provided in GitHub

# This analysis was run on the LSF - IBM HPC system, if using a different HPC system,
# like qsub, SLURM etc the user has to modify the cluster submission code chunk accordingly

#-------------------------------------------------------------------------------

ROOTDIR='/icgc/dkfzlsdf/analysis/B080/sharma/boulantLab/bulk_hypoxia/'

#-------------------------------------------------------------------------------
# FastQC
#-------------------------------------------------------------------------------

# SCRIPT=$ROOTDIR$'src/01_fastQC.sh'
# mkdir -p $ROOTDIR$'analysis/FastQC'
# mkdir -p $ROOTDIR'logs/FastQC'
#
# for i in $ROOTDIR$'data/fq/'*.fastq.gz ; do
#  bsub -n 1  -R "rusage[mem=2G]" -W 2:00 -J hyp_fastqc \
#  -o $ROOTDIR'logs/FastQC/'$i'.out' -e $ROOTDIR'logs/FastQC/'$i'.err' \
#  "$SCRIPT $ROOTDIR$'analysis/FastQC' $i"
# done

#-------------------------------------------------------------------------------
# Download genomes
#-------------------------------------------------------------------------------

# SCRIPT=$ROOTDIR$'src/03_downloadGenomes.sh'
# mkdir -p $ROOTDIR$'data/hg38'
# mkdir -p $ROOTDIR'logs/DownloadGenome'
#
# type=('genomeFasta' 'GTF')
# for files in "${type[@]}"; do
#   bsub -n 1  -R "rusage[mem=5G]" -W 5:00 -J hyp_downloadGenome \
#   -o $ROOTDIR'logs/DownloadGenome/'$files'.wget.out' -e $ROOTDIR'logs/DownloadGenome/'$files'.wget.err' \
#   "$SCRIPT $ROOTDIR$'data/hg38' $files"
# done

#-------------------------------------------------------------------------------
# Rsubread indexing
#-------------------------------------------------------------------------------

# mkdir -p $ROOTDIR'logs/Indexing'
#
# echo "source activate hypoxia; Rscript $ROOTDIR$'src/hypoxia/04_indexRsubread.R' $ROOTDIR$'data/hg38'" | \
# bsub -n 1 -R "rusage[mem=30G]" -W 10:00 -J hyp_index \
# -o $ROOTDIR'logs/Indexing/index.out' -e $ROOTDIR'logs/Indexing/index.err'

#-------------------------------------------------------------------------------
# Rsubread alignment
#-------------------------------------------------------------------------------

# mkdir -p $ROOTDIR'data/bam/alignStats/RDS'
# mkdir -p $ROOTDIR'logs/Align'
#
# for i in $ROOTDIR$'data/fq/'*.fastq.gz ; do
#
#   bn=$(basename "${i%.*.*}")
#
#   echo "source activate hypoxia; Rscript $ROOTDIR$'src/05_alignmentRsubread.R' $i $ROOTDIR$'data/hg38/index/' $ROOTDIR$'data/bam/'" | \
#   bsub -n 15 -R "rusage[mem=15G]" -W 5:00 -J hyp_align \
#   -o $ROOTDIR'logs/Align/'$bn$'.out' -e $ROOTDIR'logs/Align/'$bn$'.err'
# done

# Rerunning Failed cases

#failed=('H3CGFBGXF_D2_H_4_20s000271-1-1_Triana_lane120s000271.fastq.gz' 'H3CGFBGXF_D2_N_3_20s000265-1-1_Triana_lane120s000265.fastq.gz' 'H3CGFBGXF_D2_N_4_20s000266-1-1_Triana_lane120s000266.fastq.gz' 'H3HHMBGXF_12h_N_2_20s000252-1-1_Triana_lane120s000252.fastq.gz' 'H3HHMBGXF_D1_N_2_20s000254-1-1_Triana_lane120s000254.fastq.gz' 'H3HHMBGXF_D2_H_2_20s000261-1-1_Triana_lane120s000261.fastq.gz' 'H3HHMBGXF_D2_N_1_20s000255-1-1_Triana_lane120s000255.fastq.gz' 'H3HHMBGXF_D2_N_2_20s000256-1-1_Triana_lane120s000256.fastq.gz' 'HH7WMBGXF_6h_H_1_20s001000-1-1_Triana_lane120s001000.fastq.gz' 'HH7WMBGXF_6h_H_2_20s001001-1-1_Triana_lane120s001001.fastq.gz' 'HH7WMBGXF_6h_H_3_20s001002-1-1_Triana_lane120s001002.fastq.gz' 'HH7WMBGXF_6h_H_4_20s001003-1-1_Triana_lane120s001003.fastq.gz')

# for i in "${failed[@]}"; do
#   i=$ROOTDIR$'data/fq/'$i
#
#   bn=$(basename "${i%.*.*}")
#
#   echo "source activate hypoxia; Rscript $ROOTDIR$'src/05_alignmentRsubread.R' $i $ROOTDIR$'data/hg38/index/' $ROOTDIR$'data/bam/'" | \
#   bsub -n 15 -R "rusage[mem=30G]" -W 5:00 -J hyp_align \
#   -o $ROOTDIR'logs/Align/'$bn$'.out' -e $ROOTDIR'logs/Align/'$bn$'.err'
# done

#-------------------------------------------------------------------------------
# Rsubread featurecounts
#-------------------------------------------------------------------------------

# mkdir -p $ROOTDIR'data/counts/'
# mkdir -p $ROOTDIR'logs/Count'
#
# echo "source activate hypoxia; Rscript $ROOTDIR$'src/06_countsRsubread.R' $ROOTDIR$'data/bam' $ROOTDIR$'data/hg38/Homo_sapiens.GRCh38.102.gtf.gz' $ROOTDIR'data/counts/'" | \
# bsub -n 15 -R "rusage[mem=5G]" -W 5:00 -J hyp_count \
# -o $ROOTDIR'logs/Count/featureCount.out' -e $ROOTDIR'logs/Count/featureCount.err'

#-------------------------------------------------------------------------------
# MultiQC
#-------------------------------------------------------------------------------

SCRIPT=$ROOTDIR$'src/02_multiQC.sh'
mkdir -p $ROOTDIR$'analysis/MultiQC'
mkdir -p $ROOTDIR'logs/MultiQC'

bsub -n 1  -R "rusage[mem=5G]" -W 2:00 -J hyp_multiqc \
-o $ROOTDIR'logs/MultiQC/MultiQC.out' -e $ROOTDIR'logs/MultiQC/MultiQC.err' \
"$SCRIPT $ROOTDIR$'analysis/MultiQC' $ROOTDIR $ROOTDIR'src/miscel/renameFiles.R'"
