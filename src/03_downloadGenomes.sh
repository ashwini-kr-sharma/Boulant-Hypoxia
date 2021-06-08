#!/bin/bash
# Get arguments from main script
export OUTDIR=${1}
export TYPE=${2}

# Download genomes
if [[ $TYPE == 'genomeFasta' ]]
then
  wget http://ftp.ensembl.org/pub/release-102/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -P $OUTDIR
fi

if [[ $TYPE == 'GTF' ]]
then
  wget http://ftp.ensembl.org/pub/release-102/gtf/homo_sapiens/Homo_sapiens.GRCh38.102.gtf.gz -P $OUTDIR
fi
