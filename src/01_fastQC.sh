#!/bin/bash

# Activate environment
source activate hypoxia

# Get arguments from main script
export OUTDIR=${1}
export FQ=${2}

# Set locale values
export LANGUAGE=en_US.UTF-8
export LC_ALL=en_US.UTF-8

# FastQC
fastqc --outdir $OUTDIR $FQ
