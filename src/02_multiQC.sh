# Activate environment
source activate hypoxia

# Get arguments from main script
export OUTDIR=${1}
export FQDIR=${2}
export FCLEAN=${3}

# Set locale values
export LANGUAGE=en_US.UTF-8
export LC_ALL=en_US.UTF-8

# Rename fastQC file names to make them compatible with FeatureCounts
#Rscript $FCLEAN

# MultiQC
multiqc --title 'Bulk_Hypoxia_Data_BoulantLab' --outdir $OUTDIR $FQDIR
