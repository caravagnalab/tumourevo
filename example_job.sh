#!/bin/sh
#### Cluster specific arguments ####

#SBATCH --job-name=NF_MASTER
#SBATCH --error=err
#SBATCH --partition=THIN
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=12:00:00
#SBATCH -A cdslab
#SBATCH --output=log

module load singularity
module load java

tumourevo=''
RESULTS_PATH=$tumourevo/results/


/orfeo/LTS/LADE/LT_storage/lvaleriani/nextflow/nextflow run $tumourevo/main.nf \
    -profile singularity \
    --input $tumourevo/test_input.csv \
    --outdir $RESULTS_PATH \
    -c $tumourevo/orfeo_cdslab.config
