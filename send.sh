#!/bin/bash

#SBATCH --job-name=oscilake
#SBATCH --output=%x-%j.out
#SBATCH --mem=100G
#SBATCH --time=05-00:00:00


###SBATCH --partition=skylake

###SBATCH --nodes=1


echo EMPIEZA TODO `date`

ml R-bundle-CRAN
ml CDO
ml NCO
ml GDAL

Rscript pca_analysis_global.R

echo TERMINA TODO `date`

