#!/bin/bash

#SBATCH --job-name=oscilake
#SBATCH --output=%x-%j.out
#SBATCH --mem=100G
#SBATCH --partition=skylake
#SBATCH --time=05-00:00:00

###SBATCH --nodes=1


echo EMPIEZA TODO `date`

ml R-bundle-CRAN
ml CDO
ml NCO
ml GDAL

#rm clustering/*  
#rm figures/*
#rm pca/*
#rm pca_04/*
#rm pca_explained/* 
#rm permutation/*
#rm acf/*
#rm cluster_mdata/*
#rm pca_04_explained/*

Rscript pca_analysis_global.R

echo TERMINA TODO `date`

