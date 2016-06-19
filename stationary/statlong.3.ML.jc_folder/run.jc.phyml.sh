#!/bin/bash
#PBS -N runname
#PBS -l ncpus=1
#PBS -l walltime=5:00:00
#PBS -l mem=4GB
#PBS -P RDS-FSC-Phylogenomics-RW
#PBS -q compute
cd $PBS_O_WORKDIR
module load R
module load phyml
Rscript run.jc.phyml.Rscript
