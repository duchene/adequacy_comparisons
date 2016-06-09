#!/bin/bash
#PBS -N runname
#PBS -l ncpus=1
#PBS -l walltime=2:00:00
#PBS -l pmem=4gb
#PBS -P RDS-FSC-Phylogenomics-RW
#PBS -q compute
cd $PBS_O_WORKDIR
module load R
module load phyml
Rscript run.gtr.phyml.Rscript
