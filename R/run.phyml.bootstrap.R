library(ape)
library(phangorn)
library(doParallel)
library(foreach)


true_data <- as.DNAbin(simSeq(rtree(50)))
write.dna(true_data, file = 'test_true_data.fasta', format = 'fasta', nbcol = -1, colsep = '')



phyml_path <- '~/Desktop/phylo_programs/PhyML-3.1/PhyML-3.1_macOS-MountainLion'


