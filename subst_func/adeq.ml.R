# This is a 'wrapper' function. It takes all the information that the other functions require to perform substitution model adequacy in one step. The information required includes: The trees log file path, parameters log file path, empirical sequence data file path, sequence length, assumed tree topology, number of PP simulations to take.


adeq.ml <- function(trees.file, log.file, empdat.file, true.topo, Nsim = 100, savetrees = F){
     empdat <- as.phyDat(read.dna(empdat.file))
     seqlen <- ncol(as.matrix(as.DNAbin(empdat)))
     print(paste("sequence length is", seqlen))
     sims <- make.ml.als(trees.file, log.file, Nsim, seqlen, savetrees)
     if(length(grep("jc", trees.file, value = T)) > 0){ 
     	model <- "jc69"
     } else {
        model <- "gtr"
     }
     if(missing(true.topo)){
	sims <- make.ml.tr(sims, trees.file, log.file, empdat, model = model)
     } else {
     	sims <- make.ml.tr(sims, trees.file, empdat, truetr = true.topo, model = model)
     }     
     return(sims)
}