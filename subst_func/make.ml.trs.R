# This function processes simulated datasets created with the function make.pps.als. It takes the posterior phylograms and simulated alignments. The function creates the posterior predictive simulated distribution of phylograms. This function also requires the empirical dataset and, if desired, the true topology.

make.ml.tr <- function(sims, trees.file, log.file, empdat, truetr, phymlpath = "~/Desktop/Software/PhyML-3.1/PhyML-3.1_macOS-MountainLion", model = "jc69"){
	 require(MASS)
	 require(phangorn)
	 empmlik <- multlik(empdat)
	 empchisq <- getchisqs(empdat)
	 empconsind <- vector()
	 emplog <- readLines(log.file)
	 empuncons <- as.numeric(gsub(".*\t", "", grep("Unconstrained", emplog, value = T)))
	 emplik <- as.numeric(gsub(".*\t", "", grep("Log-likelihood", emplog, value = T)))
	 empdelta <- empuncons - emplik
	 tr <- read.tree(trees.file)
	 if(!missing(truetr)){ truelen <- sum(truetr$edge.length)
	 		       empconsind <- CI(truetr, empdat)
	 }
	 simultlik <- vector()
	 simtrlen <- vector()
	 simchisq <- vector()
	 topodists <- vector()
	 consind <- vector()
	 simdelta <- vector()
	 for(i in 1:length(sims)){
	       simultlik[i] <- multlik(sims[[i]])
	       simchisq[i] <- getchisqs(sims[[i]])
	       simtrlen[i] <- sum(tr$edge.length)
	       empconsind[i] <- CI(tr, empdat)
	       consind[i] <- CI(tr, sims[[i]])
	       if(!missing(truetr)) topodists[i] <- dist.topo(truetr, tr)
	       write.dna(as.DNAbin(sims[[i]]), file = "temp_sim.phy")
	       system(paste(phymlpath, "-m", model, "-a e --q -i temp_sim.phy"))
	       simlog <- readLines("temp_sim.phy_phyml_stats.txt")
	       simuncons <- as.numeric(gsub(".*\t", "", grep("Unconstrained", simlog, value = T)))
	       simlik <- as.numeric(gsub(".*\t", "", grep("Log-likelihood", simlog, value = T)))
	       simdelta[i] <- simuncons - simlik
	       system("rm temp_sim.phy temp_sim.phy_phyml_stats.txt temp_sim.phy_phyml_tree.txt")
	       print(paste("Simulation", i, "processed"))
	 }
	 #multlikp <- length(simultlik[which(simultlik < empmlik)]) / length(sims)
	 multlikdistr <- fitdistr(simultlik, "normal")
	 multlikp <- 1 - pnorm(empmlik, multlikdistr[[1]][1], multlikdistr[[1]][2])
	 chisqdistr <- fitdistr(simchisq, "normal")
         chisqp <- 1 - pnorm(empchisq, chisqdistr[[1]][1], chisqdistr[[1]][2])
	 consindistr <- fitdistr(consind, "normal")
	 consp <- 1 - pnorm(mean(empconsind), consindistr[[1]][1], consindistr[[1]][2])
	 deltadistr <- fitdistr(simdelta, "normal")
	 deltap <- 1 - pnorm(empdelta, deltadistr[[1]][1], deltadistr[[1]][2])

	 if(!missing(truetr)){
	       trlenp <- length(simtrlen[which(simtrlen < truelen)]) / length(sims)
	       res <- list(empirical.multinomial.likelihood = empmlik, pps.multinomial.likelihood = simultlik, multinomial.likelihood.p.value = multlikp, true.tree.length = truelen, pps.tree.lengths = simtrlen, tree.length.p.value = trlenp, topological.distances = topodists, empirical.chisq = empchisq, pps.chisq = simchisq, chisq.p.value = chisqp, pps.consistency.index = consind, consistency.index.p.value = consp, empirical.delta = empdelta, sim.deltas = simdelta, delta.p.value = deltap)
	 } else {
	       res <- list(empirical.multinomial.likelihood = empmlik, pps.multinomial.likelihood = simultlik, multinomial.likelihood.p.value = multlikp, empirical.chisq = empchisq, pps.chisq = simchisq, chisq.p.value = chisqp, pps.consistency.index = consind, consistency.index.p.value = consp, empirical.delta = empdelta, sim.deltas = simdelta, delta.p.value = deltap)
	 }

	 return(res)
}