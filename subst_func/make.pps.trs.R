# This function processes simulated datasets created with the function make.pps.als. It takes the posterior phylograms and simulated alignments. The function creates the posterior predictive simulated distribution of phylograms. This function also requires the empirical dataset and, if desired, the true topology.

make.pps.tr <- function(sims, empdat, truetr){
	 require(MASS)
	 require(phangorn)
	 empmlik <- multlik(empdat)
	 empchisq <- getchisqs(empdat)
	 empconsind <- vector()
	 if(!missing(truetr)){ truelen <- sum(truetr$edge.length)
	 		       empconsind <- CI(truetr, empdat)
	 }
	 simultlik <- vector()
	 simtrlen <- vector()
	 simchisq <- vector()
	 topodists <- vector()
	 consind <- vector()
	 for(i in 1:length(sims)){
	       simultlik[i] <- multlik(sims[[i]][[3]])
	       simchisq[i] <- getchisqs(sims[[i]][[3]])
	       simtrlen[i] <- sum(sims[[i]][[1]]$edge.length)
	       empconsind[i] <- CI(sims[[i]][[1]], empdat)
	       consind[i] <- CI(sims[[i]][[1]], sims[[i]][[3]])
	       if(!missing(truetr)) topodists[i] <- dist.topo(truetr, sims[[i]][[1]])
	       print(paste("Simulation", i, "processed"))
	 }
	 #multlikp <- length(simultlik[which(simultlik < empmlik)]) / length(sims)
	 multlikdistr <- fitdistr(simultlik, "normal")
	 multlikp <- 1 - pnorm(empmlik, multlikdistr[[1]][1], multlikdistr[[1]][2])
	 chisqdistr <- fitdistr(simchisq, "normal")
         chisqp <- 1 - pnorm(empchisq, chisqdistr[[1]][1], chisqdistr[[1]][2])
	 consindistr <- fitdistr(consind, "normal")
	 consp <- 1 - pnorm(mean(empconsind), consindistr[[1]][1], consindistr[[1]][2])

	 if(!missing(truetr)){
	       trlenp <- length(simtrlen[which(simtrlen < truelen)]) / length(sims)
	       res <- list(empirical.multinomial.likelihood = empmlik, pps.multinomial.likelihood = simultlik, multinomial.likelihood.p.value = multlikp, true.tree.length = truelen, pps.tree.lengths = simtrlen, tree.length.p.value = trlenp, topological.distances = topodists, empirical.chisq = empchisq, pps.chisq = simchisq, chisq.p.value = chisqp, pps.consistency.index = consind, consistency.index.p.value = consp)
	 } else {
	       res <- list(empirical.multinomial.likelihood = empmlik, pps.multinomial.likelihood = simultlik, multinomial.likelihood.p.value = multlikp, empirical.chisq = empchisq, pps.chisq = simchisq, chisq.p.value = chisqp, pps.consistency.index = consind, consistency.index.p.value = consp)
	 }

	 return(res)
}