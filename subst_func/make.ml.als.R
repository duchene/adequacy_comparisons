# This function takes posterior log and trees files. It simulates alignments using the model parameters from the corresponding line in the log file. This function also requires the length of the alignment.

make.pps.als <- function(trees.file, log.file, N = 100, l = 1000, savetrees = F){
	
	trees <- try(read.nexus(trees.file))
	if(class(trees) == "try-error") trees <- try(read.tree(trees.file), silent = T)
	if(class(trees) == "try-error") stop("Cannot read trees")
	print(class(trees))
	if(length(grep("[.]csv", log.file, value = T)) == 0){
		logdat <- read.table(log.file, header = T, comment = "[")
	} else {
	        logdat <- read.table(log.file, head = T, row.names = 1, sep = ",")
	}
	print(class(logdat))
	if(length(trees) == N){
		endburnin <- 1
	} else {
		endburnin <- round((length(trees) * 0.1), 0)
	}
	samp <- sample(endburnin:length(trees), N)
	trees <- trees[samp]
	logdat <- logdat[samp,]
	#if("alpha" %in% colnames(logdat)){
	#	if(length(which(logdat$alpha < 0.01)) > 0 || length(which(logdat$alpha > 1.2)) > 0){
	#		#print("Extreme values of alpha have been modified to allow simulation.")
	#		logdat$alpha[which(logdat$alpha > 1.2)] <- 1.2
	#		for(i in 5:ncol(logdat)) logdat[,i][which(logdat[,i] < 0.01)] <- 0.5
	#	}
	#}
	if(savetrees){
	      write.tree(trees, file = paste0(trees.file, N, ".tre"))
	      write.csv(logdat, file = paste0(trees.file, N, ".csv"))
	}
	sim <- list()
	for(i in 1:nrow(logdat)){
	      sim[[i]] <- list(phylogram = trees[[i]])
	      
	      if(all(c("r.A...C.", "r.A...G.", "r.A...T.", "r.C...G.", "r.C...T.", "r.G...T.") %in% colnames(logdat))){
	      	     # GENERAL TIME REVERSIBLE (GTR)
		     #print("The substitution model is GTR")
		     basef <- c(logdat$pi.A.[i], logdat$pi.C.[i], logdat$pi.G.[i], logdat$pi.T.[i])
		     qmat <- c(logdat$r.A...C[i], logdat$r.A...G.[i], logdat$r.A...T.[i], logdat$r.C...G.[i], logdat$r.C...T.[i], logdat$r.G...T.[i])
		     #print(basef)
		     #print(qmat)

		     if("alpha" %in% colnames(logdat)){
		     	   rates = phangorn:::discrete.gamma(logdat$alpha[i], k = 4)
			   rates <- rates + 0.001
        		   sim_dat_all<- lapply(rates, function(r) simSeq(sim[[i]][[1]], l = round(l/4, 0), Q = qmat, bf = basef, rate = r))
        		   sim[[i]][[3]] <- c(sim_dat_all[[1]], sim_dat_all[[2]], sim_dat_all[[3]], sim_dat_all[[4]])
		     } else {		    
	      	       	   sim[[i]][[3]] <- simSeq(sim[[i]][[1]], Q = qmat, bf = basef, l = l)
		     }		     
		     #print("DATA SIMULATION PROCESSED")
	      } else if("kappa" %in% colnames(logdat)){
		     # HASEGAWA-KISHINO-YANO (HKY)
		     #print("The substitution model is HKY")
		     basef <- c(logdat$pi.A.[i], logdat$pi.C.[i], logdat$pi.G.[i], logdat$pi.T.[i])
		     qmat <- c(1, 2*logdat$kappa[i], 1, 1, 2*logdat$kappa[i], 1)

		     if("alpha" %in% colnames(logdat)){
		     	   rates = phangorn:::discrete.gamma(logdat$alpha[i], k = 4)
			   rates <- rates + 0.001
       			   sim_dat_all<- lapply(rates, function(r) simSeq(sim[[i]][[1]], l = round(l/4, 0), Q = qmat, bf = basef, rate = r))
                           sim[[i]][[3]] <- c(sim_dat_all[[1]], sim_dat_all[[2]], sim_dat_all[[3]], sim_dat_all[[4]])
		     } else {
                           sim[[i]][[3]] <- simSeq(sim[[i]][[1]], Q = qmat, bf = basef, l = l)
	       	     }
		     
	      } else { 
	      	     # JUKES-CANTOR (JC)
		     #print("The substitution model is assumed to be JC")
	      	     sim[[i]][[3]] <- simSeq(sim[[i]][[1]], l = l)
	      }

	}
	
	return(sim)

}