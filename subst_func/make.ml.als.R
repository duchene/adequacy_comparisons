# This function takes phyml log and trees files. It simulates alignments using the model parameters from the corresponding line in the phyml file. This function also requires the length of the alignment.

make.ml.als <- function(trees.file, log.file, N = 100, l = 1000, savetrees = F){
	
	trees <- try(read.tree(trees.file))
	print("TREE:")
	print(trees)
	if(class(trees) == "try-error") stop("Cannot read trees")
	phymllog <- readLines(log.file)
	model <- phymllog[grep("Model of nucleotides substitution", phymllog)]
	#print(model)
	sim <- list()
	gammamod <- phymllog[grep("Number of classes", phymllog)]
	#print(gammamod)
	for(i in 1:N){
	if(length(grep("GTR", model)) > 0){
		qmat <- phymllog[grep("rate matrix", phymllog)+2:4]
		qmat <- as.numeric(sapply(qmat, function(x) setdiff(strsplit(x, " ")[[1]], ""))[c(2,3,4,7,8,12)])
		print(paste("Q matrix =", qmat))
		basef <- phymllog[grep("Nucleotides frequencies", phymllog)+1:4]
		basef <- as.numeric(sapply(basef, function(x) gsub(".*= ", "", x)))
		#print(paste("Base frequencies:", basef))
		if(length(grep("4", gammamod)) > 0){
			shape <- phymllog[grep("Gamma shape parameter", phymllog)]
			shape <- gsub(".*: ", "", shape)
			shape <- as.numeric(gsub("\t", "", shape))
			print(paste("Shape parameter:", shape))
			rates = phangorn::discrete.gamma(shape, 4)
			#print(paste("gamma rates", rates))
		    sim_dat_all <- lapply(rates, function(x) simSeq(trees, l = round(l/4, 0), Q = qmat, bf = basef, rate = x))
			sim[[i]] <- c(sim_dat_all[[1]], sim_dat_all[[2]], sim_dat_all[[3]], sim_dat_all[[4]])
		} else {
		        sim[[i]] <- simSeq(trees, Q = qmat, bf = basef, l = l)
		}

	} else if(length(grep("HKY", model)) > 0){



	} else if(length(grep("JC", model)) > 0){
	        if(length(grep("4", gammamod)) > 0){
                        shape <- phymllog[grep("Gamma shape parameter", phymllog)]
                        shape <- gsub(".*: ", "", shape)
			shape <- as.numeric(gsub("\t", "", shape))
			print(paste("Shape parameter:", shape))
                        rates = phangorn::discrete.gamma(shape, 4)
						sim_dat_all <- lapply(rates, function(x) simSeq(trees, l = round(l/4, 0), rate = x))
						sim[[i]] <- c(sim_dat_all[[1]], sim_dat_all[[2]], sim_dat_all[[3]], sim_dat_all[[4]])
			} else {
						sim[[i]] <- simSeq(trees, l = l)
			}

	}
	}
	return(sim)

}