run.mb <- function(sdata, format = "phyllip", temp_name, model = "GTR+G", topofix = F, ss = F, mbPath = "~/Downloads/MrBayes/mb"){
	
	if(format == "fasta"){
		d <- read.dna(sdata, format = "fasta")
		fileName <- gsub("fasta", "phy", sdata)
	} else if(format == "phyllip"){
		d <- read.dna(sdata)
		fileName <- gsub("fasta", "phy", sdata)
	} else if(format == "DNAbin"){
		fileName <- temp_name
	}
	
	### Block 1. DNA data.
	if(ss == F){
	      b1 <- write.nexus.data(as.list(d), file = paste0(file.name, ".nex"))
	} else {
	      b1 <- write.nexus.data(as.list(d), file = paste0(file.name, "_ss.nex"))
	}

	### Block 5. Topology.
	b5 <- ""
	if(topofix == T){
	tr.topo <- sim.data[[1]]
	tr.topo$edge.length <- abs(rnorm(n = length(tr.topo$edge.length)))
	b5 <- paste0("begin trees;
		tree tr = ", write.tree(tr.topo),"
		end;
		prset topologypr = fixed(tr); [ Fixes the topology ]")
	}

	### Block 2. 
	if(model == "GTR+G"){
		b2 <- "begin mrbayes;
		lset nst=6;
		lset rates=Invgamma;
		lset ngammacat=4;"
	} else if(model == "JC"){
		b2 <- "begin mrbayes;
		lset nst=1"
	}
	
	### Block 6. Final details.
	b6 <- "
	prset brlenspr=Unconstrained:Exp(10.0);
	mcmcp ngen=1000000 diagnfreq=100 samplefreq=100 printfreq=100000 Nruns=1;
	"
	
	
	b7 <- "
    	   mcmc;
	   end;
	"

	if(ss == T){
	      b7 <- "
	      ss ngen = 40000 diagnfreq = 4000;
	      end;
	      "
	}

	complete.nex <- c(b1, b2, b5, b6, b7)
	
	if(ss == F){
	      cat(complete.nex ,file= paste0(fileName, ".nex"), append=TRUE)
	      system(paste(mbPath, paste0(fileName, ".nex")))
	} else {
	      cat(complete.nex ,file= paste0(fileName, "_ss.nex"), append=TRUE)
	      system(paste(mbPath, paste0(fileName, "_ss.nex")))
	}

}
	