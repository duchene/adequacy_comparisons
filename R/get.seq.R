# This function produces sequences for six possible models of substitution, but with arbitrary values, specified in this function. When using the gamma parameter, this function assumes six rate categories.

get.seq <- function(tree, model, l = 1000, qmat = c(1.3472, 4.8145, 0.9304, 1.2491, 5.5587, 1), basef = c(0.2628, 0.2605, 0.2436, 0.2331), gamma.parameter = 1.0406, invprop = 0, rootseq = NULL){

	if(length(grep("[.]g", model)) > 0){
		if(invprop > 0){
            invlen <- round(l * invprop, 0)
            varlen <- l - invlen
		}
		rates = phangorn:::discrete.gamma(gamma.parameter, k = 4)
        sim_dat_all<- lapply(rates, function(r) simSeq(tree, l = round(l/4, 0), Q = qmat, bf = basef, rate = r, rootseq = rootseq[, 1:round(l/4, 0)]))
        al <- c(sim_dat_all[[1]], sim_dat_all[[2]], sim_dat_all[[3]], sim_dat_all[[4]])
		if(invprop > 0){
			sim_dat_all<- lapply(rates, function(r) simSeq(tree, l = round(varlen/4, 0), Q = qmat, bf = basef, rate = r, rootseq = rootseq[, 1:round(varlen/6, 0)]))
			al1 <- c(sim_dat_all[[1]], sim_dat_all[[2]], sim_dat_all[[3]], sim_dat_all[[4]])
			al2 <- simSeq(tree, l = invlen, Q = qmat, bf = basef, rate = 0.05, rootseq = rootseq)
			al <- c(al1, al2)
			return(al)
		}
    } else {
	    if(invprop > 0){
			invlen <- round(l * invprop, 0)
                	varlen <- l - invlen
			al1 <- simSeq(tree, Q = qmat, bf = basef, l = invlen, rate = 0.05, rootseq = rootseq)
			al2 <- simSeq(tree, Q = qmat, bf = basef, l = varlen, rootseq = rootseq)
			al <- c(al1, al2)
		} else {
            al <- simSeq(tree, Q = qmat, bf = basef, l = l, rootseq = rootseq)
		}
	}
	return(al)

}