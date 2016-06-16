# This function requires the function in this repository get.seq.
require(phangorn)

get.het.seq <- function(ntax = 10, exp.rate = 5, backbone.ntax = 10, trlen = 3, backbone.trlen = 1.5, basefreq.offshoot = c(0.45, 0.05, 0.05, 0.45), slen = 1000){
    tr1 <- rtree(backbone.ntax)
    tr1$edge.length <- rexp(length(tr1$edge.length), exp.rate)
    s1 <- get.seq(tr1, model = "gtr.g", l = slen)
    
    tr2 <- rtree(ntax - backbone.ntax)
    tr2$edge.length <- rexp(length(tr2$edge.length), exp.rate)
    s2 <- get.seq(tr2, model = "gtr.g", l = slen, basef = basefreq.offshoot, rootseq =  as.character(as.DNAbin(s1)[1, ]))

    tr3 <- bind.tree(tr1, tr2, where = 1)
    tr3$tip.label <- paste0('t', 1:length(tr3$tip.label))

    names(s1) <- paste0('t', 1:length(s1))
    names(s2) <- paste0('t', length(s1):(length(s1)-1 + length(s2)))

    s3 <- rbind(as.DNAbin(s1)[-length(s1), ], as.DNAbin(s2))
    
    return(list(tr3, s3))
}
