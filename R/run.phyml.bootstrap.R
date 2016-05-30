run.phyml.bootstrap <- function(phyml_path, file_name, subs_model, n_reps, n_proc){
    require(ape)
    require(phangorn)
    require(doParallel)
    require(foreach)

    nuc_data <- as.matrix(read.dna(file_name, format = 'sequential'))
    if(!is.matrix(nuc_data)) stop('data not read as a matrix. Check that the format is fasta and that it consists of nucleotides')

    run_phyml <- function(phymlPath, data_name, model='GTR+G'){
        if(model == 'GTR+G'){
            system(paste(phymlPath, '-i', data_name, '-m gtr -a e'))
            stats_raw<- readLines(grep(paste0(data_name, '_phyml_stats'), x = dir(), value = T))
            relative_rates <- gsub('[A-Z]|[-]|<|>| ', '', grep('(A <-> C)|(A <-> G)|(A <-> T)|(C <-> G)|( C <-> T)|(G <-> T)', stats_raw, value = T))
            frequencies <- gsub('[-]| |=|[A-Z]|[a-z]|[(]|[)]','',     grep('(f[(]A[)])|(f[(]C[)])|(f[(]G[)])|(f[(]T[)])', stats_raw, value = T))
            alpha <- gsub(' |[-]|[A-Z]|[a-z]|:|\t', '', grep('Gamma shape parameter', stats_raw, value = T))
            tree <- readLines(grep(paste0(data_name, '_phyml_tree'), x = dir(), value = T))
            stats <- c(relative_rates, frequencies, alpha, tree)
            return(stats)
        }else if(model == 'JC'){
            system(paste(phymlPath, '-i', data_name, '-m 000000 -a 1 -c 1'))
            tree <- readLines(grep(paste0(data_name, '_phyml_tree'), x = dir(), value = T))
            names(tree) <- 'tree'
            return(tree)
        }else{
            stop('Choose either GTR+G, or JC substitution models')
        }
    }

    get_bootstrap_replicate <- function(phymlPath, nuc_data, rep_name, subs_model){
        sample_data <- nuc_data[, sample(x = 1:ncol(nuc_data), size = ncol(nuc_data), replace = T)]
        write.dna(sample_data, file = rep_name, format = 'sequential', nbcol = -1, colsep = '')
        out <- run_phyml(phymlPath, data_name = rep_name, model = subs_model)
        system(paste0('rm ', rep_name, '*'))
        return(out)
    }

    cl <- makeCluster(n_proc)
    registerDoParallel(cl)

    reps <- foreach(x = 1:n_reps, .packages = 'ape') %dopar% get_bootstrap_replicate(phyml_path, nuc_data = nuc_data, rep_name = paste0('rep', x), subs_model = subs_model)

    rbind_list <- function(data_list){
    if(length(data_list) == 1){
        data_list
    }else if(length(data_list) == 2){
        rbind(data_list[[1]], data_list[[2]])
    }else{
        rbind(data_list[[1]], rbind_list(data_list[-1]))
    }
}
    out_log <- rbind_list(reps)
    if(subs_model == 'GTR+G'){
        colnames(out_log) <- c('ac', 'ag', 'at', 'cg', 'ct', 'gt', 'a', 'c', 'g', 't', 'alpha', 'tree')
        write.table(out_log[, -ncol(out_log)], file = paste0(file_name, '_ml_bootstrap.log'), row.names = F, quote = F)
        writeLines(text = out_log[,ncol(out_log)], con = paste0(file_name, '_ml_bootstrap.trees'), sep = '\n')
    }else{
        writeLines(text = out_log, con = paste0(file_name, '_ml_bootstrap.trees'), sep = '\n')
    }
    out_log
}


#true_data <- as.DNAbin(simSeq(rtree(10)))
#write.dna(true_data, file = 'test_true_data.phy', format = 'sequential', nbcol = -1, colsep = '')
#phyml_path <- '~/Desktop/phylo_programs/PhyML-3.1/PhyML-3.1_macOS-MountainLion'
#file_name <- 'test_true_data.phy'
#subs_model <- 'GTR+G'
#n_reps <- 50
#n_proc <- 5
#t1 <- run.phyml.bootstrap(phyml_path, 'test_true_data.phy', subs_model = 'GTR+G', n_reps = 5, n_proc = 5)



