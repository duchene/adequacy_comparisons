library(NELSI)
for(f in dir('../R/', pattern = 'R$')){
    print(f)
    source(paste0('../R/', f))
}

boot_folders <- dir(pattern = 'Lboot.+_folder')

out_matrix <- matrix(NA, length(boot_folders), 5)
colnames(out_matrix) <- c('multlik', 'delta', 'chisq', 'TL', 'node_support')


for(i in 1:length(boot_folders)){
    b <- boot_folders[i]
    sims <- dir(b, pattern = 'phy$')
    emp_data <- read.dna(gsub('_folder', '.phy', b))
    emp_chisq <- get.chisq(emp_data)
    phyml_emp_data <- run.phyml(gsub('_folder', '.phy', b))
    pvals_mat <- matrix(NA, length(sims), 5)
    for(f in 1:length(sims)){
        if(length(grep('gtr', b)) == 1){
            phyml_temp <- run.phyml(paste0(b, '/', sims[f]), model = 'GTR+G')
        }else if(length(grep('jc', b)) == 1){
            phyml_temp <- run.phyml(paste0(b, '/', sims[f]), model = 'JC')
        }else{
            stop('model not recognised')
        }
        dat_temp <- read.dna(paste0(b, '/', sims[f]))
        chisq_temp <- get.chisq(al = dat_temp)
        pvals_mat[f, ] <- c(phyml_temp$uncLikelihood, phyml_temp$uncLikelihood - phyml_temp$likelihood,
                            chisq_temp, phyml_temp$treeLength, phyml_temp$meanNodeSupport)
    }
    pvals <- c(sum(phyml_emp_data$uncLikelihood > pvals_mat[,1]),
               sum((phyml_emp_data$uncLikelihood - phyml_emp_data$likelihood) > pvals_mat[,2]),
              sum(emp_chisq > pvals_mat[,3]), sum(phyml_emp_data$treeLength > pvals_mat[,4]),
              sum(phyml_emp_data$meanNodeSupport > pvals_mat[,5]))
    out_matrix[i, ] <- pvals
}

write.table(out_matrix, file = 'out_pvals.csv', row.names = F, sep = ',')
