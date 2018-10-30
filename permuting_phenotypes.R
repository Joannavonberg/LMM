pheno <- read.table("stroke.all.consensus.AFR.EUR.QC2.phenotypes.3March2016_withheader.txt", header = T)
set.seed(seed = 584359062)
ind <- sample(1:42653, size = 42653, replace = FALSE)

pheno_permuted <- cbind(pheno[,1:14], pheno[ind,15:38])

write.table(pheno_permuted, file = "permutedlabels.stroke.all.consensus.AFR.EUR.QC2.phenotypes.3March2016_withheader.txt", quote = F, row.names = F, col.names = T, sep = " ")


options(stringsAsFactors = F)
pheno <- read.table("stroke.all.AFR.EUR.QC2.phenotypes.3Oct2014.txt", header = T)
set.seed(seed = 584359062)
ind <- pheno$FID == "0" | pheno$FID == "ID_1"
pheno <- pheno[!ind,]
ind <- sample(1:nrow(pheno), size = nrow(pheno), replace = FALSE)
pheno_permuted <- cbind(pheno[,1:14], pheno[ind,15:135])

write.table(pheno_permuted, file = "permutedlabels.stroke.all.AFR.EUR.QC2.phenotypes.3Oct2014.txt", quote = F, row.names = F, col.names = T, sep = " ")

