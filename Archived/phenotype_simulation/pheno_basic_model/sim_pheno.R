#!/bin/env Rscript
# simulate phenotype data
gscore<-read.table("sim_10Kppl_geneticscore.txt", header = F)
colnames(gscore)<-c("FID","IID","PHENO","CNT","CNT2","SCORE")
# environment effect
e = rnorm(10000, mean = 0, sd = sqrt(0.2))
# normalize genetic score
g=gscore$SCORE
g = (g - mean(g))/sd(g)
#hist(g)
y = g + e
# normalize y
y = (y - mean(y))/sd(y)
#hist(y)
# output the phenotype data
gscore$PHENO<-y
sim_pheno<-gscore[,c("FID","IID", "PHENO")]
write.table(sim_pheno, file = "sim_pheno_10Kppl_random_2.txt", quote=F, row.names = F, col.names = F)
