#!/bin/env Rscript
# simulate phenotype data
gvalue_file="sim_10Kppl_geneticscore_V1.txt"
prs=fread(gvalue_file)
colnames(prs)<-c("FID","IID","PHENO","CNT","CNT2","SCORE")
sample_size=nrow(prs)
#Add environmental component to the genetic values 
prs$environment = rnorm(n = sample_size, mean = 0,
                              sd = sqrt(1 - 0.8))
#hist(prs$environment)
#Scale the effects such that the heritability is 0.8.
#prs$environment = scale( prs$environment, scale = T) * sqrt(1 - 0.8)
# # normalize g
prs = prs %>%
  mutate(SCORE = (SCORE-mean(SCORE))/sd(SCORE))
hist(prs$SCORE)

#add prs to each of the environmental effects
prs = prs %>%
  mutate(environment = SCORE + environment)
hist(prs$environment)
# normalize prs
prs = prs %>%
  mutate(environment = (environment-mean(environment))/sd(environment))
hist(prs$environment)
#the correlation between the phenotype and the genetic values should be ~ 0.8
#print(paste("h2 (environment) :",
#            round(cor(prs$SCORE,
#                      prs$environment)^2,2)))
#0.83
# output the phenotype data
sim_pheno<-prs[,c("FID","IID", "environment")]
write.table(sim_pheno, file = "sim_pheno_10Kppl_random_V1.txt", quote=F, row.names = F, col.names = F)
