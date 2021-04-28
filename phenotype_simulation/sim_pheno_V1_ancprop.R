#!/bin/env Rscript
# simulate phenotype data
library(data.table)
library(dplyr)
library(tidyr)
library(rprojroot)
args = commandArgs(trailingOnly=TRUE)
gvalue_file=args[1] # genetic score file
Name <- sub("(.*?)\\.profile", '\\1', gvalue_file)

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
#hist(prs$SCORE)

#add prs to each of the environmental effects
prs = prs %>%
  mutate(environment = SCORE + environment)
#hist(prs$environment)
# normalize prs
prs = prs %>%
  mutate(environment = (environment-mean(environment))/sd(environment))
#hist(prs$environment)
#the correlation between the phenotype and the genetic values should be ~ 0.8
print(paste("correlation:",
            round(cor(prs$SCORE,
                      prs$environment)^2,2)))
# output the phenotype data
sim_pheno<-prs[,c("FID","IID", "environment")]
write.table(sim_pheno, file = paste0("pheno_", Name,".txt"), quote=F, row.names = F, col.names = F)
