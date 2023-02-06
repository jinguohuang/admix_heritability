#!/usr/bin/env Rscript

# simple simulation to validate hypothesis
#output1: expected value of the simulation
#output2: tped, tfam file for gcta
#output3: gvalue file for gcta

library(data.table)
library(ggplot2)
library(dplyr)
library(optparse)
option_list = list(
  make_option(c("--fstc", "-F"), type="numeric", default = 0.5, 
              help="desired proportion of total genetic variance that is between parental groups at causal loci [default= %default]", metavar="character"),
  make_option(c("--fstn", "-f"), type="numeric", default = 0.5, 
              help="desired proportion of total genetic variance that is between parental groups at noncausal loci [default= %default]", metavar="character"),
  make_option(c("--h2", "-H"), type = "numeric", default = 0.8,
              help = "simulated heritability [default = %default]"),
  make_option(c("--theta", "-t"), type="numeric", default=0.5, 
              help="admixture proportion [default= %default]", metavar="character"),
  make_option(c("--nloci", "-l"), type = "numeric", default = 1000, 
              help = "number of loci  [default = %default]", metavar="character"),
  make_option(c("--pcausal", "-p"), type = "numeric", default = 0.2, 
              help = "proportion of causal variants  [default = %default]", metavar="character"),
  make_option(c("--sample_size", "-n"), type = "numeric", default = 1e4, 
              help = "sample size of admixed people [default = %default]", metavar="character")
  #make_option(c("--seed"), type="numeric", default=120, 
             # help="random seed [default= %default]", metavar="character")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
print("setting up")

#Set up the parameters.
fbar = 0.5 # total frequency in the parental population

#vgt = opt$vgt #desired total genetic variance
vgt = opt$h2
theta = opt$theta #admixture fraction
#vgb = opt$vgb #desired genetic variance in between
vgb = opt*fstc * vgt
nloci = opt$nloci #number of loci, not necessary causal variants
pcausal = opt$pcausal # proportion of causal variants
ncausal = round(pcausal*nloci) # get number of causal variants
#seed=opt$seed #random seed, 100
venv = 1 - opt*h2 # environmental variance
sample_size = opt$sample_size # sample size of simulated admixed people

filename <- paste0("_prop", theta, "_vgb", vgb, "_vgt", vgt, "_venv", venv) #for filename
#set.seed(seed)

#vgb = vgt*Fst
#fst = vgb/vgt
#next, calculate frequency difference required to generate that Fst
#fst = (f1 - f2)^2 / 4*pbar*(1-pbar)
fdiffc = sqrt(fstc*4*fbar*(1-fbar))
#solve system of linear equation
f1_f2c = solve(a = matrix(c(0.5, 1, 0.5, -1), nrow = 2, ncol = 2), 
              b = c(fbar, fdiffc))

# 0.5*f1+0.5*f2=fbar
# f1-f2=fdiff
# given fbar and fdiff, solve the linear equation to get f1,f2
f1c = f1_f2c[1]
f2c = f1_f2c[2]

################ non causal loci ####################

#fst = (f1 - f2)^2 / 4*pbar*(1-pbar)
fdiffn = sqrt(fstn*4*fbar*(1-fbar))
#solve system of linear equation
f1_f2n = solve(a = matrix(c(0.5, 1, 0.5, -1), nrow = 2, ncol = 2), 
               b = c(fbar, fdiffn))

# 0.5*f1+0.5*f2=fbar
# f1-f2=fdiff
# given fbar and fdiff, solve the linear equation to get f1,f2
f1n = f1_f2n[1]
f2n = f1_f2n[2]

print("simulating genotypes and local ancestry")


######causal 

#simulate genotypes for pop1 (1000 loci and 1000 individuals)
#including non-causal variants
geno1 = replicate(ncausal, rbinom(1000, 2, f1))
#simulate genotypes for pop2
geno2 = replicate(ncausal, rbinom(1000, 2, f2))

#Sample effect sizes for the causal loci 
#such that the genetic variance in the total population is 1.

#choose effect sizes for causal loci, such that the total population variance is 1
bg = matrix(
  rnorm(ncausal, 0, sd= sqrt(vgt/(ncausal*2*fbar*(1-fbar)))), 
  nrow = ncausal, ncol = 1)

#standardize to remove stochasticity in realized values
bg = ((bg - mean(bg)) / sd(bg)) * sqrt(vgt/(ncausal*2*fbar*(1-fbar)))

#calculate variance between populations
# pbar = (f1 + f2)/2
# vtotal = sum(2 * bg^2 * pbar * (1 - pbar))
# vwithin = sum(2*(bg^2)*f1*(1 - f1) + 2*(bg^2)*f2*(1 - f2))/2
# vbetween = vtotal - vwithin

#Simulate local ancestry and genotypes for admixed population 
#such that the admixture fraction is 0.5.
lanc_c = replicate(ncausal, rbinom(sample_size, 2, theta))

#simulate genotypes at each locus given the ancestry at that locus.

#causal
geno.admix_c = structure(sapply(lanc_c, function(x){
  if(x == 0){
    g1 = rbinom(1, 1, f1c)
    g2 = rbinom(1, 1, f1c)
  }
  if(x == 1){
    g1 = rbinom(1, 1, f1c)
    g2 = rbinom(1, 1, f2c)
  }
  if(x == 2){
    g1 = rbinom(1, 1, f2c)
    g2 = rbinom(1, 1, f2c)
  }
  return(g1 + g2)
}), dim = dim(lanc))


#########non-causal

lanc_n = replicate(nloci - ncausal, rbinom(sample_size, 2, theta))

geno.admix_n = structure(sapply(lanc_n, function(x){
  if(x == 0){
    g1 = rbinom(1, 1, f1n)
    g2 = rbinom(1, 1, f1n)
  }
  if(x == 1){
    g1 = rbinom(1, 1, f1n)
    g2 = rbinom(1, 1, f2n)
  }
  if(x == 2){
    g1 = rbinom(1, 1, f2n)
    g2 = rbinom(1, 1, f2n)
  }
  return(g1 + g2)
}), dim = dim(lanc))

geno.admix = cbind(geno.admix_c, geno.admix_n)

bg = c(bg, rep(0, nloci - ncausal))

#Simulate genetic values in the admixed population.
gvalue.admix = t(t(bg)%*%t(geno.admix))

phenotype = as.data.table(gvalue.admix)
colnames(phenotype)=c("gvalue")
phenotype$environment = rnorm(n = sample_size, mean = 0,
                        sd = sqrt(venv))

#add prs to each of the environmental effects 
phenotype$pheno = phenotype$gvalue + phenotype$environment

# normalize gvalue
phenotype$gvalue_norm = (phenotype$gvalue - mean(phenotype$gvalue))/sd(phenotype$gvalue)

#add prs to each of the environmental effects 
phenotype$pheno_norm = phenotype$gvalue_norm + phenotype$environment

#print(paste("correlation:",
#            round(cor(phenotype$gvalue,
#                      phenotype$pheno)^2,3)))
vgeno=round(cor(phenotype$gvalue, phenotype$pheno)^2,3)

# output file with the expected value: 
# vg, admixture, vtotal, vwithin, vbetween
sink(file=paste0("sim",filename, "_expected.txt"), type = "output")
cat("Source\tExpected\n")
cat("vgb\t", vgb, "\n")
cat("admixture\t", theta, "\n")
cat("fst\t", fst, "\n")
#cat("seed\t", seed, "\n")
cat("vtotal\t", vtotal, "\n")
cat("vwithin\t", vwithin, "\n")
cat("vbetween\t", vbetween, "\n")
cat("venv\t", venv, "\n")
cat("vgeno\t", vgeno, "\n")
sink()

print("writing to plink file")
# transpose it and add map info: chrom, ID, cm (0), position
lanc.t<-t(lanc)
map = data.table(SNP = paste0(1, "_", 1:nloci, "_AG"), A1 = "A", A2 = "G")
lanc.t = cbind(map, lanc.t)
write.table(lanc.t, file = paste0("admix", filename, ".dosage"), quote=FALSE, sep = '\t' ,
            row.names = FALSE, col.names = FALSE)
# make tfam file
fam<-data.frame(FID="ADM",IID=1:sample_size,FID=0,MID=0,SEX=0,Pheno=-9 )
write.table(fam, file = paste0("admix", filename, ".tfam"), quote=FALSE, sep = '\t' ,
            row.names = FALSE, col.names = FALSE)
#phenotype = cbind(fam[,c(1,2)], c(gvalue.admix))
phenotype = cbind(fam[,c(1,2)], phenotype[,c(1,3,4,5)])
fwrite(phenotype, file = paste0("admix", filename, ".pheno"), quote=FALSE, sep = '\t' ,
                   row.names = FALSE, col.names = FALSE)





