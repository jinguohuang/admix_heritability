#!/usr/bin/env Rscript

# simple simulation to validate hypothesis
#output1: expected value of the simulation
#output2: tped, tfam file for gcta
#output3: gvalue file for gcta
#output4: expected vs estimated effect size plot

# input1: between pop genetic variance (vgb)
# input2: admixture proportion (theta)
# input3: random seed (seed)

library(data.table)
library(ggplot2)
library(dplyr)
library(optparse)
option_list = list(
  make_option(c("--vgb", "-b"), type="numeric", default=0.3, 
              help="simulated genetic variance between parental groups[default= %default]", metavar="character"),
  make_option(c("--theta", "-t"), type="numeric", default=0.5, 
              help="admixture proportion [default= %default]", metavar="character"),
  make_option(c("--venv", "-e"), type = "numeric", default = 0, 
              help = "desired environmental variance [default = %default]")
  #make_option(c("--seed"), type="numeric", default=120, 
  #help="random seed [default= %default]", metavar="character")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print("setting up")

#Set up the parameters.
fbar = 0.5 # total frequency in the parental population
nloci = 1000 #number of loci, not necessary causal variants
vgt = 1 #desired total genetic variance

theta = opt$theta #admixture fraction
vgb = opt$vgb #desired genetic variance in between
#seed=opt$seed #random seed, 100

filename <- paste0("_prop", theta, "_vgb", vgb) #for filename
#set.seed(seed)

#vgb = vgt*Fst
fst = vgb/vgt
#next, calculate frequency difference required to generate that Fst
#fst = (f1 - f2)^2 / 4*pbar*(1-pbar)
fdiff = sqrt(fst*4*fbar*(1-fbar))
#solve system of linear equation
f1_f2 = solve(a = matrix(c(0.5, 1, 0.5, -1), nrow = 2, ncol = 2), 
              b = c(fbar, fdiff))
# 0.5*f1+0.5*f2=fbar
# f1-f2=fdiff
# given fbar and fdiff, solve the linear equation to get f1,f2
f1 = f1_f2[1]
f2 = f1_f2[2]

print("simulating genotypes and local ancestry")

#simulate genotypes for pop1 (1000 loci and 1000 individuals)
#including non-causal variants
geno1 = replicate(nloci, rbinom(1000, 2, f1))
#simulate genotypes for pop2
geno2 = replicate(nloci, rbinom(1000, 2, f2))

#Sample effect sizes for the 200 loci 
# such that the genetic variance in the total population is 1.

#proportion of loci that are causal
pcausal = 0.2 #1000 loci, 200 causal variants
ncausal = round(pcausal*nloci) # get number of causal variants
#choose effect sizes for causal loci, such that the total population variance is 1
bg = matrix(
  rnorm(ncausal, 0, sd= sqrt(vgt/(ncausal*2*fbar*(1-fbar)))), 
  nrow = ncausal, ncol = 1)

#standardize to remove stochasticity in realized values
bg = ((bg - mean(bg)) / sd(bg)) * sqrt(vgt/(ncausal*2*fbar*(1-fbar)))

#add 0s for non-causal loci
bg = sample(c(bg, rep(0, nloci - ncausal)))

#calculate variance between populations
pbar = (f1 + f2)/2
vtotal = sum(2 * bg^2 * pbar * (1 - pbar))
vwithin = sum(2*(bg^2)*f1*(1 - f1) + 2*(bg^2)*f2*(1 - f2))/2
vbetween = vtotal - vwithin

#Simulate local ancestry and genotypes for admixed population 
#such that the admixture fraction is 0.5.
lanc = replicate(nloci, rbinom(10000, 2, theta))
#simulate genotypes at each locus given the ancestry at that locus.
geno.admix = structure(sapply(lanc, function(x){
  if(x == 0){
    g1 = rbinom(1, 1, f1)
    g2 = rbinom(1, 1, f1)
  }
  if(x == 1){
    g1 = rbinom(1, 1, f1)
    g2 = rbinom(1, 1, f2)
  }
  if(x == 2){
    g1 = rbinom(1, 1, f2)
    g2 = rbinom(1, 1, f2)
  }
  return(g1 + g2)
}), dim = dim(lanc))

#Simulate genetic values in the admixed population.
gvalue.admix = t(t(bg)%*%t(geno.admix))

sample_size=length(gvalue.admix)
# add environment noise to gvalue to create the pheno data
# in order to get vg/(vg+ve) = 0.8, simulate ve with var=0.25
venv=args$venv
phenotype = as.data.table(gvalue.admix)
colnames(phenotype)=c("gvalue")
phenotype$environment = rnorm(n = sample_size, mean = 0,
                              sd = sqrt(venv))
# normalize gvalue
phenotype$gvalue_norm = (phenotype$gvalue - mean(phenotype$gvalue))/sd(phenotype$gvalue)

#add prs to each of the environmental effects
phenotype$pheno = phenotype$gvalue_norm + phenotype$environment

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
#cat("seed\t", seed, "\n")
cat("vtotal\t", vtotal, "\n")
cat("vwithin\t", vwithin, "\n")
cat("vbetween\t", vbetween, "\n")
cat("venv\t", venv, "\n")
cat("vgeno\t", vgeno, "\n")
sink()

# 
# # write the id file with gvalue to make phenotype file
# #pheno<-data.frame(FID="ADM",IID=1:10000,gvalue=c(gvalue.admix))
# #fwrite(pheno, file = paste0("pheno_gvalue",filename,".txt"), 
#        #sep = '\t' ,row.names = FALSE, col.names = FALSE)
# pheno<-data.frame(FID="ADM",IID=1:10000, phenovalue=phenotype$pheno)
# fwrite(pheno, file = paste0("pheno",filename,".txt"), 
#        sep = '\t' ,row.names = FALSE, col.names = FALSE)
# 
# # plot to check the expected vs estimated effect size
# est.b.g = matrix(NA, nrow = nloci, ncol = 1)
# est.b.l = matrix(NA, nrow = nloci, ncol = 1)
# 
# lanc<-as.matrix(lanc)
# geno.admix<-as.matrix(geno.admix)
# 
# for(i in 1:nloci){
#   l1 = lm(gvalue.admix ~ geno.admix[,i])
#   s1 = summary(l1)$coefficients
#   est.b.g[i,1] = s1[2,1]
#   
#   l2 = lm(gvalue.admix ~ lanc[,i])
#   s2 = summary(l2)$coefficients
#   est.b.l[i,1] = s2[2,1]
# }
# 
# 
# #Expected effect size of local ancestry (based on analytical derivation)
# bl = bg * (f2 - f1) 
# 
# #estimated and expected effect size of genotype
# png(file=paste0("effect_geno",filename,".png"), width=6, height=6, units="in", res=300)
# plot(bg, est.b.g, main="Expected VS Estimated Effect Size of Genotype",
#      sub = paste0("Admixture: ", theta, " Vgb: ", vgb ),
#      xlab="Expected Effect Size", ylab="Estimated Effect Size")
# dev.off()
# 
# 
# #estimate and expected effect size of ancestry
# png(file=paste0("effect_lanc", filename,".png"), width=6, height=6, units="in", res=300)
# plot(bl, est.b.l, main="Expected VS Estimated Effect Size of local ancestry", 
#      sub = paste0("Admixture: ", theta, " Vgb: ", vgb ),
#      xlab="Expected Effect Size", ylab="Estimated Effect Size")
# dev.off()


# convert admix lanc to tped 
# change lanc to certain allele
lanc<-as.data.table(lanc)
lanc[lanc == 0] <- "A,A"
lanc[lanc == 1] <- "A,T"
lanc[lanc == 2] <- "T,T"
# split them into 2 rows separately
lanc.allele<-data.frame(lapply(lanc, function(x) unlist(strsplit(as.character(x), ","))))

# transpose it and add map info: chrom, ID, cm (0), position
lanc.t<-t(lanc.allele)
map<-data.frame(chrom ="1",  cm=0, pos=1:nloci)
map$ID <- paste0(map$chrom, "_", map$pos, "_AT")

# reorder
map<-map[,c(1,4,2,3)]

# bind map and lanc together
l=cbind(map, lanc.t)

# output as tped file
write.table(l, file = paste0("admix", filename, ".tped"), quote=FALSE, sep = '\t' ,
            row.names = FALSE, col.names = FALSE)

# make tfam file
fam<-data.frame(FID="ADM",IID=1:10000,FID=0,MID=0,SEX=0,Pheno=-9 )
write.table(fam, file = paste0("admix", filename, ".tfam"), quote=FALSE, sep = '\t' ,
            row.names = FALSE, col.names = FALSE)


#convert genotype data to tped for gcta estimate. (espected h2=0.8)

# change geno to certain allele
geno.admix<-as.data.table(geno.admix)
geno.admix[geno.admix == 0] <- "A,A"
geno.admix[geno.admix == 1] <- "A,T"
geno.admix[geno.admix == 2] <- "T,T"
# split them into 2 rows separately
geno.allele<-data.frame(lapply(geno.admix, function(x) unlist(strsplit(as.character(x), ","))))

# transpose it and add map info: chrom, ID, cm (0), position
geno.t<-t(geno.allele)
map<-data.frame(chrom ="1",  cm=0, pos=1:nloci)
map$ID <- paste0(map$chrom, "_", map$pos, "_AT")

# reorder
map<-map[,c(1,4,2,3)]

# bind map and lanc together
l=cbind(map, geno.t)

# output as tped file
write.table(l, file = paste0("admix_geno", filename, ".tped"), quote=FALSE, sep = '\t' ,
            row.names = FALSE, col.names = FALSE)

# make tfam file
fam<-data.frame(FID="ADM",IID=1:10000,FID=0,MID=0,SEX=0,Pheno=-9 )
write.table(fam, file = paste0("admix_geno", filename, ".tfam"), quote=FALSE, sep = '\t' ,
            row.names = FALSE, col.names = FALSE)





