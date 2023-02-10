#!/usr/bin/env Rscript
# functions to call

####################  Functions to be used ########################

# set genotypes for an admixed person from binomial sampling based on freq
LAnc2Geno <- function(index,lanc_p,lanc_m, frq_pop1, frq_pop2){
  hap_p <- rep(NA, length(lanc_p[index,]))
  hap_m <- rep(NA, length(lanc_m[index,]))
  #paternal hap
  hap_p[lanc_p[index,]==1] <- rbinom(sum(lanc_p[index,]==1), 1, frq_pop1[lanc_p[index,]==1])
  hap_p[lanc_p[index,]==0] <- rbinom(sum(lanc_p[index,]==0), 1, frq_pop2[lanc_p[index,]==0])
  #maternal hap
  hap_m[lanc_m[index,]==1] <- rbinom(sum(lanc_m[index,]==1), 1, frq_pop1[lanc_m[index,]==1])
  hap_m[lanc_m[index,]==0] <- rbinom(sum(lanc_m[index,]==0), 1, frq_pop2[lanc_m[index,]==0])
  
  hap <- hap_p + hap_m
  return(hap)
}

# prs on causal variants
indiv_prs <- function(geno, beta=beta) return (sum(geno*beta))

# function to output local ancestry to plink file
LAnc2Plink <- function(lanc, nloci, samplesize, filename, FID=FID, prs, ganc){
  # transpose and add map
  lanc.t<-t(lanc)
  map = data.frame(SNP = paste0(1, "_", 1:nloci, "_AG"), A1 = "A", A2 = "G")
  lanc.t = cbind(map, lanc.t)
  write.table(lanc.t, file = paste0(filename, ".dosage"), 
              quote=FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
  # make tfam file
  fam<-data.frame(FID=FID,IID=1:samplesize, FID=0, MID=0, SEX=0, Pheno=-9 )
  write.table(fam, file = paste0(filename, ".tfam"), 
              quote=FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
  phenotype = cbind(fam[,c(1,2)], prs)
  write.table(phenotype, file = paste0(filename, ".pheno"), 
              quote=FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
  qcovar = cbind(fam[,c(1,2)], ganc)
  write.table(qcovar, file = paste0(filename, ".qcovar"), 
              quote=FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
}

# calculate the migration rate from pop1 at each generation
CGFm = function(t, theta){
  theta2=1-theta #pop2 proportion
  m1 = 1- theta2^(1/t)
  return(m1)
}

# Calculate 4 terms of genetic variance Vg
# term 1 of Vg
vg_term1 <- function(exp.ganc, f1, f2, beta) return (sum(beta^2*(2*exp.ganc*f1*(1-f1)+2*(1-exp.ganc)*f2*(1-f2))))

# term 2 of Vg
vg_term2 <- function(exp.ganc, f1, f2, beta) return (sum(beta^2*(2*exp.ganc*(1-exp.ganc)*((f1-f2)^2))))

# term 3 of Vg
vg_term3 <- function(var.ganc, f1, f2, beta) return (sum(beta^2*(2*var.ganc*((f1-f2)^2))))

# term 4 of Vg
vg_term4 <- function(var.ganc, f1, f2, beta){
  term4=matrix(, nrow = nloci, ncol = nloci)
  for (i in 1:nloci){
    for (j in 1:nloci){
      term4[i,j]=beta[i]*beta[j]*(f1[i]-f2[i])*(f1[j]-f2[j])
    }
  }
  # make diagnal 0
  diag(term4)=0
  # sum it up
  sum.term4 = sum(term4)*4*var.ganc
  return (sum.term4)} 

# calculate variance, covariance and correlation
pvar  = function(x){
  mean((x - mean(x))^2)
}

pcov = function(x,y){
  mean((x - mean(x))*(y - mean(y)))
}

pcor = function(x, y){
  pcov(x, y)/(sqrt(pvar(x))*sqrt(pvar(y)))
}