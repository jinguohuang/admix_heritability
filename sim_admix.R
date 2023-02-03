#!/usr/bin/env Rscript
# naive simulation to get expected vs observed value in parental groups
# assortative mating 
# update in multiple generations, add time 
# update parental freq draw so rare variants are avoided 
# update effect size so it's freq dependent
# output lanc of admixed pop and metapop (Parental)
# update with fixed sample size
# update with expected value of genetic variance due to lanc
# update with both HI and CGF
# update with output genetic variation due to genotype
# update with total admixture proportion, and one direction migration
# update with output lanc and pheno file prepared for GCTA
# record mean phenotypic difference between parental groups, 
# record corr(anc, PRS), corr(anc, PRS_lanc), var(anc), mean(anc)
# partition genetic variance Vg into 4 terms, calculate and record them
# modify seeds to control the beta sign distribution to get negative/positive cov 

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
LAnc2Plink <- function(lanc, nloci, samplesize, filename, FID=FID, prs){
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

################### Set up the parameters ########################


library(data.table)
library(optparse)

option_list = list( 
  make_option(c("--fst", "-F"), type="numeric", default = 0.2, 
              help="desired fst is between parental groups at causal loci [default= %default]", metavar="character"),
  make_option(c("--nloci", "-l"), type = "numeric", default = 1e3, 
              help = "number of loci  [default = %default]", metavar="character"),
  make_option(c("--psel", "-p"), type = "numeric", default = 0, 
              help = "selection strength [default = %default]", metavar="character"),
  make_option(c("--sample_size", "-n"), type = "numeric", default = 1e3, 
              help = "sample size of admixed people [default = %default]", metavar="character"),
  make_option(c("--seed"), type="numeric", default=120, 
              help="random seed [default= %default]", metavar="character"),
  make_option(c("--pganc", "-P"), type = "numeric", default = 0.9, 
              help = "assortative mating strength [default = %default]", metavar="character"),
  make_option(c("--gen", "-t"), type = "numeric", default = 1, 
              help = "number of generations since admixture [default = %default]", metavar="character"),
  make_option(c("--theta"), type = "numeric", default = 0.5, 
              help = "total admixture proportion from pop1, range 0 to 1, [default = %default]", metavar="character"),
  make_option(c("--model", "-M"), default = "HI", 
              help = "admixture model: HI or CGF [default = %default]", metavar="character"),
  make_option(c("--cov", "-C"), default = "pos", 
              help = "covariance sign: pos or neg [default = %default]", metavar="character")
#  make_option(c("--m2"), type = "numeric", default = 0, 
#              help = "migration rate from pop2 [default = %default]", metavar="character")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
print("setting up")

fst = opt$fst # fst 
nloci = opt$nloci #number of loci
p_sel = opt$psel #selection strength
n = opt$sample_size
seed = opt$seed #seed for admix pop sample
P = opt$pganc
t = opt$gen
#m1 = opt$m1
#m2 = opt$m2
m2 = 0 # one directional CFG
theta = opt$theta
model = opt$model
cov = opt$cov

seed1= 26 #seed for sampling freq and effect size
seed2_pos=73
seed2_neg=78
seed3=seed #seed for sampling admixed pop

# if HI, m1=0
# pi is the value of nindiv_pop1, at t=0, nindiv_pop1
if( model == "HI") {
    m1=0
    pi=theta 
} else { #CGF model
    m1=CGFm(t, theta)
    pi=m1 
}

print(paste0("m1= ", m1))


# if pos, seed2
if( cov == "pos") {
    seed2=seed2_pos
} else { #neg
    seed2=seed2_neg
}

print(paste0("cov is ", cov))


nindiv_pop1 = round(n*pi)
nindiv_pop2 = n - nindiv_pop1

################### simulate parental groups genotype  ########################
print("simulate parental groups genotype")
set.seed(seed1) # keep this part same across P

# Draw allele frequency in an ancestral population prior to the divergence 
# and drawn from a uniform distribution unif(0.001,0.999)
anc.frq <- runif(nloci, min=1e-3, max=0.999)

#set freq per ancestral population at each locus, given fst and human ancestral allele maf
f1 <- rep(NA, nloci)
f2 <- rep(NA, nloci)
for (s in 1:nloci) {
  # set allele frequencies for ancestries
  p.pop1 <- 0
  p.pop2 <- 0
  anc.frq <- runif(1, min=1e-3, max=0.999)
  #Balding-Nichol model, also make sure 0.01<frq<0.99, so rare variants are avoided
  while((p.pop1 <= 0.01 || p.pop1 >= 0.99) || p.pop2 <= 0.01 || p.pop2 >= 0.99){
    p.pop1 <- rbeta(1, anc.frq*(1-fst)/fst,(1-anc.frq)*(1-fst)/fst)
    p.pop2 <- rbeta(1, anc.frq*(1-fst)/fst,(1-anc.frq)*(1-fst)/fst)
  }
  f1[s] <- p.pop1
  f2[s] <- p.pop2
}
#fbar = (f1+f2)/2
fbar = pi*f1 + (1-pi)*f2
fdiff = f1 - f2 #freq difference

# simulate parental genotype based on frequency
geno1 <- matrix(rbinom(nloci*nindiv_pop1, 2, f1), byrow=T, nrow=nindiv_pop1, ncol=nloci)
geno2 <- matrix(rbinom(nloci*nindiv_pop2, 2, f2), byrow=T, nrow=nindiv_pop2, ncol=nloci)

###################  simulate effect size  ########################
print("simulate effect size")

# Choose effect sizes from the standard normal distribution N(0, 1/2*f*(1-f))
bg = matrix(
  rnorm(nloci, 0, sd = sqrt(1/(2*nloci*fbar*(1-fbar)))), 
  nrow = nloci, ncol = 1)

#standardize to remove stochasticity in realized values
bg = ((bg - mean(bg)) / sd(bg)) * sqrt(1/(2*nloci*fbar*(1-fbar)))

set.seed(seed2)
bg_sign=sample(c(-1,1), nloci, replace=T) #for controling beta sign disribution

# For some fraction p of the loci, choose the effect direction to be such that 
# the allele that is more frequent in population 1 has a positive effect. 
# The rest are random. 
# weak selection, use prob f1/(f1+f2) decide the sign of selected loci
# strong selection, use f1-f2 decide the sign of selected loci
if (p_sel == 0){ 
  bg_sel = bg* bg_sign
  }else{
    sel_index <- sample(c(1:nloci), size=round(nloci*abs(p_sel)), replace=F) # sample selected loci
    bg_sel = matrix(NA, nrow = nloci, ncol = 1) 
    # selected index, change sign
    if (p_sel > 0){ 
      # weak selection
      bg_sel[sel_index] = abs(bg[sel_index]) * ((rbinom(nloci, 1, f1/(f1+f2))-0.5)*2)[sel_index]
      # strong selection
      #bg_sel[sel_index] = abs(bg[sel_index]) * (sign(fdiff))[sel_index]
      }else{
        bg_sel[sel_index] = abs(bg[sel_index]) * ((rbinom(nloci, 1, f1/(f1+f2))-0.5)*2)[sel_index] * (-1)
        #bg_sel[sel_index] = abs(bg[sel_index]) * (sign(fdiff))[sel_index] * (-1)
        }
      bg_sel[-sel_index] = bg[-sel_index]
}

# effect size of local ancestry
bg_lanc = bg_sel * fdiff

###################  Calculate expected genetic variance in parental ########################
print("Calculate expected genetic variance in parental")
# calculate expected value in parental groups
vtotal.genic = sum(2* bg_sel^2 * fbar * (1-fbar))
vbetween.genic = sum( bg_sel^2 * (f1-f2)^2) * 2 * pi * (1-pi)
vwithin.exp = sum(pi*2*(bg_sel^2)*f1*(1 - f1) + (1-pi)*2*(bg_sel^2)*f2*(1 - f2))

# calculate term 1 coefficient only HI
#varZ.sum1 = sum(bg_sel^2 * fdiff^2)
#varZ.coef1 = vbetween.genic

# calculate term 2 coefficient only HI
#sum2 = matrix(, nrow = nloci, ncol = nloci)
#for (i in 1:nloci){
#  for (j in 1:nloci){
#    sum2[i,j]=bg_sel[i]*bg_sel[j]*fdiff[i]*fdiff[j]
#  }
#}
# make diagnal 0
#diag(sum2) = 0
# sum it up
#varZ.coef2 = sum(sum2)

###################  Calculate observed genetic variance in parental ########################
print("Calculate observed genetic variance in parental")
# prs on causal variants
pop1_prs <- apply(geno1, MARGIN=1, FUN=indiv_prs, beta=bg_sel)
pop2_prs <- apply(geno2, MARGIN=1, FUN=indiv_prs, beta=bg_sel)

# calculate observed value
prs.mean = mean(c(pop1_prs, pop2_prs))
prs.diff = mean(pop1_prs) - mean(pop2_prs) # mean prs difference
meta_prs = c(pop1_prs, pop2_prs) #check if this equals to t=0 adm_prs
#write.table(meta_prs, file = paste0("metaprs", "_seed", seed, "_P", P, ".txt"), 
#              quote=FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
vtotal.obs = var(meta_prs)

# calculate variance within each groups
g1.var = var(pop1_prs)
g2.var = var(pop2_prs)
vwithin.obs = (g1.var + g2.var)/2
# calculate variance between populations
vbetween.obs = ((mean(pop1_prs) - prs.mean)^2 + (mean(pop2_prs) - prs.mean)^2)/2

################### simulate admixed population under assortative mating + HI/CGF ########################
print("simulate admixed population")
# output lanc of each generation
# simulate global ancestry of parental meta groups of each generation
ganc_meta = matrix(NA, nrow = t+1, ncol = n)
ganc_meta[1,] = c(rep(1, nindiv_pop1), rep(0, nindiv_pop2)) 
# get parental global ancestry combination and calculate the ganc of admix pop
#corrthresh = 0.01 #correlation threshold
corrthresh = 0.0123 # for theta=0.5 P=0.6 seed1
for (i in 1:t){ #generation of admixture
  print(paste0("finding mating pairs for generation ", i))
  #u = 100
  u = 100 # for P=0.6 seed1
  l = 0
  c = -2*corrthresh
  iter = 1
  #clist = c
  # search for mate pairs to meet the P
  while(abs(c-P) > corrthresh){
    iter = iter + 1
    #s = (u+l)/2
    s = (u+l)/2 #for P=0.6 seed1
    #print(sprintf('c = %f s = %f P = %f iter = %d', c, s, P, iter))
    #sort the score and get the index of the sorted list
    idxm = sort(ganc_meta[i,]  + rnorm(n, 0, 1)*s, index.return=TRUE)$ix 
    idxf = sort(ganc_meta[i,]  + rnorm(n, 0, 1)*s, index.return=TRUE)$ix
    # correlation of ganc of mat and fat
    c = cor(ganc_meta[i,idxf], ganc_meta[i,idxm]) # get correlation of ganc
    #clist[iter] = c # store the correlation
    if(c < P){ #adjust s to smaller value if c is smaller 
      u = s
    }else{ #adjust s to larger value if c is larger 
        l = s
        }
  }
  
  # ganc with migration
  # indidivual to be replaced
  migration1 = sample(1:n, size=round((m1+m2)*n), replace=F) 
  # CGF or HI
  if(m1+m2>0){ # confinuous geneflow
    migration2 = sample(migration1, size=round(m1*length(migration1)/(m1+m2)), replace=F)
    # get index of replaced ind of pop1
  }
  else{ #hybrid isolation
    migration2 = sample(migration1, size=0, replace=F)
    # no replacement in HI model
  }

  tmp = rep(0, n) 
  tmp[migration1] = 1 # store index of migration of pop1+pop2
  tmp1 = rep(0, n) 
  tmp1[migration2] = 1 # store index of migration of pop1
  migration1 = tmp1 # index of migrated pop1
  migration2 = tmp-migration1 # index of migrated pop2

  # update the global ancestry pool in the metapop each round with migration
  globalanc = (ganc_meta[i,idxf] + ganc_meta[i,idxm])/2
  globalanc = migration1 + (1-migration1)*globalanc # mig from pop1 ganc=1
  globalanc = (1-migration2)*globalanc # mig from pop2 ganc=0
  ganc_meta[i+1,] = globalanc
}


# output the variance of PRS due to lanc of each time
# output the expected value of genetic variance due to ancestry
# read the expected LAD
#explad = read.table("explad.txt", header = T, sep = ",")
#explad$t = seq(0,20,1) # add time
# get the column according to the P 
#Pnum=P*10
#explad.col=paste0("P",Pnum)
#explad.t=explad[explad.col]

for (i in 1:(t+1)){
  print(paste0("output for generation ", i-1))
  ganc <- ganc_meta[i,]
  # record var(theta) at each generation
  var.ganc <- var(ganc)
  set.seed(seed3) #change this seed for different replication
  # local ancestry per marker: 0/1 copies of pop1 ancestry
  # remove the byrow=T
  lanc_p <- matrix(rbinom(n*nloci, 1, ganc), nrow=n, ncol=nloci)
  lanc_m <- matrix(rbinom(n*nloci, 1, ganc), nrow=n, ncol=nloci)
  lanc_adm <- lanc_p + lanc_m
  #admixed pop genotype 
  geno_adm <- t(sapply(c(1:n), FUN=LAnc2Geno, lanc_p=lanc_p, lanc_m=lanc_m, frq_pop1=f1, frq_pop2=f2))
  # prs variance of admixed pop
  adm_prs <- apply(geno_adm, MARGIN=1, FUN=indiv_prs, beta=bg_sel)
  var.prs.geno <- var(adm_prs)
  # prs variance due to lanc of admixed pop
  prs_lanc_adm <- apply(lanc_adm, MARGIN=1, FUN=indiv_prs, beta=bg_lanc)
  var.prs.lanc <- var(prs_lanc_adm)
  # record correlation between ancestry and trait
  cor.ganc.prsgeno <- cor(ganc, adm_prs)
  cor.ganc.prslanc <- cor(ganc, prs_lanc_adm)
  # calculate dipLAD(t) 
  # calculate covariance pairwise loci 
  cov_mat=cov(lanc_adm)
  #cov_mat_out=cov_mat
  #dipLAD2=sum(cov_mat)/(nloci*nloci)
  diag(cov_mat)=0
  dipLAD=sum(cov_mat)/(nloci*(nloci-1))
  #get LAD from Zaitlen 
  #dipLAD1=4*explad.t[i,]
  # calculate the varZ(t) with the dipLAD
  #varZ.t = varZ.coef1*(1+((1+P)/2)^(i-1)) + varZ.coef2*dipLAD
  #varZ.t1 = varZ.coef1*(1+((1+P)/2)^(i-1)) + varZ.coef2*dipLAD1
  #output the lanc of each generation (including metapop, t=0)

  # Calculate 4 terms of Vg and output them with exp.ganc, var.ganc, var.PRS, avg.g1-avg.g2
  exp.ganc=mean(ganc)
  # Calculate term 1 of Vg
  vg.term1=vg_term1(exp.ganc, f1, f2, beta=bg_sel)
  # Calculate term 2 of Vg
  vg.term2=vg_term2(exp.ganc, f1, f2, beta=bg_sel)
  # Calculate term 3 of Vg
  vg.term3=vg_term3(var.ganc, f1, f2, beta=bg_sel)
  # Calculate term 4 of Vg
  vg.term4=vg_term4(var.ganc, f1, f2, beta=bg_sel)
  # Calculate Vg with 4 terms
  vg.sum=sum(vg.term1, vg.term2, vg.term3, vg.term4)

  # Output filename
  filename <- paste0("admix_", model,"_theta", theta, "_gen", t, "_P", P, "_", opt$cov, "_seed", seed,"_t", i-1) #for filename
  #LAnc2Plink(lanc_adm, nloci, n, filename, FID="ADM", adm_prs) #prs_geno as pheno output
  #LAnc2Plink(lanc_adm, nloci, n, filename, FID="ADM", prs_lanc_adm) #prs_lanc as pheno output
  # write to file to record exp. and obs. value 
  sink(file=paste0("summary_",filename, ".txt"), type = "output")
  cat("vbetween.genic \t vbetween.obs \t P \t t \t seed \t psel \t pwCov \t Var.prs.lanc \t Var.prs.geno \t m1 \t prs.diff \t exp.ganc \t var.ganc \t cor.ganc.prsgeno \t cor.ganc.prslanc \t vg.term1 \t vg.term2 \t vg.term3 \t vg.term4 \t vg.sum \n")
  cat(vbetween.genic, "\t", vbetween.obs, "\t", P, "\t", i-1, "\t", seed, "\t", p_sel, "\t", dipLAD, "\t", var.prs.lanc, "\t", var.prs.geno, "\t", m1, "\t", prs.diff, "\t", exp.ganc, "\t", var.ganc, "\t", cor.ganc.prsgeno, "\t", cor.ganc.prslanc, "\t", vg.term1, "\t", vg.term2, "\t", vg.term3, "\t", vg.term4, "\t", vg.sum, "\n")
  sink()  
}

#write.table(cov_mat_out, file = paste0(filename, ".covmat"), 
#              quote=FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)








