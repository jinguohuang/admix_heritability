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
# split the simulation files into parental simulation and admixed simulation
# output the freq and beta of parental groups and global ancestry of each generation
# update fbar with (f1+f2)/2 to keep vg same across CGF HI and different theta
# output summary file with right dir

##########initializing##########
source("main_admix.R")

suppressWarnings(suppressMessages({
  library(data.table)
  library(optparse)
}))

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

seed1= 99 #seed for sampling freq and effect size
seed2_pos=65
seed2_neg=84
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

##########simulate parental groups genotype##########

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
fbar = (f1+f2)/2
#fbar = pi*f1 + (1-pi)*f2
fdiff = f1 - f2 #freq difference

# simulate parental genotype based on frequency
geno1 <- matrix(rbinom(nloci*n, 2, f1), byrow=T, nrow=n, ncol=nloci)
geno2 <- matrix(rbinom(nloci*n, 2, f2), byrow=T, nrow=n, ncol=nloci)

##########simulate effect size##########

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

##########Calculate expected genetic variance in parental##########

print("Calculate expected genetic variance in parental")
# calculate expected value in parental groups (1/2 intead of pi)
vtotal.genic = sum(2* bg_sel^2 * fbar * (1-fbar))
vbetween.genic = sum( bg_sel^2 * (f1-f2)^2) * 2 * 0.5 * 0.5
vwithin.exp = sum(0.5*2*(bg_sel^2)*f1*(1 - f1) + 0.5*2*(bg_sel^2)*f2*(1 - f2))

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

##########Calculate observed genetic variance in parental##########
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
vtotal.obs = pvar(meta_prs)

# calculate variance within each groups
g1.var = pvar(pop1_prs)
g2.var = pvar(pop2_prs)
vwithin.obs = (g1.var + g2.var)/2
# calculate variance between populations
vbetween.obs = ((mean(pop1_prs) - prs.mean)^2 + (mean(pop2_prs) - prs.mean)^2)/2

################### simulate admixed population under assortative mating + HI/CGF ########################
print("simulate admixed population")
# output lanc of each generation

# inital number of indvidual from each parental groups
nindiv_pop1 = round(n*pi)
nindiv_pop2 = n - nindiv_pop1

# simulate global ancestry of parental meta groups of each generation
ganc_meta = matrix(NA, nrow = t+1, ncol = n)
ganc_meta[1,] = c(rep(1, nindiv_pop1), rep(0, nindiv_pop2)) 
# get parental global ancestry combination and calculate the ganc of admix pop
#corrthresh = 0.01 #correlation threshold
corrthresh = 0.01
set.seed(seed3) # so we have different var.theta
for (i in 1:t){ #generation of admixture
  print(paste0("finding mating pairs for generation ", i))
  u = 100
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
  set.seed(seed3*3)
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


# simulate admixed pop lanc from ganc
# update filename with seed
dirname <- paste0("../data/theta", theta, "_gen", t, "/summary", "/", model, "/")
filename <- paste0("admix_", model,"_theta", theta, "_gen", t, "_P", P, "_", cov, "_seed", seed)
# record summary stat of each generation
output = matrix(NA, nrow = t+1, ncol = 17)

for (i in 1:(t+1)){
  print(paste0("output for generation ", i-1))
  ganc <- ganc_meta[i,]
  
  set.seed(seed3*2) #change this seed for different replication
  # sample local ancestry given global ancestry
  # local ancestry per marker: 0/1 copies of pop1 ancestry
  # remove the byrow=T
  lanc_p <- matrix(rbinom(n*nloci, 1, ganc), nrow=n, ncol=nloci)
  lanc_m <- matrix(rbinom(n*nloci, 1, ganc), nrow=n, ncol=nloci)
  lanc_adm <- lanc_p + lanc_m 
  #admixed pop genotype 
  geno_adm <- t(sapply(c(1:n), FUN=LAnc2Geno, lanc_p=lanc_p, lanc_m=lanc_m, frq_pop1=f1, frq_pop2=f2))
  # theta: global ancestry of each individual
  #admixed pop global ancestry calculated by lanc
  ganc_adm <- apply(lanc_adm, 1, sum)/(2*nloci) #1 cause sum by row-ind
  #dim(lanc_adm) #10000  1000
  # record theta at each generation, this theta is calculated from lanc
  mean.theta <- mean(ganc_adm)
  var.theta <- pvar(ganc_adm)
  
  # trait value: PRS_geno and PRS_lanc
  # prs variance of admixed pop
  prs_geno <- apply(geno_adm, MARGIN=1, FUN=indiv_prs, beta=bg_sel)
  var.prs.geno <- pvar(prs_geno)
  
  # add noise to pheno such that h2=0.8
  noise <- rnorm(n, mean=0, sd=sqrt(var.prs.geno/4))
  pheno <- prs_geno + noise
  var.pheno <- pvar(pheno)
    
  # prs variance due to lanc of admixed pop
  prs_lanc <- apply(lanc_adm, MARGIN=1, FUN=indiv_prs, beta=bg_lanc)
  var.prs.lanc <- pvar(prs_lanc)
  
  # record covariance and correlation between ancestry and trait
  cov.ganc.prsgeno <- pcov(ganc_adm, prs_geno)
  cov.ganc.pheno <- pcov(ganc_adm, pheno)
  cov.ganc.prslanc <- pcov(ganc_adm, prs_lanc)
  cor.ganc.prsgeno <- pcor(ganc_adm, prs_geno)
  cor.ganc.pheno <- pcor(ganc_adm, pheno)
  cor.ganc.prslanc <- pcor(ganc_adm, prs_lanc)
  
  # expected Vg terms
  # Calculate term 1 of Vg
  vg.term1=vg_term1(mean.theta, f1, f2, beta=bg_sel)
  # Calculate term 2 of Vg
  vg.term2=vg_term2(mean.theta, f1, f2, beta=bg_sel)
  # Calculate term 3 of Vg
  vg.term3=vg_term3(var.theta, f1, f2, beta=bg_sel)
  # Calculate term 4 of Vg
  vg.term4=vg_term4(var.theta, f1, f2, beta=bg_sel)
  # Calculate Vg with 4 terms
  vg.sum=sum(vg.term1, vg.term2, vg.term3, vg.term4)
  
  # record to output table
  output[i, 1] <- i-1 #ith generation
  output[i, 2] <- mean.theta
  output[i, 3] <- var.theta
  output[i, 4] <- var.prs.geno # observed vg
  output[i, 5] <- cov.ganc.prsgeno 
  output[i, 6] <- cor.ganc.prsgeno
  output[i, 7] <- vg.term1
  output[i, 8] <- vg.term2
  output[i, 9] <- vg.term3
  output[i, 10] <- vg.term4
  output[i, 11] <- vg.sum # expected vg
  output[i, 12] <- var.prs.lanc # observed vgamma
  output[i, 13] <- cov.ganc.prslanc
  output[i, 14] <- cor.ganc.prslanc
  output[i, 15] <- var.pheno
  output[i, 16] <- cov.ganc.pheno 
  output[i, 17] <- cor.ganc.pheno
  
  # output PLINK files: dosage, fam, pheno, covar (global ancestry)
  #plinkname <- paste0(filename, "_t", i-1) #for filename
  #plinkname_lanc <- paste0(filename, "_t", i-1, "_lanc")
  # check the genotype output
  #Geno2Plink(geno_adm, nloci, n, plinkname, FID="ADM", prs_geno, pheno, ganc_adm) #prs_geno as pheno output
  #Geno2Plink(lanc_adm, nloci, n, plinkname_lanc, FID="ADM", prs_lanc, pheno, ganc_adm) #prs_lanc as pheno output
}


# add colnames
colnames(output) <- c("t", "mean.theta", "var.theta", "var.prs.geno", "cov.ganc.prsgeno", "cor.ganc.prsgeno", "vg.term1", "vg.term2", "vg.term3", "vg.term4", "vg.sum", "var.prs.lanc", "cov.ganc.prslanc", "cor.ganc.prslanc", "var.pheno", "cov.ganc.pheno", "cor.ganc.pheno")

# add common features: theta, P, seed
output<-as.data.table(output)
output$theta <- theta
output$P <- P
output$seed <- seed
output$vbetween.genic <- vbetween.genic
output$vbetween.obs <- vbetween.obs
output$cov <- cov
output$model <- model
output$m1 <- m1
output$prs.diff <-prs.diff

# write the output file
write.table(output,
        paste0(dirname, filename,".summary.txt",sep=""),
         col.names = TRUE,
         row.names=FALSE,
         quote=FALSE,
         sep="\t")

