#!/usr/bin/env Rscript
# naive simulation to get expected vs observed value in parental groups
# assortative mating 
# update in multiple generations, add time 
# update parental freq draw so rare variants are avoided 
# update effect size so it's freq dependent
# output lanc of admixed pop 
# output lanc of metapop (Parental)
#start.time <- Sys.time()

####################  Functions to be used ########################

# some from https://github.com/menglin44/APRICOT/blob/main/R/load_main_geno.R

# global anc proprotions from beta distribution
GAncBeta <- function(mean, var, samplesize){
  alpha <- mean^2*((1-mean)/var - 1/mean)
  beta <- (1-mean)/mean*alpha
  ganc <- rbeta(samplesize, alpha, beta)
  return(ganc)
}

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

# calc af
af <- function(genos) sum(genos)/(2*length(genos))

# calc maf
maf <- function(genos) 0.5-abs(0.5-af(genos))


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

################### Set up the parameters ########################

library(data.table)
library(optparse)

option_list = list( 
  make_option(c("--fst", "-F"), type="numeric", default = 0.2, 
              help="desired fst is between parental groups at causal loci [default= %default]", metavar="character"),
#  make_option(c("--theta", "-t"), type="numeric", default=0.5, 
#              help="admixture proportion [default= %default]", metavar="character"),
  make_option(c("--nloci", "-l"), type = "numeric", default = 1e3, 
              help = "number of loci  [default = %default]", metavar="character"),
  make_option(c("--psel", "-p"), type = "numeric", default = 0, 
              help = "selection strength [default = %default]", metavar="character"),
  make_option(c("--sample_size", "-n"), type = "numeric", default = 1e3, 
              help = "sample size of admixed people [default = %default]", metavar="character"),
  make_option(c("--seed"), type="numeric", default=120, 
              help="random seed [default= %default]", metavar="character"),
  make_option(c("--pganc", "-P"), type = "numeric", default = 0, 
              help = "assortative mating strength [default = %default]", metavar="character"),
  make_option(c("--gen", "-t"), type = "numeric", default = 1, 
              help = "number of generations since admixture [default = %default]", metavar="character")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
print("setting up")

fst = opt$fst # fst 
nloci = opt$nloci #number of loci
p_sel = opt$psel #selection strength
n1_indv = opt$sample_size #sample size in pop1
n2_indv = opt$sample_size
n = opt$sample_size
adm_indv = opt$sample_size #sample size of admixed pop
#theta = opt$theta #average global ancestry of admixed pop
#sd_ganc = 0.1 #sd of global ancestry
seed = opt$seed
P = opt$pganc
t = opt$gen
filename <- paste0("admix", "_fst", fst, "_psel", p_sel, "_seed", seed, "_P", P, "_t", t) #for filename
set.seed(seed)

################### simulate parental groups genotype  ########################
print("simulate parental groups genotype")
# Draw allele frequency in an ancestral population prior to the divergence 
#and drawn from a uniform distribution unif(0.001,0.999)
anc.frq <- runif(nloci, min=1e-3, max=0.999)

# draw allele frequencies in the two parental pops from a beta distribution 
#f1 = rbeta(nloci, anc.frq*(1-fst)/fst,(1-anc.frq)*(1-fst)/fst)
#f2 = rbeta(nloci, anc.frq*(1-fst)/fst,(1-anc.frq)*(1-fst)/fst)

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
fdiff = f1 - f2 #freq difference

# simulate parental genotype based on frequency
geno1 <- matrix(rbinom(nloci*n1_indv, 2, f1), byrow=T, nrow=n1_indv, ncol=nloci)
geno2 <- matrix(rbinom(nloci*n2_indv, 2, f2), byrow=T, nrow=n2_indv, ncol=nloci)

# simulate parental lanc, local ancestry per marker: 0/1/2 copies of pop1 ancestry
# pop1 all 2, pop2 all 0
lanc1 <- matrix(rep(2, nloci*n1_indv), byrow=T, nrow=n1_indv, ncol=nloci)
lanc2 <- matrix(rep(0, nloci*n2_indv), byrow=T, nrow=n2_indv, ncol=nloci)

# initial lanc of metapop
lanc_meta = rbind(lanc1, lanc2)

# simulate global ancestry of parental meta groups
# mother or father could from either pop
# here we could adjust the sample size to adjust theta
ganc_meta = c(rep(1, n1_indv), rep(0, n2_indv)) 

###################  simulate effect size  ########################
print("simulate effect size")
# Sample effect sizes such that beta follow normal distribution, 
# with variance of 1/2nfbar(1-fbar)

#favg = sum(f1,f2)/(2*nloci) # calculate average allele freq

# Choose effect sizes from the standard normal distribution N(0, 1/2*f*(1-f))
bg = matrix(
  rnorm(nloci, 0, sd = sqrt(1/(2*nloci*fbar*(1-fbar)))), 
  nrow = nloci, ncol = 1)


#standardize to remove stochasticity in realized values
bg = ((bg - mean(bg)) / sd(bg)) * sqrt(1/(2*nloci*fbar*(1-fbar)))

# For some fraction p of the loci, choose the effect direction to be such that 
# the allele that is more frequent in population 1 has a positive effect. 
# The rest are random. 
# weak selection, use prob f1/(f1+f2) decide the sign of selected loci
# strong selection, use f1-f2 decide the sign of selected loci

if (p_sel == 0){ 
  bg_sel = bg 
  }else{
  sel_index <- sample(c(1:nloci), size=round(nloci*p_sel), replace=F) # sample selected loci
  bg_sel = matrix(NA, nrow = nloci, ncol = 1) 
  # selected index, change sign
  # weak selection
  #bg_sel[sel_index] = abs(bg[sel_index]) * ((rbinom(nloci, 1, f1/(f1+f2))-0.5)*2)[sel_index]
  # strong selection
  bg_sel[sel_index] = abs(bg[sel_index]) * (sign(fdiff))[sel_index]
  bg_sel[-sel_index] = bg[-sel_index]
}


###################  Calculate expected genetic variance in parental ########################
print("Calculate expected genetic variance in parental")
#geno = rbind(geno1, geno2) #metapop genotype
#claculate pairwise covariance across loci in metapop
#cv = matrix(NA, nrow = nloci, ncol = nloci)
#cv_cov = matrix(NA, nrow = nloci, ncol = nloci)
#for(i in 1:nloci){
#  for(j in 1:nloci){
#    gi = geno[,i]
#    gj = geno[,j]
#    cv[i,j] = bg_sel[i]*bg_sel[j]*cov(gi, gj)
    #cv_cov[i,j] = cov(gi, gj)
#  }
#}
#diag(cv)<-0 # make diagnol 0 
#vtotal.cov=sum(cv) # covariance part

# calculate expected value in parental groups
#fbar = (f1+f2)/2 
vtotal.genic = sum(2* bg_sel^2 * fbar * (1-fbar))
vbetween.genic = sum( bg_sel^2 * (f1-f2)^2) /2
vwithin.exp = sum(2*(bg_sel^2)*f1*(1 - f1) + 2*(bg_sel^2)*f2*(1 - f2))/2
#vtotal.exp = vtotal.genic + vtotal.cov
#vbetween.exp = vbetween.genic + vtotal.cov


###################  Calculate observed genetic variance in parental ########################
print("Calculate observed genetic variance in parental")
# prs on causal variants
pop1_prs <- apply(geno1, MARGIN=1, FUN=indiv_prs, beta=bg_sel)
pop2_prs <- apply(geno2, MARGIN=1, FUN=indiv_prs, beta=bg_sel)

# calculate observed value
prs.mean = mean(c(pop1_prs, pop2_prs))
prs.diff = mean(pop1_prs) - mean(pop2_prs) # mean prs difference
meta_prs = c(pop1_prs, pop2_prs)
vtotal.obs = var(meta_prs)

# calculate variance within each groups
g1.var = var(pop1_prs)
g2.var = var(pop2_prs)
vwithin.obs = (g1.var + g2.var)/2
# calculate variance between populations
vbetween.obs = ((mean(pop1_prs) - prs.mean)^2 + (mean(pop2_prs) - prs.mean)^2)/2


################### simulate admixed population  ########################
print("simulate admixed population")
# Simulate one-pulse admixture (hybrid isolation model) between the two populations 
# and measure local ancestry genetic variance for each case.
# simulate global ancestry of admix pop
#ganc <- GAncBeta(theta, sd_ganc^2, adm_indv) #neutral

# get parental global ancestry combination 
# and calculate the ganc of admix pop
corrthresh = 0.01 #correlation threshold
for (i in 1:t){ #generation of admixture
  u = 100
  l = 0
  c = -2*corrthresh
  iter = 1
  clist = c
  # search for mate pairs to meet the P
  while(abs(c-P) > corrthresh){
      iter = iter + 1
      s = (u+l)/2
      #print(sprintf('c = %f s = %f P = %f iter = %d', c, s, P, iter))
      #sort the score and get the index of the sorted list
      idxm = sort(ganc_meta  + rnorm(n, 0, 1)*s, index.return=TRUE)$ix 
      idxf = sort(ganc_meta  + rnorm(n, 0, 1)*s, index.return=TRUE)$ix
      # correlation of ganc of mat and fat
      c = cor(ganc_meta[idxf],ganc_meta[idxm]) # get correlation of ganc
      clist[iter] = c # store the correlation
      if(c < P){ #adjust s to smaller value if c is smaller 
        u = s
      }else{ #adjust s to larger value if c is larger 
          l = s
          }
      #if(iter%%100 == 0){
      #  print(sprintf('stuck with %d iterations\n', iter))
      #  }
  }
  # calculate the global ancestry with the mate pairs
  # update the global ancestry pool in the metapop each round
  ganc_meta = (ganc_meta[idxf] + ganc_meta[idxm])/2
}


 # track the mean and var of ganc_meta each generation
ganc_mean = mean(ganc_meta)
ganc_var = var(ganc_meta)
ganc = ganc_meta

#local ancestry per marker: 0/1/2 copies of pop1 ancestry
lanc_p <- matrix(rbinom(adm_indv*nloci, 1, ganc), nrow=adm_indv, ncol=nloci)
lanc_m <- matrix(rbinom(adm_indv*nloci, 1, ganc), nrow=adm_indv, ncol=nloci)
lanc_pop1 <- lanc_p + lanc_m

#admixed pop genotype 
geno_adm <- t(sapply(c(1:adm_indv), FUN=LAnc2Geno, lanc_p=lanc_p, lanc_m=lanc_m, frq_pop1=f1, frq_pop2=f2))

# allele freq
# fadm <- apply(geno_adm, MARGIN=2, FUN=af)

# prs of admixed pop
adm_prs <- apply(geno_adm, MARGIN=1, FUN=indiv_prs, beta=bg_sel)


################### output local ancestry of admixed pop as plink input ########################

print("writing to plink file")
# output the lanc of admixed
LAnc2Plink(lanc_pop1, nloci, adm_indv, filename, FID="ADM", adm_prs)
# output the lanc of parental metapop
LAnc2Plink(lanc_meta, nloci, (n1_indv + n2_indv), paste0("meta_", filename), FID="meta", meta_prs)


################### output exp and obs value in parental groups ########################

# write to file to record exp. and obs. value 
sink(file=paste0("sim_",filename, "_exp_obs.txt"), type = "output")
cat("vtotal.genic \t vtotal.obs \t vwithin.exp \t vwithin.obs \t vbetween.genic \t vbetween.obs \t fst \t p_sel \t prs.diff \t P \t t \n")
cat(vtotal.genic, "\t", vtotal.obs, "\t", vwithin.exp, "\t", vwithin.obs, "\t", vbetween.genic, "\t", vbetween.obs, "\t",  fst, "\t", p_sel, "\t", prs.diff, "\t", P, "\t", t, "\n")
sink()

