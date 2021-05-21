#!/usr/bin/env Rscript

# skt test and plot for distribution
# adapted from Arslan's code

# function to run SKT test without replacement
# filename of the ancestry by chromosome file
# n is the number of chromosome 
# iteration is the number of iteration of SKT test
SKT_noreplace<-function(filename,n,iteration,cormethod){
library(data.table)
dat = fread(filename)
dat = as.matrix(dat)
#SKT algorithm (naive):
#1. Sample 11 chromosomes without replacement to put in one group. put the rest in second group
#2. Calculate ancestry for each individual for each group
#3. Calculate the correlation in ancestry between the two groups
skt1 = c()
for(i in 1:iteration){
  #n=22
  ix.evens = sample(1:n, n/2, replace = FALSE)
  ix.odds = setdiff(1:n, ix.evens)
  anc.evens = apply(dat[ix.evens, ], 2, mean)
  anc.odds = apply(dat[ix.odds, ], 2, mean)
  skt1[i] = cor(anc.evens, anc.odds, method=cormethod)
}
return(skt1)
}


# Technically, we need to make sure the ancestry we are calculating is across equally sized chromosomes
# The algorithm above doesn't do that. By chance, group 1 might have longer chromosomes than group 2.
# To avoid this, 1. create a set of odd and even chromosomes
# 2. resample with replacement within each set
# 3. calculate ancestry for each individual within each set
# 4. calculate the correlation between the two sets.

SKT_replace<-function(filename,n,iteration,cormethod){
  library(data.table)
  dat = fread(filename)
  dat = as.matrix(dat)
  #define indices for even and odd chromosomes
  evens = seq(2,n,2)
  odds = seq(1,n,2)
  #separate dataset for even and odd chromosomes
  dat.evens = dat[evens,]
  dat.odds = dat[odds,]
  skt2 = c()
  for(i in 1:iteration){
    ix.evens = sample(1:(n/2), (n/2), replace = TRUE)
    ix.odds = sample(1:(n/2), (n/2), replace = TRUE)
    anc.evens = apply(dat.evens[ix.evens, ], 2, mean)
    anc.odds = apply(dat.odds[ix.odds, ], 2, mean)
    skt2[i] = cor(anc.evens, anc.odds, method=cormethod)
  }
return(skt2)
}


# plot for ASW
filename="ASW_ancBYppl.txt"
#cormethod="pearson"
cormethod="spearman"


# ASW WITHOUT REPLACEMENT
ASW_no=SKT_noreplace(filename,22,1000,cormethod)
#plot
print("Plotting for ASW WITHOUT REPLACEMENT")
png(file="ASW_SKT_noreplace_s.png", width=4, height=4, units="in", res=300)
hist(ASW_no, col='blue', xlim=c(0, 1), main = "Distribution of SKT test results of ASW \n (sample without replacement)", xlab=paste0(cormethod, " correlation"))
dev.off()


# ASW WITH REPLACEMENT
ASW=SKT_replace(filename,22,1000,cormethod)
#plot
print("Plotting for ASW WITH REPLACEMENT")
png(file="ASW_SKT_replace_s.png", width=4, height=4, units="in", res=300)
hist(ASW, col='blue', xlim=c(0, 1), main = "Distribution of SKT test results of ASW \n (sample with replacement)", xlab=paste0(cormethod, " correlation"))
dev.off()


# simulated WITHOUT REPLACEMENT
file_10deme="simAA_1Mb_1Kppl_10deme_ancBYppl.txt"
file_random="simAA_1Mb_1Kppl_random_ancBYppl.txt"

no_10deme=SKT_noreplace(file_10deme,100,1000,cormethod)
no_random=SKT_noreplace(file_random,100,1000,cormethod)
#plot
print("Plotting for simulated WITHOUT REPLACEMENT")
png(file="sim_SKT_noreplace_s.png", width=4, height=4, units="in", res=300)
hist(no_10deme, col='blue', xlim=c(0, 1), main = "Distribution of SKT test results of simulated data \n (sample without replacement)", xlab=paste0(cormethod, " correlation"))
hist(no_random, col='red', add=T)
dev.off()


# simulated WITH REPLACEMENT
with_10deme=SKT_replace(file_10deme,100,1000,cormethod)
with_random=SKT_replace(file_random,100,1000,cormethod)
#plot
print("Plotting for simulated WITH REPLACEMENT")
png(file="sim_SKT_replace_s.png", width=4, height=4, units="in", res=300)
hist(with_10deme, col='blue', xlim=c(0, 1), main = "Distribution of SKT test results of simulated data \n (sample with replacement)", xlab=paste0(cormethod, " correlation"))
hist(with_random, col='red', add=T)
dev.off()



