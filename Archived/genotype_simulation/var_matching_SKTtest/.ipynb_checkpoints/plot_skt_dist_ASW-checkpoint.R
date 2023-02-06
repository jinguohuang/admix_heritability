#!/usr/bin/env Rscript

SKTtest<-function(ancfile){
anc<-read.table(ancfile, header = F)
#calculate average anc of every other 2haplotype
# transpose first
df<-t(anc)
# calculate average every other 2 rows and transpose back 
n <- 2
anc_ind<-as.data.frame(t(aggregate(df, list(rep(1:(nrow(df) %/% n + 1), each = n, len = nrow(df))), mean)[-1]))

# sample without replacement
bootfun <- function(M, reps) {
# a single bootstrap
  boot <- function(M){
    # M is the individual ancestry file
    n <- length(M)
    idx <- sample(rep(1:2, each = ceiling(n /2))[1:n], replace = FALSE)
    M1 <- M[idx == 1]
    M2 <- M[idx == 2]
    # calculate average anc across each set
    M1_anc<-colMeans(M1)
    M2_anc<-colMeans(M2)
    # return spearman correlation of 2 sets
    return(cor(M1_anc, M2_anc, method = "spearman"))
  }
  bootrep <- replicate(n=reps, boot(M))
  return(bootrep)
}

# run bootstrap
a<-bootfun(anc_ind, 10000)
return(a)
}
#load file
#deme10="simAA_1Mb_1Kppl_10deme_ancBYppl.txt"
#random="simAA_1Mb_1Kppl_random_ancBYppl.txt"
args = commandArgs(trailingOnly=TRUE)
ancfile<-args[1]  # load ancestry file
#deme10_test<-SKTtest(deme10)
#random_test<-SKTtest(random)
anc_test<-SKTtest(ancfile)
# plot distribution
#png(file=paste0(ancfile,'_anc_dist.png'), width=4, height=4, units="in", res=300)
#hist(anc_mean,main = paste("Distribution of global ancestry"), xlab = "AFR ancestry")
png(file=paste0(ancfile,'_skt_test_dist.png'), width=4, height=4, units="in", res=300)
hist(anc_test, col='blue', xlim=c(0, 1), main = "Distribution of SKT test results", xlab="Spearman correlation")
#hist(random_test, col='red', add=T)
dev.off()
