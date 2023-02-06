# estimate vgamma with greml using qgg package
args = commandArgs(trailingOnly=TRUE)
library(qgg)

seed=args[1]
t=args[2]
P=args[3]
filename=paste0("admix_seed",seed, "_P", P, "_psel0_m10_t", t)

lanc=read.table(paste0(filename, ".dosage"), header=F)
nloci=dim(lanc)[1]
#reshape to get lanc
lanc=t(lanc[4:dim(lanc)[2]])
pi=colMeans(lanc)/2 
# lanc construct GRM
W <- matrix(, nrow = nrow(lanc), ncol = ncol(lanc))
for(i in 1:nrow(W)){ 
  for (j in 1:ncol(W)){
    W[i,j]=(lanc[i,j]-2*pi[j])/sqrt(2*pi[j]*(1-pi[j]))
  }
}
# Calculate the GRM
A = W%*%t(W) / nloci

# PRS as pheno
pheno=read.table(paste0(filename, ".pheno"), header=F)
# mixed model estimate Vg
l1 = greml(y = pheno$V3, X = NULL, GRM = list(A))

#the variance due to random effects
vgamma=l1$theta[1] 
sink(file=paste0(filename, "_vgammahat.txt"), type = "output")
cat(vgamma)
sink() 
