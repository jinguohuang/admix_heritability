#HE regression without ancestry correction
library(data.table)

# load data
args <- commandArgs(trailingOnly = TRUE)
model <- args[1]
P <- args[2]
cov <- args[3]
seed <- args[4]
t <- args[5]
grm <- args[6]
grmdir <- args[7]

# model = "CGF"
# P=0.3
# cov = "pos"
# seed=1
# t=20

summarydir <- paste0("/home/aazaidi/klema030/AdjustedHE/output/", model, "/gcta/HE_", grm, "/") # directory, can be modified 

pheno <- fread(paste0("/home/aazaidi/huan2788/admix_heritability/data/theta0.5_gen20/pheno/", model, "/admix_", model, "_theta0.5_gen20_P", P , "_", cov, "_seed", seed, "_t", t, ".pheno"), header=F)$V3

ganc <- fread(paste0("/home/aazaidi/huan2788/admix_heritability/data/theta0.5_gen20/ganc/", model, "/admix_", model, "_theta0.5_gen20_P", P , "_", cov, "_seed", seed, "_t", t, ".ganc"), header=F)$V3 
ganc.ids <- fread(paste0("/home/aazaidi/huan2788/admix_heritability/data/theta0.5_gen20/ganc/", model, "/admix_", model, "_theta0.5_gen20_P", P ,"_", cov, "_seed", seed, "_t", t, ".ganc"), header=F)$V2  

prs <- fread(paste0("/home/aazaidi/huan2788/admix_heritability/data/theta0.5_gen20/PRS/", model, "/admix_", model, "_theta0.5_gen20_P", P ,"_", cov, "_seed", seed, "_t", t, ".prs"), header=F)$V3 

filename <- paste0("admix_", model,"_theta0.5_gen20_P", P ,"_", cov, "_seed", seed, "_t", t)


#Calculate variance of PRS and Pheno files
var.prs<- var(prs)
print(paste("var(prs)=",var.prs))
var.pheno<- var(pheno)
print(paste("var(pheno)=" ,var.pheno))

#h2 should be about 0.8
verify_h2<- var.prs/var.pheno
print(paste("h2 =",verify_h2))

B<- lm(prs ~ ganc)
g_star<-B$residuals
var.gstar<-var(g_star)
print(paste("var.gstar=",var.gstar))

#Part 2
ystar <- pheno

#calc y*y*' and take upper half of matrix
yy_star <- ystar %*% t(ystar)
yy_star <- yy_star[upper.tri(yy_star)]

# R script to read the GRM binary file, From GCTA
ReadGRMBin=function(prefix, AllN=F, size=4){
  sum_i=function(i){
    return(sum(1:i))
  }
  BinFileName=paste(prefix,".grm.bin",sep="")
  NFileName=paste(prefix,".grm.N.bin",sep="")
  IDFileName=paste(prefix,".grm.id",sep="")
  id = read.table(IDFileName)
  n=dim(id)[1]
  BinFile=file(BinFileName, "rb");
  grm=readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
  NFile=file(NFileName, "rb");
  if(AllN==T){
    N=readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
  }
  else N=readBin(NFile, n=1, what=numeric(0), size=size)
  i=sapply(1:n, sum_i)
  return(list(diag=grm[i], off=grm[-i], id=id, N=N))
}
# G<-ReadGRMBin(paste0("/home/aazaidi/huan2788/admix_heritability/data/theta0.5_gen20/gcta/grm/admix_", model, "_theta0.5_gen20_P", P ,"_zero_seed", seed, "_t", t))

G<-ReadGRMBin(paste0(grmdir,"/admix_", model, "_theta0.5_gen20_P", P ,"_", cov,"_seed", seed, "_t", t))

G.ids = G$id$V2
G.diag = G$diag
G.off = G$off

length(G.off)

#regress y*y*' ~ A 

y<- lm(yy_star ~ G.off)

summary(y)

#t, var.gstar, coefficient from G.off
G.off.coeff <- summary(y)$coefficients[2, 1]
print(paste("Vg= ",G.off.coeff))

#Output
output = data.table(t=t, P=P, seed = seed, cov=cov, var.gstar = var.gstar, G.off.coeff = G.off.coeff)
fwrite(output,
       paste0(summarydir, filename,"noa.txt",sep=""),
       col.names = TRUE,
       row.names=FALSE,
       quote=FALSE,
       sep="\t")





