
#Part 1: Verify h2
#Step 1 - calculate var(prs)
#Step 2 - calculate var(pheno)
#Step 3 - Verify that var(prs)\var(pheno) is about 0.8 = h2
#Step 4 - g* = g - B0, var(g*) = expectation
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

# model = "HI"
# P=0
# cov = "pos"
# seed=9
# t=12

summarydir <- paste0("/home/aazaidi/klema030/AdjustedHE/output/", model, "/gcta/HE_", grm, "/") # directory, can be modified 

pheno <- fread(paste0("/home/aazaidi/huan2788/admix_heritability/data/theta0.5_gen20/pheno/", model, "/admix_", model, "_theta0.5_gen20_P", P , "_", cov, "_seed", seed, "_t", t, ".pheno"), header=F)$V3

ganc <- fread(paste0("/home/aazaidi/huan2788/admix_heritability/data/theta0.5_gen20/ganc/", model, "/admix_", model, "_theta0.5_gen20_P", P , "_", cov, "_seed", seed, "_t", t, ".ganc"), header=F)$V3 
ganc.ids <- fread(paste0("/home/aazaidi/huan2788/admix_heritability/data/theta0.5_gen20/ganc/", model, "/admix_", model, "_theta0.5_gen20_P", P ,"_", cov, "_seed", seed, "_t", t, ".ganc"), header=F)$V2  

prs <- fread(paste0("/home/aazaidi/huan2788/admix_heritability/data/theta0.5_gen20/PRS/", model, "/admix_", model, "_theta0.5_gen20_P", P ,"_", cov, "_seed", seed, "_t", t, ".prs"), header=F)$V3 

filename <- paste0("admix_", model,"_theta0.5_gen20_P", P ,"_", cov, "_seed", seed, "_t", t)

#head(ganc.ids)

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
print(paste("Vg=",var.gstar))

#Part 2
print("Calculating y*y*' and 00'")
#Step 1 - Calculate g* (residualized g, regress g~0 and take residuals) 
s_pheno<-scale(pheno, scale=FALSE)
ly<- lm(s_pheno~ganc)
ystar<- ly$residuals
#print(g_resid)

#Step 2 - calc g*g*' and take lower half of matrix
yy_star <- ystar%*%t(ystar)
yy_star <- yy_star[upper.tri(yy_star)]
#print(gg_star)

#Step 3 - Calculate 00' and take lower half of matrix
tt <- ganc %*% t(ganc)
#print(theta)
tt<- tt[upper.tri(tt)]
#print(theta)

#head(ganc)

print("reading GRM")
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
G<-ReadGRMBin(paste0(grmdir,"/admix_", model, "_theta0.5_gen20_P", P ,"_", cov,"_seed", seed, "_t", t))

G.ids = G$id$V2
G.diag = G$diag
G.off = G$off

length(G.off)

#Step 5: regress g*g*' ~ A + 00'

y<- lm(yy_star ~ (G.off + tt))

summary(y)

#t, var.gstar, coefficient from G.off
G.off.coeff <- summary(y)$coefficients[2, 1]
print(paste("Adj Vg= ",G.off.coeff))

print("Writing output")
#Output
output = data.table(t=t, P=P, seed = seed, cov=cov, var.gstar = var.gstar, G.off.coeff = G.off.coeff)
fwrite(output,
       paste0(summarydir, filename,".txt",sep=""),
       col.names = TRUE,
       row.names=FALSE,
       quote=FALSE,
       sep="\t")
print("Done")


