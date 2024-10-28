#!/usr/bin/env Rscript

suppressWarnings(suppressMessages({
  library(data.table)
  library(rprojroot)
}))

F = is_git_root$make_fix_file()

args=commandArgs(TRUE)

p=args[1]
rep=args[2]

set.seed(rep)

filename = paste("data/ldsc/plink/admix_chrall.cm_thin", p , "_" , rep , ".raw", sep = "")
dat = fread(F(filename))
mat = as.matrix(dat[,-c(1:6)])
m = ncol(mat)
n = nrow(mat)

snps = data.table(snps = colnames(mat))
snps = snps[, tstrsplit(snps, split = "_",names = c("chr","pos","a1","a2","a3"))]
snps[, id := paste(chr, "_",pos,"_",a1,"_",a2, sep = "")]
snps[, a3freq := apply(mat, 2, mean)/2]

h_l = 0.8/m

#sample maf-dependent effects using the model above
snps$beta = sapply(snps$a3freq, function(x){rnorm( 1 , mean = 0, sd = sqrt(h_l * (2*x*(1-x))^-1 ))})

# let's calculate sigma2_g to confirm that the total genetic variance is indeed 0.8
# sigma2_g = sum( mapply(function(b,p){ b^2* 2*p*(1-p) }, snps2$beta, snps2$ALT_FREQS))

effects_file = paste("data/ldsc/effects/admix_chrall.cm_thin", p , "_" , rep , ".effects", sep = "")
#save the effect sizes to file and use plink2 to generate PRS
fwrite(snps[, .(id,a3,beta)], F(effects_file), row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")


#simulate phenotype
g = mat%*%snps$beta
y = g + rnorm(n, 0, sqrt(0.2))
pgs = dat[, .(FID, IID)]
pgs[, g := g ]
pgs[, y := y]
pgs_file = paste("data/ldsc/phenotypes/admix_chrall_p",p,"_",rep,".pgs", sep = "")
fwrite(pgs, F(pgs_file), row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")




