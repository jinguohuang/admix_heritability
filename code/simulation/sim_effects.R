#!/usr/bin/env Rscript

args = commandArgs(TRUE)

F = rprojroot::is_git_root$make_fix_file()

print("reading in raw genotypes")
infile = args[1]
iter = args[2]
pop=args[3]
outfile=args[4]
library(data.table)
#infile = F(paste("data/",pop,"/1kg.",pop,".rmdup.thinned.", iter, ".raw", sep = ""))

dat= fread(infile, header = TRUE)
ninds = dim(dat)[1]
nvariants = 1000
mat = as.matrix(dat[,c(7:1006)])

freq = apply(mat, 2, mean)/2
pvar = apply(mat, 2, var)

print("assigning effect sizes")

b = matrix(sapply(pvar, function(x){
	rnorm(1, 0, sd = sqrt(1/(1000*x)))
	}), nvariants, 1)

# sum(b^2 * pvar)
# sum(b^2* 2*freq*(1-freq))

print("computing genetic variance components")

ldmat = cov(mat)

g = mat%*%b
vg = var(g)
vg.mat = b%*%t(b) * ldmat
vgenic = sum(diag(vg.mat))

vld = 2*sum(vg.mat[lower.tri(vg.mat, diag = FALSE)])

print("writing output")

# output = data.table(vg = vg, vgenic = vgenic, vld = vld)
# fwrite(output, "~/admixture/asw/1kg.asw.rmdup.thinned.1.vg", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

#outfile = F(paste("data/",pop,"/1kg.",pop,".rmdup.fst0.05.thinned.",iter,".vg", sep = ""))
output=formatC(c(vg, vgenic, vld), digits = 3)
#output = c(output, "\n")
zz=file(outfile, "wb")
writeBin(paste(output, collapse = "\t"), zz)
close(zz)


