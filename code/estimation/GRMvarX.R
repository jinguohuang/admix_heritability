#!/usr/bin/env Rscript
# calculate GRM scaled by Var(X) instead of 2P(1-P) as GCTA did
# output GRM to feed GCTA to see if it give us 2.1+2.2+2.3


library(data.table)
args = commandArgs(trailingOnly=TRUE)

#filename="admix_HI_theta0.5_gen20_P0.9_pos_seed1_t5"
filename=args[1]
plinkdir=args[2] 
grmdir=args[3] #GRM output directory

# load the genotype data
print(paste0("Reading file: ", filename))
geno=fread(paste0(plinkdir, "/", filename, ".dosage"))
# deal with the 3 col map, now we have 1000 row as SNP, 10000 col as ind
X=geno[, 4:10003]

# given the geno matrix, try to calculate the pi for each loci.
# transpose geno 
M=t(X)
#dim(M) #10000ppl(row)   1000loci(col)
pi=colMeans(M)/2  

# variance of genotype of each SNP
VarMj=apply(M, 2, var)

print("Standardizing ...")
W <- matrix(, nrow = nrow(M), ncol = ncol(M))
for(i in 1:nrow(W)){ 
  for (j in 1:ncol(W)){
    W[i,j]=(M[i,j]-2*pi[j])/sqrt(VarMj[j])
  }
}

#dim(W) #10000  1000
#sum(is.na(W)) #0

#We can calculate the GRM with A=WW'/N.
print("Constructing GRM ...")
A=W %*% t(W)/1000 # Does 1/N matter? Yes
#dim(A) #10000 10000
#sum(is.na(A))

# save the lower triangle elements
A[upper.tri(A)]<-NA

#Save the lower triangle elements to a compressed text file (e.g. test.grm.gz) and save the IDs in a plain text file (e.g. test.grm.id).
#Reformat the GRM to GCTA format: test.grm.gz (no header line; columns are indices of pairs of individuals (row numbers of the test.grm.id),
#number of non-missing SNPs and the estimate of enetic relatedness)

# convert triangle matrix to pairwise list
rownames(A)<-1:10000
colnames(A)<-1:10000
A_grm<-data.frame(ID1=rownames(A)[row(A)],
                  ID2=colnames(A)[col(A)], 
                  nomiss=1000, 
                  relatedness=c(A))
# DELETE na
print("Getting lower triangle of GRM ...")
A_grm<-A_grm[complete.cases(A_grm), ]

print("Saving GRM ...")
# output this as GCTA GRM format
fwrite(A_grm, 
       file = paste0(grmdir, "/", filename, ".grm.gz"),
       #file = "admix_geno.grm.gz", 
       sep = '\t',
       row.names = FALSE, 
       col.names = FALSE, 
       compress="gzip")

# write the id file for gcta test.grm.id
id_grm<-data.frame(FID="ADM",IID=1:10000)
fwrite(id_grm, 
       file = paste0(grmdir, "/", filename, ".grm.id"),
       #file = "admix_geno.grm.id", 
       sep = '\t',
       row.names = FALSE, 
       col.names = FALSE)

# output this to gcta
