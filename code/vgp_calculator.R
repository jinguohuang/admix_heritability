# calculate vg proportion for each trait
# load trait file
args <- commandArgs(trailingOnly = TRUE)
trait <- args[1]
# trait="Height"
betafile=paste0(trait, "_allele_beta.txt")
f1file=paste0("CEU_", trait, ".frq")
f2file=paste0("YRI_", trait, ".frq")
betai=read.table(betafile, header = T)
#colClasses=c("character","character","numeric"))
f1=read.table(f1file, header = T,comment.char = "" )
f2=read.table(f2file, header = T,comment.char = "" )
# merge 2 freq files
f1f2=merge(f1,f2, by=c("CHR", "SNP", "A1", "A2"), suffixes = c("_EUR", "_AFR"))
# merge freq and beta
betafreq=merge(f1f2, betai, by=c("SNP", "A1"))

# output
write.table(betafreq, file=paste0(trait, "_betafreq.txt"), sep = '\t',
            quote = F, col.names = T, row.names = F)

# Calculate 4 terms of genetic variance Vg
# term 1 of Vg
vg_term1 <- function(exp.ganc, f1, f2, beta) return (sum(beta^2*(2*exp.ganc*f1*(1-f1)+2*(1-exp.ganc)*f2*(1-f2))))

# term 2 of Vg
vg_term2 <- function(exp.ganc, f1, f2, beta) return (sum(beta^2*(2*exp.ganc*(1-exp.ganc)*((f1-f2)^2))))

# term 3 of Vg
vg_term3 <- function(var.ganc, f1, f2, beta) return (sum(beta^2*(2*var.ganc*((f1-f2)^2))))

nloci=nrow(betafreq)
# term 4 of Vg
vg_term4 <- function(var.ganc, f1, f2, beta){
  nloci=nloci
  term4=matrix(, nrow = nloci, ncol = nloci)
  for (i in 1:nloci){
    for (j in 1:nloci){
      term4[i,j]=beta[i]*beta[j]*(f1[i]-f2[i])*(f1[j]-f2[j])
    }
  }
  # make diagnal 0
  diag(term4)=0
  # sum it up
  sum.term4 = sum(term4)*4*var.ganc
  return (sum.term4)} 


# calculate each term and total
vartheta=0.018
exptheta=0.767
vg1=vg_term1(exp.ganc=exptheta, 
             f1=betafreq$MAF_EUR, 
             f2=betafreq$MAF_AFR, 
             beta=betafreq$BETA)
vg2=vg_term2(exp.ganc=exptheta, 
             f1=betafreq$MAF_EUR, 
             f2=betafreq$MAF_AFR, 
             beta=betafreq$BETA)
vg3=vg_term3(var.ganc=vartheta, 
             f1=betafreq$MAF_EUR, 
             f2=betafreq$MAF_AFR, 
             beta=betafreq$BETA)
vg4=vg_term4(var.ganc=vartheta, 
             f1=betafreq$MAF_EUR, 
             f2=betafreq$MAF_AFR, 
             beta=betafreq$BETA)
vg=vg1+vg2+vg3+vg4

# add nloci and vg
vgp=data.frame(vg1p=vg1/vg,
               vg2p=vg2/vg,
               vg3p=vg3/vg,
               vg4p=vg4/vg, 
               vg=vg,
               nloci=nloci,
               Traits=trait)

# output the vgp for each trait
write.table(vgp, file=paste0("../",trait, "_vgp.txt"), sep = '\t',
            quote = F, col.names = T, row.names = F)

