#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# reformat the result to compare expected vs estimated vgb
# add geno gcta results in to the table
library(data.table)
theta=args[1]
vgb=args[2]

simfile<-paste0("sim_prop",theta, "_vgb",vgb,"_expected.txt")
estfile<-paste0("admix_prop",theta, "_vgb",vgb, ".hsq")
genofile<-paste0("admix_geno_prop",theta, "_vgb",vgb, ".hsq")
  
sim<-fread(simfile)
est<-fread(estfile, fill=TRUE)
geno<-fread(genofile, fill=TRUE)
# function to reformat hsq file
reformat <- function(x){
  x<-x[complete.cases(x), ] # delete NA
  x<-as.data.frame(melt(x, id.vars = "Source"))
  x$Source <- paste(x$Source, x$variable, sep=".")
  x<-x[c("Source", "value")]
}
# reformat
geno<-reformat(geno)
est<-reformat(est)
#add extra for geno
geno$Source=paste0(geno$Source, ".geno")
# bind all together
colnames(sim)<-c("Source", "value")
est_sim<-rbind(sim, est, geno)
est_sim_t<-t(est_sim)
# output the reformat info
write.table(est_sim_t, file = paste0("simVSest_prop",theta, "_vgb",vgb,".txt"), 
            quote = F, col.names = F, row.names = F)