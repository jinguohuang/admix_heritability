#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# reformat the result to compare expected vs estimated vgb
# add geno gcta results in to the table
library(data.table)
theta=args[1]
fstc=args[2]
fstn=args[3]
pcausal=args[4]

filename <- paste0("_prop", theta, "_fstc", fstc, "_fstn", fstn, "_pcausal", pcausal) #for filename

simfile<-paste0("sim", filename ,"_expected.txt")
estfile<-paste0("admix",filename, ".reml.hsq")


sim<-fread(simfile)
est<-fread(estfile, fill=TRUE)
# function to reformat hsq file
reformat <- function(x){
  x<-x[complete.cases(x), ] # delete NA
  x<-as.data.frame(melt(x, id.vars = "Source"))
  x$Source <- paste(x$Source, x$variable, sep=".")
  x<-x[c("Source", "value")]
}
# reformat
est<-reformat(est)
# bind all together
colnames(sim)<-c("Source", "value")
est_sim<-rbind(est, sim)
est_sim_t<-t(est_sim)
# output the reformat info
write.table(est_sim_t, file = paste0("simVSest", filename, ".txt"), 
            quote = F, col.names = F, row.names = F)