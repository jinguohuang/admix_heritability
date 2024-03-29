#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# reformat the result to compare expected vs estimated vgb
filename = args[1]

library(data.table)
estfile<-paste0(filename, ".reml.hsq")
metafile<-paste0("meta_", filename, ".reml.hsq")
simfile<-paste0("sim_",filename, "_exp_obs.txt")

sim<-fread(simfile, header = FALSE)
est<-fread(estfile, fill=TRUE)
meta<-fread(metafile, fill=TRUE)
#rownames(est)
# function to reformat hsq file
reformat <- function(x){
  x<-x[complete.cases(x), ] # delete NA row
  # change data type to numeric
  x$Variance<-as.numeric(x$Variance)
  x<-as.data.frame(melt(x, id.vars = "Source"))
  x$Source <- paste(x$Source, x$variable, sep=".")
  x<-t(x[c("Source", "value")]) # transpose
}
# reformat
est<-reformat(est)
rownames(est)<-NULL #delete rownames

meta<-reformat(meta)
rownames(meta)<-NULL 

# change rownames of meta
meta[1,] = paste0(meta[1,], "_meta")

# bind all together
est_sim<-cbind(est, sim, meta)

# output
write.table(est_sim, paste0("summary_",filename, ".txt"), quote=FALSE, sep = '\t' ,
            row.names = FALSE, col.names = FALSE)

