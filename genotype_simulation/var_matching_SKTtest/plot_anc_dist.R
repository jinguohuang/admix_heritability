#!/usr/bin/env Rscript

# for plotting ancestry at each locus and calculate mean and variance across chromosomes

#Usage:
#-Place the following files in the same folder as this script
# 	1) calculated global ancestry for each person


# load file
args = commandArgs(trailingOnly=TRUE)
ancfile<-args[1]  # load ancestry file

#load the data
anc<-read.table(ancfile, header = F)
anc_mean<-colMeans(anc)
png(file=paste0(ancfile,'_anc_dist.png'), width=4, height=4, units="in", res=300)
hist(anc_mean,main = paste("Distribution of global ancestry"), xlab = "AFR ancestry")
dev.off()
avganc<-mean(anc_mean)

#output the mean result to a file
sink(file=paste0(ancfile,'_mean'), type = "output")
cat("/* File created on", date(), "*/\n")
cat("/* Mean of AFR ancestry:", avganc, "*/\n")
sink()
