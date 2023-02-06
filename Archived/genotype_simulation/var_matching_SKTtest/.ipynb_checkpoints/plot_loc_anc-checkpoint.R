#!/usr/bin/env Rscript

# for plotting ancestry at each locus and calculate mean and variance across chromosomes

#Usage:
#-Place the following files in the same folder as this script
# 	1) calculated local ancestry at each locus file


# load file
args = commandArgs(trailingOnly=TRUE)
ancfile<-args[1]  # load ancestry file

library("CMplot")
#load the data
anc<-read.table(ancfile, header = F)
colnames(anc)<-c("SNP","CHR","POS","ANC")
avganc<-mean(anc$ANC)  
varanc<-var(anc$ANC) 
CMplot(anc, plot.type="m", band=0.2, LOG10=FALSE, ylab="AFR anc per locus",threshold=avganc,
        threshold.lty=2, threshold.lwd=3, threshold.col="red", amplify=FALSE, width=14,height=6,
        signal.col=NULL, chr.den.col=NULL, file="jpg",memo="",dpi=300,file.output=TRUE,
        verbose=TRUE,cex=0.8)

#output the mean and variance result to a file
sink(file=paste0(ancfile,'mean_var.txt'), type = "output")
cat("/* File created on", date(), "*/\n")
cat("/* Mean of AFR ancestry:", avganc, "*/\n")
cat("/* Variance of AFR ancestry:", varanc, "*/\n")
sink()
