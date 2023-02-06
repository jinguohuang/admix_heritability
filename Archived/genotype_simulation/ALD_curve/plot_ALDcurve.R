#!/usr/bin/env Rscript

ald_plot<-function(mapfile, lancfile){
	library(data.table)
	map<-fread(mapfile)
	colnames(map)<-c("pos","g_dist")
	lanc<-fread(lancfile)
	a<-cbind(map, lanc)
	# create dataframe to put data
	b<-data.frame()
	k=1
	for (i in 1:10){
	  for (j in 1:500){
	    b[k,1]=abs(a[i,2]-a[j,2])
	    b[k,2]=cor(unlist(a[i,3:402]),unlist(a[j,3:402]))
	    k=k+1
	  }
	} 
	colnames(b)<-c("distance","cor")
	png(file=paste0(lancfile,'_ALDcurve.png'), width=4, height=4, units="in", res=300)
	plot(b$distance, b$cor, main="Genetic distance vs SNP correlation",
	   xlab="Genetic distance (cM) ", ylab="SNP correlation ", pch=1)
	dev.off()
}

args = commandArgs(trailingOnly=TRUE)
mapfile=args[1]
lancfile=args[2]
ald_plot(mapfile, lancfile)
