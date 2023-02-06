library("CMplot")
#load the data
filename="adm_400hap_10deme_local_anc.AFR.annotated"
anc<-read.table(filename, header = F)
colnames(anc)<-c("SNP","CHR","POS","ANC")
avganc<-mean(anc$ANC)  
varanc<-var(anc$ANC) 
CMplot(anc, plot.type="m", band=0.2, LOG10=FALSE, ylab="AFR anc per locus",threshold=avganc,
        threshold.lty=2, threshold.lwd=3, threshold.col="red", amplify=FALSE, width=14,height=6,
        signal.col=NULL, chr.den.col=NULL, file="jpg",memo="",dpi=300,file.output=TRUE,
        verbose=TRUE,cex=0.8)

#output the mean and variance result to a file
sink(file=paste0(filename,'mean_var.txt'), type = "output")
cat("/* File created on", date(), "*/\n")
cat("/* Mean of AFR ancestry:", avganc, "*/\n")
cat("/* Variance of AFR ancestry:", varanc, "*/\n")
sink()
