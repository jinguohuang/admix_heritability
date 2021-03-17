#!/usr/bin/env Rscript

# for plotting structure of admixture results for each K

#Usage:
#-Place the following files in the same folder as this script
# 	1) your plink binary files (.fam), 
#	2) ADMIXTURE results file (.Q), 
# 3) ASW 66 diploid AFR anc file 


# load file
args = commandArgs(trailingOnly=TRUE)
Qfile<-args[1]  # load .Q file from ADMIXUTRE result
famfile<-args[2] # load fam file from dataset

# get filename and K 
#Name <- sub("\\_.*", '', Qfile)
Name <- sub("\\_AT.*", '', Qfile) #Name before _AT, to include mig_rate
K <- as.numeric(sub('.*?\\.(.*?)\\.Q', '\\1', Qfile))

require(reshape2)
require(ggplot2)

#load the .Q file
qfile=read.table(Qfile,header=F,stringsAsFactors = F)
colnames(qfile)=c(paste0("Anc", 1:K))
#add individual ID labels from .fam file
fam<-read.table(famfile,header=F,stringsAsFactors = F)
colnames(fam)=c("FID","IID","Mot_ID","Fat_ID","Sex","Pheno")
qfile$IID<-fam$IID
#melt the data.frame - i.e. convert to long format
mqfile<-melt(qfile,id.vars=c("IID"))
colnames(mqfile)<-c("IID","Ancestry_component","Ancestry")
#change the order to IID
mqfile$IID<-factor(mqfile$IID, levels = qfile$IID)
#mqfile$Ancestry_component <- factor(mqfile$Ancestry_component, levels =c("Anc1","Anc2"))

#plt<-ggplot(mqfile[order(mqfile$IID),],aes(x=factor(IID), y=Ancestry, fill=Ancestry_component))+
plt<-ggplot(mqfile,aes(x=IID, y=Ancestry, fill=Ancestry_component))+
  geom_bar(stat="identity",width=1,position="fill")
#remove individual IDs from x axis as it can get crowded if you have a lot of samples
plt<-plt+
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line=element_blank(),
        panel.background=element_blank())+
  labs(y="Ancestry",fill="Ancestral group") #label axes and legends

#plt
#one final tweak - change the colors as they are horrendous
#create new pallete to choose from
plt<-plt+
  scale_fill_manual(values=c(
   "skyblue2","#1b9e77","#7570b3","#d95f02", "khaki2","#FB9A99", # lt pink
  "gray70", "yellow3",
  "maroon", "steelblue4",
   "yellow4", "darkturquoise",
  "darkorange4",
  "palegreen2","#6A3D9A", # purple
  "#FF7F00", # orange
  "gold1", "deeppink1", "blue1",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
"dodgerblue2", "#E31A1C", # red
  "green4"))

# save plot as pdf
ggsave(plt,file=paste0(Name,"_K",K,"_admixture.pdf"))


# plot for histogram of 10 demes global anc 
deme<-read.table(Qfile,header=F,stringsAsFactors = F)
# check the anc order in ADMIXTURE result before assignment
colnames(deme)<-c("EUR","AFR")
# check the number of individuals in ref and adm samples first
deme_adm<-deme[101:10100,2]
pdf(file = paste0(Name,"_AFR_dist.pdf"), width = 4,height = 4)
hist(deme_adm, main = "10 deme AFR anc distribution", xlab = "AFR")
dev.off()

#plot the density plot to compare ASW and AA
# density plot
deme_adm<-deme[101:10100,][2]
deme_adm$pop <-"Simulated_AA"
ASW_file<-args[3]
ASW_hap <- read.table(ASW_file, header= F)
colnames(ASW_hap)<-c("AFR")
ASW_hap$pop <-"ASW"
ASW_AA<-rbind(ASW_hap,deme_adm)
#Calculate the mean of each group :
require(plyr)
mu <- ddply(ASW_AA, "pop", summarise, grp.mean=mean(AFR))
#plot
# Change density plot line colors by groups
# Add mean lines
p<-ggplot(ASW_AA, aes(x=AFR, fill=pop)) +
  geom_density(alpha=0.4)+
  ggtitle(paste0(Name,"_ASW_density_plot"))+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=pop),
             linetype="dashed") +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
ggsave(p,file=paste0(Name,"_ASW_density.pdf"))


#output the ks test result
# Start writing to an output file
sink(file=paste0(Name,'kstest.txt'))
# ks test 10 deme vs ASW
ASW<-as.numeric(unlist(ASW_hap[1]))
AA<-as.numeric(unlist(deme_adm[1]))
cat("=============================\n")
cat("ks-test between ASW and AA\n")
cat("=============================\n")
# Do x and y come from the same distribution?
ks.test(ASW, AA)
# Stop writing to the file
sink()






