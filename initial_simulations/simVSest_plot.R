#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(data.table)
#filename<-"simVSest_ALL.txt"
filename<-args[1]
simVSest<-fread(filename)
library(ggplot2)

# calculate the upper and lower for errorbar 95% confidence interval
simVSest$adjustVG = simVSest$`V(G).Variance` /( 4 * simVSest$admixture * (1 - simVSest$admixture))
#Adjust SE or not?
simVSest$VG.SE=1.96*simVSest$`V(G).SE`/( 4 * simVSest$admixture * (1 - simVSest$admixture))
#simVSest$VG.SE=1.96*simVSest$`V(G).SE`

simVSest$adjustVG.lower=simVSest$adjustVG-simVSest$VG.SE 
simVSest$adjustVG.upper=simVSest$adjustVG+simVSest$VG.SE

# plot to compare vgb
png(file=paste0("simVSest_ALL_vgb_addenv.png"), width=12, height=6, units="in", res=300)
ggplot(simVSest, aes(x=admixture, y=adjustVG)) +
#ggplot(simVSest, aes(x=admixture, y=`V(G).Variance`)) +
  geom_errorbar(width=.005,aes(ymin = adjustVG.lower, ymax = adjustVG.upper)) +
  geom_point(shape=21, size=3, fill="white") +
  #add facet for different admixture
  facet_grid(cols = vars(vgb)) +
  #facet_grid(vars(vgb), vars(admixture))+
  ggtitle("Expected VS Estimated Between Population Genetic Variance Under different Admixture Proportion (CI=95%) Ve=0.2" ) + 
  xlab("Admixture Proportion") + ylab("Estimated Genetic Variance (Adjusted)") +
  #geom_abline(slope=0, intercept = vbetween, size = .5, color = "red") + 
  geom_point(aes(x=admixture, y=vbetween), colour="red") +
  #ylim(0,1.2) + 
  theme_bw()
dev.off()