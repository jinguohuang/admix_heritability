# vgexplained validation
# plot Vg (Analytical) vs Variance of polygenic scores

vgsum=read.table("../data/vgprop_multi", header = F)
colnames(vgsum)=c("vg1p",	"vg2p",	"vg3p",	"vg4p",	"vg", "nloci", "trait")
vgsum$trait=c("SkinPigmentation","Height","BMI","BirthWeight","SBP","DBP","PP","HDL","LDL","TG","TC","RBC","HGB","HCT","MCH","MCV","MCHC","RDW","WBC","NEU","LYM","MON","BAS","EOS","PLT","MPV")
vgvalue=read.table("../data/vgvalueCI_multi", header = F)
colnames(vgvalue)=c("vgvalue", "CI95l", "CI95r", "trait")
vgvalue$trait=c("SkinPigmentation","Height","BMI","BirthWeight","SBP","DBP","PP","HDL","LDL","TG","TC","RBC","HGB","HCT","MCH","MCV","MCHC","RDW","WBC","NEU","LYM","MON","BAS","EOS","PLT","MPV")

library(plyr)
vgall=plyr::join(vgsum, vgvalue, by='trait') # keep order
vgall$group=c(rep("Anthropometric", 4), rep("Blood Pressure", 3), rep("Lipids", 4), rep("Blood Cell", 15))

library(ggplot2)
library(ggpubr)
# colorbline-friendly palette
cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
p=ggplot(vgall, aes(x=vg, y=vgvalue)) + 
  theme_classic() +
  xlab("Vg (Analytical)")+
  ylab("Variance of polygenic scores")+
  geom_errorbar(aes(ymin = CI95l, ymax = CI95r),
                colour="grey", width=.01) +
  geom_abline(slope=1, intercept=0, color="darkred", linewidth=1) + 
  geom_text(x=0.42, y=0.35, label="y=x", color="darkred", size=5) +
  geom_smooth(method='lm', formula= y~x, color="darkblue", linewidth=1, se = FALSE)+
  stat_regline_equation(label.x = 0.15, 
                        label.y = 0.50, 
                        aes(label = after_stat(eq.label))) +
  stat_regline_equation(label.x = 0.15, 
                        label.y = 0.43, 
                        aes(label = ..rr.label..)) +
  geom_point(size=2, alpha=0.5, aes(color=group)) +
  
  labs(color=NULL) + # remove legend title
  scale_colour_manual(values=cbp2)

ggsave("../figs/traits_vgpropVSvgvalueCI_multi.png", plot=p,
       width = 5, height = 5, dpi = 300, units = "in", device='png')
