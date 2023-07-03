# plot decomposing variance explained in 26 traits

alltraits=read.table("../data/vgprop_multi", header = F)
colnames(alltraits)=c("vg1p",	"vg2p",	"vg3p",	"vg4p",	"vg", "nloci", "trait")

# update trait name
alltraits$trait=c("SkinPigmentation","Height","BMI","BirthWeight","SBP","DBP","PP","HDL","LDL","TG","TC","RBC","HGB","HCT","MCH","MCV","MCHC","RDW","WBC","NEU","LYM","MON","BAS","EOS","PLT","MPV")

# reformat data for plot
# Load the package
library(tidyr)
library(dplyr)

# transform the format, get the order straight
data_long <- gather(alltraits, Term, Percentage, vg1p:vg4p) %>%
  arrange(factor(trait, levels = c("SkinPigmentation","Height","BMI","BirthWeight","SBP","DBP","PP","HDL","LDL","TG","TC","RBC","HGB","HCT","MCH","MCV","MCHC","RDW","WBC","NEU","LYM","MON","BAS","EOS","PLT","MPV"))) %>% 
  mutate(trait=factor(trait, levels=unique(trait)))

# add group
data_long$group=c(rep("Anthropometric", 4*4), rep("Blood Pressure", 3*4), rep("Lipids", 4*4), rep("Blood Cell", 15*4))

# change color
# add nloci and vg info

library(ggplot2)
library(hrbrthemes)
library(viridis)
# update variance explained percent format
pct_format = scales::percent_format(accuracy = .1)

plt=ggplot(data_long, aes(fill=Term, y=Percentage, x=trait)) + 
  geom_bar(position="stack", stat="identity") + 
  theme_classic()+
  ylab("Proportion of variance explained") + 
  xlab("Trait") + 
  scale_y_continuous(limits = c(-0.1, 1.1), breaks = seq(0, 1, .25)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    geom_text(
    aes(label = sprintf( '%d \n %s', nloci, pct_format(vg))),
    y=-0.08,  #angle = 50, 
    colour = 'black', size = 2.5)+
    facet_grid(~ factor(group, 
              levels=c("Anthropometric", "Blood Pressure", "Lipids", "Blood Cell")),
             scales = "free_x", space = "free_x", switch = "x") +
    theme(panel.spacing = unit(0, units = "cm"), # removes space between panels
        #strip.placement = "outside", # moves the states down
        strip.text.x = element_text(size = 8),
        strip.background = element_rect(fill = NA) ) +# removes the background from the names
guides(fill=guide_legend(title="Variance\nComponent")) +
 scale_fill_viridis(discrete=TRUE, labels=c('(1.1)', '(1.2)', '(1.3)', '(1.4)')) 

ggsave("../figs/traits_vgprop_multi.png", plot=plt,
       width = 10, height = 5.5, dpi = 300, units = "in", 
       device='png')

