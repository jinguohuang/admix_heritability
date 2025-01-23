library(tidyverse)
library(data.table)

#load in data
vgCI <- fread("vgpCI.txt", header = TRUE)
head(vgCI)
dim(vgCI)

#melt df
df_long <- vgCI %>%
  pivot_longer(cols = starts_with("vg"), 
               names_to = c("vg", "stat"), 
               names_pattern = "(vg\\d)_(.*)") %>%
  pivot_wider(names_from = stat, values_from = value)

head(df_long)

#keep order
df_long <- 
  arrange(df_long, factor(trait, levels = c("SkinPigmentation","Height","BMI","BirthWeight","SBP","DBP","PP","HDL","LDL","TG","TC","RBC","HGB","HCT","MCH","MCV","MCHC","RDW","WBC","NEU","LYM","MON","BAS","EOS","PLT","MPV"))) %>% 
  mutate(trait=factor(trait, levels=unique(trait)))

# add group
df_long$group=c(rep("Anthropometric", 16), rep("Blood Pressure", 12), rep("Lipids", 16), rep("Blood Cell", 15*4))
head(df_long)

#plot 
plt_bar <- ggplot(df_long, aes(x = trait, y = avg, fill = factor(vg, levels = c("vg1", "vg2",'vg3','vg4')))) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) + 
  geom_errorbar(aes(ymin = CI95l, ymax = CI95r), 
                position = position_dodge(width = 0.9), width = 0.25) +
  geom_text(aes(label = nloci), y=-0.12, colour = 'black', size = 3)+
  scale_y_continuous(limits = c(-0.1, 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_classic() +
  ylab("Proportion of genetic variance") +    
  xlab("Trait") +
  scale_fill_manual(values = c("vg1" = "steelblue", "vg2" = "cornflowerblue", 
                               "vg3" = "steelblue2", "vg4" = "firebrick"), 
                    name = "Variance\nComponent", 
                    labels = c("(1.1)", "(1.2)",'(1.3)','(1.4)')) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=12),
        axis.text.y = element_text(angle = 0, hjust = 1, size =12),  
        axis.title.x = element_text(hjust = 0.5, size = 16), 
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 12)) +
  facet_grid(~ factor(group, 
                      levels=c("Anthropometric", "Blood Pressure", "Lipids", "Blood Cell")),
             scales = "free_x", space = "free_x", switch = "x") +
  theme(panel.spacing = unit(0, units = "cm"), 
        strip.text.x = element_text(size = 12),
        strip.background = element_rect(fill = NA))

plt_bar

#save plot
ggsave("vgp_terms_CI_multi.png", plot=plt_bar,
       width = 15, height = 8, dpi = 300, units = "in", 
       device='png')
