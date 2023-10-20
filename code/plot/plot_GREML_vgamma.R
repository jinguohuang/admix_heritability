# plot GREML estimated vgamma: genetic variance due to local ancestry

# load files
library(data.table)
df_HI=fread("../data/admix_HI_GREML_vgamma.txt", header = T)
df_CGF=fread("../data/admix_CGF_GREML_vgamma.txt", header = T)

# function to plot 
# expectation VS standard scaled GCTA results
fig4A=function(data, yobs, legend.position="right"){
  library(ggplot2)
  ggplot() +
    geom_line(data=data, linewidth=0.9, alpha = 0.65, #transparent this line
              aes(x=t, y = yobs,                               
                  linetype = "obs",
                  group=interaction(P, cov), 
                  color=interaction(P, cov))) +
    geom_line(data=data, linewidth=0.9, 
              aes(x=t, y = va.term2+va.term3+va.term4, 
                  linetype = "exp",
                  group=interaction(P, cov),
                  color=interaction(P, cov))) +
    ylim(c(-0.05, 1.20))+
    scale_linetype_manual("", 
                          breaks = c("obs",   "exp"),
                          values = c("solid",  "11"),
                          labels = c("Estimated\n(standard)",
                                     expression(V[gamma]))) +
    scale_colour_manual("", 
                        values = c('#92c5de','#053061',
                                   '#f4a582','#67001f') ,
                        labels = c("P=0", "P=0.9", 
                                   "P=0", "P=0.9"
                        )) +
    theme_classic() +
    xlab("t") +
    ylab(expression(hat(sigma)[v]^2))  + 
    #ggtitle(title) +
    theme(aspect.ratio = 1, 
          plot.title = element_text(hjust = 0.5), # to center the title
          legend.position = legend.position, 
          legend.text.align = 0, #left align legend
          text = element_text(size = 12),
          plot.margin = unit(c(0, 0, 0, 0), 'cm')
    ) + 
    guides(color = "none", # no show color legend
           linetype = guide_legend(order = 2, reverse = T)
    )
}


# standard scale gcta vgamma vs expected (1.2)+(1.3)
fig4B=function(data, yobs, legend.position="right"){
  library(ggplot2)
  ggplot() +
    geom_line(data=data, linewidth=0.9, color="black", 
              aes(x=t, y = va.term2,
                  linetype = "exp",
                  group=interaction(P, cov)
              )) + 
    geom_line(data=data, linewidth=0.9, alpha = 0.65, #transparent this line
              aes(x=t, y = yobs,
                  linetype = "obs",
                  group=interaction(P, cov),
                  color=interaction(P, cov))) +
    ylim(c(-0.05, 0.15))+
    scale_linetype_manual("", 
                          breaks = c("obs",   "exp"),
                          values = c("solid",  "11"),
                          labels = c("Estimated\n(standard)",
                                     "(1.2)")) + 
    scale_colour_manual("", 
                        values = c('#92c5de','#053061',
                                   '#f4a582','#67001f') ,
                        labels = c("P=0", "P=0.9", 
                                   "P=0", "P=0.9"
                        )) +
    theme_classic() +
    xlab("t") +
    ylab(expression(hat(sigma)[v]^2))  + 
    theme(aspect.ratio = 1, 
          plot.title = element_text(hjust = 0.5), # to center the title
          legend.position = legend.position, 
          legend.text.align = 0, #left align legend
          text = element_text(size = 12),
          plot.margin = unit(c(0, 0, 0, 0), 'cm')
    ) + 
    guides(color = "none", 
           linetype = guide_legend(order = 1, reverse = T
           ) )}


# (1.2)+(1.3) vs GRMvarX scaled vgamma gcta estimates

fig4C=function(data, yobs, legend.position="right"){
  library(ggplot2)
  ggplot() +
    geom_line(data=data, linewidth=0.9, color="black", 
              aes(x=t, y = va.term2+va.term3, 
                  linetype = "exp",
                  group=interaction(P, cov)
              )) +
    geom_line(data=data, linewidth=0.9, alpha = 0.65, #transparent this line
              aes(x=t, y = yobs,                               
                  linetype = "obs",
                  group=interaction(P, cov), 
                  color=interaction(P, cov))) +
    # label P
    annotate(geom = "text", x=16, y=0.12, label="P=0.9") +
    annotate(geom = "text", x=16, y=0.05, label="P=0") +
    ylim(c(0, 0.25))+
    scale_linetype_manual("", 
                          breaks = c("obs",   "exp"),
                          values = c("solid",  "11"),
                          labels = c(#expression(paste("Estimated\n(V(", gamma, ") scaled)", sep ="")),
                                     paste("Estimated\n(V(\u03B3) scaled)"),
                                     "(1.2)+(1.3)")) +
    scale_colour_manual("", 
                        values = c('#92c5de','#053061',
                                   '#f4a582','#67001f') ,
                        labels = c("P=0", "P=0.9", 
                                   "P=0", "P=0.9"
                        )) +
    theme_classic() +
    xlab("t") +
    ylab(expression(hat(sigma)[v]^2))  + 
    theme(aspect.ratio = 1, 
          plot.title = element_text(hjust = 0.5), # to center the title
          legend.position = legend.position, 
          legend.text.align = 0, #left align legend
          text = element_text(size = 12),
          plot.margin = unit(c(0, 0, 0, 0), 'cm')
    ) + 
    guides(color = "none", # no show color legend
           linetype = guide_legend(order = 2, reverse = T
           ) )
}


fig4D=function(data, yobs, legend.position="right"){
  library(ggplot2)
  ggplot() +
    geom_line(data=data, linewidth=0.9, 
            aes(x=t, y = va.term2+va.term3+va.term4, 
                linetype = "exp",
                group=interaction(P, cov),
                color=interaction(P, cov))) +
  geom_line(data=data, linewidth=0.9, alpha = 0.65, #transparent this line
            aes(x=t, y = yobs,                               
                linetype = "obs",
                group=interaction(P, cov), 
                color=interaction(P, cov))) +
  
  ylim(c(-0.05, 1.20))+
  scale_linetype_manual("", 
                       breaks = c("obs",   "exp"),
                       values = c("solid",  "11"),
                       labels = c("Estimated\n(LD scaled)",
                                  expression(V[gamma]))) +
  scale_colour_manual("", 
                     values = c('#92c5de','#053061',
                                '#f4a582','#67001f'),
  ) +
  theme_classic() +
  xlab("t") +
  ylab(expression(hat(sigma)[v]^2))  + 
  theme(aspect.ratio = 1, 
        plot.title = element_text(hjust = 0.5), # to center the title
        legend.position = legend.position, 
        legend.text.align = 0, #left align legend
        text = element_text(size = 12),
        plot.margin = unit(c(0, 0, 0, 0), 'cm')
        ) + 
  guides(color = "none", # no show color legend
          linetype = guide_legend(order = 2, reverse = T)
         )
}


fig4D_w=function(data, yobs, legend.position="right"){
  library(ggplot2)
  ggplot() +
    geom_line(data=data, linewidth=0.9, #color="black",
            aes(x=t, y = va.term2-va.term3, 
                linetype = "exp",
                group=interaction(P, cov),
                color=interaction(P, cov) #in order to keep this order
                )) +
  geom_line(data=data, linewidth=0.9, alpha = 0.65, #transparent this line
            aes(x=t, y = yobs,                               
                linetype = "obs",
                group=interaction(P, cov), 
                color=interaction(P, cov))) +
    
    geom_line(data=data, linewidth=0.9, color="black", # cover it
            aes(x=t, y = va.term2-va.term3, 
                linetype = "exp",
                group=interaction(P, cov),
               # color=interaction(P, cov)
                )) +
     # label P
  annotate(geom = "text", x=16, y=0.027, label="P=0.9") +
  annotate(geom = "text", x=16, y=0.07, label="P=0") +
  ylim(c(-0.01, 0.1))+
  scale_linetype_manual("", 
                       breaks = c("obs",   "exp"),
                       values = c("solid",  "11"),
                       labels = c("Estimated\n(LD scaled)",
                                  "(1.2)-(1.3)")) +
  scale_colour_manual("", 
                     values = c('#92c5de','#053061',
                                '#f4a582','#67001f'),
  ) +
  theme_classic() +
  xlab("t") +
  ylab(expression(hat(sigma)[v]^2))  + 
  theme(aspect.ratio = 1, 
        plot.title = element_text(hjust = 0.5), # to center the title
        legend.position = legend.position, 
        legend.text.align = 0, #left align legend
        text = element_text(size = 12),
        plot.margin = unit(c(0, 0, 0, 0), 'cm')
        ) + 
  guides(color = "none", 
          linetype = guide_legend(order = 2, reverse = T)
         )
}


#color legend 
fig4color=function(data){
  library(ggplot2)
  ggplot() +
    geom_line(data=data, linewidth=0.9, alpha = 0.65, #transparent this line
              aes(x=t, y = vgamma_GRMstd,                               
                  linetype = "obs",
                  group=interaction(P, cov), 
                  color=interaction(P, cov))) +
    geom_line(data=data, linewidth=0.9, 
              aes(x=t, y = va.term2+va.term3+va.term4, 
                  linetype = "exp",
                  group=interaction(P, cov),
                  color=interaction(P, cov))) +
    scale_linetype_manual("", 
                          breaks = c("obs",   "exp"),
                          values = c("solid",  "11"),
                          labels = c("Estimated\n(standard)",
                                     "Expected")) +
    scale_colour_manual("", 
                        values = c('#92c5de','#053061',
                                   '#f4a582','#67001f') ,
                        labels = c("P=0", "P=0.9", 
                                   "P=0", "P=0.9"
                        )) +
    theme_classic() +
    xlab("t") +
    ylab(expression(hat(sigma)[v]^2))  + 
    theme(aspect.ratio = 1, 
          legend.position = "bottom", #hide legend for HI
          text = element_text(size = 12)
    ) + 
    guides(color =  guide_legend(order = 1, nrow=2, 
                                 byrow = T, reverse = T,
                                 override.aes = list(linewidth = 2)),  # thicken the line in legend
           linetype = "none"
    )}

# plot
HI_color=fig4color(data=df_HI)


HIa=fig4A(data=df_HI, 
          yobs=df_HI$vgamma_GRMstd, 
          #title="Hybrid Isolation", 
          legend.position = "none")
CGFa=fig4A(data=df_CGF, 
           yobs=df_CGF$vgamma_GRMstd) 
           #title="Continuous Gene Flow")

HIb=fig4B(data=df_HI,
          yobs=df_HI$vgamma_GRMstd,
          legend.position = "none")

CGFb=fig4B(data=df_CGF,
           yobs=df_CGF$vgamma_GRMstd)

HIc=fig4C(data=df_HI, 
          yobs=df_HI$vgamma_GRMvarX,
          legend.position = "none")
CGFc=fig4C(data=df_CGF,
           yobs=df_CGF$vgamma_GRMvarX)

HId=fig4D(data=df_HI, 
          yobs=df_HI$vgamma_GRMld,
          legend.position = "none")
CGFd=fig4D(data=df_CGF,
           yobs=df_CGF$vgamma_GRMld)

library(ggpubr)
plt_wo=ggarrange(HIa, CGFa, HIb, CGFb, HIc, CGFc, HId, CGFd,
              ncol = 2, nrow = 4, 
              labels = c("A", "", "B", "","C", "","D", ""),
              align = c("h")) 

# plot for gcta_wganc
HIa_w=fig4A(data=df_HI, 
            yobs = df_HI$vgamma_GRMstd_ganc,
            #title="Hybrid Isolation", 
            legend.position = "none")
CGFa_w=fig4A(data=df_CGF, 
             yobs=df_CGF$vgamma_GRMstd_ganc)
             #title="Continuous Gene Flow")

HIb_w=fig4B(data=df_HI, 
            yobs=df_HI$vgamma_GRMstd_ganc,
            legend.position = "none")
CGFb_w=fig4B(data=df_CGF,
             yobs=df_CGF$vgamma_GRMstd_ganc)

HIc_w=fig4C(data=df_HI,
            yobs=df_HI$vgamma_GRMvarX_ganc,
            legend.position = "none")
CGFc_w=fig4C(data=df_CGF,
             yobs=df_CGF$vgamma_GRMvarX_ganc)

HId_w=fig4D_w(data=df_HI,
            yobs=df_HI$vgamma_GRMld_ganc,
            legend.position = "none")
CGFd_w=fig4D_w(data=df_CGF,
             yobs=df_CGF$vgamma_GRMld_ganc)


plt_wganc=ggarrange(HIa_w, CGFa_w, HIb_w, CGFb_w, HIc_w, CGFc_w, HId_w, CGFd_w, 
              ncol = 2, nrow = 4, 
              labels = c("E", "", "F", "","G", "","H", ""),
              align = c("h")) 

# add space between these two
plt=ggarrange(plt_wo, '', plt_wganc, 
              ncol = 3, nrow = 1, 
             # labels = c("A", "B",  "C", "D"),
              widths = c(8,1,8),
              align = c("v")) %>% # to move the legend closer 
  gridExtra::grid.arrange(get_legend(HI_color), 
                          #heights = unit(c(250, 10), "mm")
                          heights = unit(c(9, 0.8), "in")
                          ) 

ggsave("../figs/GREML_vgamma_wowganc.png", plot=plt,
       width = 16, height = 10, dpi = 300, units = "in", device='png')
