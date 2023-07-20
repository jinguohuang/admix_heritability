# plot GREML estimated vgamma: genetic variance due to local ancestry

# load files
library(data.table)
df_HI=fread("../data/admix_HI_GREML_vgamma.txt", header = T)
df_CGF=fread("../data/admix_CGF_GREML_vgamma.txt", header = T)

# function to plot 
# expectation VS standard scaled GCTA results
fig4A=function(data, yobs, title, legend.position="right"){
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
                                     "Expected")) +
    scale_colour_manual("", 
                        values = c('#92c5de','#053061',
                                   '#f4a582','#67001f') ,
                        labels = c("P=0", "P=0.9", 
                                   "P=0", "P=0.9"
                        )) +
    theme_classic() +
    xlab("t") +
    ylab(expression(hat(V)[gamma]))  + 
    ggtitle(title) +
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
              aes(x=t, y = va.term2+va.term3,
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
                                     "(1.2)+(1.3)")) + 
    scale_colour_manual("", 
                        values = c('#92c5de','#053061',
                                   '#f4a582','#67001f') ,
                        labels = c("P=0", "P=0.9", 
                                   "P=0", "P=0.9"
                        )) +
    theme_classic() +
    xlab("t") +
    ylab(expression(hat(V)[gamma]))  + 
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
                          labels = c("Estimated\n(adjusted)",
                                     "(1.2)+(1.3)")) +
    scale_colour_manual("", 
                        values = c('#92c5de','#053061',
                                   '#f4a582','#67001f') ,
                        labels = c("P=0", "P=0.9", 
                                   "P=0", "P=0.9"
                        )) +
    theme_classic() +
    xlab("t") +
    ylab(expression(hat(V)[gamma]))  + 
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



#color legend 
fig4color=function(data){
  library(ggplot2)
  ggplot() +
    geom_line(data=data, linewidth=0.9, alpha = 0.65, #transparent this line
              aes(x=t, y = vgamma_gcta,                               
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
    ylab(expression(V[gamma]))  + 
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
          yobs=df_HI$vgamma_gcta, 
          title="Hybrid Isolation", 
          legend.position = "none")
CGFa=fig4A(data=df_CGF, 
           yobs=df_CGF$vgamma_gcta, 
           title="Continuous Gene Flow")

HIb=fig4B(data=df_HI,
          yobs=df_HI$vgamma_gcta,
          legend.position = "none")

CGFb=fig4B(data=df_CGF,
           yobs=df_CGF$vgamma_gcta)

HIc=fig4C(data=df_HI, 
          yobs=df_HI$vgamma_GRMvarX,
          legend.position = "none")
CGFc=fig4C(data=df_CGF,
           yobs=df_CGF$vgamma_GRMvarX)


library(ggpubr)
plt=ggarrange(HIa, CGFa, HIb, CGFb, HIc, CGFc, 
              ncol = 2, nrow = 3, 
              labels = c("A", "", "B", "","C", ""),
              align = c("h")) %>% # to move the legend closer 
  gridExtra::grid.arrange(get_legend(HI_color), 
                          #heights = unit(c(250, 10), "mm")
                          heights = unit(c(8, 0.8), "in")
  ) 

ggsave("../figs/GREML_vgamma.png", plot=plt,
       width = 9, height = 9, dpi = 300, units = "in", device='png')



# plot for gcta_wganc
HIa_w=fig4A(data=df_HI, 
            yobs = df_HI$vgamma_GCTA_wganc,
            title="Hybrid Isolation", 
            legend.position = "none")
CGFa_w=fig4A(data=df_CGF, 
             yobs=df_CGF$vgamma_GCTA_wganc,
             title="Continuous Gene Flow")

HIb_w=fig4B(data=df_HI, 
            yobs=df_HI$vgamma_GCTA_wganc,
            legend.position = "none")
CGFb_w=fig4B(data=df_CGF,
             yobs=df_CGF$vgamma_GCTA_wganc)


HIc_w=fig4C(data=df_HI,
            yobs=df_HI$vgamma_GRMvarX_ganc,
            legend.position = "none")
CGFc_w=fig4C(data=df_CGF,
             yobs=df_CGF$vgamma_GRMvarX_ganc)


plt=ggarrange(HIa_w, CGFa_w, HIb_w, CGFb_w, HIc_w, CGFc_w, 
              ncol = 2, nrow = 3, 
              labels = c("A", "", "B", "","C", ""),
              align = c("h")) %>% # to move the legend closer 
  gridExtra::grid.arrange(get_legend(HI_color), 
                          #heights = unit(c(250, 10), "mm")
                          heights = unit(c(8, 0.8), "in")
  ) 

ggsave("../figs/GREML_vgamma_wganc.png", plot=plt,
       width = 9, height = 9, dpi = 300, units = "in", device='png')

