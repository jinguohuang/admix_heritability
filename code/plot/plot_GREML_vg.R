#Plot GREML results

# CGF
filename="../data/admix_CGF_GREML_vg.txt"
df_CGF=read.table(filename, header=T)
# HI
filename="../data/admix_HI_GREML_vg.txt"
df_HI=read.table(filename, header=T)

# A function factory for getting integer y-axis values.
integer_breaks <- function(n = 5, ...) {
  fxn <- function(x) {
    breaks <- floor(pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
  return(fxn)
}


# plot gcta estimate and expected 
fig4A=function(data, yobs, yexp, ylab, legend.position="right"){
  library(ggplot2)
  ggplot() +
  geom_line(data=data, linewidth=0.9, alpha = 0.65, #transparent this line
            aes(x=t, y = yobs,                               
                linetype = "obs",
                group=interaction(P, cov), 
                color=interaction(P, cov))) +
  geom_line(data=data, linewidth=0.9, 
            aes(x=t, y = yexp, 
                linetype = "exp",
                group=interaction(P, cov),
                color=interaction(P, cov))) +
  scale_y_log10(limits=c(0.9, 2.1)) +
  scale_linetype_manual("", 
                       breaks = c("obs",   "exp"),
                       values = c("solid",  "11"),
                       labels = c("Estimated\n(standard)",
                                  expression(V[g])
                                  )) +
  scale_colour_manual("", 
                      values = c('#92c5de','#053061',
                                 '#f4a582','#67001f') ,
                      labels = c("P=0", "P=0.9", 
                                 "P=0", "P=0.9"
                                 )) +
  theme_classic() +
  xlab("t") +
  ylab(ylab)  + 
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


HIa_wo=fig4A(data=df_HI, 
          yobs = df_HI$vg_gcta, 
          yexp = df_HI$exp.vg,
          ylab = expression(hat(sigma)[u]^2),
          legend.position = "none") 

CGFa_wo=fig4A(data=df_CGF, 
           yobs = df_CGF$vg_gcta, 
           yexp = df_CGF$exp.vg,
           ylab = expression(hat(sigma)[u]^2),
          legend.position = "right") 

HIa_wganc=fig4A(data=df_HI, 
          yobs = df_HI$vg_gcta_wganc, 
          yexp = df_HI$exp.vg,
          ylab = expression(hat(sigma)[u]^2),
          legend.position = "none") 
 
CGFa_wganc=fig4A(data=df_CGF, 
           yobs = df_CGF$vg_gcta_wganc, 
           yexp = df_CGF$exp.vg,
           ylab = expression(hat(sigma)[u]^2),
          legend.position = "right") 

fig4B=function(data, yobs, legend.position="right"){
library(ggplot2)
ggplot() +
  geom_line(data=data, linewidth=0.9, color="black", 
            aes(x=t, y = va.term1+va.term2,
                linetype = "exp",
                group=interaction(P, cov)
                                     )) + 
  geom_line(data=data, linewidth=0.9, alpha = 0.65, #transparent this line
            aes(x=t, y = yobs,
                linetype = "obs",
                group=interaction(P, cov),
                color=interaction(P, cov))) +
  
   ylim(c(0.93, 1.01))+
  scale_linetype_manual("", 
                       breaks = c("obs",   "exp"),
                       values = c("solid",  "11"),
                       labels = c("Estimated\n(standard)",
                                  "(1.1)+(1.2)")) + 
  scale_colour_manual("", 
                      values = c('#92c5de','#053061',
                                 '#f4a582','#67001f') ,
                      labels = c("P=0", "P=0.9", 
                                 "P=0", "P=0.9"
                                 )) +
  theme_classic() +
  xlab("t") +
  ylab(expression(hat(sigma)[u]^2))  + 
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
HIb_wo=fig4B(data=df_HI, 
             yobs = df_HI$vg_gcta,
             legend.position = "none")

CGFb_wo=fig4B(data=df_CGF, 
              yobs=df_CGF$vg_gcta,
          legend.position = "right")

HIb_wganc=fig4B(data=df_HI,
                yobs = df_HI$vg_gcta_wganc, 
                legend.position = "none")
CGFb_wganc=fig4B(data=df_CGF, 
                 yobs = df_CGF$vg_gcta_wganc,
          legend.position = "right")


fig4C=function(data, yobs, legend.position="right", y0.9=1.05){
  library(ggplot2)
  ggplot() +
     geom_line(data=data, linewidth=0.9, color="black", 
            aes(x=t, y = va.term1+va.term2+va.term3, 
                linetype = "exp",
                group=interaction(P, cov)#,
                #color=interaction(P, cov)
                )) +
  geom_line(data=data, linewidth=0.9, alpha = 0.65, #transparent this line
            aes(x=t, y = yobs,                               
                linetype = "obs",
                group=interaction(P, cov), 
                color=interaction(P, cov))) +
 
    annotate(geom = "text", x=16, y=y0.9, label="P=0.9") +
    annotate(geom = "text", x=16, y=0.985, label="P=0") +
    ylim(c(0.93, 1.08))+
  scale_linetype_manual("", 
                       breaks = c("obs",   "exp"),
                       values = c("solid",  "11"),
                       labels = c("Estimated\n(V(x) scaled)",
                                  "(1.1)+(1.2)\n+(1.3)")) +
  scale_colour_manual("", 
                      values = c('#92c5de','#053061',
                                 '#f4a582','#67001f') ,
                     labels = c("P=0", "P=0.9", 
                                 "P=0", "P=0.9"
                                 )) +
  theme_classic() +
  xlab("t") +
  ylab(expression(hat(sigma)[u]^2))  + 
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


HIc_wo=fig4C(data=df_HI, 
          yobs = df_HI$vg_GRMvarX,
          legend.position = "none")
CGFc_wo=fig4C(data=df_CGF,
           yobs = df_CGF$vg_GRMvarX,
           y0.9 = 1.06,
          legend.position = "right")

HIc_wganc=fig4C(data=df_HI, 
          yobs = df_HI$vg_GRMvarX_ganc,
          legend.position = "none")
CGFc_wganc=fig4C(data=df_CGF,
           yobs = df_CGF$vg_GRMvarX_ganc,
           y0.9 = 1.06,
          legend.position = "right")


# LD matrix scaled results panel D
fig4D=function(data, yobs, yexp, ylab, legend.position="right"){
  library(ggplot2)
  ggplot() +
  geom_line(data=data, linewidth=0.9, alpha = 0.65, #transparent this line
            aes(x=t, y = yobs,                               
                linetype = "obs",
                group=interaction(P, cov), 
                color=interaction(P, cov))) +
  geom_line(data=data, linewidth=0.9, 
            aes(x=t, y = yexp, 
                linetype = "exp",
                group=interaction(P, cov),
                color=interaction(P, cov))) +
  scale_y_log10(limits=c(0.9, 2.1)) +
  scale_linetype_manual("", 
                       breaks = c("obs",   "exp"),
                       values = c("solid",  "11"),
                       labels = c("Estimated\n(LD scaled)",
                                  #"Expected"
                                  expression(V[g])
                                  )) +
  scale_colour_manual("", 
                      values = c('#92c5de','#053061',
                                 '#f4a582','#67001f') ,
                      labels = c("P=0", "P=0.9", 
                                 "P=0", "P=0.9"
                                 )) +
  theme_classic() +
  xlab("t") +
  ylab(ylab)  + 
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



HId_wo=fig4D(data=df_HI, 
          yobs = df_HI$vg_GRMld, 
          yexp = df_HI$exp.vg,
          ylab = expression(hat(sigma)[u]^2),
          legend.position = "none") 

CGFd_wo=fig4D(data=df_CGF, 
           yobs = df_CGF$vg_GRMld, 
           yexp = df_CGF$exp.vg,
           ylab = expression(hat(sigma)[u]^2),
          legend.position = "right") 

fig4D_w=function(data, yobs, yexp, ylab, legend.position="right", y0.9=0.93){
  library(ggplot2)
  ggplot() +
  geom_line(data=data, linewidth=0.9, color="black", 
            aes(x=t, y = yexp, 
                linetype = "exp",
                group=interaction(P, cov))) +
    
  geom_line(data=data, linewidth=0.9, alpha = 0.65, #transparent this line
            aes(x=t, y = yobs,                               
                linetype = "obs",
                group=interaction(P, cov), 
                color=interaction(P, cov))) +
  
    annotate(geom = "text", x=16, y=y0.9, label="P=0.9") +
    annotate(geom = "text", x=16, y=1.03, label="P=0") +
  scale_y_log10(limits=c(0.9, 1.1)) +
  scale_linetype_manual("", 
                       breaks = c("obs",   "exp"),
                       values = c("solid",  "11"),
                       labels = c("Estimated\n(LD scaled)",
                                  "(1.1)+(1.2)\n-(1.3)")) +
  scale_colour_manual("", 
                      values = c('#92c5de','#053061',
                                 '#f4a582','#67001f') ,
                     labels = c("P=0", "P=0.9", 
                                "P=0", "P=0.9"
                                 )) +
  theme_classic() +
  xlab("t") +
  ylab(ylab)  + 
  #ggtitle(title) +
  theme(aspect.ratio = 1, 
        plot.title = element_text(hjust = 0.5), # to center the title
        legend.position = legend.position, 
        legend.text.align = 0, #left align legend
        text = element_text(size = 12),
        plot.margin = unit(c(0, 0, 0, 0), 'cm')
        ) + 
  guides(color = "none", # no show color legend
          linetype = guide_legend(order = 2, reverse = T)#, 
         )
}
HId_wganc=fig4D_w(data=df_HI, 
          yobs = df_HI$vg_GRMld_ganc, 
          yexp = df_HI$va.term1 + df_HI$va.term2 - df_HI$va.term3, # not this
          #title="Hybrid Isolation wganc",
          ylab = expression(hat(sigma)[u]^2),
          legend.position = "none", y0.9=0.95) 
 
CGFd_wganc=fig4D_w(data=df_CGF, 
           yobs = df_CGF$vg_GRMld_ganc, 
           yexp = df_CGF$va.term1 + df_CGF$va.term2 - df_CGF$va.term3, 
           ylab = expression(hat(sigma)[u]^2),
          legend.position = "right") 


#color legend only
fig4color=function(data){
  library(ggplot2)
ggplot() +
  geom_line(data=data, linewidth=0.9, alpha = 0.65, #transparent this line
            aes(x=t, y = vg_gcta,                               
                linetype = "obs",
                group=interaction(P, cov), 
                color=interaction(P, cov))) +
  geom_line(data=data, linewidth=0.9, 
            aes(x=t, y = exp.vg, 
                linetype = "exp",
                group=interaction(P, cov),
                color=interaction(P, cov))) +
  scale_y_log10(limits=c(0.9, 2.1)) +
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
  ylab(expression(hat(V)[g]))  + 
  theme(aspect.ratio = 1, 
        legend.position = "bottom", #hide legend for HI
        text = element_text(size = 12)
        ) + 
  guides(color =  guide_legend(order = 1, nrow=2, 
                              byrow = T, reverse = T,
                              override.aes = list(linewidth = 2)),  # thicken the line in legend
          linetype = "none" 
         )}
HI_color=fig4color(df_HI)



# GREML Vg wo and wganc
library(ggpubr)
plt_wo=ggarrange(HIa_wo, CGFa_wo, 
              HIb_wo, CGFb_wo, 
              HIc_wo, CGFc_wo, 
              HId_wo, CGFd_wo,
              ncol = 2, nrow = 4, 
              labels = c("A", "", "B", "","C", "","D", ""),
              align = c("h"))

plt_wganc=ggarrange(HIa_wganc, CGFa_wganc, 
              HIb_wganc, CGFb_wganc, 
              HIc_wganc, CGFc_wganc, 
              HId_wganc, CGFd_wganc,
              ncol = 2, nrow = 4, 
              labels = c("E", "", "F", "","G", "","H", ""),
              align = c("h"))

# add space between these two
plt=ggarrange(plt_wo, '', plt_wganc, 
              ncol = 3, nrow = 1, 
              widths = c(8,1,8),
              align = c("v")) %>% # to move the legend closer 
  gridExtra::grid.arrange(get_legend(HI_color), 
                          heights = unit(c(9, 0.8), "in")
                          ) 

ggsave("../figs/GREML_vg_wowganc.png", plot=plt,
       width = 16, height = 10, dpi = 300, units = "in", device='png')
