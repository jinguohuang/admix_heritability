#Plot GREML results
library(data.table)

# CGF
filename="admix_CGF_HE_vg_CI_P0_P9_posneg_grms.txt"
df_CGF_CI=read.table(filename, header=T)
# HI
filename="admix_HI_HE_vg_CI_P0_P9_posneg_grms.txt"
df_HI_CI=read.table(filename, header=T)


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
fig4Awo=function(data, yobs, yexp, ylab, 
               CIl, CIr,
               legend.position="right"){
  library(ggplot2)
  ggplot() +
    geom_ribbon(data=data, alpha=0.2, linetype = 0, #remove the boarder
              aes(x=t, 
                  ymin = CIl, ymax = CIr,
                  group=interaction(P, cov),
                fill=interaction(P, cov))) +
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

  scale_y_continuous(limits=c(0, 6.1)) +
  scale_linetype_manual("", 
                       breaks = c("obs",   "exp"),
                       values = c("solid",  "11"),
                       labels = c("Estimated\n(standard)",
                                  expression(V[g])
                                  )) +
  scale_colour_manual("", 
                      values = c('#92c5de','#053061',
                                 '#f4a582','#67001f',
                                 '#79bc5c','#2a6112') ,
                      labels = c("P=0", "P=0.9", 
                                 "P=0", "P=0.9",
                                 "P=0", "P=0.9"
                                 )) +
    scale_fill_manual("", 
                      values = c('#92c5de','#053061',
                                 '#f4a582','#67001f',
                                 '#79bc5c','#2a6112') ,
                      labels = c("P=0", "P=0.9", 
                                 "P=0", "P=0.9",
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
  guides(color = "none", fill = "none",# no show color legend
          linetype = guide_legend(order = 2, reverse = T)
         )
}



HIa_wo=fig4Awo(data=df_HI_CI, 
          yobs = df_HI_CI$vg_HEreg_noa.mean, 
          yexp = df_HI_CI$exp.vg,
          ylab = expression(hat(sigma)[u]^2),
          CIl = df_HI_CI$vg_HEreg_noa.CI95l,
          CIr = df_HI_CI$vg_HEreg_noa.CI95r,
          legend.position = "none") 

print(HIa_wo)

CGFa_wo=fig4Awo(data=df_CGF_CI, 
          yobs = df_CGF_CI$vg_HEreg_noa.mean, 
          yexp = df_CGF_CI$exp.vg,
          ylab = expression(hat(sigma)[u]^2),
          CIl = df_CGF_CI$vg_HEreg_noa.CI95l,
          CIr = df_CGF_CI$vg_HEreg_noa.CI95r,
          legend.position = "right") 

#fig 4A with ganc
fig4Awganc=function(data, yobs, yexp, ylab, 
                 CIl, CIr, y0.9 = 0.3, y0 = 0.87,
                 legend.position="right"){
  library(ggplot2)
  ggplot() +
    geom_ribbon(data=data, alpha=0.2, linetype = 0, #remove the boarder
                aes(x=t, 
                    ymin = CIl, ymax = CIr,
                    group=interaction(P, cov),
                    fill=interaction(P, cov))) +
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
    
    scale_y_continuous(limits=c(0, 2.1)) +
    annotate(geom = "text", x=17, y=y0.9, label="P=0.9") +
    annotate(geom = "text", x=17, y=y0, label="P=0") +
    scale_linetype_manual("", 
                          breaks = c("obs",   "exp"),
                          values = c("solid",  "11"),
                          labels = c("Estimated\n(standard)",
                                     expression(V[g])
                          )) +
    scale_colour_manual("", 
                        values = c('#92c5de','#053061',
                                   '#f4a582','#67001f',
                                   '#79bc5c','#2a6112') ,
                        labels = c("P=0", "P=0.9", 
                                   "P=0", "P=0.9",
                                   "P=0", "P=0.9"
                        )) +
    scale_fill_manual("", 
                      values = c('#92c5de','#053061',
                                 '#f4a582','#67001f',
                                 '#79bc5c','#2a6112') ,
                      labels = c("P=0", "P=0.9", 
                                 "P=0", "P=0.9",
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
    guides(color = "none", fill = "none",# no show color legend
           linetype = guide_legend(order = 2, reverse = T)
    )
}


HIa_wganc=fig4Awganc(data=df_HI_CI, 
          yobs = df_HI_CI$vg_HEreg_adj.mean, 
          yexp = df_HI_CI$exp.vg,
          ylab = expression(hat(sigma)[u]^2),
          CIl = df_HI_CI$vg_HEreg_adj.CI95l,
          CIr = df_HI_CI$vg_HEreg_adj.CI95r,
          legend.position = "none") 
print(HIa_wganc)
 
CGFa_wganc=fig4Awganc(data=df_CGF_CI, 
          yobs = df_CGF_CI$vg_HEreg_adj.mean, 
          yexp = df_CGF_CI$exp.vg,
          ylab = expression(hat(sigma)[u]^2),
          CIl = df_CGF_CI$vg_HEreg_adj.CI95l,
          CIr = df_CGF_CI$vg_HEreg_adj.CI95r,
          y0.9 = 0.2, y0 = 0.65,
          legend.position = "right") 
print(CGFa_wganc)


fig4Bwo=function(data, yobs, CIl, CIr, legend.position="right"){
library(ggplot2)
ggplot() +
    geom_ribbon(data=data, alpha=0.2, linetype = 0,
              aes(x=t, 
                  ymin = CIl, ymax = CIr,
                  group=interaction(P, cov),
                fill=interaction(P, cov))) +
  geom_line(data=data, linewidth=0.9, 
            aes(x=t, y = exp.vg,
                linetype = "exp",
                group=interaction(P, cov),
                color=interaction(P, cov)
                                     )) + 
  geom_line(data=data, linewidth=0.9, alpha = 0.65, #transparent this line
            aes(x=t, y = yobs,
                linetype = "obs",
                group=interaction(P, cov),
                color=interaction(P, cov))) +

  ylim(c(0, 6.1))+
  scale_linetype_manual("", 
                       breaks = c("obs",   "exp"),
                       values = c("solid",  "11"),
                       labels = c("Estimated\n(standard)",
                                  expression(V[g]))) + 
  scale_colour_manual("", 
                      values = c('#92c5de','#053061',
                                 '#f4a582','#67001f',
                                 '#79bc5c','#2a6112') ,
                      labels = c("P=0", "P=0.9", 
                                 "P=0", "P=0.9",
                                 "P=0", "P=0.9"
                                 )) +
  scale_fill_manual("", 
                      values = c('#92c5de','#053061',
                                 '#f4a582','#67001f',
                                 '#79bc5c','#2a6112') ,
                      labels = c("P=0", "P=0.9", 
                                 "P=0", "P=0.9",
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
  guides(color = "none", fill = "none",
          linetype = guide_legend(order = 1, reverse = T
          ) )}


HIb_wo=fig4Bwo(data=df_HI_CI, 
             yobs = df_HI_CI$vg_HEreg_noa.mean,
             CIl = df_HI_CI$vg_HEreg_noa.CI95l,
             CIr = df_HI_CI$vg_HEreg_noa.CI95r,
             legend.position = "none")
print(HIb_wo)

CGFb_wo=fig4Bwo(data=df_CGF_CI, 
             yobs = df_CGF_CI$vg_HEreg_noa.mean,
             CIl = df_CGF_CI$vg_HEreg_noa.CI95l,
             CIr = df_CGF_CI$vg_HEreg_noa.CI95r,
             legend.position = "right")
print(CGFb_wo)


fig4Bwganc=function(data, yobs, CIl, CIr, legend.position="right", y0.9=0.3, y0 = 0.87){
  library(ggplot2)
  ggplot() +
    geom_ribbon(data=data, alpha=0.2, linetype = 0,
                aes(x=t, 
                    ymin = CIl, ymax = CIr,
                    group=interaction(P, cov),
                    fill=interaction(P, cov))) +
    geom_line(data=data, linewidth=0.9, 
              aes(x=t, y = exp.vg,
                  linetype = "exp",
                  group=interaction(P, cov),
                  color=interaction(P, cov)
              )) + 
    geom_line(data=data, linewidth=0.9, alpha = 0.65, #transparent this line
              aes(x=t, y = yobs,
                  linetype = "obs",
                  group=interaction(P, cov),
                  color=interaction(P, cov))) +
    
    ylim(c(0, 2.1))+
    annotate(geom = "text", x=17, y=y0.9, label="P=0.9") +
    annotate(geom = "text", x=17, y=y0, label="P=0") +
    scale_linetype_manual("", 
                          breaks = c("obs",   "exp"),
                          values = c("solid",  "11"),
                          labels = c("Estimated\n(standard)",
                                     expression(V[g]))) + 
    scale_colour_manual("", 
                        values = c('#92c5de','#053061',
                                   '#f4a582','#67001f',
                                   '#79bc5c','#2a6112') ,
                        labels = c("P=0", "P=0.9", 
                                   "P=0", "P=0.9",
                                   "P=0", "P=0.9"
                        )) +
    scale_fill_manual("", 
                      values = c('#92c5de','#053061',
                                 '#f4a582','#67001f',
                                 '#79bc5c','#2a6112') ,
                      labels = c("P=0", "P=0.9", 
                                 "P=0", "P=0.9",
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
    guides(color = "none", fill = "none",
           linetype = guide_legend(order = 1, reverse = T
           ) )}


HIb_wganc=fig4Bwganc(data=df_HI_CI, 
             yobs = df_HI_CI$vg_HEreg_adj.mean,
             CIl = df_HI_CI$vg_HEreg_adj.CI95l,
             CIr = df_HI_CI$vg_HEreg_adj.CI95r,
             legend.position = "none")
print(HIb_wganc)

CGFb_wganc=fig4Bwganc(data=df_CGF_CI, 
             yobs = df_CGF_CI$vg_HEreg_adj.mean,
             CIl = df_CGF_CI$vg_HEreg_adj.CI95l,
             CIr = df_CGF_CI$vg_HEreg_adj.CI95r,
             y0.9 = 0.2, y0 = 0.65,
             legend.position = "right")
print(CGFb_wganc)


fig4Cwo=function(data, yobs, CIl, CIr, 
               legend.position="right", y0.9=5){
  library(ggplot2)
  ggplot() +
    geom_ribbon(data=data, alpha=0.2, linetype = 0,
              aes(x=t, 
                  ymin = CIl, ymax = CIr,
                  group=interaction(P, cov),
                fill=interaction(P, cov))) +
     geom_line(data=data, linewidth=0.9, 
            aes(x=t, y = exp.vg, 
                linetype = "exp",
                group=interaction(P, cov),
                color=interaction(P, cov)
                )) +
  geom_line(data=data, linewidth=0.9, alpha = 0.65, #transparent this line
            aes(x=t, y = yobs,                               
                linetype = "obs",
                group=interaction(P, cov), 
                color=interaction(P, cov))) +
 
    # annotate(geom = "text", x=17, y=y0.9, label="P=0.9") +
    # annotate(geom = "text", x=17, y=1.5, label="P=0") +
    ylim(c(0, 8))+
  scale_linetype_manual("", 
                       breaks = c("obs",   "exp"),
                       values = c("solid",  "11"),
                       labels = c("Estimated\n(V(x) scaled)",
                                  expression(V[g]))) +
  scale_colour_manual("", 
                      values = c('#92c5de','#053061',
                                 '#f4a582','#67001f',
                                 '#79bc5c','#2a6112') ,
                     labels = c("P=0", "P=0.9", 
                                 "P=0", "P=0.9",
                                "P=0", "P=0.9"
                                 )) +
    scale_fill_manual("", 
                      values = c('#92c5de','#053061',
                                 '#f4a582','#67001f',
                                 '#79bc5c','#2a6112') ,
                      labels = c("P=0", "P=0.9", 
                                 "P=0", "P=0.9",
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
  guides(color = "none", fill = "none", # no show color legend
          linetype = guide_legend(order = 2, reverse = T
                              ) )
}


HIc_wo=fig4Cwo(data=df_HI_CI, 
          yobs = df_HI_CI$vg_GRMvarX_noa.mean,
          CIl = df_HI_CI$vg_GRMvarX_noa.CI95l,
          CIr = df_HI_CI$vg_GRMvarX_noa.CI95r,
          legend.position = "none")
print(HIc_wo)

CGFc_wo=fig4Cwo(data=df_CGF_CI, 
          yobs = df_CGF_CI$vg_GRMvarX_noa.mean,
          CIl = df_CGF_CI$vg_GRMvarX_noa.CI95l,
          CIr = df_CGF_CI$vg_GRMvarX_noa.CI95r,
          y0.9 = 6,
          legend.position = "right")
print(CGFc_wo)


fig4Cwganc=function(data, yobs, CIl, CIr, 
                 legend.position="right", y0.9=0.34, y0 = 0.87){
  library(ggplot2)
  ggplot() +
    geom_ribbon(data=data, alpha=0.2, linetype = 0,
                aes(x=t, 
                    ymin = CIl, ymax = CIr,
                    group=interaction(P, cov),
                    fill=interaction(P, cov))) +
    geom_line(data=data, linewidth=0.9,  
              aes(x=t, y = exp.vg, 
                  linetype = "exp",
                  group=interaction(P, cov),
                  color=interaction(P, cov)
              )) +
    geom_line(data=data, linewidth=0.9, alpha = 0.65, #transparent this line
              aes(x=t, y = yobs,                               
                  linetype = "obs",
                  group=interaction(P, cov), 
                  color=interaction(P, cov))) +
    
    annotate(geom = "text", x=16, y=y0.9, label="P=0.9") +
    annotate(geom = "text", x=16, y=y0, label="P=0") +
    ylim(c(0, 2.1))+
    scale_linetype_manual("", 
                          breaks = c("obs",   "exp"),
                          values = c("solid",  "11"),
                          labels = c("Estimated\n(V(x) scaled)",
                                     expression(V[g]))) +
    scale_colour_manual("", 
                        values = c('#92c5de','#053061',
                                   '#f4a582','#67001f',
                                   '#79bc5c','#2a6112') ,
                        labels = c("P=0", "P=0.9", 
                                   "P=0", "P=0.9",
                                   "P=0", "P=0.9"
                        )) +
    scale_fill_manual("", 
                      values = c('#92c5de','#053061',
                                 '#f4a582','#67001f',
                                 '#79bc5c','#2a6112') ,
                      labels = c("P=0", "P=0.9", 
                                 "P=0", "P=0.9",
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
    guides(color = "none", fill = "none", # no show color legend
           linetype = guide_legend(order = 2, reverse = T
           ) )
}


HIc_wganc=fig4Cwganc(data=df_HI_CI, 
          yobs = df_HI_CI$vg_GRMvarX_ganc.mean,
          CIl = df_HI_CI$vg_GRMvarX_ganc.CI95l,
          CIr = df_HI_CI$vg_GRMvarX_ganc.CI95r,
          legend.position = "none")
print(HIc_wganc)

CGFc_wganc=fig4Cwganc(data=df_CGF_CI, 
          yobs = df_CGF_CI$vg_GRMvarX_ganc.mean,
          CIl = df_CGF_CI$vg_GRMvarX_ganc.CI95l,
          CIr = df_CGF_CI$vg_GRMvarX_ganc.CI95r,
          y0.9 = 0.27, y0 = 0.6,
          legend.position = "right")
print(CGFc_wganc)


# LD matrix scaled results panel D
fig4D=function(data, yobs, CIl, CIr, yexp, ylab, legend.position="right"){
  library(ggplot2)
  ggplot() +
  geom_ribbon(data=data, alpha=0.2, linetype = 0,
            aes(x=t, 
                ymin = CIl, ymax = CIr,
                group=interaction(P, cov),
              fill=interaction(P, cov))) +
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
  scale_fill_manual("", 
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
  guides(color = "none", fill = "none",# no show color legend
          linetype = guide_legend(order = 2, reverse = T) 
         )
}



HId_wo=fig4D(data=df_HI_CI, 
          yobs = df_HI_CI$vg_GRMld_noa.mean, 
          yexp = df_HI_CI$exp.vg,
          CIl = df_HI_CI$vg_GRMld_noa.CI95l,
          CIr = df_HI_CI$vg_GRMld_noa.CI95r,
          ylab = expression(hat(sigma)[u]^2),
          legend.position = "none")
print(HId_wo)

CGFd_wo=fig4D(data=df_CGF_CI, 
          yobs = df_CGF_CI$vg_GRMld_noa.mean, 
          yexp = df_CGF_CI$exp.vg,
          CIl = df_CGF_CI$vg_GRMld_noa.CI95l,
          CIr = df_CGF_CI$vg_GRMld_noa.CI95r,
          ylab = expression(hat(sigma)[u]^2),
          legend.position = "right")
print(CGFd_wo)


fig4D_w=function(data, yobs, CIl, CIr, ylab, legend.position="right", y0.9=0.91){
  library(ggplot2)
  ggplot() +
  geom_ribbon(data=data, alpha=0.2, linetype = 0,
            aes(x=t, 
                ymin = CIl, ymax = CIr,
                group=interaction(P, cov),
              fill=interaction(P, cov))) +
  geom_line(data=data, linewidth=0.9, color="black", 
            aes(x=t, y = va.term1+va.term2-va.term3, 
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
  scale_fill_manual("", 
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
  guides(color = "none", fill = "none",# no show color legend
          linetype = guide_legend(order = 2, reverse = T)#, 
         )
}
HId_wganc=fig4D_w(data=df_HI_CI, 
          yobs = df_HI_CI$vg_GRMld_ganc.mean, 
          CIl = df_HI_CI$vg_GRMld_ganc.CI95l,
          CIr = df_HI_CI$vg_GRMld_ganc.CI95r,
          ylab = expression(hat(sigma)[u]^2),
          legend.position = "none", y0.9=0.92) 
print(HId_wganc)

CGFd_wganc=fig4D_w(data=df_CGF_CI, 
          yobs = df_CGF_CI$vg_GRMld_ganc.mean, 
          CIl = df_CGF_CI$vg_GRMld_ganc.CI95l,
          CIr = df_CGF_CI$vg_GRMld_ganc.CI95r,
          ylab = expression(hat(sigma)[u]^2),
          legend.position = "right") 
print(CGFd_wganc)



#color legend only
fig4color=function(data){
  library(ggplot2)
ggplot() +
  geom_line(data=data, linewidth=0.9, alpha = 0.65, #transparent this line
            aes(x=t, y = vg_HEreg_noa.mean,                               
                linetype = "obs",
                group=interaction(P, cov), 
                color=interaction(P, cov))) +
  geom_line(data=data, linewidth=0.9, 
            aes(x=t, y = exp.vg, 
                linetype = "exp",
                group=interaction(P, cov),
                color=interaction(P, cov))) +
  #scale_y_log10(limits=c(0.9, 2.1)) +
  scale_linetype_manual("", 
                       breaks = c("obs",   "exp"),
                       values = c("solid",  "11"),
                       labels = c("Estimated\n(standard)",
                                  "Expected")) +
  scale_colour_manual("", 
                      values = c('#92c5de','#053061',
                                 '#f4a582','#67001f',
                                 '#79bc5c','#2a6112') ,
                      labels = c("P=0", "P=0.9", 
                                 "P=0", "P=0.9",
                                 "P=0", "P=0.9"
                                 )) +
 theme_classic() +
  xlab("t") +
  ylab(expression(hat(sigma)[u]^2))  + 
  theme(aspect.ratio = 1, 
        legend.position = "bottom", #hide legend for HI
        text = element_text(size = 12)
        ) + 
  guides(color =  guide_legend(order = 1, nrow=3, 
                              byrow = T, reverse = T,
                              override.aes = list(linewidth = 2)),  # thicken the line in legend
          linetype = "none" 
         )}
HI_color=fig4color(df_HI_CI)



# HE Vg wo and wganc
library(ggpubr)
plt_wo=ggarrange(plotlist = list(HIa_wo, CGFa_wo, 
                 HIb_wo, CGFb_wo, 
                 HIc_wo, CGFc_wo, 
                 HId_wo, CGFd_wo),
                 ncol = 2, nrow = 4, 
                 labels = c("A", "", "B", "","C", "","D", ""),
                 align = c("h"))

plt_wganc=ggarrange(plotlist = list(HIa_wganc, CGFa_wganc, 
                    HIb_wganc, CGFb_wganc, 
                    HIc_wganc, CGFc_wganc, 
                    HId_wganc, CGFd_wganc),
                    ncol = 2, nrow = 4, 
                    labels = c("E", "", "F", "","G", "","H", ""),
                    align = c("h"))

# add space between these two
plt=ggarrange(plotlist = list(plt_wo, '', plt_wganc), 
              ncol = 3, nrow = 1, 
              widths = c(8,1,8),
              align = c("v")) %>% # to move the legend closer 
  gridExtra::grid.arrange(get_legend(HI_color), 
                          heights = unit(c(9, 0.8), "in")
  ) 

ggsave("HE_vg_wowganc_CI.png", plot=plt,
       width = 16, height = 10, dpi = 300, units = "in", device='png')
