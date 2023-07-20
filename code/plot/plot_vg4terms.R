# Plot behavior of vg 4 terms

# load data
filename1="../data/admix_CGF_vg_vgamma.txt"
df_CGF=read.table(filename1, header=T)
filename2="../data/admix_HI_vg_vgamma.txt"
df_HI=read.table(filename2, header=T)


# A function factory for getting integer y-axis values.
integer_breaks <- function(n = 5, ...) {
  fxn <- function(x) {
    breaks <- floor(pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
  return(fxn)
}


# plot just with theta=0.5, filter for theta=0.5
CGF0.5=df_CGF[which(df_CGF$theta==0.5),]
HI0.5=df_HI[which(df_HI$theta==0.5),]

#plot function
# middle plot, without title without x 
myplot_mid = function(data, yobs, yexp, ylab, #title, 
                   ylim1=0.9, ylim2=2.1) {
library(ggplot2)
p = ggplot() +
  geom_line(data=data, alpha = 0.65, linewidth=0.5, 
            aes(x=t, y = yobs, linetype = "obs",
                group=interaction(P, cov), 
                color=interaction(P, cov))) +
  geom_line(data=data, linewidth=0.5, 
            aes(x=t, y = yexp, linetype = "exp",
                group=interaction(P, cov),
                color=interaction(P, cov))) + 
  scale_linetype_manual("", 
                       breaks = c("obs", "exp"),
                       values = c("solid", "dotted"),
                       labels = c("Observed", "Expected")) +
   scale_colour_manual("", 
                      values = c('#92c5de','#4393c3','#2166ac','#053061',
                                 '#f4a582','#d6604d','#b2182b','#67001f'),
                      labels = c("P=0", "P=0.3", "P=0.6", "P=0.9",
                                 "P=0", "P=0.3", "P=0.6", "P=0.9"
                     )) +
  ylab(ylab)  +  
  theme_bw() +
  scale_x_continuous(breaks = integer_breaks()) +
  ylim(ylim1, ylim2) +
  theme(#aspect.ratio = 1, 
        plot.title = element_text(hjust = 0.5), # to center the title
        legend.position = "bottom", 
        text = element_text(size = 10),
        axis.title.x=element_blank(), axis.text.x=element_blank()) 

# New facet label names
theta.labs <- c("\u03B8 = 0.1", "\u03B8 = 0.2", "\u03B8 = 0.5")
names(theta.labs) <- c("0.1", "0.2", "0.5")

gen.labs <- c("10 generations", "20 generations","50 generations","100 generations")
names(gen.labs) <- c("10", "20", "50", "100")

vgCGF = p + facet_grid(rows = vars(theta), cols = vars(gen), #row as theta col as gen
             scale="free_x",  #adjust the xlim for each
  ) + 
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=7),
        axis.ticks.x=element_blank(), axis.text.x=element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)
        ) + 
   guides(color = guide_legend(order = 1, ncol=4, 
                              byrow = T, reverse = T,
                              override.aes = list(linewidth = 2)),  # thicken the line in legend
          linetype = guide_legend(order = 2))}



#plot top: with title and the generation
myplot_top = function(data, yobs, yexp, ylab, #title, 
                   ylim1=0.9, ylim2=2.1) {
library(ggplot2)
p = ggplot() +
  geom_line(data=data, alpha = 0.65, linewidth=0.5,
            aes(x=t, y = yobs, linetype = "obs",
                group=interaction(P, cov),
                color=interaction(P, cov))) +
  geom_line(data=data, linewidth=0.5, 
            aes(x=t, y = yexp, linetype = "exp", 
                group=interaction(P, cov),
                color=interaction(P, cov))) + 
  scale_linetype_manual("", 
                       breaks = c("obs", "exp"),
                       values = c("solid", "dotted"),
                       labels = c("Observed", "Expected")) +
   scale_colour_manual("", 
                      values = c('#92c5de','#4393c3','#2166ac','#053061',
                                 '#f4a582','#d6604d','#b2182b','#67001f'),
                      labels = c("P=0", "P=0.3", "P=0.6", "P=0.9",
                                 "P=0", "P=0.3", "P=0.6", "P=0.9"
                     )) +
  ylab(ylab)  +  
  theme_bw() +
  #ggtitle(title) +
  scale_x_continuous(breaks = integer_breaks()) +
  ylim(ylim1, ylim2) +
  theme(#aspect.ratio = 1, 
        plot.title = element_text(hjust = 0.5), # to center the title
        legend.position = "bottom", 
        text = element_text(size = 10),
        axis.title.x=element_blank(), axis.text.x=element_blank()) 

# New facet label names
theta.labs <- c("\u03B8 = 0.1", "\u03B8 = 0.2", "\u03B8 = 0.5")
names(theta.labs) <- c("0.1", "0.2", "0.5")

gen.labs <- c("10 generations", "20 generations","50 generations","100 generations")
names(gen.labs) <- c("10", "20", "50", "100")

vgCGF = p + facet_grid(rows = vars(theta), cols = vars(gen), #row as theta col as gen
             scale="free_x",  #adjust the xlim for each
  labeller = labeller(theta = theta.labs, gen = gen.labs)
  ) + 
  theme(
    strip.background = element_rect(color="black", fill="white", size=0.8, linetype="solid"),
    strip.text.x = element_text(size = 6, color = "black", face = "bold"),
    strip.text.y = element_blank(),
    panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=7),
        axis.ticks.x=element_blank(), axis.text.x=element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)
        ) + 
   guides(color = guide_legend(order = 1, ncol=4, 
                              byrow = T, reverse = T,
                              override.aes = list(linewidth = 2)),  # thicken the line in legend
          linetype = guide_legend(order = 2))}


#plot bottom: with x axis and ticks 
myplot_bot = function(data, yobs, yexp, ylab, #title, 
                   ylim1=0.9, ylim2=2.1) {
library(ggplot2)
p = ggplot() +
  geom_line(data=data, alpha = 0.65, linewidth=0.5,
            aes(x=t, y = yobs, linetype = "obs",
                group=interaction(P, cov),
                color=interaction(P, cov))) +
  geom_line(data=data, linewidth=0.5, 
            aes(x=t, y = yexp, linetype = "exp", 
                group=interaction(P, cov),
                color=interaction(P, cov))) + 
  scale_linetype_manual("", 
                       breaks = c("obs", "exp"),
                       values = c("solid", "dotted"),
                       labels = c("Observed", "Expected")) +
   scale_colour_manual("", 
                      values = c('#92c5de','#4393c3','#2166ac','#053061',
                                 '#f4a582','#d6604d','#b2182b','#67001f'),
                      labels = c("P=0", "P=0.3", "P=0.6", "P=0.9",
                                 "P=0", "P=0.3", "P=0.6", "P=0.9"
                     )) +
  ylab(ylab)  +
  theme_bw() +
  scale_x_continuous(breaks = integer_breaks()) +
  #scale_y_log10(limits=c(0.9, 2.1)) +
  ylim(ylim1, ylim2) +
  theme(#aspect.ratio = 1, 
        plot.title = element_text(hjust = 0.5), # to center the title
        legend.position = "bottom", 
        text = element_text(size = 10))#,
        #axis.title.x=element_blank(), axis.text.x=element_blank()) 

# New facet label names
theta.labs <- c("\u03B8 = 0.1", "\u03B8 = 0.2", "\u03B8 = 0.5")
names(theta.labs) <- c("0.1", "0.2", "0.5")

gen.labs <- c("10 generations", "20 generations","50 generations","100 generations")
names(gen.labs) <- c("10", "20", "50", "100")

vgCGF = p + facet_grid(rows = vars(theta), cols = vars(gen), #row as theta col as gen
             scale="free_x",  #adjust the xlim for each
  labeller = labeller(theta = theta.labs, gen = gen.labs)
  ) + 
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),#
    strip.text.y = element_blank(),
    panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=7),
        panel.border = element_rect(colour = "black", fill = NA)
        ) + 
   guides(color = guide_legend(order = 1, ncol=4, 
                              byrow = T, reverse = T,
                              override.aes = list(linewidth = 2)),  # thicken the line in legend
          linetype = guide_legend(order = 2))}

vg1CGF=myplot_top(data=CGF0.5, 
         yobs=CGF0.5$vg.term1, 
         yexp=CGF0.5$va.term1,
         ylab="(1.1)", 
         #title="Continuous Gene Flow",
         ylim1=0.92, ylim2 = 0.95 )


vg1HI=myplot_top(data=HI0.5, 
         yobs=HI0.5$vg.term1, 
         yexp=HI0.5$va.term1,
         ylab="(1.1)", 
         #title="Hybrid Isolation",
         ylim1=0.92, ylim2 = 0.95 )

vg2CGF=myplot_mid(data=CGF0.5, 
         yobs=CGF0.5$vg.term2, 
         yexp=CGF0.5$va.term2,
         ylab="(1.2)", 
         #title="Continuous Gene Flow",
         ylim1=0, ylim2 = 0.1 )


vg2HI=myplot_mid(data=HI0.5, 
         yobs=HI0.5$vg.term2, 
         yexp=HI0.5$va.term2,
         ylab="(1.2)", 
         #title="Hybrid Isolation",
         ylim1=0, ylim2 = 0.1 )

vg3CGF=myplot_mid(data=CGF0.5, 
         yobs=CGF0.5$vg.term3, 
         yexp=CGF0.5$va.term3,
         ylab="(1.3)", 
         #title="Continuous Gene Flow",
         ylim1=0, ylim2 = 0.065)

vg3HI=myplot_mid(data=HI0.5, 
         yobs=HI0.5$vg.term3, 
         yexp=HI0.5$va.term3,
         ylab="(1.3)", 
         #title="Hybrid Isolation",
         ylim1=0, ylim2 = 0.065 )

vg4CGF=myplot_bot(data=CGF0.5, 
         yobs=CGF0.5$vg.term4, 
         yexp=CGF0.5$va.term4,
         ylab="(1.4)", 
         #title="Continuous Gene Flow",
         ylim1=-0.25, ylim2 = 1.1)


vg4HI=myplot_bot(data=HI0.5, 
         yobs=HI0.5$vg.term4, 
         yexp=HI0.5$va.term4,
         ylab="(1.4)", 
         #title="Hybrid Isolation",
         ylim1=-0.25, ylim2 = 1.1 )

library(ggpubr)
require(grid)
plt=
  ggarrange(vg1HI, 
            vg2HI, 
            vg3HI, 
            vg4HI + rremove("xlab"), 
            ncol = 1, nrow = 4, 
            align = "v",
            legend = "none") 
 
plt1=annotate_figure(plt, top = textGrob("Hybrid Isolation", vjust = 1, gp = gpar(cex = 1.1)),
                    bottom = textGrob("t", gp = gpar(cex = 1)))
  
plt=
  ggarrange(vg1CGF,
            vg2CGF,
            vg3CGF,
            vg4CGF + rremove("xlab"),
            ncol = 1, nrow = 4, 
            align = "v",
            legend = "none")  
plt2=annotate_figure(plt, top = textGrob("Continuous Gene Flow", vjust = 1, gp = gpar(cex = 1.1)),
                    bottom = textGrob("t", gp = gpar(cex = 1)))

plt=
  ggarrange(plt1, plt2,
            ncol = 2, nrow = 1, 
           labels = c("A", "B"),
            align = "v",
            legend = "none") %>% # to move the legend closer 
  gridExtra::grid.arrange(get_legend(vg1HI), 
                          heights = unit(c(5, 0.6), "in")
                          )

ggsave("../figs/Fig_vgterms.png", plot=plt,
       width = 8, height = 6,
       dpi = 300, units = "in", device='png')

