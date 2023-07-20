# Plot behavior of vg

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

plot_vg = function(data, yobs, yexp, ylab, title) {
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
  xlab("t") +
  ylab(ylab)  +  
  ggtitle(title) +
  theme_bw() +
  scale_x_continuous(breaks = integer_breaks()) +
  scale_y_log10(limits=c(0.9, 2.1)) +
  theme(aspect.ratio = 1, 
        plot.title = element_text(hjust = 0.5), # to center the title
        legend.position = "bottom", 
        text = element_text(size = 10)) 

# New facet label names
theta.labs <- c("\u03B8 = 0.1", "\u03B8 = 0.2", "\u03B8 = 0.5")
names(theta.labs) <- c("0.1", "0.2", "0.5")

gen.labs <- c("10 generations", "20 generations","50 generations","100 generations")
names(gen.labs) <- c("10", "20", "50", "100")

vgCGF = p + facet_grid(rows = vars(theta), cols = vars(gen), #row as theta col as gen
             scale="free_x",  #adjust the xlim for each
  labeller = labeller(theta = theta.labs, gen = gen.labs)) + 
  theme(
    strip.background = element_rect(
      color="black", fill="white", size=0.8, linetype="solid"),
    strip.text.x = element_text(
        size = 6, color = "black", face = "bold"),
    strip.text.y = element_text(
        size = 7, color = "black", face = "bold"),
    panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=7),
        #strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)
        ) + 
   guides(color = guide_legend(order = 1, ncol=4, 
                              byrow = T, reverse = T,
                              override.aes = list(linewidth = 2)
                             ),  # thicken the line in legend
          linetype = guide_legend(order = 2))
          }

vgHI=plot_vg(data=df_HI, 
         yobs=df_HI$var.prs.geno, 
         yexp=df_HI$exp.vg,
         ylab=expression(V[g]), 
         title="Hybrid Isolation")

vgCGF=plot_vg(data=df_CGF, 
         yobs=df_CGF$var.prs.geno, 
         yexp=df_CGF$exp.vg,
         ylab=expression(V[g]), 
         title="Continuous Gene Flow")

library(ggpubr)
plt=ggarrange(vgHI, vgCGF, ncol = 2, nrow = 1, 
          labels = c("A", "B"),
          align = c("hv"),
          legend = "none") %>% # to move the legend closer 
  gridExtra::grid.arrange(get_legend(vgHI), heights = unit(c(100, 5), "mm"))

ggsave("../figs/Fig_Vg.png", plot=plt,
       width = 8, height = 5, dpi = 300, units = "in", device='png')
