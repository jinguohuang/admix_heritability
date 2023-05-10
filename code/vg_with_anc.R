

library(data.table)
library(ggplot2)

F = rprojroot::is_rstudio_project$make_fix_file()

prs = fread(F("data/admix_HI_theta0.5_gen20_P0.9_pos_sall_tall.prs"))
colnames(prs) = c("rep","time","iid","gvalue")

ganc = fread(F("data/admix_HI_theta0.5_gen20_P0.9_pos_sall_tall.ganc"))
colnames(ganc) = c("rep","time","iid","glanc")

prs = merge(prs, ganc, by = c("rep","time","iid"))

arch.pos = fread("data/admix_FreqBeta_pos.txt")
arch.neg = fread("data/admix_FreqBeta_neg.txt")

a = with(arch.pos, sum( 2*0.5*beta_pos^2 *f1*(1 - f1) + 2*0.5*beta_pos^2 *f2*(1-f2) ))
b = with(arch.pos, sum(2*0.25*beta_pos^2 *(f1 - f2)^2))
c = with(arch.pos, sum(beta_pos^2 *(f1 - f2)^2))
d1 = as.matrix(with(arch.pos, beta_pos*(f1 - f2)))
d1 = d1%*%t(d1)
d1 = sum(d1[lower.tri(d1)])

d2 = as.matrix(with(arch.neg, beta_neg*(f1 - f2)))
d2 = d2%*%t(d2)
d2 = sum(d2[lower.tri(d2)])

prs$gvalue.res = 0
# atrait.r2 = matrix(NA, 10, 20)
atrait.vgres = matrix(NA, 10, 20)
exp.vgall = matrix(NA, 10, 20)
exp.vgres = matrix(NA, 10, 20)
exp1 = matrix(NA, 10, 20)
exp2 = matrix(NA, 10, 20)
exp3 = matrix(NA, 10, 20)
exp4.pos = matrix(NA, 10, 20)
exp4.neg = matrix(NA, 10, 20)
for(i in 1:10){
  for(j in 1:20){
    l1 = lm(data = prs[rep==i & time == j], 
            gvalue ~ glanc)
    vtheta = var(prs[rep==i & time == j, glanc])
    vres = summary(l1)$r.squared * var(prs[rep==i & time == j, gvalue])
    atrait.vgres[i,j] = var(prs[rep==i & time == j, gvalue]) - vres 
    exp1[i,j] = a
    exp2[i,j] = a + 2*0.5*(1-0.5)*b
    exp3[i,j] = a + 2*0.5*(1-0.5)*b + 2*vtheta*c
    exp4.pos[i,j] = a + 2*vtheta*c +  8*vtheta*d1
    exp4.neg[i,j] = a + 2*vtheta*c +  8*vtheta*d2
    
  }
}

prs.summary = prs[, .(vg1 = var(gvalue)), by = c("rep","time")]
prs.summary = prs.summary[, .(vg1 = mean(vg1)), by = "time"]
#prs.summary$obs.c = apply(atrait.c, 2, mean)

# prs.summary$exp.c = apply(exp.c, 2, mean)
# prs.summary$exp.b = apply(exp.b, 2, mean)
# prs.summary$obs.c = apply(atrait.c, 2, mean)
# prs.summary$obs.b = apply(atrait.b, 2, mean)
prs.summary$res.v = apply(atrait.vgres, 2, mean)
# prs.summary$r2 = apply(atrait.r2, 2, mean)
prs.summary$exp1 = apply(exp1, 2, mean)
prs.summary$exp2 = apply(exp2, 2, mean)
prs.summary$exp3 = apply(exp3, 2, mean)
prs.summary$exp4.pos = apply(exp4.pos, 2, mean)
prs.summary$exp4.neg = apply(exp4.neg, 2, mean)


prs.summary = melt(prs.summary, id.vars = c("time","vg1","res.v"),
                    variable.name = "exp",
                    value.name = "value")

plt1 = ggplot(prs.summary)+
  geom_line(aes(time, value, color = exp),
            linetype = "dashed")+
  geom_line(aes(time, res.v), color = "black", 
            size = 0.4, alpha = 0.5)+
  theme_classic()+
  scale_color_manual(values = c("#fdbe85","#fd8d3c","#e6550d","#a63603","#045a8d"))+
  theme(legend.position = "none")+
  labs(x = "Time",
       y = "Vg")

plt1


ggsave(F("plots/plt_vg_anc_corrected.pdf"),
       plt1,
       height = 3, width = 3)





