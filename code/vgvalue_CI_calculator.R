# calculate var(gvalue) for each trait
# output the value and the trait
# bootstrap get CI 
args <- commandArgs(trailingOnly = TRUE)
trait <- args[1]

pvar  = function(x){
  mean((x - mean(x))^2)
}

#trait="Height"
gvaluefile=paste0("ASW_",trait, ".profile")
gvalue=read.table(gvaluefile, header = T)
#vg=pvar(gvalue$SCORE)
library(dplyr)
# Resample 10000 times, and find the var of each
boot=tibble(num = 1:10000) %>% 
  group_by(num) %>% 
  mutate(vars = pvar(sample(gvalue$SCORE, replace = TRUE))) 
# hist(boot$vars)
# report median
boot_median=median(boot$vars)
# Bootstrap 95% percentile confidence interval
boot_CI=quantile(boot$vars, c(0.025,0.975))
# report the median and CI95%

vg_output=data.frame(
  vg=boot_median,
  CI95l=boot_CI[1],
  CI95r=boot_CI[2],
  Traits=trait)

# output the vgp for each trait
write.table(vg_output, file=paste0("../",trait, "_vgvalueCI.txt"), sep = '\t',
            quote = F, col.names = T, row.names = F)
