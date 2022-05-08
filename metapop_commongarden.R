
library(data.table)
library(ggplot2)

#checking if Vgb is equal to local ancestry genetic variance in a common garden
nloci = 100
ninds = 2000
theta = 0.5
f1 = 0.1
f2 = 0.7
fbar = (0.1+0.7)/2

#bg = effect size of locus - constant across loci to keep things simple
bg = 1
by = bg*(f2-f1)

#since this is a common garden, global ancestry is basically either 0 or 1
global.anc = c(rep(0, theta*ninds), rep(1, (1-theta)*ninds))
#local ancestry at each locus is the same as global ancestry
local.anc = t(sapply(global.anc, function(x){rbinom(nloci, 2, x)}))

#function to generate genetic value 
gen.gvalue = function(ninds, nloci){
  genos = matrix(NA, nrow = ninds, ncol = nloci)
  for(i in 1:ninds){
    if(global.anc[i]==0){
      genos[i,] = rbinom(nloci, 2, f1)
    }
    if(global.anc[i]==1){
      genos[i,] = rbinom(nloci, 2, f2)
    }
  }
  
  gvalue1 = apply(bg*genos, 1, sum)
  return(gvalue1)
}

#covariance in genetic values in the meta-population
cov.genos = cov(genos)
diag(cov.genos) = 0 # set diagonal to zero (this is the variance)
exp.vt = 2*nloci*bg^2 *fbar*(1-fbar) + sum(cov.genos) # expected total variance
exp.vw = nloci*((bg^2)*f1*(1-f1) + (bg^2)*f2*(1-f2)) # expected variance within
exp.vb = exp.vt - exp.vw #expected variance between (formula1)
exp.vb2 = 0.5*nloci*bg^2*(f1-f2)^2 + sum(cov.genos) #expected variance between (formula 2)


# simulate multiple iterations of genetic values and estimate the variance components in each
nreps = 100
obs.vt = rep(NA, nreps)
obs.vw = rep(NA, nreps)
obs.vb = rep(NA, nreps)
npop1 = (theta*ninds)
for(i in 1:nreps){
  gvalue1 = gen.gvalue(ninds, nloci)
  obs.vt[i] = sum((gvalue1 - mean(gvalue1))^2) /ninds
  v1 = sum((gvalue1[c(1:npop1)] - mean(gvalue1[c(1:npop1)]))^2) /npop1
  v2 = sum((gvalue1[c((npop1+1):ninds)] - mean(gvalue1[c((npop1+1):ninds)]))^2) /(ninds - npop1)
  obs.vw[i] = (v1+v2)/2
  obs.vb[i] = obs.vt[i] - obs.vw[i]
}


#plot the obvserved vs expected
obs.dat = data.table(vw = obs.vw, vt = obs.vt, vb = obs.vb)
obs.dat = melt(obs.dat)
exp.dat = data.table(variable = c("vw","vt","vb"),
                     value = c(exp.vw, exp.vt, exp.vb))

ggplot()+
  geom_point(data = obs.dat, aes("a", value), position = "jitter")+
  geom_point(data = exp.dat, aes("a", value), color = "red")+
  facet_wrap(~variable, scales = "free")

# now, calculate genetic values using local ancestry and local ancestry effect size
gvalue2 = apply(by*local.anc, 1,sum)
#calculate genetic variance due to local ancestry
obs.vl =sum((gvalue2 - mean(gvalue2))^2) /ninds


