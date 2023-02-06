
library(MASS)


#generate one phenotype

#the frequency of 1000 causal variants
freqs = runif(1000,0.1,0.9)

#sample genotypes for 100 individual for each variant
x1 = sapply(freqs, function(f){rbinom(100, 2, f)})
#normalize x
x1 = apply(x1, 2, function(f){(f - mean(f))/sd(f)})

#sample the effect sizes of each variant
#the desired heritability is 0.8
#the sd is therefore sqrt(0.8)
beta1 = rnorm(1000, mean = 0, sd = sqrt(0.8))

#sample random error
#sd = sqrt(remaining variance = 0.2)
#we want the phenotype variance to be 1
#var(y) = var(x*beta) + var(e) => 1 = 0.8 + 0.2
e1 = rnorm(100, mean = 0, sd = sqrt(0.2))

#phenotype of all individuals
y1 = x1 %*% beta1 + e1

#normalize y 
y1 = (y1 - mean(y1))/sd(y1)

h2_1 = 0.8 # known heritability of the first trait
h2_2 = 0.5 # known heritability of the second trait
cov_12 = 0.5 # known genetic covariance

beta_12 = mvrnorm(n = 1000, 
                  mu = c(0,0),
                  Sigma = matrix(c(h2_1, cov_12, cov_12, h2_2), 2,2))


