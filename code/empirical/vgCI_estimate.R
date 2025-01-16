#estimate 4 terms for each trait by parametric bootstrap of 100 betas from a normal distribution, 
#calculate 4 terms, bootstrap these results for a mean and CI's for each term

library(readxl)
library(tidyverse)
library(data.table)

#load in trait data
sd_names <- excel_sheets("TableS1_Traits_AlleleFreq_EffectSize_sd.xlsx")
print(sd_names)

#read all sheets into a list
sd_sheets <- lapply(sd_names, function(x) read_excel("TableS1_Traits_AlleleFreq_EffectSize_sd.xlsx", sheet = x))

#name each element of the list after the sheet name
#replace "-*" with ""
fixed_sd_names <- gsub("-.*", "", sd_names)
# print(fixed_sheet_names)

#assign names
names(sd_sheets) <- fixed_sd_names
names(sd_sheets)
length(names(sd_sheets))

#subset just trait sheets
sd_sheets <- sd_sheets[4:29]
names(sd_sheets)
length(sd_sheets)

#calculate 4 terms of genetic variance Vg

# term 1 of Vg
vg_term1 <- function(exp.ganc, f1, f2, beta) return (sum(beta^2*(2*exp.ganc*f1*(1-f1)+2*(1-exp.ganc)*f2*(1-f2))))

# term 2 of Vg
vg_term2 <- function(exp.ganc, f1, f2, beta) return (sum(beta^2*(2*exp.ganc*(1-exp.ganc)*((f1-f2)^2))))

# term 3 of Vg
vg_term3 <- function(var.ganc, f1, f2, beta) return (sum(beta^2*(2*var.ganc*((f1-f2)^2))))

# term 4 of Vg
vg_term4 <- function(var.ganc, f1, f2, beta, nloci){
  nloci=nloci
  term4=matrix(, nrow = nloci, ncol = nloci)
  for (i in 1:nloci){
    for (j in 1:nloci){
      term4[i,j]=beta[i]*beta[j]*(f1[i]-f2[i])*(f1[j]-f2[j])
    }
  }
  # make diagnal 0
  diag(term4)=0
  # sum it up
  sum.term4 = sum(term4)*4*var.ganc
  return (sum.term4)
  }

#function to bootstrap 4 terms using 100 pseudo betas for each SNP for each trait
boot_beta <- function(betafreq, n = 100) {
  nloci = nrow(betafreq)
  
  #empty mat for pseudo betas
  beta_mat <- matrix(NA, nrow = nloci, ncol = 100)
  
  #for each beta get 100 pseudo betas
  for (i in 1:nloci) {
    b <- betafreq$BETA[i]
    b_sd <- betafreq$SD[i]
    
    #randomly sample 100 beta values from a normal distribution
    pseudo_betas <- rnorm(n, mean = b, sd = b_sd)
    
    #combine all pseudo betas into df
    beta_mat[i, ] <- pseudo_betas
  }
  #print(beta_mat)
  
  boot_mat <- matrix(NA, nrow = 100, ncol = 4) 
  
  #loop through each col in beta_mat (i.e., each set of 100 pseudo betas for each variant)
  for (j in 1:100) { 
    pseudo_beta <- beta_mat[, j]
    
    #get 4 terms for pseudo-betas
    vg1 <- vg_term1(exp.ganc = 0.767,
                    f1 = betafreq$AF_EUR,
                    f2 = betafreq$AF_AFR,
                    beta = pseudo_beta)
    
    vg2 <- vg_term2(exp.ganc = 0.767,
                    f1 = betafreq$AF_EUR,
                    f2 = betafreq$AF_AFR,
                    beta = pseudo_beta)
    
    vg3 <- vg_term3(var.ganc = 0.018,
                    f1 = betafreq$AF_EUR,
                    f2 = betafreq$AF_AFR,
                    beta = pseudo_beta)
    
    vg4 <- vg_term4(var.ganc = 0.018,
                    f1 = betafreq$AF_EUR,
                    f2 = betafreq$AF_AFR,
                    beta = pseudo_beta,
                    nloci = nloci)
    
    boot_mat[j, 1] <- vg1
    boot_mat[j, 2] <- vg2
    boot_mat[j, 3] <- vg3
    boot_mat[j, 4] <- vg4
  }
  
  #print(boot_mat)
  
  ####output df
  mat <- matrix(NA, nrow = 1, ncol = 12)
  
  boot_df <- as.data.table(boot_mat)
  colnames(boot_df) <- c("vg1", "vg2", "vg3", "vg4")
  #print(boot_df)
  
  #get bootstrapped results
  mat[1,1] <- mean(boot_df$vg1)
  mat[1,2] <- mean(boot_df$vg2)
  mat[1,3] <- mean(boot_df$vg3)
  mat[1,4] <- mean(boot_df$vg4)
  
  #get CI's
  mat[1,5] <- quantile(boot_df$vg1, 0.025)  #lower CI
  mat[1,6] <- quantile(boot_df$vg1, 0.975)  #upper CI
  
  mat[1,7] <- quantile(boot_df$vg2, 0.025)  #lower CI
  mat[1,8] <- quantile(boot_df$vg2, 0.975)  #upper CI
  
  mat[1,9] <- quantile(boot_df$vg3, 0.025)  #lower CI
  mat[1,10] <- quantile(boot_df$vg3, 0.975)  #upper CI
  
  mat[1,11] <- quantile(boot_df$vg4, 0.025)  #lower CI
  mat[1,12] <- quantile(boot_df$vg4, 0.975)  #upper CI
  
  #final output
  df <- as.data.table(mat)
  colnames(df) <- c("vg1_avg", "vg2_avg", "vg3_avg", "vg4_avg",
                    "vg1_CI95l", "vg1_CI95r", "vg2_CI95l", "vg2_CI95r",
                    "vg3_CI95l", "vg3_CI95r", "vg4_CI95l", "vg4_CI95r")
  
  return(df)
}

#run boot_beta function for traits
#this will take a while to run as there are over 15,000 variants
vgpCI_list <- list()  

for (i in seq_along(sd_sheets)) {
  betafreq <- sd_sheets[[i]]
  
  #get bootstrap results
  boot_CIs <- boot_beta(betafreq)
  
  #output
  ntrait <- data.table(nloci = nrow(betafreq), trait = names(sd_sheets)[i])
  vgp <- cbind(boot_CIs, ntrait)
  
  vgpCI_list[[i]] <- vgp
}

head(vgpCI_list)

#combine all terms
boot_terms <- do.call(rbind, vgpCI_list)
head(boot_terms)

#save output
write.table(boot_terms, file = "vgCI.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

