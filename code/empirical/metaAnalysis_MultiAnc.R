# meta analysis for trans ancestry summary stat to get beta weighting by inverse variance
args <- commandArgs(trailingOnly = TRUE)
trait <- args[1]
df=read.table(paste0(trait, ".txt"), header=T )

# filter : D I
df=df[!(df$SNP==":" | df$A1=="D" | df$A1=="I") ,]
df=df[order(df$SNP),] # reorder
# if only one, keep the beta
df_nodup=df[!(duplicated(df$SNP) | duplicated(df$SNP, fromLast=TRUE)),]
df_dup=df[(duplicated(df$SNP) | duplicated(df$SNP, fromLast=TRUE)),]

# if different effect allele, flip the allele and the beta
# get unique SNP list
SNPLIST=df_dup[duplicated(df_dup$SNP),]$SNP
df_dup_flip=df_dup[0,] # empty table to put new data
# empty table to put new data
df_dup_flip_weight = data.frame(matrix(nrow = length(SNPLIST), 
                                       ncol = ncol(df_dup)))
colnames(df_dup_flip_weight)=colnames(df_dup)

if (length(SNPLIST)!=0){
  for (n in 1:length(SNPLIST)){
    a=df_dup[(df_dup$SNP==SNPLIST[n]) ,]
    for (i in 2:nrow(a)){
      if (a$A1[1]!=a$A1[i]) { #Loop to flip not matched A1
        a$A1[i]=a$A1[1] #flip to the first one
        a$BETA[i]=a$BETA[i]*(-1) #update the beta sign
        df_dup_flip=rbind(df_dup_flip, a)
      }
      #update a 
      else{df_dup_flip=rbind(df_dup_flip, a)}
    }
  } 
  
  library(data.table)
  # after SNP flip, if more than one, calculated SE with pvalue
  df_dup_flip= as.data.table(df_dup_flip) 
  df_dup_flip[, ':='(se=BETA/qnorm(P, lower.tail = FALSE))]
  df_dup_flip[, ':='(inv_var= 1/(se^2))] # inverse variance
  #  calculate inverse-variance weighting beta

  # function to calculate Inverse-variance weighting
  IVW=function(y, inv_var){
    yhat=sum(y*inv_var)/sum(inv_var)
    return(yhat)
  }
  
  for (n in 1:length(SNPLIST)){
    df_dup_flip_weight$SNP[n]=SNPLIST[n]
    # for each duplicated SNP, get each of them
    a=df_dup_flip[(df_dup_flip$SNP==SNPLIST[n]) ,]
    # calculated new beta
    # count how many Inf
    InfCount=nrow(a[a$inv_var==Inf,])
    if (InfCount==0) {# if InfCount=0, calculate the weighted beta
      # calculate the weighted beta
      df_dup_flip_weight$A1[n]=a$A1[1]
      df_dup_flip_weight$A1_FREQ[n]=a$A1_FREQ[1]
      df_dup_flip_weight$P[n]=a$P[1]
      df_dup_flip_weight$BETA[n]=IVW(y=a$BETA, inv_var=a$inv_var)  
    }
    
    else {# if InfCount>= 1, just use that average of that beta
      a_inf=a[a$inv_var==Inf,]
      df_dup_flip_weight$A1[n]=a_inf$A1[1]
      df_dup_flip_weight$A1_FREQ[n]=a_inf$A1_FREQ[1]
      df_dup_flip_weight$P[n]=a_inf$P[1]
      df_dup_flip_weight$BETA[n]=mean(a_inf$BETA)
    }
  }
}
# put the nodup and updated dup together
df_new=rbind(df_nodup, df_dup_flip_weight)

# output the updated summary stat
write.table(df_new, file = paste0(trait, "_meta.txt"), 
            quote=FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
