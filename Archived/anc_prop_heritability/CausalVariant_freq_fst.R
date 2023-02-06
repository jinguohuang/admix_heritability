# extract causal variants and 
# output frequency dependent effect size file
library(data.table)
library(dplyr)
library(tidyr)
library(rprojroot)
library(ggplot2)
set.seed(120)
# load file
args = commandArgs(trailingOnly=TRUE)
#frequency file path
freq_file=args[1] # frequency file
#freq_file="ADM_0.2prop.frq"
# load fst file 
fst_file=args[2] # frequency file
#fst_file="ref_0.2prop_updateid.fst"

Name <- sub("(.*?)\\.frq", '\\1', freq_file) # extract name 

# load variant frequency file
p = fread(freq_file)
colnames(p)=c("CHROM","ID","REF","ALT","ALT_FREQS","COUNT")
p=p[,c("CHROM","ID","REF","ALT_FREQS")]
#extract the ID to get chr+position+ref+alt
p[, c("chr", "position","ref","alt") := tstrsplit(ID, "_", fixed=TRUE)]
p = p[,c("CHROM","ID","position","REF", "ALT_FREQS")]
p$position = as.numeric(p$position)

# load fst file
df = fread(fst_file)
colnames(df)=c("CHR","ID","POS","NMISS","FST")
df=df[,c("ID","FST")]
#dist(fst_file$FST) #Error: vector memory exhausted (limit reached?)
# filter for 0.15<fst<0.25
df<-subset(df, FST > 0.15 & FST < 0.25)

# merge freq and fst file so we filter out for fst
p_df<-merge(p,df, by = "ID")

#write function to do this.
sample.variant=function(df1){
  #sample first variant
  position1 = as.numeric( sample_n( df1[ position < 1e5, 'position' ], 1 ))
  #select all other variants to be at least 100kb apart
  #minimum positions for each window 
  positions = position1 + seq(0,199)*1e5 #200 per chromosome
  #pick variants that are further than these
  positions.adj = lapply( positions, function(x){
    ix = min( df1[position > x, which =TRUE ] )
    return(df1[ix])
  })
  #return datatable
  positions.adj = bind_rows(positions.adj)
  return(positions.adj)
}
#carry this out grouped by chromosome
causal.variants <- p[, sample.variant(.SD), by=CHROM]
#let's remove NAs if there are any.
causal.variants = causal.variants%>%drop_na(ID)

# simulate effect sizes for these variants
# calculate the independent component of variance required
sigma2_l = 0.8 / sum( sapply( causal.variants$ALT_FREQS,function(x){
  beta= ( 2*x*(1-x)) ^ (1-0.4)
  return(beta)
}))

#  sample the effect sizes for each variant.
# sample maf-dependent effects using the model above
causal.variants$beta = sapply( causal.variants$ALT_FREQS , function(x){
  beta = rnorm( 1 , mean = 0, sd = sqrt(sigma2_l * (2*x*(1-x))^-0.4 ))
})
#let's calculate sigma2_g to confirm that the total genetic variance is indeed 0.8
sigma2_g = sum( mapply(function(b,p){ b^2* 2*p*(1-p) }, causal.variants$beta, causal.variants$ALT_FREQS))
#print(paste("sigma2_g : ",round(sigma2_g,3)))
# "sigma2_g :  0.800"

# output file for plink to process
score=causal.variants[,c("ID","REF", "beta")]
# output the file

write.table(score, file = paste0("CausalVariantScore_", Name,".txt"), quote=F, row.names = F, col.names = F)
# output the ID list
ID=causal.variants[,c("ID")]
write.table(ID, file = paste0("CausalVariantID_", Name,".txt"), quote=F, row.names = F, col.names = F)
