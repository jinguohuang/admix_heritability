#!/usr/bin/env Rscript

# getting causal variants 1/100kb

library(data.table)
library(dplyr)
library(tidyr)
library(rprojroot)
#library(ggplot2)
set.seed(123)
# load file
args = commandArgs(trailingOnly=TRUE)
#frequency file path
#freq_file="AA_10Mb_10deme_0.05_chr1_10_AT_maf0.05.frq"
freq_file=args[1] # frequency file
# load variant frequency file
p = fread(freq_file)
colnames(p)=c("CHROM","ID","REF","ALT","ALT_FREQS","COUNT")
p=p[,c("CHROM","ID","ALT_FREQS")]
#extract the ID to get chr+position+ref+alt
p[, c("chr", "position","ref","alt") := tstrsplit(ID, "_", fixed=TRUE)]
p = p[,c("CHROM","ID","position","ALT_FREQS")]
p$position = as.numeric(p$position)

#Now sample 1,000 variants at random such that they are uniformly spaced 
# across the genome with ~ 100Kb between them. We do this for each chromosome separately.
#write function to do this.
sample.variant=function(df1){
  #sample first variant
  position1 = as.numeric( sample_n( df1[ position < 1e5, 'position' ], 1 ))
  #select all other variants to be at least 100kb apart
  #minimum positions for each window
  positions = position1 + seq(0,99)*1e5
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
# Output causal variants ID list to extract with plink and get dosage (additive) file.
id_list<-causal.variants$ID
write.table(id_list, file = "CausalVariantExtract.txt", quote=F, row.names = F, col.names = F)

