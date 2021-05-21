#!/usr/bin/env Rscript

library(data.table)
# load file
args = commandArgs(trailingOnly=TRUE)
mapname<-args[1]  # load ancestry file

dat<-fread(mapname)
colnames(dat)<-c("pos","g_dist")
#sample 10 positions per cM, get the index
s_list=c()
#sample 50cM
for(i in 1:50){
#for(i in 1:5){
  s=100 # start sample position 100cM
  a <- which(dat$g_dist>(s+i-1) & dat$g_dist<(s+i))
  n=length(a) # sample from the range
  ix = sample(1:n, 10, replace = FALSE) # sample 10 positions
  s_list[(1+(i-1)*10):(i*10)]=a[ix] #put in the list
}
# get the map
map<-dat[s_list,]
write.table(s_list, file = "index_sample500.txt", quote = FALSE, sep = "\n",
            row.names = FALSE, col.names = FALSE)
write.table(map, file = "map_sample500.txt", quote = FALSE, sep = " ",
            row.names = FALSE, col.names = FALSE)
# extract to get the local ancestry
lanc_file<-args[2]
lanc<-fread(lanc_file)
lanc_sample<-lanc[s_list,]
write.table(lanc_sample, file = "lanc_sample500.txt", quote = FALSE, sep = " ",
            row.names = FALSE, col.names = FALSE)


# join map and lanc_sample
# calculate genetic distance and correlation 


