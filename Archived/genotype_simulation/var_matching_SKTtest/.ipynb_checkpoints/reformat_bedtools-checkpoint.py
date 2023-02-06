# This script is for reformating the output from 
# bedtools intersect output: chr  left   right   loc_anc   individual
# to a local ancestry file with row as variants column as individuals
# and another file for storing the variants position
import sys
import pandas as pd
import numpy as np
# load the file without duplicate
bedfile=sys.argv[1]
# get filename
filename=bedfile.rsplit('.', 3)[0]
# read file
df=pd.read_csv(bedfile, sep="\t", header=None)
df.columns =['CHR', 'START', "END", "ANC", "IND"]
# get end, anc and ind column only
df=df[['END', 'ANC', 'IND']]
# unmelt!
df_unmelt=df.pivot(index='END', columns= 'IND')
# output the unmelt dataframe in a file, without index
print("writing the reformated local ancestry ...")
df_unmelt.to_csv(f'{filename}.reformat', sep='\t', encoding='utf-8', header= False, index=False)
# output the index only
print("wrting the variants to file ...")
a=df_unmelt.index.tolist()
np.savetxt(f'{filename}.map', a, delimiter ="\n", fmt ='% s')

