# Convert the local ancestry to Alicia's format and plot for ancestry percentage at each locus

# processing my compact local ancestry result
import pandas as pd
import numpy as np
# load the local ancestry file
df=pd.read_csv("loc_anc_AA_100Mb_random_400hap_chr1.vcf.txt", delimiter='\t') 

# delete child 0-199 and extract poulation 1:EUR only 
df=df.loc[df['child'] > 199].loc[df['populations'] == 1].sort_values('child')

# convert left and right to integer
df['left'] = df['left'].fillna(0.0).astype(int)
df['right'] = df['right'].fillna(0.0).astype(int)

# reindex the df
df=df.reset_index(drop=True)

# create a datafame with all 0
newdf = pd.DataFrame(np.zeros((100000000,400), dtype=np.uint8))
#change all EUR position to 1
for i in range(len(df)):
    l=df.loc[i].left
    r=df.loc[i].right
    c=df.loc[i].child-200
    newdf.loc[l:r,c]=1

# output data as a file
np.savetxt('adm_400hap_random_local_anc.txt', newdf, fmt="%d") 
