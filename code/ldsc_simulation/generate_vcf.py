#script generates simulated genotype data
#the primary motviation to write this is to generate genotype data (a time-consuming step), which can then be used for downstream analyses
import argparse

parser=argparse.ArgumentParser()
req_grp=parser.add_argument_group(title="Required arguments")

req_grp.add_argument("--sample_size","-s",dest="ss",help="diploid sample size within each deme",type=int,default = 2000, nargs = "?")
req_grp.add_argument("--outpre","-o",dest="outpre",help="output file prefix",type=str,required=True)
req_grp.add_argument("--chr","-c",dest="chr",help="chromosome number",type=str,required=True)
parser.add_argument("--length","-L",dest="length",help="length of chromosome (bp) (def:1e6)",type=int,default=1000000,nargs="?")
parser.add_argument("--rho","-r",dest="rho",help="recombination rate (def:1e-08)",type=float,default=1e-08,nargs="?")
parser.add_argument("--mu","-u",dest="mu",help="mutation rate (def:1e-08)",type=float,default=1e-08,nargs="?")
parser.add_argument("--migrate","-m",dest="migrate",help="migration rate (def:0.05)",type=float,default=0.05,nargs="?")
args=parser.parse_args()

print(args)

import msprime
import numpy as np

#random_seed = int(sys.argv[1]) # input the seed
#random_seed = 1

mu=args.mu # mutation rate per bp
rho=args.rho # recombination rate per bp
L=args.length # sequence length
n = args.ss
chr=args.chr
outpre=args.outpre

ABsplit=2040 # number of generations back to Out of Africa
Na = 1e4
Nb = 1e4 #gives an Fst of 0.2
N0 = 1e4

Tadmix=20 # time of admixture
Nadmix=1e4 # initial size of admixed population

mdeme=args.migrate # migration rate between neighboring demes, one direction

demography = msprime.Demography()
demography.add_population(name = "A", initial_size= Na)
demography.add_population(name = "B", initial_size= Nb)
demography.add_population(name = "Anc", initial_size = N0)

for i in range(10):
    demography.add_population(name = "adm"+str(i), initial_size=Nadmix)


for i in np.arange(4,13):
    demography.migration_matrix[i,i-1] = mdeme
    demography.migration_matrix[i-1,i] = mdeme


anc = np.arange(0.05, 1, 0.1)
for i in range(10):
    demography.add_admixture(time=Tadmix, derived="adm"+str(i), ancestral=["A", "B"], proportions=[anc[i], 1 - anc[i]])

demography.add_population_split(time=ABsplit, derived=["A", "B"], ancestral="Anc")


# Initializing an empty dictionary
samples = {}

# Using a for loop to populate the dictionary
for i in range(10):
    samples["adm"+str(i)] = n/10


ts = msprime.sim_ancestry(samples=samples, sequence_length = L, recombination_rate=rho, demography=demography)

ts = msprime.sim_mutations(ts, rate = mu)


with open(outpre + "_chr" + chr + ".vcf","w") as vcf_file:
    ts.write_vcf(vcf_file, contig_id=chr)