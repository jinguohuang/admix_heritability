# This script is for simulating genotype of admixed population under random mating
# update: update the local ancestry inference function: loop for samples so it won't be killed 

# variabiles: 
import argparse

parser=argparse.ArgumentParser()
req_grp=parser.add_argument_group(title="Required arguments")

# required input: output filename prefix, chromosome number(seed number) 
req_grp.add_argument("--outpre","-o",dest="outpre",help="output file prefix",type=str,required=True)
req_grp.add_argument("--chr","-c",dest="chr",help="chromosome number",type=str,required=True)

#optional input, with default setting

parser.add_argument("--mu","-u",dest="mu",help="mutation rate (def:1.25e-8)",type=float,default=1.25e-8,nargs="?")
parser.add_argument("--rho","-r",dest="rho",help="recombination rate (def:1e-08)",type=float,default=1e-08,nargs="?")
parser.add_argument("--length","-L",dest="length",help="length of chromosome (bp) (def:1e7)",type=int,default=10000000,nargs="?")
#parser.add_argument("--model",dest="model",help="population structure (either random or 10deme)",type=str,default="uniform",nargs="?")
parser.add_argument("--tadmix","-t",dest="tadmix",help="time of admixture",type=int,default=16,nargs="?")
parser.add_argument("--tcensus",dest="tcensus",help="time of census event",type=int,default=17,nargs="?")
parser.add_argument("--nadmix",dest="nadmix",help="initial size of admixed population",type=int,default=30000,nargs="?")
parser.add_argument("--radmix",dest="radmix",help="growth rate of admixed population",type=float,default=.05,nargs="?")
parser.add_argument("--padmix","-p",dest="padmix",help="proportion of admixture",type=float,default=.76517,nargs="?")
#parser.add_argument("--mdeme","-m",dest="mdeme",help="one directional migration rate between demes(def:0.025)",type=float,default=0.025,nargs="?")
parser.add_argument("--refsamplesize","-ref",dest="refsamplesize",help="number of reference sample",type=int,default=100,nargs="?")
parser.add_argument("--admsamplesize","-adm",dest="admsamplesize",help="number of admixed sample",type=int,default=400,nargs="?")

args=parser.parse_args()

print(args)

# simulate genotype

import msprime, sys
from math import log
from math import exp

#seed=1
#mu=1.25e-8 # mutation rate per bp
#rho=1e-8 # recombination rate per bp
#nbp=1e8 # generate 100Mb
#nbp=1e7 # generate 10Mb
N0=7310 # initial population size
Thum=5920 # time (gens) of advent of modern humans
Naf=14474 # size of african population
Tooa=2040 # number of generations back to Out of Africa
Nb=1861 # size of out of africa population
mafb=1.5e-4 # migration rate Africa and Out-of-Afica
Teu=920 # number generations back to Asia-Euroupe split
Neu=1032; Nas=554 # bottleneck population sizes
mafeu=2.5e-5;mafas=7.8e-6;meuas=3.11e-5 #mig.rates
reu=0.0038 #growth rate per generation in Europe
ras=0.0048 #growth rate per generation in Asia
#Tadmix=16 # time of admixture
#Nadmix=30000 # initial size of admixed population
#radmix=.05 # growth rate of admixed population
#Tcensus=17

# pop0 is Africa, pop1 is Europe,  pop2 is admixed
#refsamplesize=100
#admsamplesize=1000 #sample size of admixted pop
#admsamplesize=400
pop_config = [
    msprime.PopulationConfiguration(sample_size=args.refsamplesize, initial_size=Naf), #growth_rate=0.0),
    msprime.PopulationConfiguration(sample_size=args.refsamplesize, initial_size=Neu*exp(reu*Teu)), #growth_rate=reu),
    msprime.PopulationConfiguration(sample_size=args.admsamplesize, initial_size=args.nadmix*exp(args.radmix*args.tadmix), )]#growth_rate=radmix)]
mig_mat=[
    [0,mafeu,0],
    [mafeu,0,0],
    [0,0,0]]

# Admixture event, 0.76517 Africa, rest Europe
# change the AFR prop according to the 10deme model avg anc
admixture_event = [
    msprime.MigrationRateChange(time=args.tadmix, rate=0.0), #change the migration rate to 0
    msprime.MassMigration(time=args.tadmix, source=2, destination=0, proportion=args.padmix), # African came to America        
    msprime.MassMigration(time=args.tadmix+0.0001, source=2, destination=1, proportion=1.0)] # European came to America

# census events
census_event = [
    msprime.CensusEvent(time=args.tcensus)
]

# Out of Africa event
ooa_event = [
    msprime.MigrationRateChange(time=Tooa, rate=0.0),
    msprime.MassMigration(time=Tooa+0.0001, source=1, destination=0, proportion=1.0)
]

# initial population size
init_event=[msprime.PopulationParametersChange(time=Thum, initial_size=N0, population_id=0)]

#events = admixture_event + eu_event + ooa_event + init_event
events = admixture_event + census_event + ooa_event + init_event

ts = msprime.simulate(
    population_configurations=pop_config,
    migration_matrix=mig_mat,
    demographic_events=events,
    length=args.length,
    recombination_rate=args.rho, 
    mutation_rate=args.mu,
    random_seed=args.chr)

#print("writing genotype to vcf file")
#export to diploid to validate simulated global anc
filename=f'{args.outpre}_{args.padmix}prop_chr{args.chr}.vcf'
#n_dip_indv = int(ts.num_samples / 2)
#indv_names = [f"AA_{str(i)}indv" for i in range(n_dip_indv)]
#with open(filename, "w") as vcf_file:
    #ts.write_vcf(vcf_file, ploidy=2, contig_id=args.chr, individual_names=indv_names)


# output the local ancestry chunk of each person
import tskit
import numpy as np
import msprime 

# replace the parent to population label
def get_population_id(node, ts):
    return ts.tables.nodes.population[node]

def local_ancestry(Tcensus, ts, filename):
    # select nodes that are Tcensus generations ago
    ancestors_Tcensus_gens = np.where(ts.tables.nodes.time == Tcensus)[0]
    # which of your samples descend from which ancestors
    #ancestrytable_Tcensus_gens = ts.tables.link_ancestors(
    #samples=ts.samples(), ancestors=ancestors_Tcensus_gens)
    ### loop on a small set 
    a_list = ts.samples() # samples list
    length = len(a_list) #total number of sample size
    lu=100 #loop unit:loop 100 at a time
    ln=length//lu +1 #loop number
    x=0  #initial number for the start index
    for i in range(1,ln):
        mysample = a_list[x:(x+lu)]
        x=lu*i #new start index
        ancestrytable_Tcensus_gens = ts.tables.link_ancestors(
        samples=mysample, ancestors=ancestors_Tcensus_gens)
        import pandas as pd
        ancestry_table = pd.DataFrame(
            data = {
                'left': ancestrytable_Tcensus_gens.left,
                'right': ancestrytable_Tcensus_gens.right,
                'populations': [get_population_id(u, ts) for u in ancestrytable_Tcensus_gens.parent],
                'child': ancestrytable_Tcensus_gens.child
            }
        )
        #return ancestry_table
        # output the ancestry table to a file, append mode so it will accumulate for each loop
        ancestry_table.to_csv(f'loc_anc_{filename}.txt', sep='\t', encoding='utf-8', index=False, header=False, mode='a')

# run the function
local_ancestry(args.tcensus, ts, filename)



