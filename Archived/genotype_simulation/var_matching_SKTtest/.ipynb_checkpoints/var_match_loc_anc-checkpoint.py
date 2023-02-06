# This script is for simulating genotype of admixed population under different model
# changes: incorporate 2 models in one script, one can call by either random or 10deme after -model
## changed the random mating admixture to 0.77
## changed the 10 deme model to calculate the ppl in each deme automatically according to input adm sample size
## changed the initial pop size for 10deme to 3k

# variables: 
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
parser.add_argument("--model", "-model", dest="model",help="population structure (either random or 10deme)",type=str,default="random",nargs="?")
parser.add_argument("--tadmix","-t",dest="tadmix",help="time of admixture",type=int,default=16,nargs="?")
parser.add_argument("--tcensus",dest="tcensus",help="time of census event",type=int,default=17,nargs="?")
parser.add_argument("--nadmix",dest="nadmix",help="initial size of admixed population",type=int,default=30000,nargs="?")
parser.add_argument("--radmix",dest="radmix",help="growth rate of admixed population",type=float,default=.05,nargs="?")
parser.add_argument("--padmix","-p",dest="padmix",help="proportion of admixture",type=float,default=.77,nargs="?")
parser.add_argument("--mdeme","-m",dest="mdeme",help="one directional migration rate between demes(def:0.025)",type=float,default=0.025,nargs="?")
parser.add_argument("--refsamplesize","-ref",dest="refsamplesize",help="number of reference sample",type=int,default=100,nargs="?")
parser.add_argument("--admsamplesize","-adm",dest="admsamplesize",help="number of admixed sample",type=int,default=400,nargs="?")

args=parser.parse_args()

print(args)

# simulate genotype

import msprime, sys
from math import log
from math import exp

### define some parameters
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


##### define function to simulate genotypes under a random mating migration model
def random_mating():

    # pop0 is Africa, pop1 is Europe,  pop2 is admixed
    pop_config = [
        msprime.PopulationConfiguration(sample_size=args.refsamplesize, initial_size=Naf), #growth_rate=0.0),
        msprime.PopulationConfiguration(sample_size=args.refsamplesize, initial_size=Neu*exp(reu*Teu)), #growth_rate=reu),
        msprime.PopulationConfiguration(sample_size=args.admsamplesize, initial_size=args.nadmix*exp(args.radmix*args.tadmix), )]#growth_rate=radmix)]
    mig_mat=[
        [0,mafeu,0],
        [mafeu,0,0],
        [0,0,0]]

    # Admixture event, 0.77 Africa, rest Europe
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
    return(ts)

##### define function to simulate genotypes under a 10deme migration model
def deme10():
    mdeme=args.mdeme # migration rate between neighboring demes
    # pop0 is Africa, pop1 is Europe,  pop2-10 is admixed
    # sample size of each admixed pop
    n = args.admsamplesize
    adm1 = int(n*0.01)
    adm2 = int(n*0.01)
    adm3 = int(n*0.01)
    adm4 = int(n*0.01)
    adm5 = int(n*0.01)
    adm6 = int(n*0.03)
    adm7 = int(n*0.13)
    adm8 = int(n*0.26)
    adm9 = int(n*0.37)
    adm10 = int(n*0.16)
    Nadmix = int(args.nadmix/10) #initial population size for each deme
    #print("inital size for each deme is:" + str(Nadmix))

    # population configuration 
    pop_config = [
        msprime.PopulationConfiguration(sample_size=args.refsamplesize, initial_size=Naf), #AFR
        msprime.PopulationConfiguration(sample_size=args.refsamplesize, initial_size=Neu*exp(reu*Teu)), #EUR
        msprime.PopulationConfiguration(sample_size=adm1, initial_size=Nadmix*exp(args.radmix*args.tadmix)), #adm1
        msprime.PopulationConfiguration(sample_size=adm2, initial_size=Nadmix*exp(args.radmix*args.tadmix)), #adm2
        msprime.PopulationConfiguration(sample_size=adm3, initial_size=Nadmix*exp(args.radmix*args.tadmix)), #adm3
        msprime.PopulationConfiguration(sample_size=adm4, initial_size=Nadmix*exp(args.radmix*args.tadmix)), #adm4
        msprime.PopulationConfiguration(sample_size=adm5, initial_size=Nadmix*exp(args.radmix*args.tadmix)), #adm5
        msprime.PopulationConfiguration(sample_size=adm6, initial_size=Nadmix*exp(args.radmix*args.tadmix)), #adm6
        msprime.PopulationConfiguration(sample_size=adm7, initial_size=Nadmix*exp(args.radmix*args.tadmix)), #adm7
        msprime.PopulationConfiguration(sample_size=adm8, initial_size=Nadmix*exp(args.radmix*args.tadmix)), #adm8
        msprime.PopulationConfiguration(sample_size=adm9, initial_size=Nadmix*exp(args.radmix*args.tadmix)), #adm9
        msprime.PopulationConfiguration(sample_size=adm10, initial_size=Nadmix*exp(args.radmix*args.tadmix)), #adm10 
    ]
    
    # adjust the edge mdeme in the migration matrix
    mig_mat=[
        [    0,mafeu,0,0,0,0,0,0,0,0,0,0], #AFR
        [mafeu,0,    0,0,0,0,0,0,0,0,0,0], #EUR
        [    0,0,0,2*mdeme,0,0,0,0,0,0,0,0], #ADM1
        [0,0,mdeme,0,mdeme,0,0,0,0,0,0,0], #ADM2
        [0,0,0,mdeme,0,mdeme,0,0,0,0,0,0], #ADM3
        [0,0,0,0,mdeme,0,mdeme,0,0,0,0,0], #ADM4
        [0,0,0,0,0,mdeme,0,mdeme,0,0,0,0], #ADM5
        [0,0,0,0,0,0,mdeme,0,mdeme,0,0,0], #ADM6
        [0,0,0,0,0,0,0,mdeme,0,mdeme,0,0], #ADM7
        [0,0,0,0,0,0,0,0,mdeme,0,mdeme,0], #ADM8
        [0,0,0,0,0,0,0,0,0,mdeme,0,mdeme], #ADM9
        [0,0,0,0,0,0,0,0,0,0,2*mdeme,0], #ADM10
    ]

    # Admixture event, 0.8 Africa, 0.2 Europe in total
    # 0 Afr, 1 Eur, source 2-12 admixed 
    Tadmix = args.tadmix # assign the time of admixure
    admixture_event = [
        msprime.MigrationRateChange(time=Tadmix, rate=0.0), #newly added, change the migration rate to 0
        msprime.MassMigration(time=Tadmix, source=2, destination=0, proportion=0.05), # Afr 0.05
        msprime.MassMigration(time=Tadmix+0.0001, source=2, destination=1, proportion=1.0), # Eur
        msprime.MassMigration(time=Tadmix, source=3, destination=0, proportion=0.15), # Afr 0.15
        msprime.MassMigration(time=Tadmix+0.0001, source=3, destination=1, proportion=1.0), # Eur
        msprime.MassMigration(time=Tadmix, source=4, destination=0, proportion=0.25), # Afr 0.25
        msprime.MassMigration(time=Tadmix+0.0001, source=4, destination=1, proportion=1.0), # Eur
        msprime.MassMigration(time=Tadmix, source=5, destination=0, proportion=0.35), # Afr 0.35
        msprime.MassMigration(time=Tadmix+0.0001, source=5, destination=1, proportion=1.0), # Eur
        msprime.MassMigration(time=Tadmix, source=6, destination=0, proportion=0.45), # Afr 0.45
        msprime.MassMigration(time=Tadmix+0.0001, source=6, destination=1, proportion=1.0), # Eur
        msprime.MassMigration(time=Tadmix, source=7, destination=0, proportion=0.55), # Afr 0.55
        msprime.MassMigration(time=Tadmix+0.0001, source=7, destination=1, proportion=1.0), # Eur
        msprime.MassMigration(time=Tadmix, source=8, destination=0, proportion=0.65), # Afr 0.65
        msprime.MassMigration(time=Tadmix+0.0001, source=8, destination=1, proportion=1.0), # Eur
        msprime.MassMigration(time=Tadmix, source=9, destination=0, proportion=0.75), # Afr 0.75
        msprime.MassMigration(time=Tadmix+0.0001, source=9, destination=1, proportion=1.0), # Eur
        msprime.MassMigration(time=Tadmix, source=10, destination=0, proportion=0.85), # Afr 0.85
        msprime.MassMigration(time=Tadmix+0.0001, source=10, destination=1, proportion=1.0), # Eur
        msprime.MassMigration(time=Tadmix, source=11, destination=0, proportion=0.95), # Afr 0.95
        msprime.MassMigration(time=Tadmix+0.0001, source=11, destination=1, proportion=1.0), # Eur                   
    ] 

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
    return(ts)


###simulate according to the model
model=args.model
print("model is:" + model)

if model=="random":
    ts=random_mating()

elif model=="10deme":
    ts=deme10()

else:
    print("Please type either 'random' or '10deme' after -model")
    sys.exit() #exit if not right model


print("writing genotype to vcf file")
#export to diploid to validate simulated global anc
filename=f'{args.outpre}_{args.model}_chr{args.chr}.vcf'
n_dip_indv = int(ts.num_samples / 2)
indv_names = [f"AA_{str(i)}indv" for i in range(n_dip_indv)]
with open(filename, "w") as vcf_file:
    ts.write_vcf(vcf_file, ploidy=2, contig_id=args.chr, individual_names=indv_names)


##### output the local ancestry chunk of each person
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
    ancestrytable_Tcensus_gens = ts.tables.link_ancestors(
    samples=ts.samples(), ancestors=ancestors_Tcensus_gens)
    # replace the parent to population label
    #nodeTable = ts.tables.nodes
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
    # output the ancestry table to a file
    ancestry_table.to_csv(f'loc_anc_{filename}.txt', sep='\t', encoding='utf-8', index=False)

# run the function
print("writing local ancestry file")
local_ancestry(args.tcensus, ts, filename)


