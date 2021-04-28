#simulate 400hap 100MB genotype and get local ancestry for each person
import msprime, sys
from math import log
from math import exp

# simulated genotype model
#seed = int(sys.argv[1]) # input the seed
seed=1
mu=1.25e-8 # mutation rate per bp
rho=1e-8 # recombination rate per bp
nbp=1e8 # generate 100Mb
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
Tadmix=16 # time of admixture
Nadmix=30000 # initial size of admixed population
radmix=.05 # growth rate of admixed population
Tcencus=17
mdeme=0.05/2 # migration rate between neighboring demes

# pop0 is Africa, pop1 is Europe,  pop2-10 is admixed
refsamplesize=100 #change to 100 haplotype
admsamplesize=400 #sample size of admixted pop
# sample size of each admixed pop, 400 in total
adm1=6
adm2=0
adm3=3
adm4=6
adm5=3
adm6=18
adm7=52
adm8=118
adm9=149
adm10=45
pop_config = [
    msprime.PopulationConfiguration(sample_size=refsamplesize, initial_size=Naf), #AFR
    msprime.PopulationConfiguration(sample_size=refsamplesize, initial_size=Neu*exp(reu*Teu)), #EUR
    msprime.PopulationConfiguration(sample_size=adm1, initial_size=Nadmix*exp(radmix*Tadmix)), #adm1
    msprime.PopulationConfiguration(sample_size=adm2, initial_size=Nadmix*exp(radmix*Tadmix)), #adm2
    msprime.PopulationConfiguration(sample_size=adm3, initial_size=Nadmix*exp(radmix*Tadmix)), #adm3
    msprime.PopulationConfiguration(sample_size=adm4, initial_size=Nadmix*exp(radmix*Tadmix)), #adm4
    msprime.PopulationConfiguration(sample_size=adm5, initial_size=Nadmix*exp(radmix*Tadmix)), #adm5
    msprime.PopulationConfiguration(sample_size=adm6, initial_size=Nadmix*exp(radmix*Tadmix)), #adm6
    msprime.PopulationConfiguration(sample_size=adm7, initial_size=Nadmix*exp(radmix*Tadmix)), #adm7
    msprime.PopulationConfiguration(sample_size=adm8, initial_size=Nadmix*exp(radmix*Tadmix)), #adm8
    msprime.PopulationConfiguration(sample_size=adm9, initial_size=Nadmix*exp(radmix*Tadmix)), #adm9
    msprime.PopulationConfiguration(sample_size=adm10, initial_size=Nadmix*exp(radmix*Tadmix)), #adm10 
]
# adjust the edge mdeme
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

# cencus events
census_event = [
    msprime.CensusEvent(time=Tcencus)
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
    length=nbp,
    recombination_rate=rho, 
    mutation_rate=mu,
    #record_migrations=True, # newly added, Needed for tracking segments.
    random_seed=seed)

#export to diploid to validate simulated global anc
filename=f'AA_100Mb_10deme_{admsamplesize}hap_chr{seed}.vcf'
n_dip_indv = int(ts.num_samples / 2)
indv_names = [f"AA_{str(i)}indv" for i in range(n_dip_indv)]
with open(filename, "w") as vcf_file:
    ts.write_vcf(vcf_file, ploidy=2, contig_id=str(seed), individual_names=indv_names)

# get local ancestry 
import tskit
import numpy as np
import msprime 

# replace the parent to population label
def get_population_id(node, ts):
    return ts.tables.nodes.population[node]

def local_ancestry(Tcencus, ts):
    # select nodes that are Tcencus generations ago
    ancestors_Tcencus_gens = np.where(ts.tables.nodes.time == Tcencus)[0]
    # which of your samples descend from which ancestors
    ancestrytable_Tcencus_gens = ts.tables.link_ancestors(
    samples=ts.samples(), ancestors=ancestors_Tcencus_gens)
    # replace the parent to population label
    #nodeTable = ts.tables.nodes
    import pandas as pd
    ancestry_table = pd.DataFrame(
        data = {
            'left': ancestrytable_Tcencus_gens.left,
            'right': ancestrytable_Tcencus_gens.right,
            'populations': [get_population_id(u, ts) for u in ancestrytable_Tcencus_gens.parent],
            'child': ancestrytable_Tcencus_gens.child
        }
    )
    #return ancestry_table
    # output the ancestry table to a file
    ancestry_table.to_csv(f'loc_anc_{filename}.txt', sep='\t', encoding='utf-8', index=False)

# run the function
local_ancestry(Tcencus, ts)

